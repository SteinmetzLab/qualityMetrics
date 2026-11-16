"""HTML pages that show a folder of QC figures in reading order.

One page per shank (``index.html`` beside its figures) and one per session
(``index.html`` beside the session figures, linking to each shank's page).
Both are plain files with relative links, so they work opened straight from
the lab share, with no server (Nick, 2026-11-16).

The figures stay the record. A page only arranges what is in the folder when
it is written; a figure that is missing is listed as missing, with the reason
where the report recorded one, rather than silently left out.

The shank page also redraws the template panels in the browser: the
templates of every unit are embedded in the page, so a button can show six
more units, or units typed by number, without Python. That is the limit of
interactivity in a file opened from the share: whatever a control shows must
already be in the page. Anything that needs the raw recording would need a
server (SortingManager's dashboard is the place for that).
"""
from __future__ import annotations

import base64
import html
import json
import os
import re
from pathlib import Path

import numpy as np

PAGE_NAME = "index.html"

#: The shank figures in reading order, grouped, with what each shows. A
#: figure in the folder that is not listed here goes under "Other figures";
#: one listed here and not in the folder is named under "Not drawn".
SHANK_SECTIONS: tuple[tuple[str, tuple[tuple[str, str, str], ...]], ...] = (
    ("Units", (
        ("amp_depth_scatter", "Unit amplitude against depth",
         "Each unit's amplitude and depth, colored by firing rate and by "
         "waveform duration. Large markers passed QC."),
        ("templates_grid", "Templates",
         "Two units at each extreme (rate, amplitude, footprint, width), then "
         "six at random, each on its 32 nearest sites drawn at their positions."),
        ("example_neurons", "Example neurons",
         "Three QC-passing units across the amplitude range: raw snippets, "
         "template over depth, autocorrelogram."),
        ("noise_cutoff_diagnostic", "Noise cutoff",
         "How the noise cutoff metric reaches its number: amplitude over the "
         "session and the amplitude histogram it works on."),
        ("depth_profiles", "Profiles down the probe",
         "Firing rate, spike amplitude and unit count against depth."),
        ("firing_rate_image", "Firing rate over time and depth",
         "Binned rate image; shows rate changes inside dense bands."),
        ("depth_correlation", "Rate correlation between depths",
         "Correlated blocks along the diagonal suggest one structure; their "
         "edges suggest boundaries."),
    )),
    ("Drift", (
        ("dredge_driftmap", "DREDGE motion estimate",
         "The peaks DREDGE detected, with its displacement estimate on top."),
        ("drift_map", "Drift map, every unit",
         "Every spike's depth over time, darker for larger amplitude."),
        ("drift_map_qc-passing", "Drift map, QC-passing units",
         "As above, only units that passed QC, same axes and shading."),
        ("unit_drift", "Spike depth by unit, every unit",
         "A depth band that changes color partway is one neuron split by "
         "the sorter."),
        ("unit_drift_qc-passing", "Spike depth by unit, QC-passing units",
         "As above, only units that passed QC."),
        ("tracking_continuity", "Firing continuity",
         "Where firing stops in one unit and starts in another at the same "
         "depth: a sign of a drift split."),
        ("amplitude_cdf_over_depth", "Amplitude distribution over depth",
         "The rate of spikes at each amplitude and depth."),
    )),
    ("Signal", (
        ("raw_traces", "Traces, every channel",
         "One second of every channel as the sorting pipeline's destripe "
         "left it, shallow at the bottom."),
        ("raw_with_spikes", "Traces with the sorter's spikes",
         "Raw traces with each unit's spikes marked at its depth."),
        ("wall_heatmap", "Voltage across the shank",
         "Depth by time at three moments, as stored and as sorted."),
        ("filter_state", "Median spectrum per processing stage",
         "The session noise figure's stages, for this shank alone."),
        ("depth_power_wide", "Power by depth, whole band",
         "As stored, on a logarithmic frequency axis."),
        ("depth_power_wide_destriped", "Power by depth, after the destripe",
         "The same seconds and channels as stored and after the sorting "
         "pipeline's destripe, one color scale."),
        ("channel_health", "Channel health",
         "RMS per channel against depth, with bad sites flagged."),
        ("band_rms", "RMS over time and depth",
         "AP and LFP band RMS in blocks of time and depth."),
        ("saturation", "Saturation",
         "Samples at the ADC's rails, per channel and over the session."),
    )),
    ("LFP", (
        ("depth_power", "LFP power by depth",
         "Depth by frequency below 300 Hz, normalized per frequency."),
        ("lfp_band_profiles", "LFP band profiles",
         "Power in five named bands down the probe."),
    )),
    ("Sorter", (
        ("kilosort_diagnostics", "Kilosort diagnostics",
         "The diagnostics figure Kilosort drew itself."),
    )),
)

SESSION_FIGURES: tuple[tuple[str, str, str], ...] = (
    ("session_raster_all_units", "Every sorted unit",
     "Rates over the session for every unit the sorter produced."),
    ("session_raster_passing", "Units passing QC",
     "The same, for the units that passed QC."),
    ("session_rms_across_shanks", "Noise on every shank",
     "AP band RMS per site, every shank."),
    ("session_noise_rms_summary", "Median RMS per shank",
     "With a 95% bootstrap confidence interval."),
    ("session_noise_rms_over_time", "RMS over the recording",
     "After demux CAR, median across channels."),
    ("session_noise_spectra", "Median spectrum per processing stage",
     "Every shank, every stage."),
    ("session_noise_depth_power", "Power by depth per processing stage",
     "Every shank, every stage."),
    ("session_noise_rms_by_channel", "RMS per channel per stage",
     "Every shank, every stage."),
    ("session_noise_snippets", "Voltage snippets per stage",
     "33 ms of voltage from the middle window."),
)

#: Why each figure was not drawn, kept beside the figures when a report
#: writes its page, so a page written again later (``sm redo-qc
#: --pages-only``) still gives the reasons (review, 2026-11-16).
REASONS_NAME = "not_drawn.json"
#: Files a page never lists.
_IGNORED = {PAGE_NAME, REASONS_NAME, "Thumbs.db", "desktop.ini", ".DS_Store"}
#: "p1_g0_imec0_shank2", or "imec0_shank2" from before the run was in the name.
SHANK_FOLDER = re.compile(r"^(?:(?P<run>.+?)_)?imec(?P<probe>\d+)_shank(?P<shank>\d+)$")


# ------------------------------------------------------------------ writing

def write_shank_page(out_dir: str | Path, *, title: str, ks=None,
                     skipped: dict | None = None, exclude_units=(),
                     session_link: bool = True) -> Path:
    """Write ``index.html`` for one shank's figure folder and return its path.

    ``ks`` (:class:`~qualitymetrics.ksdata.KilosortResults`) embeds the
    templates for the redraw control beside templates_grid; without it the
    page has no control. ``exclude_units`` are the units the static grid
    already shows, left out of the first redraws. ``skipped`` is the report's
    ``{figure: reason}``, kept in ``not_drawn.json``; None reads the reasons
    kept there by an earlier report.
    """
    out_dir = Path(out_dir)
    skipped = _reasons(out_dir, skipped)
    files = _files(out_dir)
    pngs = {f[:-4] for f in files if f.lower().endswith(".png")}
    listed = {name for _s, figs in SHANK_SECTIONS for name, *_ in figs}

    body, toc = [], []
    for section, figures in SHANK_SECTIONS:
        present = [f for f in figures if f[0] in pngs]
        if not present:
            continue
        toc.append(f'<a href="#{_slug(section)}">{_e(section)}</a>')
        body.append(f'<h2 id="{_slug(section)}">{_e(section)}</h2>')
        for name, heading, what in present:
            extra = ""
            if name == "templates_grid" and ks is not None:
                try:
                    extra = _template_control(ks, exclude_units)
                except Exception as exc:  # noqa: BLE001 - the figure still shows
                    extra = (f'<p class="note">No redraw control: '
                             f'{_e(type(exc).__name__)}: {_e(exc)}</p>')
            body.append(_figure(name, heading, what, extra, out_dir))
    others = sorted(p for p in pngs if p not in listed)
    if others:
        toc.append('<a href="#other-figures">Other figures</a>')
        body.append('<h2 id="other-figures">Other figures</h2>')
        body.extend(_figure(name, name.replace("_", " ").capitalize(), "",
                            out_dir=out_dir) for name in others)
    body.append(_not_drawn(_missing(listed, pngs, skipped)))
    body.append(_other_files(files))

    # Not checked: the shank's figures are drawn on the worker, away from the
    # session folder, and the session page is written only once every shank
    # of the session is sorted. The tooltip says so.
    nav = (['<a href="../index.html" title="The session page, written once every '
            'shank of the session is sorted">Session page</a>'] if session_link else [])
    return _write(out_dir, title, nav + toc, body,
                  script=TEMPLATE_SCRIPT if "tpl-data" in "".join(body) else "")


def write_session_page(out_dir: str | Path, *, title: str,
                       shank_folders=(), skipped: dict | None = None) -> Path:
    """Write ``index.html`` for a session's figure folder and return its path.

    Links every shank folder named in ``shank_folders`` and every subfolder
    named like one or holding a shank page, by probe and shank. A shank
    folder with no page yet (figures drawn by older code) links to the folder
    itself. ``skipped`` as for :func:`write_shank_page`.
    """
    out_dir = Path(out_dir)
    skipped = _reasons(out_dir, skipped)
    files = _files(out_dir)
    pngs = {f[:-4] for f in files if f.lower().endswith(".png")}
    folders = set(shank_folders)
    try:
        folders.update(p.name for p in out_dir.iterdir()
                       if p.is_dir() and (SHANK_FOLDER.match(p.name)
                                          or (p / PAGE_NAME).exists()))
    except OSError:
        pass
    folders = {f for f in folders if not f.startswith((".", "_"))}

    body, toc = [], []
    if folders:
        toc.append('<a href="#shanks">Shanks</a>')
        body.append('<h2 id="shanks">Shanks</h2>')
        body.append(_shank_table(out_dir, folders))
    listed = {name for name, *_ in SESSION_FIGURES}
    present = [f for f in SESSION_FIGURES if f[0] in pngs]
    if present:
        toc.append('<a href="#figures">Session figures</a>')
        body.append('<h2 id="figures">Session figures</h2>')
        body.extend(_figure(*f, out_dir=out_dir) for f in present)
    others = sorted(p for p in pngs if p not in listed)
    if others:
        toc.append('<a href="#other-figures">Other figures</a>')
        body.append('<h2 id="other-figures">Other figures</h2>')
        body.extend(_figure(n, n.replace("_", " ").capitalize(), "", out_dir=out_dir)
                    for n in others)
    body.append(_not_drawn(_missing(listed, pngs, skipped)))
    body.append(_other_files(files))
    return _write(out_dir, title, toc, body)


# ------------------------------------------------------------ the templates

def template_data(ks, n_channels: int = 32) -> dict:
    """Every unit's template on its nearest sites, for the page's redraw.

    Each unit's waveforms are stored as int8 scaled to that unit's largest
    excursion, which is how the panels draw them anyway (each scaled to its
    own peak): 32 sites x 61 samples is about 200 kB per 100 units, 270 kB
    after base64. A template with no finite value draws as flat lines and is
    named in the panel (``bad``).
    """
    from .plots.rawviews import _site_spacing
    from .plots.units import template_sites

    unw = ks.templates_unwhitened
    pos = np.asarray(ks.channel_positions, dtype=float)
    units = [int(u) for u in ks.unit_ids if int(u) < len(unw)]
    if not units:
        raise ValueError("no units with templates")
    nc = min(n_channels, len(pos))
    chans = np.zeros((len(units), nc), dtype=np.int32)
    waves = np.zeros((len(units), unw.shape[1], nc), dtype=np.int8)
    bad = []
    for i, u in enumerate(units):
        c, w = template_sites(ks, u, nc)
        chans[i] = c
        finite = np.isfinite(w)
        if not finite.all():
            bad.append(u)
        w = np.where(finite, w, 0.0)
        peak = float(np.abs(w).max()) or 1.0
        waves[i] = np.round(w / peak * 127.0).astype(np.int8)
    dx, dy = _site_spacing(pos)
    amps = ks.unit_amplitude_uv if ks.uv_per_bit else {}
    return {
        "units": units,
        "passing": [bool(ks.qc_pass.get(u, False)) for u in units],
        "amp_uv": [_finite(amps.get(u)) for u in units],
        "peak_site": [int(ks.peak_channel[u]) for u in units],
        "chans": chans.ravel().tolist(),
        "waves": base64.b64encode(waves.tobytes()).decode("ascii"),
        "n_sites": nc,
        "n_samples": int(unw.shape[1]),
        "fs": float(ks.fs),
        "positions": np.round(pos, 3).tolist(),
        "dx": float(dx),
        "dy": float(dy),
        "bad": bad,
    }


def _template_control(ks, exclude_units) -> str:
    data = template_data(ks)
    data["exclude"] = sorted(int(u) for u in exclude_units)
    n_pool = sum(data["passing"]) or len(data["units"])
    what = "QC-passing units" if any(data["passing"]) else "units"
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    return (
        '<div class="tpl">'
        '<div class="tpl-controls">'
        f'<button type="button" id="tpl-next">Draw 6 more {what}</button>'
        '<label>Units <input id="tpl-ids" type="text" size="18" '
        'placeholder="e.g. 12, 40, 117"></label>'
        '<button type="button" id="tpl-show">Show</button>'
        f'<span id="tpl-status" class="note">Drawn in the browser from the '
        f'templates of all {len(data["units"])} units, from {n_pool} {what} at '
        f'random, without repeats until all have been shown.</span>'
        '</div>'
        '<div id="tpl-panels"></div>'
        f'<script type="application/json" id="tpl-data">{payload}</script>'
        '</div>')


# ------------------------------------------------------------------ pieces

def _figure(name: str, heading: str, what: str, extra: str = "",
            out_dir: Path | None = None) -> str:
    src = _href(f"{name}.png")
    # The image's own size, so the page reserves its space before it loads:
    # without it, lazily loaded figures above a link's target push the
    # target down after the jump.
    size = _png_size(out_dir / f"{name}.png") if out_dir is not None else None
    dims = f' width="{size[0]}" height="{size[1]}"' if size else ""
    return (f'<section class="fig" id="fig-{_slug(name)}">'
            f'<h3>{_e(heading)} <span class="name">{_e(name)}</span></h3>'
            + (f'<p>{_e(what)}</p>' if what else "")
            + f'<a href="{src}"><img src="{src}"{dims} loading="lazy" '
              f'alt="{_e(heading)}"></a>'
            + extra + '</section>')


def _png_size(path: Path) -> tuple[int, int] | None:
    """Width and height from a PNG's header, or None."""
    try:
        with open(path, "rb") as f:
            head = f.read(24)
    except OSError:
        return None
    if len(head) < 24 or head[:8] != b"\x89PNG\r\n\x1a\n" or head[12:16] != b"IHDR":
        return None
    return int.from_bytes(head[16:20], "big"), int.from_bytes(head[20:24], "big")


def _reasons(out_dir: Path, skipped: dict | None) -> dict:
    """The report's reasons, kept beside the figures; or, given None, those
    kept by an earlier report."""
    path = out_dir / REASONS_NAME
    if skipped is None:
        try:
            kept = json.loads(path.read_text(encoding="utf-8"))
            return {str(k): str(v) for k, v in kept.items()} if isinstance(kept, dict) else {}
        except (OSError, ValueError):
            return {}
    skipped = {str(k): str(v) for k, v in skipped.items()}
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(skipped, indent=1, sort_keys=True),
                        encoding="utf-8", newline="\n")
    except OSError:
        pass                    # the page still says why; only a rewrite would not
    return skipped


def _missing(listed: set, pngs: set, skipped: dict) -> dict:
    """Every expected figure not in the folder, and every recorded reason for
    something not in it."""
    missing = {name: skipped.get(name, "no reason recorded")
               for name in sorted(listed - pngs)}
    missing.update({k: v for k, v in skipped.items()
                    if k not in pngs and k not in missing and k != PAGE_NAME})
    return missing


def _not_drawn(missing: dict) -> str:
    if not missing:
        return ""
    rows = "".join(f"<tr><td>{_e(k)}</td><td>{_e(v)}</td></tr>"
                   for k, v in sorted(missing.items()))
    return ('<h2 id="not-drawn">Not drawn</h2>'
            '<table class="plain"><tr><th>Figure</th><th>Why</th></tr>'
            f'{rows}</table>')


def _other_files(files: list[str]) -> str:
    rest = [f for f in files if not f.lower().endswith(".png")]
    if not rest:
        return ""
    items = "".join(f'<li><a href="{_href(f)}">{_e(f)}</a></li>' for f in rest)
    return f'<h2 id="files">Other files</h2><ul>{items}</ul>'


def _shank_table(out_dir: Path, folders: set[str]) -> str:
    parsed, loose = {}, []
    for f in folders:
        m = SHANK_FOLDER.match(f)
        if m:
            probe = f"{m['run']}_imec{m['probe']}" if m["run"] else f"imec{m['probe']}"
            parsed.setdefault(probe, {})[int(m["shank"])] = f
        else:
            loose.append(f)
    shanks = sorted({s for row in parsed.values() for s in row})
    head = "".join(f"<th>Shank {s}</th>" for s in shanks)
    rows = []
    for probe in sorted(parsed, key=_natural):
        cells = []
        for s in shanks:
            f = parsed[probe].get(s)
            cells.append(f"<td>{_shank_link(out_dir, f)}</td>" if f else "<td></td>")
        rows.append(f"<tr><th>{_e(probe)}</th>{''.join(cells)}</tr>")
    table = (f'<table class="grid"><tr><th>Probe</th>{head}</tr>{"".join(rows)}</table>'
             if parsed else "")
    if loose:
        table += "<ul>" + "".join(f"<li>{_e(f)}: {_shank_link(out_dir, f)}</li>"
                                  for f in sorted(loose)) + "</ul>"
    return table


def _shank_link(out_dir: Path, folder: str) -> str:
    if (out_dir / folder / PAGE_NAME).exists():
        return f'<a href="{_href(folder)}/{PAGE_NAME}" title="{_e(folder)}">Figures</a>'
    if (out_dir / folder).is_dir():
        return (f'<a href="{_href(folder)}/" title="{_e(folder)}">Folder</a> '
                f'<span class="note">(no page yet)</span>')
    return '<span class="note">No figures</span>'


def _files(out_dir: Path) -> list[str]:
    try:
        return sorted(p.name for p in out_dir.iterdir()
                      if p.is_file() and p.name not in _IGNORED
                      and not p.name.startswith(".") and not p.name.endswith(".tmp"))
    except OSError:
        return []


def _write(out_dir: Path, title: str, nav: list[str], body: list[str],
           script: str = "") -> Path:
    from . import __version__

    page = PAGE_TEMPLATE.format(
        title=_e(title), nav=" ".join(nav), body="\n".join(b for b in body if b),
        version=_e(__version__), script=script)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / PAGE_NAME
    tmp = out_dir / f".{PAGE_NAME}.tmp"
    # Written with \n line ends on every platform, so a page is byte for byte
    # the same whichever machine wrote it.
    tmp.write_text(page, encoding="utf-8", newline="\n")
    os.replace(tmp, path)
    return path


def _e(text) -> str:
    return html.escape(str(text), quote=True)


def _href(name: str) -> str:
    from urllib.parse import quote
    return quote(name)


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def _natural(text: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", text)]


def _finite(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


PAGE_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
:root {{ --ink: #1d1d1f; --muted: #5f6368; --rule: #d8d8dc; --bg: #ffffff;
         --link: #0b57d0; --peak: crimson; }}
body {{ font-family: Arial, Helvetica, sans-serif; color: var(--ink);
        background: var(--bg); margin: 0; }}
header {{ position: sticky; top: 0; background: var(--bg); z-index: 2;
          border-bottom: 1px solid var(--rule); padding: 10px 16px; }}
header h1 {{ font-size: 18px; margin: 0 0 6px 0; }}
nav a {{ margin-right: 14px; font-size: 14px; }}
main {{ padding: 0 16px 48px 16px; max-width: 1800px; }}
a {{ color: var(--link); }}
h2 {{ font-size: 17px; margin: 28px 0 4px 0; padding-top: 8px;
      border-top: 2px solid var(--rule); scroll-margin-top: 80px; }}
h3 {{ font-size: 15px; margin: 18px 0 2px 0; }}
.fig {{ scroll-margin-top: 80px; }}
.fig p {{ margin: 2px 0 6px 0; color: var(--muted); font-size: 13px; }}
.fig img {{ max-width: 100%; height: auto; border: 1px solid var(--rule); }}
.name {{ font-weight: normal; color: var(--muted); font-size: 12px;
         margin-left: 6px; }}
.note {{ color: var(--muted); font-size: 12px; }}
table.plain td, table.plain th {{ text-align: left; padding: 3px 12px 3px 0;
                                  font-size: 13px; vertical-align: top; }}
table.grid {{ border-collapse: collapse; }}
table.grid td, table.grid th {{ border: 1px solid var(--rule); padding: 6px 10px;
                               font-size: 13px; text-align: left; }}
.tpl {{ margin-top: 10px; }}
.tpl-controls {{ display: flex; flex-wrap: wrap; gap: 10px; align-items: center;
                 margin: 6px 0; font-size: 13px; }}
.tpl-controls button {{ font: inherit; padding: 4px 10px; cursor: pointer; }}
.tpl-controls input {{ font: inherit; }}
#tpl-panels {{ display: grid; grid-template-columns: repeat(6, minmax(0, 1fr));
               gap: 4px; }}
#tpl-panels svg {{ width: 100%; height: auto; }}
@media (max-width: 900px) {{ #tpl-panels {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }} }}
footer {{ color: var(--muted); font-size: 12px; padding: 16px; }}
</style>
</head>
<body>
<header><h1>{title}</h1><nav>{nav}</nav></header>
<main>
{body}
</main>
<footer>Written by qualitymetrics {version}. The figures are the record; this
page only arranges them, and lists what was in the folder when it was written.</footer>
{script}
</body>
</html>
"""

#: Draws template panels the way plots.units.templates_grid does: each site's
#: waveform at its position, the unit scaled to its own largest excursion over
#: three rows, the peak site in red, depth on the y axis, a 1 ms bar.
TEMPLATE_SCRIPT = r"""<script>
(function () {
  "use strict";
  var el = document.getElementById("tpl-data");
  if (!el) return;
  var D = JSON.parse(el.textContent);
  var raw = atob(D.waves), bytes = new Int8Array(raw.length);
  for (var i = 0; i < raw.length; i++) {
    var b = raw.charCodeAt(i); bytes[i] = b > 127 ? b - 256 : b;
  }
  var nS = D.n_sites, nT = D.n_samples, index = {};
  D.units.forEach(function (u, i) { index[u] = i; });
  var anyPass = D.passing.some(Boolean);
  var exclude = {}; (D.exclude || []).forEach(function (u) { exclude[u] = true; });
  var pool = D.units.filter(function (u, i) { return anyPass ? D.passing[i] : true; });
  var queue = [], firstRound = true;
  function shuffle(a) {
    for (var i = a.length - 1; i > 0; i--) {
      var j = Math.floor(Math.random() * (i + 1)), t = a[i]; a[i] = a[j]; a[j] = t;
    }
    return a;
  }
  function refill() {
    // The first round leaves out the units the static grid already shows.
    var fresh = pool.filter(function (u) { return !(firstRound && exclude[u]); });
    if (!fresh.length) fresh = pool.slice();
    queue = shuffle(fresh); firstRound = false;
  }
  var NS = "http://www.w3.org/2000/svg";
  function node(tag, attrs, parent) {
    var n = document.createElementNS(NS, tag);
    for (var k in attrs) n.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(n);
    return n;
  }
  var spanMs = (nT - 1) / D.fs * 1000 || 1, xScale = 0.85 * D.dx / spanMs;
  var bad = {}; (D.bad || []).forEach(function (u) { bad[u] = true; });
  function panel(u, why, first) {
    var i = index[u], chans = D.chans.slice(i * nS, i * nS + nS);
    var xs = chans.map(function (c) { return D.positions[c][0]; });
    var ys = chans.map(function (c) { return D.positions[c][1]; });
    var x0 = Math.min.apply(null, xs) - 0.15 * D.dx, x1 = Math.max.apply(null, xs) + D.dx;
    var y0 = Math.min.apply(null, ys) - 3.2 * D.dy, y1 = Math.max.apply(null, ys) + 3.2 * D.dy;
    var W = 200, H = 380, L = 40, T = 34, B = 18, pw = W - L - 6, ph = H - T - B;
    function px(x) { return L + (x - x0) / (x1 - x0) * pw; }
    function py(y) { return T + (y1 - y) / (y1 - y0) * ph; }
    var svg = node("svg", {viewBox: "0 0 " + W + " " + H, role: "img",
                           "font-family": "Arial, Helvetica, sans-serif"});
    var amp = D.amp_uv[i];
    node("text", {x: L + pw / 2, y: 12, "text-anchor": "middle", "font-size": 10},
         svg).textContent = why;
    node("text", {x: L + pw / 2, y: 25, "text-anchor": "middle", "font-size": 10},
         svg).textContent = "Unit " + u + (amp === null ? "" : ", " + Math.round(amp) + " µV") +
                            (bad[u] ? " (template not finite)" : "");
    node("line", {x1: L, x2: L, y1: T, y2: T + ph, stroke: "#222", "stroke-width": 1}, svg);
    var step = 50, t0 = Math.ceil(y0 / step) * step;
    for (var y = t0; y <= y1; y += step) {
      node("line", {x1: L - 3, x2: L, y1: py(y), y2: py(y), stroke: "#222"}, svg);
      node("text", {x: L - 5, y: py(y) + 3, "text-anchor": "end", "font-size": 8},
           svg).textContent = y;
    }
    if (first) {
      node("text", {x: 9, y: T + ph / 2, "font-size": 9, "text-anchor": "middle",
                    transform: "rotate(-90 9 " + (T + ph / 2) + ")"}, svg).textContent = "Depth (µm)";
      var bx = px(x0), by = T + ph + 4;
      node("line", {x1: bx, x2: px(x0 + xScale), y1: by, y2: by, stroke: "#222",
                    "stroke-width": 2}, svg);
      node("text", {x: (bx + px(x0 + xScale)) / 2, y: by + 11, "font-size": 8,
                    "text-anchor": "middle"}, svg).textContent = "1 ms";
    }
    var yScale = 3.0 * D.dy / 127.0, base = i * nT * nS;
    for (var j = 0; j < nS; j++) {
      var pts = [];
      for (var t = 0; t < nT; t++) {
        var v = bytes[base + t * nS + j];
        pts.push(px(xs[j] + (t / D.fs * 1000) * xScale).toFixed(1) + "," +
                 py(ys[j] + v * yScale).toFixed(1));
      }
      var isPeak = chans[j] === D.peak_site[i];
      node("polyline", {points: pts.join(" "), fill: "none",
                        stroke: isPeak ? "crimson" : "#333",
                        "stroke-width": isPeak ? 1.1 : 0.8}, svg);
    }
    return svg;
  }
  var box = document.getElementById("tpl-panels");
  var status = document.getElementById("tpl-status");
  function draw(units, why) {
    box.innerHTML = "";
    units.forEach(function (u, k) {
      // The depth label on the first panel of each row of six, as in the grid.
      box.appendChild(panel(u, typeof why === "function" ? why(u) : why, k % 6 === 0));
    });
  }
  document.getElementById("tpl-next").addEventListener("click", function () {
    var picked = [];
    while (picked.length < 6 && picked.length < pool.length) {
      if (!queue.length) refill();
      var u = queue.pop();
      if (picked.indexOf(u) < 0) picked.push(u);
    }
    draw(picked.sort(function (a, b) { return a - b; }), "Random");
    status.textContent = "Random units " + picked.join(", ") + ". " +
      queue.length + " not yet shown in this round.";
  });
  function showTyped() {
    var asked = (document.getElementById("tpl-ids").value.match(/\d+/g) || [])
      .map(Number), known = [], unknown = [];
    asked.forEach(function (u) { (u in index ? known : unknown).push(u); });
    known = known.filter(function (u, k) { return known.indexOf(u) === k; });
    var dropped = known.length > 18 ? known.length - 18 : 0;
    known = known.slice(0, 18);
    if (known.length) draw(known, function (u) {
      return D.passing[index[u]] ? "Passed QC" : "Did not pass QC";
    });
    status.textContent = known.length
      ? "Units " + known.join(", ") + "." + (unknown.length ? " No template for " + unknown.join(", ") + "." : "") +
        (dropped ? " Only the first 18 are drawn; " + dropped + " more left out." : "")
      : "No template for " + (unknown.join(", ") || "those units") + ".";
  }
  document.getElementById("tpl-show").addEventListener("click", showTyped);
  document.getElementById("tpl-ids").addEventListener("keydown", function (e) {
    if (e.key === "Enter") showTyped();
  });
})();
</script>"""
