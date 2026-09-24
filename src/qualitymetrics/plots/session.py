"""Figures that describe a whole session rather than one shank.

Every other module here takes one :class:`~qualitymetrics.ksdata.KilosortResults`
and draws what that shank shows. Two questions are not answerable that way:

* **How noisy is this recording?** Per-site AP band RMS for every shank at once,
  drawn as distributions. A shank with a dozen bad channels and a shank that is
  uniformly noisy have the same median and look nothing alike, so a bar of
  medians hides the two things a noise measurement is read for.
* **What did the whole probe see?** Every unit on one time axis, hue by probe
  and lightness by shank, so drift, a dead shank or a session-wide event is one
  glance rather than sixteen.

The unit figures come in two forms, every sorted unit and only those that passed
quality control, and the second needs nothing computed: ``quality_metrics.json``
is written beside each sort and already carries the per-unit verdicts. Nothing
here recomputes a metric. A figure that disagreed with the pipeline's own
numbers would be worse than no figure.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..ksdata import KilosortResults
from ..style import despine, save, use_lab_style

#: ``imec<P>_shank<S>``, the directory one sort writes.
SHANK_DIR_RE = re.compile(r"^imec(?P<probe>\d+)_shank(?P<shank>\d+)$")

#: Hue identifies the probe, lightness identifies the shank.
#:
#: Generated in OKLCH and checked all-pairs for separation under protanopia and
#: deuteranopia (Machado, Oliveira and Fernandes 2009, severity 1.0), and within
#: each hue for monotone lightness so shank order is visible in the color.
#: Sixteen colors can never be pairwise distinct, so color is never the only
#: cue here: every figure also labels the probe and shank on an axis.
PALETTE: dict[tuple[int, int], str] = {
    (0, 0): "#006eae", (0, 1): "#0090e2", (0, 2): "#50b3ff", (0, 3): "#a2d4ff",
    (1, 0): "#9f5000", (1, 1): "#ce6a00", (1, 2): "#f28e42", (1, 3): "#ffbe91",
    (2, 0): "#007d51", (2, 1): "#00a36a", (2, 2): "#38c789", (2, 3): "#8ae4b5",
    (3, 0): "#9b418b", (3, 1): "#c25faf", (3, 2): "#e484d0", (3, 3): "#fbb2ea",
}
#: Grey, for a layout beyond the sixteen the palette names. A wrong hue would
#: read as a probe identity that does not exist.
UNKNOWN_COLOR = "#8a8a8a"
INK = "#222222"
MUTED = "#666666"
GRID = "#dddddd"

#: The AP band: above the LFP, below where a spike has no energy left.
SPIKE_BAND_HZ = 300.0
#: Four one-second windows per shank measures a noise floor perfectly well and
#: keeps a 16-shank session to about a minute of reading over a network share.
NOISE_WINDOWS = 4
NOISE_WINDOW_S = 1.0
#: Skipped at each end: the first moments of a recording carry settling
#: transients that are not the noise floor.
NOISE_EDGE_SKIP_S = 20.0
#: Point size of the per-shank labels down the left of the raster.
SHANK_LABEL_PT = 7.5
#: Height of the raster axes in inches, near enough, after tight_layout on the
#: 8 inch figure below. Only used to convert a label height into rows.
RASTER_AXES_IN = 6.2


def min_block_rows(max_rows: int, *, label_pt: float = SHANK_LABEL_PT,
                   axes_in: float = RASTER_AXES_IN, spacing: float = 1.6) -> int:
    """Rows a shank's block needs before its label stops hitting its neighbour.

    Computed rather than picked, because the two quantities are in different
    units and guessing converts between them badly. A block's height is a share
    of ``max_rows`` spread over the axes; a label's height is a share of an
    inch. Fixed guesses of 8 and then 20 rows both looked reasonable and both
    collided, the second one by about a pixel.

    This bites only when quality control leaves a shank with almost no units,
    which is exactly when the figure most needs to say which shank that is.
    """
    import math

    return int(math.ceil(max_rows * (label_pt / 72.0) / axes_in * spacing))


def color_for(probe: int, shank: int) -> str:
    return PALETTE.get((probe, shank), UNKNOWN_COLOR)


# --------------------------------------------------------------------------
# Finding the sorts
# --------------------------------------------------------------------------
@dataclass
class SessionShank:
    """One sorted shank of a session, and which shank of which probe it is."""

    probe: int
    shank: int
    directory: Path
    ks: KilosortResults
    duration_s: float

    @property
    def label(self) -> str:
        return f"imec{self.probe}_shank{self.shank}"

    @property
    def short_label(self) -> str:
        """``0.2`` for imec0 shank 2, which is what fits on a crowded axis."""
        return f"{self.probe}.{self.shank}"

    @property
    def color(self) -> str:
        return color_for(self.probe, self.shank)

    def quality_metrics(self) -> dict:
        """The per-unit verdicts the pipeline recorded, or ``{}`` if absent.

        Read rather than recomputed. The sliding refractory test, the noise
        cutoff and the amplitude are all decided at sort time and written here;
        deciding them again in a plotting library invites two answers to one
        question.
        """
        path = self.directory / "quality_metrics.json"
        if not path.exists():
            return {}
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {}

    def passing_units(self, criterion: str = "pass_rescued") -> set[int] | None:
        """Unit ids that passed, or None when this sort recorded no verdicts.

        None and an empty set mean different things: nothing was recorded, as
        against nothing passed. A figure must not present the first as the
        second.
        """
        metrics = self.quality_metrics()
        units = metrics.get("units")
        if not units:
            return None
        if not any(criterion in unit for unit in units):
            return None
        return {int(unit["unit_id"]) for unit in units if unit.get(criterion)}


@dataclass
class SessionSorts:
    """What was found under a session, including what was not usable."""

    session_dir: Path
    shanks: list[SessionShank] = field(default_factory=list)
    #: Directories that look sorted but cannot be read, with the reason.
    incomplete: list[str] = field(default_factory=list)

    @property
    def probes(self) -> list[int]:
        return sorted({s.probe for s in self.shanks})

    @property
    def describe(self) -> str:
        return (f"{len(self.shanks)} shanks across "
                f"{len(self.probes)} probe{'s' if len(self.probes) != 1 else ''}")


def find_session_shanks(session_dir: str | Path) -> SessionSorts:
    """Every sorted shank under a session, ordered by probe then shank.

    ``session_dir`` may be the session folder or its ``sorting`` directory.

    Unusable directories are collected rather than raised on, so the figures
    are still drawn for the shanks that are there. They are reported so the
    caller can say so *on the figure*: a session quietly becoming a subset of
    itself is the failure worth avoiding, and a subtitle naming what is missing
    avoids it without throwing the picture away.
    """
    session_dir = Path(session_dir)
    sorting_dir = (session_dir if session_dir.name == "sorting"
                   else session_dir / "sorting")
    found = SessionSorts(session_dir=session_dir)
    if not sorting_dir.is_dir():
        found.incomplete.append(f"no sorting directory under {session_dir}")
        return found

    for entry in sorted(sorting_dir.iterdir()):
        match = SHANK_DIR_RE.match(entry.name)
        if not (match and entry.is_dir()):
            continue
        provenance_path = entry / "provenance.json"
        if not provenance_path.exists():
            found.incomplete.append(f"{entry.name} (no provenance.json)")
            continue
        try:
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            found.incomplete.append(f"{entry.name} (provenance unreadable: {exc})")
            continue
        recording = provenance.get("recording", {})
        gain = provenance.get("gain_to_uv", recording.get("uv_per_bit"))
        if gain is None:
            # Every voltage in the noise figure would otherwise be ADC counts
            # wearing a microvolt label.
            found.incomplete.append(f"{entry.name} (no gain_to_uv or uv_per_bit)")
            continue
        label = f"imec{match.group('probe')}_shank{match.group('shank')}"
        try:
            # Imported here rather than at module scope: report.py imports the
            # plots package, so a top-level import would close a cycle.
            from ..report import find_sorter_output

            # A SortingManager sort nests its output under kilosort4/, one
            # level deeper than KilosortResults.load looks by itself.
            ks = KilosortResults.load(find_sorter_output(entry),
                                      uv_per_bit=float(gain), label=label)
        except Exception as exc:  # noqa: BLE001 - reported, never raised
            found.incomplete.append(f"{entry.name} ({type(exc).__name__}: {exc})")
            continue
        found.shanks.append(SessionShank(
            probe=int(match.group("probe")),
            shank=int(match.group("shank")),
            directory=entry,
            ks=ks,
            duration_s=float(provenance.get("duration_s", 0.0)),
        ))
    found.shanks.sort(key=lambda s: (s.probe, s.shank))
    return found


# --------------------------------------------------------------------------
# Noise
# --------------------------------------------------------------------------
def _ap_paths(shank: SessionShank) -> tuple[Path, Path]:
    """The archived ``.cbin`` and ``.ch`` this shank was sorted from."""
    raw_root = shank.directory.parent.parent / "raw_ephys_data"
    pattern = f"*imec{shank.probe}/*imec{shank.probe}.sh{shank.shank}.ap.cbin"
    matches = sorted(raw_root.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"no AP file matching {pattern} under {raw_root}")
    if len(matches) > 1:
        raise ValueError(f"{len(matches)} AP files match {pattern}; refusing to guess")
    return matches[0], matches[0].with_suffix(".ch")


def channel_rms(shank: SessionShank, *, n_windows: int = NOISE_WINDOWS,
                window_s: float = NOISE_WINDOW_S) -> dict:
    """Per-site AP band RMS for one shank, in microvolts.

    Two numbers per site, because they answer different questions. The
    unreferenced RMS is the conventional figure quoted for Neuropixels and
    includes whatever the whole shank picks up in common. The residual, after
    subtracting the across-site median sample by sample, is what sets the spike
    detection floor. A grounding problem moves the first and leaves the second
    alone, which is how it is told apart from a genuinely noisy probe.

    Median across windows, per site, so one transient in one window cannot set
    a site's level.
    """
    import mtscomp
    from scipy import signal

    cbin, ch = _ap_paths(shank)
    reader = mtscomp.Reader()
    reader.open(str(cbin), str(ch))
    try:
        total, n_stored = reader.shape
        positions = np.load(shank.ks.path / "channel_positions.npy")
        n_ap = positions.shape[0]
        if n_stored < n_ap:
            raise ValueError(
                f"{shank.label}: the AP file stores {n_stored} channels but the "
                f"geometry describes {n_ap}. Refusing to guess which is which.")

        fs = float(shank.ks.fs)
        per_window = int(round(window_s * fs))
        skip = int(round(NOISE_EDGE_SKIP_S * fs))
        start, stop = skip, total - skip - per_window
        if stop <= start:                       # a recording shorter than 40 s
            start, stop = 0, max(total - per_window, 1)
        starts = np.linspace(start, stop, n_windows).astype(np.int64)

        sos = signal.butter(3, SPIKE_BAND_HZ, "hp", fs=fs, output="sos")
        common_rms, residual_rms, raw_rms = [], [], []
        for begin in starts:
            # The trailing stored channel is the SY sync word, a bit field
            # rather than a signal, so it is never part of a noise estimate.
            block = np.asarray(reader[begin:begin + per_window],
                               dtype=np.float32)[:, :n_ap] * shank.ks.uv_per_bit
            block = signal.sosfiltfilt(sos, block, axis=0)
            common = np.median(block, axis=1)
            common_rms.append(float(np.sqrt(np.mean(common ** 2))))
            residual_rms.append(np.std(block - common[:, None], axis=0))
            raw_rms.append(np.std(block, axis=0))
    finally:
        reader.close()

    residual = np.median(np.stack(residual_rms), axis=0)
    raw = np.median(np.stack(raw_rms), axis=0)
    return {
        "label": shank.label,
        "probe": shank.probe,
        "shank": shank.shank,
        "n_sites": int(n_ap),
        "common_mode_uv": float(np.median(common_rms)),
        "residual_uv": float(np.median(residual)),
        "unreferenced_uv": float(np.median(raw)),
        "per_channel": {
            "depth_um": positions[:, 1].astype(float),
            "unreferenced_uv": raw,
            "residual_uv": residual,
        },
    }


def rms_across_shanks(measurements: list[dict], *, title: str = "",
                      subtitle: str = ""):
    """Per-site AP band RMS for every shank, as distributions.

    One violin per shank rather than a bar of medians: how tightly the sites
    agree, and whether a handful are far out, are both invisible in a median
    and are the two things this figure is read for.
    """
    import matplotlib.pyplot as plt

    use_lab_style()
    rows = [
        ("unreferenced_uv", "AP band RMS, unreferenced\n(the conventional number)"),
        ("residual_uv", "AP band RMS, common-median referenced\n(the detection floor)"),
    ]
    labels = [f"{m['probe']}.{m['shank']}" for m in measurements]
    figure, axes = plt.subplots(
        len(rows), 1, figsize=(max(7.0, 0.55 * len(labels) + 3.0), 8.0),
        sharex=True, squeeze=False)

    for row, (measure, ylabel) in enumerate(rows):
        axis = axes[row][0]
        groups = [np.asarray(m["per_channel"][measure], dtype=float)
                  for m in measurements]
        parts = axis.violinplot(groups, positions=np.arange(len(labels)),
                                widths=0.85, showextrema=False, showmedians=True)
        for body, measurement in zip(parts["bodies"], measurements, strict=False):
            body.set_facecolor(color_for(measurement["probe"], measurement["shank"]))
            body.set_alpha(0.85)
            body.set_edgecolor("none")
        parts["cmedians"].set_color(INK)
        parts["cmedians"].set_linewidth(1.2)

        overall = float(np.median(np.concatenate(groups)))
        axis.axhline(overall, color=INK, linewidth=0.9, linestyle="--", alpha=0.6)
        axis.text(len(labels) - 0.4, overall, f" {overall:.1f}", va="center",
                  ha="left", fontsize=8, color=INK, clip_on=False)
        axis.set_ylabel(f"{ylabel}\nRMS (µV)", fontsize=9)
        axis.grid(axis="y", color=GRID, linewidth=0.6)
        axis.set_axisbelow(True)
        # A violin tail is a kernel estimate and runs past the data, so an
        # automatic limit leaves the panel mostly empty.
        axis.set_ylim(0, float(np.percentile(np.concatenate(groups), 99.9)) * 1.12)
        despine(axis)

    axes[-1][0].set_xticks(np.arange(len(labels)))
    axes[-1][0].set_xticklabels(labels, fontsize=8)
    axes[-1][0].set_xlabel("probe.shank", fontsize=9, color=MUTED)
    _caption(figure, title, subtitle)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    return figure


# --------------------------------------------------------------------------
# The raster
# --------------------------------------------------------------------------
def unit_rates(shank: SessionShank, edges: np.ndarray,
               keep: set[int] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Binned spike counts per unit, ordered by depth from the shank tip.

    Depth order puts neighbouring rows near each other in the brain rather than
    near each other in cluster numbering. ``keep`` selects unit ids.
    """
    times = shank.ks.spike_times_s
    clusters = shank.ks.spike_clusters
    if times.size == 0:
        return np.zeros((0, len(edges) - 1), dtype=np.float32), np.zeros(0)

    unit_ids = np.unique(clusters)
    if keep is not None:
        unit_ids = np.array([u for u in unit_ids if int(u) in keep], dtype=unit_ids.dtype)
    if unit_ids.size == 0:
        return np.zeros((0, len(edges) - 1), dtype=np.float32), np.zeros(0)

    try:
        spike_depths = shank.ks.spike_depths_um
    except Exception:  # noqa: BLE001 - depth is an ordering nicety, not the data
        spike_depths = None

    counts = np.zeros((unit_ids.size, len(edges) - 1), dtype=np.float32)
    depths = np.zeros(unit_ids.size, dtype=float)
    for row, unit in enumerate(unit_ids):
        mine = clusters == unit
        counts[row], _ = np.histogram(times[mine], bins=edges)
        if spike_depths is not None and mine.any():
            depths[row] = float(np.median(spike_depths[mine]))

    if spike_depths is not None:
        order = np.argsort(depths, kind="stable")
        counts, depths = counts[order], depths[order]
    return counts, depths


def _fit_rows(counts: np.ndarray, target: int) -> np.ndarray:
    """Make a block exactly ``target`` rows tall, both directions.

    Shrinking averages neighbouring units: handing ``imshow`` more rows than it
    has pixels lets it drop rows silently, which on a five thousand unit
    session means units missing from the picture with nothing saying so.

    Growing repeats rows, which is what gives a shank with one surviving unit
    enough height to carry its own label. A row is therefore not a unit in
    either direction, which is why the axis says "units, by probe.shank then
    depth" rather than labelling row counts.
    """
    n = counts.shape[0]
    if n == target or n == 0:
        return counts
    if n > target:
        edges = np.linspace(0, n, target + 1).astype(int)
        return np.stack([counts[a:b].mean(axis=0) if b > a else counts[a]
                         for a, b in zip(edges[:-1], edges[1:], strict=False)])
    return counts[np.linspace(0, n - 1, target).round().astype(int)]


def session_raster(shanks: list[SessionShank], *, duration_s: float | None = None,
                   criterion: str | None = None, bin_s: float = 1.0,
                   max_rows: int = 1400, title: str = "", subtitle: str = ""):
    """Every unit of every shank on one time axis.

    ``criterion`` names a column of ``quality_metrics.json`` (``pass_strict``
    or ``pass_rescued``) to filter on. None draws every sorted unit, which is
    the right input to "is this recording worth curating".

    Each unit is a row shaded by its rate normalised to its own maximum, so a
    quiet unit is as visible as a loud one and the figure shows *when* things
    fired rather than how hard.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb
    from matplotlib.lines import Line2D

    use_lab_style()
    duration_s = duration_s or max((s.duration_s for s in shanks), default=0.0)
    if duration_s <= 0:
        raise ValueError("the session has no recorded duration to plot against")
    edges = np.linspace(0.0, duration_s, max(int(round(duration_s / bin_s)), 1) + 1)

    blocks, unverdicted = [], []
    for shank in shanks:
        keep = None
        if criterion is not None:
            keep = shank.passing_units(criterion)
            if keep is None:
                unverdicted.append(shank.label)
                continue
        counts, _ = unit_rates(shank, edges, keep=keep)
        if counts.shape[0]:
            blocks.append((shank, counts))
    if not blocks:
        raise ValueError("no units to draw for this session")

    total_units = sum(counts.shape[0] for _, counts in blocks)
    floor_rows = min_block_rows(max_rows)
    image_rows, boundaries, row_colors = [], [], []
    for shank, counts in blocks:
        # Proportional to unit count, so a shank that produced few units looks
        # like it, but floored so every shank still has room for its label.
        share = max(int(round(max_rows * counts.shape[0] / total_units)),
                    floor_rows)
        reduced = _fit_rows(counts, share)
        peak = reduced.max(axis=1, keepdims=True)
        peak[peak == 0] = 1.0
        image_rows.append(reduced / peak)
        row_colors.extend([shank.color] * reduced.shape[0])
        boundaries.append((shank, reduced.shape[0], counts.shape[0]))

    intensity = np.vstack(image_rows)
    rgb = np.ones((*intensity.shape, 3), dtype=float)
    for row in range(intensity.shape[0]):
        color = np.asarray(to_rgb(row_colors[row]))
        # Composite onto white: a silent bin is the page, a busy one is the
        # shank's color at full strength.
        rgb[row] = 1.0 - intensity[row][:, None] * (1.0 - color)[None, :]

    figure, axis = plt.subplots(figsize=(13.0, 8.0))
    axis.imshow(rgb, aspect="auto", interpolation="nearest", origin="upper",
                extent=(0.0, duration_s / 60.0, intensity.shape[0], 0))

    position = 0
    for index, (shank, height, _n) in enumerate(boundaries):
        axis.text(-0.008, position + height / 2, shank.short_label,
                  transform=axis.get_yaxis_transform(), ha="right",
                  va="center", fontsize=SHANK_LABEL_PT, color=shank.color)
        position += height
        if position < intensity.shape[0]:
            following = (boundaries[index + 1][0]
                         if index + 1 < len(boundaries) else None)
            # A heavier rule between probes, where the shade ramp restarts.
            if following is not None and following.probe != shank.probe:
                axis.axhline(position, color=INK, linewidth=1.0, alpha=0.5)
            else:
                axis.axhline(position, color="white", linewidth=0.8)

    axis.set_xlabel("Time in the recording (minutes)", fontsize=9)
    # labelpad clears the per-shank labels, drawn in the same margin with the
    # y-axis transform. Without it the two are painted over each other.
    axis.set_ylabel(f"{total_units} units, by probe.shank then depth",
                    fontsize=9, labelpad=26)
    axis.set_yticks([])
    axis.spines["left"].set_visible(False)
    despine(axis)
    axis.legend(handles=[
        Line2D([0], [0], color=color_for(probe, 1), linewidth=6,
               label=f"imec{probe}")
        for probe in sorted({s.probe for s, _h, _n in boundaries})],
        loc="upper right", frameon=False, fontsize=8, ncol=4,
        bbox_to_anchor=(1.0, 1.06))

    if unverdicted:
        subtitle += (f"   No recorded verdicts, so left out: "
                     f"{', '.join(unverdicted)}")
    _caption(figure, title, subtitle)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    return figure


def _caption(figure, title: str, subtitle: str) -> None:
    if title:
        figure.text(0.006, 0.992, title, fontsize=13, fontweight="bold", va="top")
    if subtitle:
        figure.text(0.006, 0.960, subtitle, fontsize=9, color=MUTED, va="top")


def noise_table(measurements: list[dict], path: str | Path) -> str:
    """The medians behind the violins, so two sessions compare as numbers."""
    lines = ["label,probe,shank,n_sites,unreferenced_uv,residual_uv,common_mode_uv"]
    for m in measurements:
        lines.append(f"{m['label']},{m['probe']},{m['shank']},{m['n_sites']},"
                     f"{m['unreferenced_uv']:.3f},{m['residual_uv']:.3f},"
                     f"{m['common_mode_uv']:.3f}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


__all__ = ["SessionShank", "SessionSorts", "find_session_shanks", "channel_rms",
           "rms_across_shanks", "session_raster", "unit_rates", "noise_table",
           "color_for", "PALETTE", "save"]
