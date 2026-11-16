"""Build the whole figure set for one sorted shank.

Every figure is attempted independently and a failure is recorded rather than
raised, because a QC report that produces nothing when one panel cannot be drawn
is a QC report nobody runs. The return value says what was made and what was
not, so a missing figure is visible rather than silently absent.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from . import plots
from .ksdata import KilosortResults, label_from_path
from .metrics import phy_metrics_from_results
from .pipelinefiles import (
    LIGHT_EVENTS_FILE,
    PIPELINE_FILE,
    load_destriped_windows,
    load_light_events,
)
from .raw import RawError, RawRecording, find_band, parse_meta, resolve_uv_per_bit
from .style import save


@dataclass
class ReportResult:
    """What the report managed to build."""

    out_dir: Path
    made: dict[str, str] = field(default_factory=dict)
    skipped: dict[str, str] = field(default_factory=dict)

    def summary(self) -> str:
        lines = [f"{len(self.made)} figures in {self.out_dir}"]
        for name in sorted(self.made):
            lines.append(f"  made    {name}")
        for name, why in sorted(self.skipped.items()):
            lines.append(f"  skipped {name}: {why}")
        return "\n".join(lines)


def _attempt(result: ReportResult, name: str, fn):
    """Run one figure builder, recording success or the reason it failed."""
    try:
        fig = fn()
    except Exception as exc:  # noqa: BLE001 - one bad panel must not stop the rest
        result.skipped[name] = f"{type(exc).__name__}: {exc}"
        return
    if fig is None:
        result.skipped[name] = "produced no figure"
        return
    result.made[name] = save(fig, result.out_dir / f"{name}.png")


def find_sorter_output(sorting_dir: str | Path) -> Path:
    """Locate the sorter_output inside a shank directory."""
    sorting_dir = Path(sorting_dir)
    if (sorting_dir / "spike_times.npy").exists():
        return sorting_dir
    for cand in (sorting_dir / "sorter_output",
                 sorting_dir / "kilosort4" / "sorter_output"):
        if (cand / "spike_times.npy").exists():
            return cand
    hits = sorted(sorting_dir.glob("**/spike_times.npy"))
    if hits:
        return hits[0].parent
    raise FileNotFoundError(f"no Kilosort output under {sorting_dir}")


def find_recording(sorting_dir: str | Path, band: str = "ap") -> Path | None:
    """Find the archived recording that goes with a shank directory.

    Walks up to the session root and looks in raw_ephys_data. Returns None
    rather than raising, because the sorter-only figures are still worth having
    when the raw file is on another machine.
    """
    sorting_dir = Path(sorting_dir).resolve()
    run, trigger, probe, shank = _shank_dir_key(sorting_dir)
    if probe is None:
        return None
    for parent in list(sorting_dir.parents)[:4]:
        session_root = (parent / "raw_ephys_data").is_dir()
        try:
            return find_band(parent, probe, shank, band, run=run, trigger=trigger)
        except RawError:
            pass
        if session_root:
            # The session's own raw_ephys_data has no match. Folders higher
            # up hold other sessions of the subject. Another session's file
            # for this probe and shank would be wrong.
            return None
    return None


def _shank_dir_key(sorting_dir: Path) -> tuple[str | None, str | None,
                                               str | None, int]:
    """Returns ``(run, trigger, "imec<P>", shank)`` for a shank directory.

    The name is ``p1_g0_imec0_shank2``. Sorts made before SortingManager
    added the run are named ``imec0_shank2``. The run is read from
    ``provenance.json`` if it names a source. For details, see
    ``plots.session.sort_run``. Returns ``(None, None, None, 0)`` for any
    other name.
    """
    from .plots.session import SHANK_DIR_RE, sort_run, sort_trigger

    match = SHANK_DIR_RE.match(sorting_dir.name)
    if match is None:
        return None, None, None, 0
    try:
        recording = json.loads((sorting_dir / "provenance.json").read_text(
            encoding="utf-8")).get("recording", {})
        if not isinstance(recording, dict):
            recording = {}
    except (OSError, ValueError, AttributeError):
        recording = {}
    run = sort_run(match.group("run"), recording)
    trigger = sort_trigger(match.group("run"), recording)
    return run, trigger, f"imec{match.group('probe')}", int(match.group("shank"))


def build_report(sorting_dir: str | Path, out_dir: str | Path,
                 raw_path: str | Path | None = None,
                 lf_path: str | Path | None = None,
                 artifact_z_um: float | None = None,
                 t_start_s: float | None = None,
                 write_phy: bool = True,
                 example_units: int = 12,
                 uv_per_bit: float | None = None) -> ReportResult:
    """Make every figure that this shank's available data supports.

    raw_path and lf_path are found automatically from the session layout when
    not given. Figures that need the raw file are skipped, with a reason, when
    it cannot be found.

    ``uv_per_bit`` is the AP band's microvolts per bit, for a recording whose
    meta states no gain (SpikeGLX 20201024 writes no ``imChan0apGain``), or
    with no raw file at all. Where the meta states one, the two must agree
    (``raw.resolve_uv_per_bit``); when they do not, every figure in
    microvolts and every figure reading the raw file is skipped with the
    reason, because one of the two numbers is wrong. It applies to the AP
    file only: an LF file's gain can differ (Neuropixels 1.0), so its meta
    must state its own.

    Every figure that puts a unit or a Kilosort channel beside raw data reads
    the recording column the sort's ``channel_map.npy`` names
    (``KilosortResults.channel_map``), and is skipped with the reason when
    that map does not fit the recording.
    """
    sorting_dir = Path(sorting_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    result = ReportResult(out_dir=out_dir)

    # Independent of raw traces, amplitude calibration, and template loading.
    # A failed diagnostic must never stop the rest of this shank's report.
    from .tracking import write as write_tracking
    try:
        tracking = write_tracking(sorting_dir, out_dir, label=label_from_path(sorting_dir))
        result.made["tracking_continuity"] = str(out_dir / tracking["made"][0])
        for name in tracking["made"][1:]:
            result.made[name] = str(out_dir / name)
    except Exception as exc:  # noqa: BLE001 - keep other figures available
        result.skipped["tracking_continuity"] = f"{type(exc).__name__}: {exc}"

    sorter_output = find_sorter_output(sorting_dir)
    if raw_path is None:
        raw_path = find_recording(sorting_dir, "ap")
    if lf_path is None:
        lf_path = find_recording(sorting_dir, "lf")

    # The gain lives in the recording's meta, and every microvolt axis needs
    # it; or the caller's, where the meta states none or there is no raw file.
    scale, scale_problem = None, None
    meta_path = Path(raw_path).with_suffix(".meta") if raw_path is not None else None
    try:
        if meta_path is not None and meta_path.exists():
            scale = resolve_uv_per_bit(parse_meta(meta_path), "ap", uv_per_bit)
        elif uv_per_bit is not None:
            scale = resolve_uv_per_bit({}, "ap", uv_per_bit)
    except RawError as exc:
        scale_problem = str(exc)

    ks = KilosortResults.load(sorter_output, uv_per_bit=scale,
                              label=label_from_path(sorting_dir))

    rec = lf_rec = None
    raw_problem = scale_problem
    try:
        if raw_path is not None and scale_problem is None:
            try:
                rec = RawRecording.open(raw_path, uv_per_bit=uv_per_bit)
            except RawError as exc:
                raw_problem = str(exc)
        if lf_path is not None:
            lf_rec = RawRecording.open(lf_path)

        duration = rec.duration_s if rec is not None else None
        if t_start_s is None:
            # A window a third of the way in: past any settling transient, and
            # not so late that a short recording has nothing there.
            t_start_s = (duration or ks.duration_s) / 3.0

        # ---- saturation: counts made by SortingManager's worker, if present.
        # Positions travel in the counts file, not from the sort, which drops
        # the channels above a brain edge.
        from .plots.saturation import SATURATION_FILE, load_counts
        if not (Path(sorting_dir) / SATURATION_FILE).exists():
            result.skipped["saturation"] = ("no saturation.npz in the sort folder (SortingManager's "
                                            "worker counts saturation during its LFP pass; sorts "
                                            "made before 2026-10 have none)")
        else:
            # Loading inside the attempt, so a damaged file costs this figure only.
            _attempt(result, "saturation", lambda: plots.saturation_report(
                load_counts(sorting_dir),
                title=f"{label_from_path(sorting_dir)}: samples at the ADC's rails"))

        # ---- sorter-only figures
        if scale is None:
            reason = (f"microvolts are unknown: {scale_problem}" if scale_problem
                      else "no recording meta found, so microvolts are unknown; "
                           "amplitude figures need it")
            for name in ("amp_depth_scatter", "unit_drift", "drift_map",
                         "unit_drift_qc-passing", "drift_map_qc-passing",
                         "amplitude_cdf_over_depth", "firing_rate_image",
                         "depth_profiles", "depth_correlation",
                         "noise_cutoff_diagnostic"):
                result.skipped[name] = reason
        else:
            _attempt(result, "amp_depth_scatter",
                     lambda: plots.amp_depth_scatter(
                         ks, duration_s=duration, artifact_z_um=artifact_z_um))
            _attempt(result, "unit_drift",
                     lambda: plots.unit_drift(ks, artifact_z_um=artifact_z_um))
            _attempt(result, "drift_map", lambda: plots.drift_map(ks))
            # The same two with only the units that passed QC, on the same
            # axes, so they lie over the every-unit versions (Nick, 2026-11-15).
            # Read inside a guard: a malformed quality_metrics.json must cost
            # these two figures, not the rest of the report (review, 2026-11-15).
            try:
                passing, passing_problem = plots.passing_units(ks), None
            except Exception as exc:  # noqa: BLE001 - recorded as the reason
                passing, passing_problem = [], f"{type(exc).__name__}: {exc}"
            if passing_problem:
                for name in ("unit_drift_qc-passing", "drift_map_qc-passing"):
                    result.skipped[name] = f"QC verdicts unreadable: {passing_problem}"
            elif passing:
                _attempt(result, "unit_drift_qc-passing",
                         lambda: plots.unit_drift(ks, artifact_z_um=artifact_z_um,
                                                  unit_ids=passing))
                _attempt(result, "drift_map_qc-passing",
                         lambda: plots.drift_map(ks, unit_ids=passing))
            else:
                for name in ("unit_drift_qc-passing", "drift_map_qc-passing"):
                    result.skipped[name] = f"no unit passed QC by {ks.qc_source}"
            _attempt(result, "amplitude_cdf_over_depth",
                     lambda: plots.amplitude_cdf_over_depth(ks))
            _attempt(result, "firing_rate_image",
                     lambda: plots.firing_rate_image(ks, duration_s=duration))
            _attempt(result, "depth_profiles",
                     lambda: plots.depth_profiles(ks, duration_s=duration))
            _attempt(result, "depth_correlation",
                     lambda: plots.depth_correlation(ks, duration_s=duration))

        _attempt(result, "templates_grid", lambda: _templates_grid(
            ks, duration, example_units))
        if scale is not None:
            _attempt(result, "noise_cutoff_diagnostic",
                     lambda: plots.noise_cutoff_diagnostic(ks))

        # ---- figures that need the raw file
        if rec is None:
            reason = (f"the archived recording was not opened: {raw_problem}"
                      if raw_problem else "no archived recording found for this shank")
            for name in ("wall_heatmap", "raw_traces", "raw_with_spikes",
                         "channel_health", "filter_state", "band_rms",
                         "example_neurons", "depth_power_wide_destriped"):
                result.skipped[name] = reason
        else:
            # The pipeline's own destriped seconds, for the figures that show
            # what Kilosort was given. A damaged file costs those figures
            # their destriped half, never the rest of the report.
            cleaned, cleaned_problem = None, None
            try:
                cleaned = load_destriped_windows(sorting_dir)
                if cleaned is None:
                    cleaned_problem = (f"no {PIPELINE_FILE} in the sort folder "
                                       "(re-running the sort's quality metrics "
                                       "in SortingManager makes it)")
            except Exception as exc:  # noqa: BLE001 - recorded as the reason
                cleaned_problem = f"{PIPELINE_FILE} unreadable ({type(exc).__name__}: {exc})"
            # Units and Kilosort channels beside the raw traces: checked once,
            # so a sort and a file that do not fit cost those two figures with
            # the reason, rather than drawing "No snippets" panels.
            try:
                plots.sorted_columns(rec, ks)
                pairing = None
            except Exception as exc:  # noqa: BLE001 - recorded as the reason
                pairing = f"{type(exc).__name__}: {exc}"
            _attempt(result, "wall_heatmap", lambda: _wall_heatmap(
                rec, sorting_dir, artifact_z_um))
            # Every channel, as the pipeline's destripe left it; the archive,
            # high-passed and referenced here, only when the sort kept none.
            # A probe-wide sort names its kept channels AP<n>, whose depths are
            # not mapped here: its traces come from the archive, as before.
            traced = cleaned if (cleaned is not None
                                 and cleaned.channels is not None) else None
            _attempt(result, "raw_traces",
                     lambda: plots.raw_traces(rec, t_start_s=t_start_s,
                                              cleaned=traced))
            if cleaned is None:
                result.skipped["depth_power_wide_destriped"] = cleaned_problem
            else:
                _attempt(result, "depth_power_wide_destriped",
                         lambda: plots.depth_power_wide_destriped(rec, cleaned))
            if pairing:
                result.skipped["raw_with_spikes"] = pairing
            else:
                _attempt(result, "raw_with_spikes", lambda: plots.raw_with_spikes(
                    rec, ks, t_start_s=t_start_s))
            _attempt(result, "channel_health",
                     lambda: plots.plot_channel_health(rec))
            # The session's noise stages for this shank alone, at the same
            # seconds (plots.noise.shank_measurement); it replaced a single
            # raw spectrum (Nick, 2026-11-15).
            _attempt(result, "filter_state", lambda: _filter_state(
                rec, ks, sorting_dir))
            _attempt(result, "band_rms",
                     lambda: plots.band_rms(rec, lf_rec=lf_rec))
            if pairing:
                result.skipped["example_neurons"] = pairing
            else:
                _attempt(result, "example_neurons",
                         lambda: plots.example_neurons(rec, ks))

        power_rec = lf_rec or rec
        if power_rec is None:
            why = (f"the archived recording was not opened: {raw_problem}"
                   if raw_problem else None)
            result.skipped["depth_power"] = why or "no recording found for this shank"
            result.skipped["lfp_band_profiles"] = why or "no recording for this shank"
            result.skipped["depth_power_wide"] = why or "no recording for this shank"
        else:
            safe_start = min(t_start_s, power_rec.duration_s - 11)
            _attempt(result, "depth_power", lambda: plots.depth_power(
                power_rec, t_start_s=safe_start))
            _attempt(result, "lfp_band_profiles",
                     lambda: plots.lfp_band_profiles(power_rec,
                                                     t_start_s=safe_start))
            # The wideband version goes on the AP recording, not the LFP one:
            # the point is the noise above 300 Hz, which the LF band does not
            # contain.
            if rec is not None:
                _attempt(result, "depth_power_wide",
                         lambda: plots.depth_power_wide(
                             rec, t_start_s=min(t_start_s,
                                                rec.duration_s - 5)))
            else:
                result.skipped["depth_power_wide"] = (
                    f"the archived AP recording was not opened: {raw_problem}"
                    if raw_problem else "no AP recording for this shank")
        # No destriped depth_power: the pipeline high-passes at 300 Hz first,
        # and depth_power stops at 300 Hz, so it would show the filter.
        # depth_power_wide_destriped shows the whole band.

        # ---- phy columns
        if write_phy and scale is not None:
            try:
                phy_metrics_from_results(ks, write=True)
                result.made["phy_tsvs"] = str(sorter_output)
            except Exception as exc:  # noqa: BLE001
                result.skipped["phy_tsvs"] = f"{type(exc).__name__}: {exc}"
    finally:
        for r in (rec, lf_rec):
            if r is not None:
                r.close()

    write_page(result, ks, duration_s=duration)
    return result


def write_page(result: ReportResult, ks: KilosortResults, *,
               duration_s: float | None = None, title: str | None = None,
               kept_reasons: bool = False) -> None:
    """Write the shank's ``index.html`` (``pages.write_shank_page``) over
    whatever figures are in ``result.out_dir``, recording it in ``result``.
    ``kept_reasons``: give the reasons an earlier report kept beside the
    figures, for a page written again without drawing. Never raises: a page
    is not worth the figures."""
    from .pages import PAGE_NAME, write_shank_page

    try:
        try:
            # The units the static grid shows, left out of the page's first
            # random redraws. The picks are seeded, so these are the same.
            from .plots.units import diverse_template_units
            shown = ([u for u, _ in diverse_template_units(ks, duration_s)]
                     if ks.uv_per_bit else _biggest_units(ks, 12))
        except Exception:  # noqa: BLE001 - only the first redraw differs
            shown = []
        path = write_shank_page(result.out_dir, title=title or ks.title("QC figures"),
                                ks=ks, exclude_units=shown,
                                skipped=None if kept_reasons else result.skipped)
        result.made[PAGE_NAME] = str(path)
    except Exception as exc:  # noqa: BLE001 - recorded, never raised
        result.skipped[PAGE_NAME] = f"{type(exc).__name__}: {exc}"


def _wall_heatmap(rec, sorting_dir: Path, artifact_z_um: float | None):
    """The wall heatmap, with the as-sorted row when the sort kept its windows.

    The pipeline's files are read here rather than inside the attempt's figure
    call, so a damaged one costs the second row, not the raw figure: the
    figure then says why the row is missing instead of being absent itself.
    """
    cleaned, light, notes = None, None, []
    try:
        cleaned = load_destriped_windows(sorting_dir)
        if cleaned is None:
            notes.append(f"No as-sorted row: no {PIPELINE_FILE} in the sort "
                         "folder (re-running the sort's quality metrics in "
                         "SortingManager makes it).")
    except Exception as exc:  # noqa: BLE001 - the raw row is still worth having
        notes.append(f"No as-sorted row: {PIPELINE_FILE} unreadable "
                     f"({type(exc).__name__}: {exc}).")
    try:
        light = load_light_events(sorting_dir)
    except Exception as exc:  # noqa: BLE001 - the marks are not worth the figure
        notes.append(f"Laser edges not marked: {LIGHT_EVENTS_FILE} unreadable "
                     f"({type(exc).__name__}: {exc}).")
    return plots.wall_heatmap(rec, artifact_z_um=artifact_z_um, cleaned=cleaned,
                              light_events=light, note="  ".join(notes) or None)


def _templates_grid(ks: KilosortResults, duration_s: float | None, fallback_n: int):
    """Templates chosen for variety (plots.units.diverse_template_units): two
    of each extreme, then six at random (Nick, 2026-11-15). Footprints need
    microvolts; without them, the units with the most spikes, as before."""
    from .plots.units import PICKS_EACH, RANDOM_PICKS, diverse_template_units

    if not ks.uv_per_bit:
        return plots.templates_grid(ks, unit_ids=_biggest_units(ks, fallback_n),
                                    title=ks.title(
                                        f"Templates of the {fallback_n} units with "
                                        f"the most spikes (microvolts unknown)"))
    picks = diverse_template_units(ks, duration_s)
    with_template = {int(u) for u in ks.unit_ids
                     if int(u) < len(ks.templates_unwhitened)}
    passing = {int(u) for u, ok in ks.qc_pass.items() if ok} & with_template
    pool = (f"the {len(passing)} units passing QC" if passing
            else f"all {len(with_template)} units (none passed QC)")
    return plots.templates_grid(
        ks, unit_ids=[u for u, _ in picks], reasons=dict(picks),
        title=ks.title(f"Templates from {pool}: {PICKS_EACH} at each extreme, "
                       f"then {RANDOM_PICKS} at random; 32 nearest sites at "
                       f"their positions, each scaled to its own peak, peak "
                       f"site in red"))


def _filter_state(rec, ks: KilosortResults, sorting_dir: Path):
    """The session's median-spectrum figure, for this shank alone."""
    from .plots import noise

    _run, _trigger, probe, shank = _shank_dir_key(Path(sorting_dir))
    number = int(probe[4:]) if probe else 0
    measured = noise.shank_measurement(rec, ks, sorting_dir, probe=number,
                                       shank=shank)
    note = measured.get("pipeline_note")
    missing = f" No pipeline destripe line: {note}." if note else ""
    return noise.median_spectra(
        [measured], panel_w=7.0, panel_h=4.6,
        title=f"{label_from_path(sorting_dir)}: median spectrum across "
              f"channels, per processing stage",
        subtitle=f"{_windows_text(noise)}. {NOISE_STAGES_TEXT}{missing}")


def _windows_text(noise) -> str:
    return (f"1 s each at {', '.join(f'{int(f * 100)}%' for f in noise.WINDOW_FRACTIONS)} "
            f"of the recording")


#: What the noise figures' stages are, said once for the session figures
#: and the per-shank filter_state.
NOISE_STAGES_TEXT = (
    "Stages: raw (channel medians removed); band-pass 0.5-10 kHz; "
    "the shank's across-channel median removed (global CAR), then "
    "band-passed; then each simultaneously sampled group's median "
    "also removed (demux CAR), then band-passed; last, the sorting "
    "pipeline's own destripe of the same seconds (ADC phase shift, "
    "from sort v3 the ADC-group median, spatial high-pass), band-passed "
    "likewise, where the sort kept it.")


def _biggest_units(ks: KilosortResults, n: int) -> list[int]:
    """The n units with the most spikes, which are the ones worth showing."""
    counts = ks.n_spikes
    return [u for u, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:n]]


def build_session_report(session_dir: str | Path, out_dir: str | Path,
                         noise: bool = True) -> ReportResult:
    """Make the figures that describe a whole session rather than one shank,
    and its ``index.html`` (:func:`write_session_page`).

    ``session_dir`` is the numbered session folder, or its ``sorting``
    directory. Every shank under it is read; unusable ones are named on the
    figures rather than quietly dropped, so a partial session cannot be
    mistaken for a whole one.

    Two rasters are drawn, every sorted unit and only those that passed quality
    control, the second from the verdicts already in ``quality_metrics.json``.
    Nothing is recomputed. ``noise`` also measures the AP band from the
    archived files, which needs ``mtscomp`` and reads about a minute of data
    for a sixteen shank session; it is skipped with a reason when the files or
    the package are not there.
    """
    from .plots import session as sessionplots

    session_dir = Path(session_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    result = ReportResult(out_dir=out_dir)
    found = sessionplots.find_session_shanks(session_dir)
    try:
        _session_figures(result, session_dir, found, noise)
    finally:
        # No page for a session with nothing sorted: it would link nothing.
        if found.shanks:
            write_session_page(result, session_dir,
                               [s.directory.name for s in found.shanks])
    return result


def write_session_page(result: ReportResult, session_dir: str | Path,
                       shank_folders=(), kept_reasons: bool = False) -> None:
    """Write the session's ``index.html`` over the figures in
    ``result.out_dir``, linking each shank's page. ``kept_reasons`` as for
    :func:`write_page`. Never raises."""
    from .pages import PAGE_NAME
    from .pages import write_session_page as write

    session_dir = Path(session_dir)
    if session_dir.name == "sorting":
        session_dir = session_dir.parent
    label = " ".join(p for p in (session_dir.parent.parent.name,
                                 session_dir.parent.name, session_dir.name) if p)
    try:
        path = write(result.out_dir, title=f"{label} - QC figures",
                     shank_folders=shank_folders,
                     skipped=None if kept_reasons else result.skipped)
        result.made[PAGE_NAME] = str(path)
    except Exception as exc:  # noqa: BLE001 - recorded, never raised
        result.skipped[PAGE_NAME] = f"{type(exc).__name__}: {exc}"


def _session_figures(result: ReportResult, session_dir: Path, found,
                     noise: bool) -> None:
    from .plots import session as sessionplots

    out_dir = result.out_dir
    if not found.shanks:
        reason = "; ".join(found.incomplete) or "no sorted shanks found"
        result.skipped["session_drift_traces"] = reason
        result.skipped["session_raster_all_units"] = reason
        result.skipped["session_raster_passing"] = reason
        result.skipped["session_rms_across_shanks"] = reason
        return

    session_root = session_dir.parent if session_dir.name == "sorting" else session_dir
    label = f"{session_root.parent.parent.name} {session_root.parent.name}".strip()
    session_label = f"{label} {session_root.name}".strip()
    # Named on the figure, not only in a log beside it.
    caveat = (f"   Unusable and left out: {'; '.join(found.incomplete)}"
              if found.incomplete else "")
    totals = sum(int(s.quality_metrics().get("n_units_total", 0))
                 for s in found.shanks)

    _attempt(result, "session_drift_traces",
             lambda: sessionplots.session_drift_traces(
                 found.shanks,
                 title=f"{session_label}: DREDGE displacement, every shank".strip(),
                 subtitle=(f"{found.describe}. Each panel is one shank; colored "
                           "lines are the depth-window estimates. Time and "
                           f"displacement axes are shared across panels.{caveat}")))

    _attempt(result, "session_raster_all_units", lambda: sessionplots.session_raster(
        found.shanks,
        title=f"{label}: every sorted unit".strip(),
        subtitle=(f"{found.describe}. Every unit the sorter produced, with no "
                  f"quality filtering. Row shading is each unit's rate "
                  f"normalised to its own maximum.{caveat}")))

    passing = sum(len(s.passing_units("pass_rescued") or ()) for s in found.shanks)
    _attempt(result, "session_raster_passing", lambda: sessionplots.session_raster(
        found.shanks, criterion="pass_rescued",
        title=f"{label}: units passing quality control".strip(),
        subtitle=(f"{found.describe}. {passing} of {totals} units pass, by the "
                  f"verdicts recorded in quality_metrics.json at sort time; "
                  f"nothing is recomputed here.{caveat}")))

    if not noise:
        result.skipped["session_rms_across_shanks"] = "not requested"
        for name in NOISE_FIGURES:
            result.skipped[name] = "not requested"
        return

    _noise_detail(result, found, label, caveat)

    measured, failures = [], {}
    for shank in found.shanks:
        try:
            measured.append(sessionplots.channel_rms(shank))
        except Exception as exc:  # noqa: BLE001 - one shank must not stop the rest
            failures[shank.label] = f"{type(exc).__name__}: {exc}"
    if not measured:
        result.skipped["session_rms_across_shanks"] = (
            "; ".join(f"{k}: {v}" for k, v in failures.items())
            or "no shank yielded a measurement")
        return

    missing = f"   No measurement: {', '.join(failures)}" if failures else ""
    _attempt(result, "session_rms_across_shanks",
             lambda: sessionplots.rms_across_shanks(
                 measured,
                 title=f"{label}: recording noise, every shank".strip(),
                 subtitle=(f"{found.describe}. Per-site AP band RMS from "
                           f"{sessionplots.NOISE_WINDOWS} windows of "
                           f"{sessionplots.NOISE_WINDOW_S:g} s per shank. "
                           f"Dashed line is the session median."
                           f"{missing}{caveat}")))
    result.made["session_rms_per_shank.csv"] = sessionplots.noise_table(
        measured, out_dir / "session_rms_per_shank.csv")


#: The noise figures ported from the quad-base noise testing (plots/noise.py).
NOISE_FIGURES = ("session_noise_depth_power", "session_noise_snippets",
                 "session_noise_rms_by_channel", "session_noise_spectra",
                 "session_noise_rms_over_time", "session_noise_rms_summary")


def _noise_detail(result: ReportResult, found, label: str, caveat: str) -> None:
    """Measure every shank once, then draw the noise figures from it."""
    from .plots import noise

    from concurrent.futures import ThreadPoolExecutor

    def attempt(shank):
        try:
            return noise.measure_shank(shank), None
        except Exception as exc:  # noqa: BLE001 - one shank must not stop the rest
            return None, f"{type(exc).__name__}: {exc}"

    # Four shanks at once: the decompression, medians and filters do their
    # work outside the interpreter lock, and a 16-shank session took 23 min
    # one shank at a time.
    with ThreadPoolExecutor(max_workers=4) as pool:
        outcomes = list(pool.map(attempt, found.shanks))
    measured = [m for m, _err in outcomes if m is not None]
    failures = {s.label: err for s, (_m, err) in zip(found.shanks, outcomes) if err}
    if not measured:
        reason = ("; ".join(f"{k}: {v}" for k, v in failures.items())
                  or "no shank yielded a measurement")
        for name in NOISE_FIGURES:
            result.skipped[name] = reason
        return
    missing = f"   No measurement: {', '.join(failures)}" if failures else ""
    windows = _windows_text(noise)
    stages = NOISE_STAGES_TEXT
    spec = [
        ("session_noise_depth_power", noise.depth_power_grid,
         "power by depth and frequency, per stage and shank",
         f"{windows}. {stages} Log-frequency bins keep each bin's maximum."),
        ("session_noise_snippets", noise.snippet_grid,
         "33 ms of voltage, per stage and shank",
         f"From the middle window. {stages}"),
        ("session_noise_rms_by_channel", noise.rms_by_channel,
         "RMS per channel after each stage",
         f"{windows}. {stages} Stages overlaid within a shank only."),
        ("session_noise_spectra", noise.median_spectra,
         "median spectrum across channels, per stage",
         f"{windows}. {stages}"),
        ("session_noise_rms_over_time", noise.rms_over_time,
         "RMS over the recording, after demux CAR",
         f"{noise.N_TIME_WINDOWS} windows of 1 s spread across the recording; "
         "median across channels in each."),
        ("session_noise_rms_summary", noise.rms_summary,
         "median RMS per shank, 95% bootstrap CI",
         f"Median over {noise.N_TIME_WINDOWS} windows of 1 s after demux CAR; "
         "the CI resamples those windows."),
    ]
    for name, draw, what, how in spec:
        _attempt(result, name, lambda draw=draw, what=what, how=how: draw(
            measured, title=f"{label}: {what}".strip(),
            subtitle=f"{how}{missing}{caveat}"))
