"""Build the whole figure set for one sorted shank.

Every figure is attempted independently and a failure is recorded rather than
raised, because a QC report that produces nothing when one panel cannot be drawn
is a QC report nobody runs. The return value says what was made and what was
not, so a missing figure is visible rather than silently absent.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from . import plots
from .ksdata import KilosortResults, label_from_path
from .metrics import phy_metrics_from_results
from .raw import RawError, RawRecording, find_band, parse_meta, uv_per_bit
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
    name = sorting_dir.name                          # e.g. imec0_shank0
    probe, _, shank_part = name.partition("_")
    try:
        shank = int("".join(c for c in shank_part if c.isdigit()))
    except ValueError:
        return None
    for parent in list(sorting_dir.parents)[:4]:
        try:
            return find_band(parent, probe, shank, band)
        except RawError:
            continue
    return None


def build_report(sorting_dir: str | Path, out_dir: str | Path,
                 raw_path: str | Path | None = None,
                 lf_path: str | Path | None = None,
                 artifact_z_um: float | None = None,
                 t_start_s: float | None = None,
                 write_phy: bool = True,
                 example_units: int = 12) -> ReportResult:
    """Make every figure that this shank's available data supports.

    raw_path and lf_path are found automatically from the session layout when
    not given. Figures that need the raw file are skipped, with a reason, when
    it cannot be found.
    """
    sorting_dir = Path(sorting_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    result = ReportResult(out_dir=out_dir)

    sorter_output = find_sorter_output(sorting_dir)
    if raw_path is None:
        raw_path = find_recording(sorting_dir, "ap")
    if lf_path is None:
        lf_path = find_recording(sorting_dir, "lf")

    # The gain lives in the recording's meta, and every microvolt axis needs it.
    scale = None
    if raw_path is not None:
        meta_path = Path(raw_path).with_suffix(".meta")
        if meta_path.exists():
            scale = uv_per_bit(parse_meta(meta_path), "ap")

    ks = KilosortResults.load(sorter_output, uv_per_bit=scale,
                              label=label_from_path(sorting_dir))

    rec = lf_rec = None
    try:
        if raw_path is not None:
            rec = RawRecording.open(raw_path)
        if lf_path is not None:
            lf_rec = RawRecording.open(lf_path)

        duration = rec.duration_s if rec is not None else None
        if t_start_s is None:
            # A window a third of the way in: past any settling transient, and
            # not so late that a short recording has nothing there.
            t_start_s = (duration or ks.duration_s) / 3.0

        # ---- sorter-only figures
        if scale is None:
            reason = ("no recording meta found, so microvolts are unknown; "
                      "amplitude figures need it")
            for name in ("amp_depth_scatter", "unit_drift", "drift_map",
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
            _attempt(result, "amplitude_cdf_over_depth",
                     lambda: plots.amplitude_cdf_over_depth(ks))
            _attempt(result, "firing_rate_image",
                     lambda: plots.firing_rate_image(ks, duration_s=duration))
            _attempt(result, "depth_profiles",
                     lambda: plots.depth_profiles(ks, duration_s=duration))
            _attempt(result, "depth_correlation",
                     lambda: plots.depth_correlation(ks, duration_s=duration))

        _attempt(result, "templates_grid", lambda: plots.templates_grid(
            ks, unit_ids=_biggest_units(ks, example_units)))
        if scale is not None:
            _attempt(result, "noise_cutoff_diagnostic",
                     lambda: plots.noise_cutoff_diagnostic(ks))

        # ---- figures that need the raw file
        if rec is None:
            reason = "no archived recording found for this shank"
            for name in ("wall_heatmap", "raw_traces", "raw_with_spikes",
                         "channel_health", "filter_state", "band_rms",
                         "example_neurons"):
                result.skipped[name] = reason
        else:
            _attempt(result, "wall_heatmap", lambda: plots.wall_heatmap(
                rec, artifact_z_um=artifact_z_um))
            _attempt(result, "raw_traces",
                     lambda: plots.raw_traces(rec, t_start_s=t_start_s))
            _attempt(result, "raw_with_spikes", lambda: plots.raw_with_spikes(
                rec, ks, t_start_s=t_start_s))
            _attempt(result, "channel_health",
                     lambda: plots.plot_channel_health(rec))
            _attempt(result, "filter_state",
                     lambda: plots.filter_state(rec)[0])
            _attempt(result, "band_rms",
                     lambda: plots.band_rms(rec, lf_rec=lf_rec))
            _attempt(result, "example_neurons",
                     lambda: plots.example_neurons(rec, ks))

        power_rec = lf_rec or rec
        if power_rec is None:
            result.skipped["depth_power"] = "no recording found for this shank"
            result.skipped["lfp_band_profiles"] = "no recording for this shank"
            result.skipped["depth_power_wide"] = "no recording for this shank"
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

    return result


def _biggest_units(ks: KilosortResults, n: int) -> list[int]:
    """The n units with the most spikes, which are the ones worth showing."""
    counts = ks.n_spikes
    return [u for u, _ in sorted(counts.items(), key=lambda kv: -kv[1])[:n]]


def build_session_report(session_dir: str | Path, out_dir: str | Path,
                         noise: bool = True) -> ReportResult:
    """Make the figures that describe a whole session rather than one shank.

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
    if not found.shanks:
        reason = "; ".join(found.incomplete) or "no sorted shanks found"
        result.skipped["session_raster_all_units"] = reason
        result.skipped["session_raster_passing"] = reason
        result.skipped["session_rms_across_shanks"] = reason
        return result

    label = f"{session_dir.parent.parent.name} {session_dir.parent.name}".strip()
    # Named on the figure, not only in a log beside it.
    caveat = (f"   Unusable and left out: {'; '.join(found.incomplete)}"
              if found.incomplete else "")
    totals = sum(int(s.quality_metrics().get("n_units_total", 0))
                 for s in found.shanks)

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
        return result

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
        return result

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
    return result


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
    windows = (f"1 s each at {', '.join(f'{int(f * 100)}%' for f in noise.WINDOW_FRACTIONS)} "
               f"of the recording")
    stages = ("Stages: raw (channel medians removed); band-pass 0.5-10 kHz; "
              "the shank's across-channel median removed (global CAR), then "
              "band-passed; then each simultaneously sampled group's median "
              "also removed (demux CAR), then band-passed.")
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
