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
                         "depth_profiles", "depth_correlation"):
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
        biggest = _biggest_units(ks, 1)
        if biggest:
            _attempt(result, "template_waterfall",
                     lambda: plots.template_waterfall(ks, biggest[0]))

        # ---- figures that need the raw file
        if rec is None:
            reason = "no archived recording found for this shank"
            for name in ("wall_heatmap", "raw_traces", "raw_with_spikes",
                         "channel_health", "filter_state", "band_rms",
                         "sample_waveforms", "example_neurons"):
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
            if biggest:
                _attempt(result, "sample_waveforms",
                         lambda: plots.sample_waveforms(rec, ks, biggest[0]))
            _attempt(result, "example_neurons",
                     lambda: plots.example_neurons(rec, ks))

        power_rec = lf_rec or rec
        if power_rec is None:
            result.skipped["depth_power"] = "no recording found for this shank"
            result.skipped["lfp_band_profiles"] = "no recording for this shank"
        else:
            safe_start = min(t_start_s, power_rec.duration_s - 11)
            _attempt(result, "depth_power", lambda: plots.depth_power(
                power_rec, t_start_s=safe_start))
            _attempt(result, "lfp_band_profiles",
                     lambda: plots.lfp_band_profiles(power_rec,
                                                     t_start_s=safe_start))

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
