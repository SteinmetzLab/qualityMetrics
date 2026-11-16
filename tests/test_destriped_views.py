"""The 2.4 figure changes (Nick, 2026-11-15): QC-passing drift figures,
destriped traces of every channel, destriped depth power, per-stage spectra
for one shank, templates on 32 sites at their positions, log amplitudes.

What is checked is state, not wording: which spikes and channels a figure
drew, from which data, at which depth.
"""
from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from qualitymetrics import KilosortResults, RawRecording, plots  # noqa: E402

from .conftest import UV_PER_BIT  # noqa: E402
from .test_wall_heatmap import UV_PER_BIT as WALL_UV_PER_BIT  # noqa: E402
from .test_wall_heatmap import _keep, session  # noqa: E402,F401  (fixture)


@pytest.fixture(autouse=True)
def _close_figures():
    import matplotlib.pyplot as plt
    yield
    plt.close("all")


@pytest.fixture
def ks(sorter_output):
    return KilosortResults.load(sorter_output, uv_per_bit=UV_PER_BIT,
                                label="TEST 2026-01-01 imec0_shank0")


def _scatter_points(ax):
    from matplotlib.collections import PathCollection
    return np.vstack([c.get_offsets() for c in ax.collections
                      if isinstance(c, PathCollection)])


# ---------------------------------------------------------------- drift

def test_qc_passing_drift_figures_draw_only_those_units_spikes(ks):
    chosen = [1, 3]
    wanted = np.isin(ks.spike_clusters, chosen)
    for build in (plots.drift_map, plots.unit_drift):
        fig = build(ks, unit_ids=chosen, max_spikes=None)
        points = _scatter_points(fig.axes[0])
        assert len(points) == int(wanted.sum())
        # Every drawn depth is one of the chosen units' spike depths.
        assert np.isin(np.round(points[:, 1], 3),
                       np.round(ks.spike_depths_um[wanted], 3)).all()


def test_qc_passing_drift_figures_keep_the_whole_shanks_depth_range(ks):
    """So the QC-passing version lies over the every-unit one."""
    full = plots.drift_map(ks).axes[0].get_ylim()
    only = plots.drift_map(ks, unit_ids=[0]).axes[0].get_ylim()
    assert only == pytest.approx(full)


def test_an_empty_selection_is_refused_not_drawn_empty(ks):
    with pytest.raises(ValueError):
        plots.drift_map(ks, unit_ids=[])


def test_report_draws_qc_passing_versions_of_only_passing_units(
        sorter_output, recording_dir, tmp_path):
    from qualitymetrics.report import build_report

    # Units 0 and 1 fail by Kilosort's label; 2 to 4 pass.
    (sorter_output / "cluster_KSLabel.tsv").write_text(
        "cluster_id\tKSLabel\n0\tmua\n1\tmua\n2\tgood\n3\tgood\n4\tgood\n",
        encoding="utf-8")
    ks = KilosortResults.load(sorter_output, uv_per_bit=UV_PER_BIT)
    assert plots.passing_units(ks) == [2, 3, 4]
    result = build_report(sorter_output, tmp_path / "qc", raw_path=recording_dir)
    assert "drift_map_qc-passing" in result.made
    assert "unit_drift_qc-passing" in result.made


# ---------------------------------------------------------------- amplitudes

def test_amplitude_axis_is_logarithmic(ks):
    fig = plots.amp_depth_scatter(ks, duration_s=4.0)
    scatter_axes = [a for a in fig.axes if a.get_xlabel()]
    assert scatter_axes and all(a.get_xscale() == "log" for a in scatter_axes)


# ---------------------------------------------------------------- templates

def test_templates_are_drawn_at_their_sites_positions(ks):
    """Each line starts at its site's x and sits around its site's depth."""
    pos = np.asarray(ks.channel_positions, dtype=float)
    fig = plots.templates_grid(ks, unit_ids=[2])
    lines = [ln for ln in fig.axes[0].get_lines() if len(ln.get_xdata()) > 2]
    # 32 nearest sites, or every site when the probe has fewer (16 here).
    assert len(lines) == min(32, len(pos))
    starts = {(round(float(ln.get_xdata()[0]), 3), round(float(np.median(ln.get_ydata())), 1))
              for ln in lines}
    sites = {(round(float(x), 3), round(float(y), 1)) for x, y in pos}
    assert starts == sites


# ---------------------------------------------------------------- destriped

def test_raw_traces_show_every_kept_channel_from_the_destriped_data(session):  # noqa: F811 (fixture)
    """Each trace is the kept destriped samples of one channel, at its depth."""
    path, sort_dir, data = session
    kept = _keep(sort_dir, data)
    with RawRecording.open(path) as rec:
        fig = plots.raw_traces(rec, cleaned=kept, win_ms=20.0)
        ax = fig.axes[0]
        traces = [ln for ln in ax.get_lines() if len(ln.get_xdata()) > 2]
        assert len(traces) == kept.channels.size
        # The bottom trace is the shallowest channel; its samples are the
        # kept windows', not the archive's.
        geom = rec.geometry
        shallowest = int(np.argmin(geom.y[kept.channels]))
        middle = kept.n_windows // 2
        lo, hi = kept.core(middle)
        width = int(round(0.020 * rec.fs))
        i0 = lo + (hi - lo - width) // 2
        expected = kept.samples_uv(i0, i0 + width)[:, shallowest]
        assert np.allclose(traces[0].get_ydata(), expected, atol=1e-3)


def test_raw_traces_without_kept_windows_draw_every_archive_channel(session):  # noqa: F811 (fixture)
    path, _sort_dir, _data = session
    with RawRecording.open(path) as rec:
        fig = plots.raw_traces(rec, t_start_s=1.0, win_ms=20.0)
        traces = [ln for ln in fig.axes[0].get_lines() if len(ln.get_xdata()) > 2]
        assert len(traces) == rec.n_data_channels


def test_destriped_depth_power_puts_both_panels_on_one_color_scale(session):  # noqa: F811 (fixture)
    path, sort_dir, data = session
    kept = _keep(sort_dir, data)
    with RawRecording.open(path) as rec:
        fig = plots.depth_power_wide_destriped(rec, kept)
    from matplotlib.collections import QuadMesh
    meshes = [c for ax in fig.axes if ax.get_label() != "<colorbar>"
              for c in ax.collections if isinstance(c, QuadMesh)]
    assert len(meshes) == 2
    assert meshes[0].get_clim() == meshes[1].get_clim()


def test_destriped_depth_power_refuses_windows_without_channel_numbers(session):  # noqa: F811 (fixture)
    path, sort_dir, data = session
    kept = _keep(sort_dir, data, acquired_names=True)
    with RawRecording.open(path) as rec:
        with pytest.raises(ValueError):
            plots.depth_power_wide_destriped(rec, kept)


def test_filter_state_has_every_stage_for_one_shank(session, tmp_path):  # noqa: F811 (fixture)
    """The session figure's five stages, the destripe's from the kept file."""
    from qualitymetrics.plots import noise

    path, sort_dir, data = session
    _keep(sort_dir, data)
    out = tmp_path / "ks"
    out.mkdir()
    n = data.shape[1]
    np.save(out / "channel_positions.npy",
            np.column_stack([np.zeros(n), np.arange(n) * 15.0]))
    np.save(out / "channel_map.npy", np.arange(n, dtype=np.int32))
    np.save(out / "spike_times.npy", np.array([100], dtype=np.int64))
    np.save(out / "spike_clusters.npy", np.array([0], dtype=np.int32))
    np.save(out / "amplitudes.npy", np.array([10.0], dtype=np.float32))
    np.save(out / "templates.npy", np.zeros((1, 61, n), dtype=np.float32))
    (out / "params.py").write_text(f"n_channels_dat = {n}\nsample_rate = 30000.0\n"
                                   "dtype = 'int16'\n", encoding="utf-8")
    ks = KilosortResults.load(out, uv_per_bit=WALL_UV_PER_BIT)
    with RawRecording.open(path) as rec:
        measured = noise.shank_measurement(rec, ks, sort_dir)
        assert rec._reader is not None or rec._memmap is not None  # left open
    assert all(spectrum is not None for spectrum in measured["median_psd_db"])
    assert len(measured["median_psd_db"]) == len(noise.STAGES)
    fig = noise.median_spectra([measured], panel_w=7.0, panel_h=4.6)
    lines = [ln for ln in fig.axes[0].get_lines() if len(ln.get_xdata()) > 2]
    assert len(lines) == len(noise.STAGES)


# ------------------------------------------------- review findings, 2026-11-15

def test_a_probe_wide_sorts_traces_still_come_from_the_archive(session, tmp_path):  # noqa: F811 (fixture)
    """Kept channels named AP<n> have no depths here: raw_traces is drawn
    from the archive, not skipped."""
    import shutil

    from qualitymetrics.report import build_report

    path, sort_dir, data = session
    _keep(sort_dir, data, acquired_names=True)
    shutil.copytree(_minimal_sorter_output(tmp_path, data.shape[1]),
                    sort_dir / "sorter_output")
    result = build_report(sort_dir, tmp_path / "qc", raw_path=path,
                          write_phy=False, example_units=1)
    assert "raw_traces" in result.made, result.skipped.get("raw_traces")


def test_unreadable_qc_verdicts_cost_only_the_qc_passing_figures(sorter_output,
                                                                 recording_dir, tmp_path):
    import json

    from qualitymetrics.report import build_report

    (sorter_output.parent / "quality_metrics.json").write_text(
        json.dumps({"units": [{"pass_strict": True}]}), encoding="utf-8")
    shutil_target = tmp_path / "sort"
    shutil_target.mkdir()
    import shutil
    shutil.copytree(sorter_output, shutil_target / "sorter_output")
    shutil.copy(sorter_output.parent / "quality_metrics.json", shutil_target)
    try:
        result = build_report(shutil_target, tmp_path / "qc", raw_path=recording_dir,
                              write_phy=False, example_units=1)
    except Exception as exc:  # noqa: BLE001
        pytest.fail(f"build_report raised {type(exc).__name__}: {exc}")
    assert "raw_traces" in result.made
    assert "drift_map_qc-passing" not in result.made


def test_a_passing_unit_with_no_spikes_is_left_out(ks):
    from qualitymetrics.plots.drift import _only_units

    mask, units = _only_units(ks, [2, 99])
    assert units == [2]
    assert mask.sum() == int((ks.spike_clusters == 2).sum())
    fig = plots.unit_drift(ks, unit_ids=[2, 99])
    assert "only the 1 of" in fig._suptitle.get_text()


def test_the_qc_passing_drift_map_shares_the_full_maps_shading(ks):
    from matplotlib.collections import PathCollection

    def clim(fig):
        return [c.get_clim() for c in fig.axes[0].collections
                if isinstance(c, PathCollection)][0]

    assert clim(plots.drift_map(ks, unit_ids=[0])) == pytest.approx(
        clim(plots.drift_map(ks)))


def _minimal_sorter_output(tmp_path, n):
    out = tmp_path / "minimal_ks"
    out.mkdir()
    np.save(out / "channel_positions.npy",
            np.column_stack([np.zeros(n), np.arange(n) * 15.0]))
    np.save(out / "channel_map.npy", np.arange(n, dtype=np.int32))
    np.save(out / "spike_times.npy", np.array([100, 5000], dtype=np.int64))
    np.save(out / "spike_clusters.npy", np.array([0, 0], dtype=np.int32))
    np.save(out / "spike_templates.npy", np.array([0, 0], dtype=np.int64))
    np.save(out / "amplitudes.npy", np.array([10.0, 12.0], dtype=np.float32))
    np.save(out / "templates.npy", np.ones((1, 61, n), dtype=np.float32))
    np.save(out / "whitening_mat_inv.npy", np.eye(n, dtype=np.float32))
    (out / "params.py").write_text(f"n_channels_dat = {n}\nsample_rate = 30000.0\n"
                                   "dtype = 'int16'\n", encoding="utf-8")
    return out
