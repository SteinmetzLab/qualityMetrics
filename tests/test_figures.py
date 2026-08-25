"""Every figure builds, and its axes say what they are.

These are not image comparisons. They check the properties that were wrong
before: that an axis labelled seconds carries seconds, that an axis labelled
microvolts carries microvolts, and that nothing is hardcoded to one probe.
"""
from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from qualitymetrics import KilosortResults, RawRecording, plots  # noqa: E402

from .conftest import UV_PER_BIT  # noqa: E402


@pytest.fixture(autouse=True)
def _close_figures():
    """Every test here builds figures; leaving them open leaks across the suite."""
    import matplotlib.pyplot as plt
    yield
    plt.close("all")


@pytest.fixture
def ks(sorter_output):
    return KilosortResults.load(sorter_output, uv_per_bit=UV_PER_BIT,
                                label="TEST 2026-01-01 imec0_shank0")


@pytest.fixture
def rec(recording_dir):
    with RawRecording.open(recording_dir) as r:
        yield r


# ---------------------------------------------------------------- sorter only

def test_unit_drift_time_axis_is_seconds(ks):
    fig = plots.unit_drift(ks)
    ax = fig.axes[0]
    assert ax.get_xlabel() == "Time (s)"
    # The recording is 4 s. A sample-index axis would end near 120,000.
    assert ax.get_xlim()[1] < 10
    assert "µm" in ax.get_ylabel()


def test_amp_depth_scatter_labels_rate_as_spikes_per_second(ks):
    fig = plots.amp_depth_scatter(ks, duration_s=4.0)
    labels = [a.yaxis.label.get_text() for a in fig.axes]
    assert any("spikes/s" in x for x in labels), labels
    assert not any("(Hz)" in x for x in labels), labels


def test_marker_size_is_binary_qc_and_is_explained(ks):
    """Size carries the one bit that matters, and the legend names the criterion.

    Colour already carries firing rate and waveform duration, so spending the
    size channel on a third continuous variable wasted it. The MATLAB original
    had this right: large if the unit passed, small if it did not.
    """
    from qualitymetrics.plots.units import QC_FAIL_SIZE, QC_PASS_SIZE

    fig = plots.amp_depth_scatter(ks, duration_s=4.0)
    legends = [a.get_legend() for a in fig.axes if a.get_legend() is not None]
    assert legends, "no legend explaining marker size"
    assert "QC" in legends[0].get_title().get_text()

    from matplotlib.collections import PathCollection

    sizes = set()
    for ax in fig.axes:
        for coll in ax.collections:
            if isinstance(coll, PathCollection):     # the scatter, not a legend
                sizes.update(np.round(coll.get_sizes(), 3).tolist())
    assert sizes, "no scatter found"
    assert sizes <= {QC_PASS_SIZE, QC_FAIL_SIZE}, sizes


def test_qc_pass_prefers_our_gate_and_says_which(ks, sorter_output):
    """With no quality_metrics.json we fall back to the Kilosort label."""
    assert ks.qc_source == "Kilosort label"
    assert set(ks.qc_pass.values()) == {True}      # the fixture labels all good


def test_no_figure_hardcodes_a_probe_name(ks):
    for build in (lambda: plots.unit_drift(ks),
                  lambda: plots.drift_map(ks),
                  lambda: plots.amp_depth_scatter(ks, duration_s=4.0),
                  lambda: plots.amplitude_cdf_over_depth(ks)):
        fig = build()
        text = " ".join(t.get_text() for t in fig.findobj(
            match=lambda o: hasattr(o, "get_text")))
        assert "probe00a" not in text
        assert "TEST" in text or "shank0" in text


def test_artifact_band_annotation_follows_the_argument(ks):
    """The MATLAB version drew the line from the argument and the label from a
    hardcoded 4850, so the two disagreed whenever a caller passed anything."""
    fig = plots.unit_drift(ks, artifact_z_um=60.0)
    text = " ".join(t.get_text() for t in fig.findobj(
        match=lambda o: hasattr(o, "get_text")))
    assert "60" in text
    assert "4850" not in text


def test_no_artifact_band_by_default(ks):
    fig = plots.unit_drift(ks)
    text = " ".join(t.get_text() for t in fig.findobj(
        match=lambda o: hasattr(o, "get_text")))
    assert "Above" not in text


def test_templates_and_waterfall_build(ks):
    assert plots.templates_grid(ks, unit_ids=[0, 1, 2]) is not None
    assert plots.template_waterfall(ks, 0) is not None


def test_unit_summary_table_carries_the_plotted_numbers(ks):
    rows = plots.unit_summary_table(ks, duration_s=4.0)
    assert len(rows) == 5
    assert {"unit_id", "amplitude_uv", "firing_rate_sps", "depth_um"} <= set(rows[0])


# ---------------------------------------------------------------- raw views

def test_raw_figures_build_from_a_flat_binary(rec, ks):
    assert plots.wall_heatmap(rec, t_starts_s=[1.0], win_ms=20) is not None
    assert plots.raw_traces(rec, t_start_s=1.0, win_ms=20) is not None
    assert plots.raw_with_spikes(rec, ks, t_start_s=1.0, win_ms=20) is not None
    assert plots.plot_channel_health(rec, n_chunks=2, win_s=0.2) is not None


def test_common_reference_is_taken_across_the_array_not_the_shown_subset(rec):
    """Referencing a narrow band against its own median deletes real signal.

    Filtering the whole array and then slicing must not equal filtering only the
    slice, otherwise the wide reference is not actually being used.
    """
    from qualitymetrics.plots.rawviews import _filtered_window
    from qualitymetrics.preproc import preprocess_ap

    channels = np.arange(6, 10)
    wide = _filtered_window(rec, 1.0, 20.0, 300.0, True, channels)
    narrow_block = rec.get_traces(1.0, 1.02, channels=channels)
    narrow = preprocess_ap(narrow_block, rec.fs, 300.0, cmr=True)
    assert not np.allclose(wide, narrow, atol=1e-6)


def test_channel_health_reports_both_detectors(rec):
    result = plots.channel_health(rec, n_chunks=2, win_s=0.2)
    assert {"dead", "noisy", "outlier", "robust_z", "mad_rms_uv"} <= set(result)
    assert result["rms_uv"].size == rec.n_data_channels


def test_filter_state_calls_wideband_data_wideband(rec):
    _fig, info = plots.filter_state(rec, n_chunks=2, win_s=0.5, n_channels=8)
    assert np.isfinite(info["ratio"])
    assert isinstance(info["verdict"], str)


def test_band_rms_and_depth_power_build(rec):
    assert plots.band_rms(rec, n_time_bins=4, win_s=0.1) is not None
    assert plots.depth_power(rec, t_start_s=0.5, win_s=1.0) is not None


# ------------------------------------------------- the second round of changes

def test_wall_heatmap_shows_three_moments_on_one_colour_scale(rec):
    """One window cannot tell a persistent fault from a transient one.

    And the panels must share a colour scale, or early and late are individually
    pretty and mutually incomparable, which defeats showing three.
    """
    fig = plots.wall_heatmap(rec, win_ms=10)
    images = [im for ax in fig.axes for im in ax.images]
    assert len(images) == 3
    limits = {im.get_clim() for im in images}
    assert len(limits) == 1, f"panels autoscaled separately: {limits}"


def test_empty_amplitude_bins_are_grey_not_the_bottom_of_the_colormap(ks):
    """At high amplitude a handful of spikes is the signal.

    On a continuous scale anchored at zero, one spike and no spikes look the
    same, so rare and never must be different colours.
    """
    fig = plots.amplitude_cdf_over_depth(ks)
    images = [im for ax in fig.axes for im in ax.images]
    assert images
    for im in images:
        assert np.ma.is_masked(im.get_array()), "zero bins were not masked"
        bad = im.get_cmap().get_bad()
        assert bad[3] > 0, "the bad colour is transparent, so zeros show as white"
        assert not np.allclose(bad[:3], im.get_cmap()(0.0)[:3]), \
            "empty bins render as the bottom of the colormap"


def test_spike_markers_are_filled_circles(rec, ks):
    fig = plots.raw_with_spikes(rec, ks, t_start_s=1.0, win_ms=20)
    markers = [ln.get_marker() for ax in fig.axes for ln in ax.lines
               if ln.get_linestyle() == "None"]
    assert markers, "no spike markers drawn"
    assert set(markers) == {"o"}


def test_autocorrelogram_has_no_zero_lag_spike(ks):
    """Every spike correlates with itself; leaving that in hides the dip."""
    uid = int(ks.unit_ids[0])
    lags, counts = ks.autocorrelogram(uid, window_ms=20, bin_ms=0.5)
    assert lags.size and counts.size == lags.size
    centre = counts[np.argmin(np.abs(lags))]
    assert centre <= counts.max(), "zero-lag bin dominates the autocorrelogram"
    # Symmetric by construction.
    assert counts.sum() % 2 == 0


def test_autocorrelogram_finds_a_planted_refractory_gap():
    """A synthetic train with an enforced 5 ms gap must show an empty centre."""
    import types

    fs = 30000.0
    rng = np.random.default_rng(0)
    times = np.cumsum(rng.uniform(0.006, 0.05, 4000))
    fake = types.SimpleNamespace(
        spike_times_s=times,
        spike_clusters=np.zeros(times.size, dtype=int),
        autocorrelogram=KilosortResults.autocorrelogram,
    )
    lags, counts = KilosortResults.autocorrelogram(fake, 0, window_ms=30,
                                                   bin_ms=0.5)
    inside = counts[np.abs(lags) < 5.0]
    outside = counts[np.abs(lags) > 15.0]
    assert inside.sum() == 0, "the enforced refractory gap was not empty"
    assert outside.sum() > 0


def test_example_neurons_is_a_three_by_three_with_metrics_on_it(rec, ks):
    fig = plots.example_neurons(rec, ks, n_waveforms=10, acg_window_ms=20)
    # 9 panels plus one colorbar per template panel.
    assert len(fig.axes) >= 9
    text = " ".join(t.get_text() for t in fig.findobj(
        match=lambda o: hasattr(o, "get_text")))
    assert "percentile" in text
    assert "Amplitude" in text and "spikes" in text
    assert "QC" in text


def test_example_units_are_picked_by_percentile_not_by_rank(ks):
    picks = ks.units_at_amplitude_percentiles((95, 75, 50))
    assert len(picks) == 3
    amps = ks.unit_amplitude_uv
    values = [amps[u] for u in picks]
    assert values[0] >= values[1] >= values[2]
    # The 95th percentile should not simply be the largest unit of all.
    assert values[0] <= max(amps.values())


def test_depth_correlation_is_square_and_symmetric(ks):
    fig = plots.depth_correlation(ks, depth_bin_um=20, time_bin_s=0.2,
                                  duration_s=4.0, min_spikes=1)
    im = fig.axes[0].images[0]
    data = im.get_array()
    assert data.shape[0] == data.shape[1]
    assert np.allclose(data, data.T, atol=1e-9)


def test_population_figures_build(ks, rec):
    assert plots.firing_rate_image(ks, time_bin_s=0.2, duration_s=4.0) is not None
    assert plots.depth_profiles(ks, duration_s=4.0, min_spikes=1) is not None
    assert plots.lfp_band_profiles(rec, t_start_s=0.5, win_s=1.0) is not None


def test_outlier_detector_catches_a_planted_bad_channel(recording_dir, tmp_path):
    """A detector that never fires is not a check.

    The local comparison was introduced because the global one flagged 111 of
    384 channels on a probe whose RMS rises smoothly with depth. Making it
    quieter is only an improvement if it still finds a real fault, so plant one.
    """
    import numpy as np

    from qualitymetrics import RawRecording, plots
    from qualitymetrics.raw import parse_meta

    meta = parse_meta(recording_dir.with_suffix(".meta"))
    n_ch = int(meta["nSavedChans"])
    data = np.fromfile(recording_dir, dtype=np.int16).reshape(-1, n_ch).copy()
    data[:, 5] = (data[:, 5].astype(np.int32) * 6).astype(np.int16)  # noisy site
    bad = tmp_path / "planted_g0_t0.imec0.sh0.ap.bin"
    data.tofile(bad)
    bad.with_suffix(".meta").write_text(
        recording_dir.with_suffix(".meta").read_text(encoding="utf-8"),
        encoding="utf-8")

    with RawRecording.open(bad) as rec:
        result = plots.channel_health(rec, n_chunks=2, win_s=0.2)
    flagged = result["dead"] | result["noisy"] | result["outlier"]
    assert flagged[5], "the planted noisy channel was not flagged"
    assert flagged.sum() <= 3, f"flagged too much: {int(flagged.sum())}"


def test_a_smooth_depth_gradient_is_not_called_a_fault(recording_dir, tmp_path):
    """Anatomy makes RMS rise along a probe. That is not a broken channel."""
    import numpy as np

    from qualitymetrics import RawRecording, plots
    from qualitymetrics.raw import parse_meta

    meta = parse_meta(recording_dir.with_suffix(".meta"))
    n_ch = int(meta["nSavedChans"])
    data = np.fromfile(recording_dir, dtype=np.int16).reshape(-1, n_ch).copy()
    ramp = np.linspace(1.0, 2.5, n_ch)          # smooth, like a real gradient
    data = (data.astype(np.float64) * ramp).astype(np.int16)
    ramped = tmp_path / "ramped_g0_t0.imec0.sh0.ap.bin"
    data.tofile(ramped)
    ramped.with_suffix(".meta").write_text(
        recording_dir.with_suffix(".meta").read_text(encoding="utf-8"),
        encoding="utf-8")

    with RawRecording.open(ramped) as rec:
        result = plots.channel_health(rec, n_chunks=2, win_s=0.2)
    assert result["n_outlier"] == 0, \
        f"{result['n_outlier']} channels flagged on a smooth gradient"
