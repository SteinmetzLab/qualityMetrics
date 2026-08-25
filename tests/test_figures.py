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


def test_amp_depth_scatter_has_a_size_legend(ks):
    """Marker area encodes spike count; without a key it is undecodable."""
    fig = plots.amp_depth_scatter(ks, duration_s=4.0)
    legends = [a.get_legend() for a in fig.axes if a.get_legend() is not None]
    assert legends, "no legend explaining marker size"
    assert "Spikes" in legends[0].get_title().get_text()


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
    assert plots.wall_heatmap(rec, t_start_s=1.0, win_ms=20) is not None
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
