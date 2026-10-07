"""The saturation figure: drawn from the worker's counts, and only when they exist."""
from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from qualitymetrics import plots  # noqa: E402
from qualitymetrics.plots.saturation import SATURATION_FILE, load_counts  # noqa: E402


@pytest.fixture(autouse=True)
def _close_figures():
    import matplotlib.pyplot as plt
    yield
    plt.close("all")


def _counts(n_ch=16, seconds=60, saturated=True):
    fs = 30000.0
    neg = np.zeros(n_ch, np.int64)
    pos = np.zeros(n_ch, np.int64)
    per_bin = np.zeros((seconds, n_ch), np.uint16)
    if saturated:
        neg[3] = 900
        pos[5] = 40
        neg[7] = seconds * int(fs)                       # pinned throughout
        per_bin[10, 3] = 900
        per_bin[20, 5] = 40
        per_bin[:, 7] = int(fs)
    return dict(negative=neg, positive=pos, events=(neg + pos > 0).astype(np.int64), per_bin=per_bin,
                bin_s=np.float64(1.0), threshold_bits=np.float64(1986.56), fs=np.float64(fs),
                n_samples=np.int64(seconds * fs))


def test_the_figure_shows_counts_per_channel_on_a_log_axis_with_zero_visible():
    fig = plots.saturation_report(_counts(), title="TEST")
    ax = fig.axes[0]
    assert ax.get_yscale() == "symlog"
    assert ax.get_ylim()[0] < 0.5, "channels that never saturated must stay on the axis"
    assert "Channel" in ax.get_xlabel() and "count" in ax.get_ylabel()
    assert "pinned" in ax.get_title() and "[7]" in ax.get_title()


def test_a_session_with_no_saturation_still_gets_its_figure():
    fig = plots.saturation_report(_counts(saturated=False))
    assert "No sample reached" in fig.axes[0].get_title()


def test_the_counts_are_read_from_the_sort_folder_and_absent_is_none(tmp_path):
    assert load_counts(tmp_path) is None
    np.savez(tmp_path / SATURATION_FILE, **_counts())
    got = load_counts(tmp_path)
    assert int(got["negative"][3]) == 900


def test_channels_are_drawn_in_depth_order_when_the_counts_carry_positions():
    from qualitymetrics.plots.saturation import depth_order

    counts = _counts(n_ch=8)
    # NP2 Quad order: channel number is not depth (0, 2880, 15, 2895, ...).
    counts["z_um"] = np.array([0.0, 2880.0, 15.0, 2895.0, 30.0, 2910.0, 45.0, 2925.0])
    counts["x_um"] = np.zeros(8)
    expected = [0, 2, 4, 6, 1, 3, 5, 7]                # tip first
    np.testing.assert_array_equal(depth_order(counts), expected)
    fig = plots.saturation_report(counts)
    line = fig.axes[0].lines[1]                      # the negative-rail points, as drawn
    neg = counts["negative"]
    np.testing.assert_array_equal(line.get_ydata(), neg[expected])
    assert "depth order" in fig.axes[0].get_xlabel()


def test_without_positions_the_channels_are_drawn_by_number():
    from qualitymetrics.plots.saturation import depth_order

    counts = _counts()
    assert depth_order(counts) is None
    fig = plots.saturation_report(counts)
    np.testing.assert_array_equal(fig.axes[0].lines[1].get_ydata(), counts["negative"])


def test_the_report_draws_the_figure_from_the_sort_folder(tmp_path, sorter_output):
    from qualitymetrics.report import build_report

    np.savez(tmp_path / SATURATION_FILE, **_counts(n_ch=16))
    result = build_report(tmp_path, tmp_path / "qc", write_phy=False, example_units=2)
    assert "saturation" in result.made, result.skipped.get("saturation")


def test_a_damaged_counts_file_costs_only_the_saturation_figure(tmp_path, sorter_output):
    from qualitymetrics.report import build_report

    (tmp_path / SATURATION_FILE).write_bytes(b"not a zip file")
    result = build_report(tmp_path, tmp_path / "qc", write_phy=False, example_units=2)
    assert "saturation" in result.skipped and "saturation" not in result.made
    assert result.made, "the other figures must still be made"


def test_a_detected_rail_is_named_and_an_assumed_threshold_is_flagged():
    detected = dict(_counts(), rail_bits=np.float64(1955.0), disagreeing=np.array([4]),
                    held_below_band=np.array([], int), minimum=np.full(16, -1955),
                    maximum=np.full(16, 1955))
    title = plots.saturation_report(detected).axes[0].get_title()
    assert "Rail +/-1955 ADC counts, found in the data" in title and "channels 4 hold a different rail" in title
    old = dict(_counts(), minimum=np.full(16, -1955), maximum=np.full(16, 1955))
    assert "Assumed threshold" in plots.saturation_report(old).axes[0].get_title()


def test_no_rail_found_says_so():
    none = dict(_counts(saturated=False), rail_bits=np.float64(np.nan))
    title = plots.saturation_report(none).axes[0].get_title()
    assert "No rail found: no channel held an extreme value" in title and "brief touch" in title


def test_long_channel_lists_are_shortened_in_the_title():
    many = dict(_counts(), rail_bits=np.float64(1999.0), disagreeing=np.arange(20))
    title = plots.saturation_report(many).axes[0].get_title()
    assert "and 12 more hold a different rail" in title
