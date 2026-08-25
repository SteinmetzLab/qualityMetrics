"""Round 4: waveform spacing, the noise cutoff diagnostic, and wideband power."""
from __future__ import annotations

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from qualitymetrics import KilosortResults, RawRecording, plots  # noqa: E402

from .conftest import UV_PER_BIT  # noqa: E402


@pytest.fixture(autouse=True)
def _close_figures():
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


def test_waveform_rows_keep_whitespace_between_channels(rec, ks):
    """The trace cloud must not run into the row above it.

    Scaling only by the peak made a quiet unit's noise cloud span three rows.
    Scaling only by the noise threw a loud unit's mean many rows clear of its
    own channel. The rule takes whichever of the two is tighter.
    """
    import matplotlib.pyplot as plt

    from qualitymetrics.plots.rawviews import _site_spacing

    uid = int(ks.unit_ids[0])
    t_ms, stack, chans = plots.spike_snippets(rec, ks, uid, n_waveforms=30,
                                              n_channels=12)
    pos = np.asarray(ks.channel_positions)[np.asarray(chans, dtype=int)]
    _dx, dy = _site_spacing(pos)

    _fig, ax = plt.subplots()
    y_scale = plots.draw_waveforms_on_probe(ax, ks, t_ms, stack, chans)

    noise_uv = float(np.median(np.std(stack, axis=0)))
    # The cloud spans roughly plus or minus 2.5 standard deviations, and that
    # has to fit inside one row for the channels to read separately.
    assert 5.0 * noise_uv * y_scale < dy, "trace clouds would touch"


def test_noise_cutoff_diagnostic_redraws_the_metric_faithfully(ks):
    """The figure recomputes the intermediates; they must match the metric.

    metrics.noise_cutoff is vendored verbatim from ibllib and must stay that
    way, so the drawing keeps its own copy of the arithmetic. If the two drift
    apart, the picture is of a different calculation than the number beside it.
    """
    from qualitymetrics.metrics import noise_cutoff
    from qualitymetrics.plots.diagnostics import _noise_cutoff_parts

    amps = ks.spike_amplitudes_uv()
    checked = 0
    for uid in ks.unit_ids:
        mine = amps[ks.spike_clusters == uid]
        if mine.size < 50:
            continue
        _passed, cutoff, first_low = noise_cutoff(mine)
        parts = _noise_cutoff_parts(mine)
        if np.isfinite(cutoff):
            assert parts["cutoff"] == pytest.approx(cutoff, rel=1e-9)
            assert parts["first_low"] == pytest.approx(first_low, rel=1e-9)
            checked += 1
    assert checked > 0, "no unit was actually compared"


def test_noise_cutoff_diagnostic_builds_and_labels_its_pieces(ks):
    fig = plots.noise_cutoff_diagnostic(ks, n_units=2, time_bins=20,
                                        amp_bins=20)
    text = " ".join(t.get_text() for t in fig.findobj(
        match=lambda o: hasattr(o, "get_text")))
    assert "noise cutoff" in text.lower()
    assert "first low bin" in text.lower()
    assert "µV" in text


def test_units_spanning_noise_cutoff_returns_distinct_units(ks):
    picks = plots.units_spanning_noise_cutoff(ks, n_units=3, min_spikes=1)
    assert 1 <= len(picks) <= 3
    assert len(set(picks)) == len(picks)


def test_wideband_depth_power_runs_to_nyquist_on_a_log_axis(rec):
    """The point of this one is the noise above 300 Hz, so it has to reach it."""
    fig = plots.depth_power_wide(rec, t_start_s=0.5, win_s=2.0)
    ax = fig.axes[0]
    assert ax.get_xscale() == "log"
    assert ax.get_xlim()[1] > rec.fs / 2 * 0.9, "does not reach Nyquist"
    assert "µm" in ax.get_ylabel()


def test_the_narrow_depth_power_is_still_there(rec):
    """Keeping both was the request; the wideband one does not replace it."""
    fig = plots.depth_power(rec, t_start_s=0.5, win_s=1.0)
    assert fig.axes[0].get_xscale() == "linear"
