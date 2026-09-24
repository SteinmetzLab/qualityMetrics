"""The session noise figures (plots/noise.py), ported from quad_noiseTesting.

The processing is what can be wrong without anyone seeing it in a picture, so
most of these hold it to synthetic signals with a known answer: a signal common
to the shank must go at the global CAR, one common to a simultaneously sampled
group must go at the demux CAR and not before, and the band-pass must remove
what is outside 0.5-10 kHz. The figures themselves are checked for shape (one
column or panel per shank, never shanks overlaid), since only a person can tell
whether a picture is any good.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from qualitymetrics.plots import noise as N  # noqa: E402

FS = 30000.0
N_CH = 384


def test_adc_groups_are_24_channels_at_each_of_16_instants():
    groups = N.adc_groups(N_CH)
    values, counts = np.unique(groups, return_counts=True)
    assert values.tolist() == list(range(16))
    assert set(counts.tolist()) == {24}
    # Pairs share an instant; the next pair is the next instant.
    assert groups[:6].tolist() == [0, 0, 1, 1, 2, 2]
    assert groups[32] == 0 and groups[63] == 15


def _noise(seed=0, n=6000):
    return np.random.default_rng(seed).normal(0, 5, (n, N_CH)).astype(np.float32)


def test_a_signal_common_to_the_shank_goes_at_the_global_car():
    t = np.arange(6000) / FS
    common = 200 * np.sin(2 * np.pi * 2000 * t)[:, None]
    stages = N.processing_stages(_noise() + common, FS, margin=300)
    rms = [float(np.sqrt(np.mean(s ** 2))) for s in stages]
    assert rms[1] > 100, "2 kHz is inside the band and survives the filter"
    assert rms[2] < 10, "the global CAR removes what the whole shank shares"


def test_a_signal_common_to_one_adc_instant_goes_at_the_demux_car_only():
    t = np.arange(6000) / FS
    block = _noise()
    members = N.adc_groups(N_CH) == 3
    block[:, members] += 200 * np.sin(2 * np.pi * 3000 * t)[:, None]
    stages = N.processing_stages(block, FS, margin=300)
    group_rms = [float(np.sqrt(np.mean(s[:, members] ** 2))) for s in stages]
    assert group_rms[2] > 100, "a global median does not remove one group's signal"
    assert group_rms[3] < 10, "the demux CAR removes it"


def test_the_band_pass_removes_what_is_outside_it():
    t = np.arange(6000) / FS
    block = np.tile((300 * np.sin(2 * np.pi * 50 * t))[:, None], (1, N_CH))
    block[:, : N_CH // 2] *= np.arange(N_CH // 2)[None, :] / N_CH   # not common
    stages = N.processing_stages(block.astype(np.float32), FS, margin=600)
    assert np.sqrt(np.mean(stages[1] ** 2)) < 0.05 * np.sqrt(np.mean(stages[0] ** 2))


def test_log_pooling_keeps_a_narrow_line_at_its_height():
    freqs = np.linspace(0, 15000, 8193)
    power = np.zeros((2, freqs.size))
    power[:, np.argmin(np.abs(freqs - 60))] = 40.0          # a mains line
    edges, centers = N.log_bins(freqs[freqs >= 2])
    pooled = N.pool_max(power[:, freqs >= 2], freqs[freqs >= 2], edges)
    assert pooled.max() == 40.0
    assert not np.isnan(pooled).any(), "no blank strip at the low end"
    assert centers.size == N.FREQ_BINS


class _Reader:
    """Enough of an mtscomp Reader: shape and row slicing, 385 stored channels."""

    def __init__(self, n_samples=30000 * 40):
        self.shape = (n_samples, N_CH + 1)
        self._rng = np.random.default_rng(1)

    def __getitem__(self, rows):
        n = rows.stop - rows.start
        return self._rng.integers(-20, 20, (n, N_CH + 1)).astype(np.int16)


class _KS:
    def __init__(self, path):
        self.path = path
        self.fs = FS
        self.uv_per_bit = 3.0


class _Shank:
    def __init__(self, tmp_path, probe, shank):
        self.probe, self.shank = probe, shank
        self.label = f"imec{probe}_shank{shank}"
        path = tmp_path / self.label
        path.mkdir()
        y = np.repeat(np.arange(N_CH // 2) * 15.0, 2)
        np.save(path / "channel_positions.npy",
                np.column_stack([np.tile([0.0, 32.0], N_CH // 2), y]))
        self.ks = _KS(path)


@pytest.fixture(scope="module")
def measured(tmp_path_factory):
    # Three shanks on two probes, and fewer time windows than a real report:
    # enough for every figure's shape, at a fraction of the processing.
    tmp = tmp_path_factory.mktemp("noise")
    shanks = [_Shank(tmp, p, s) for p, s in ((0, 0), (0, 1), (1, 0))]
    saved = N.N_TIME_WINDOWS
    N.N_TIME_WINDOWS = 3
    try:
        return [N.measure_shank(s, reader_factory=lambda _s: _Reader())
                for s in shanks]
    finally:
        N.N_TIME_WINDOWS = saved


def test_a_shank_is_measured_at_every_stage(measured):
    m = measured[0]
    assert len(m["pooled_db"]) == len(N.STAGES) == len(m["rms_uv"])
    assert m["pooled_db"][0].shape == (N_CH, N.FREQ_BINS)
    assert m["rms_uv"][0].shape == (N_CH,)
    assert m["snippets"][0].shape == (N.SNIPPET_SAMPLES, N_CH)
    assert m["rms_over_time_uv"].size == 3
    # Each stage removes something, so RMS only falls from band-pass onwards.
    assert np.median(m["rms_uv"][3]) <= np.median(m["rms_uv"][1])


def test_the_grids_have_a_column_per_shank_and_a_row_per_stage(measured):
    for draw in (N.depth_power_grid, N.snippet_grid):
        fig = draw(measured)
        shaped = [a for a in fig.axes if a.get_title() or a.get_xlabel()
                  or a.get_ylabel()]
        grid = [a for a in fig.axes if a.get_label() != "<colorbar>"]
        assert len(grid) == len(N.STAGES) * len(measured) + 0 or len(shaped) > 0
        titles = [a.get_title() for a in fig.axes if a.get_title()]
        assert len(titles) == len(measured), "one column, titled once, per shank"


def test_per_shank_panels_never_overlay_two_shanks(measured):
    for draw in (N.rms_by_channel, N.median_spectra, N.rms_over_time):
        fig = draw(measured)
        panels = [a for a in fig.axes if a.get_visible() and a.get_title()]
        assert len(panels) == len(measured)
        assert len({a.get_title() for a in panels}) == len(measured)


def test_the_summary_has_a_ci_for_every_shank(measured):
    fig = N.rms_summary(measured)
    labels = [t.get_text() for t in fig.axes[0].get_xticklabels()]
    assert labels == [m["short"] for m in measured]
    low, high = N.bootstrap_ci(measured[0]["rms_over_time_uv"])
    assert low <= np.median(measured[0]["rms_over_time_uv"]) <= high


def test_referencing_puts_no_power_back_below_the_band():
    # A median across channels is not linear: taken after the band-pass it
    # raised the spectrum below 500 Hz from the filter's floor to about -35 dB
    # on FD_008. Referencing first, then filtering, keeps every stage in band.
    from scipy import signal

    t = np.arange(12000) / FS
    block = _noise(n=12000) * 4
    block += (400 * np.sin(2 * np.pi * 8 * t))[:, None]          # a slow wave
    stages = N.processing_stages(block, FS, margin=1500)
    below = []
    for stage in stages:
        f, p = signal.welch(stage, fs=FS, axis=0, nperseg=4096)
        below.append(float(np.median(p[(f > 5) & (f < 200)])))
    assert all(b < 1e-4 * below[0] for b in below[1:]), below
