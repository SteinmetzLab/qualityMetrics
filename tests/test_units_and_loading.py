"""The things that were wrong in the MATLAB version, pinned so they stay right.

Every test here corresponds to a specific defect: a time axis in samples, an
amplitude axis in raw integers, a hardcoded title, a metric silently dropped.
"""
from __future__ import annotations

import numpy as np
import pytest

from qualitymetrics import KilosortResults, RawRecording, parse_meta
from qualitymetrics.metrics import (
    CONTAMINATION_CENSORED,
    noise_cutoff,
    write_cluster_tsv,
)
from qualitymetrics.raw import geometry_from_meta, uv_per_bit

from .conftest import FS, N_CHANNELS, UV_PER_BIT

# ---------------------------------------------------------------- metadata

def test_uv_per_bit_matches_the_hand_calculation(recording_dir):
    meta = parse_meta(recording_dir.with_suffix(".meta"))
    assert uv_per_bit(meta, "ap") == pytest.approx(0.62 / 2048 / 100 * 1e6)
    assert uv_per_bit(meta, "ap") == pytest.approx(UV_PER_BIT, abs=1e-3)


def test_lf_band_falls_back_to_the_ap_gain(recording_dir):
    """Our derived LF meta carries only imChan0apGain, and must still open.

    Neuropixels 2.0 has no LF stream from the probe, so a pipeline that derives
    one writes the AP gain, which is the gain that actually applied. A reader
    that demands imChan0lfGain cannot open those files at all.
    """
    meta = parse_meta(recording_dir.with_suffix(".meta"))
    assert "imChan0lfGain" not in meta
    assert uv_per_bit(meta, "lf") == pytest.approx(uv_per_bit(meta, "ap"))


def test_uv_per_bit_refuses_to_guess():
    from qualitymetrics.raw import RawError

    with pytest.raises(RawError):
        uv_per_bit({"imAiRangeMax": "0.62"}, "ap")


def test_geometry_comes_out_aligned_to_channels(recording_dir):
    geom = geometry_from_meta(parse_meta(recording_dir.with_suffix(".meta")))
    assert geom.n_channels == N_CHANNELS
    assert set(np.unique(geom.x)) == {27.0, 59.0}
    assert geom.pitch() == pytest.approx(15.0)


# ---------------------------------------------------------------- recording

def test_traces_are_microvolts_not_integers(recording_dir):
    with RawRecording.open(recording_dir) as rec:
        uv = rec.get_traces(0.0, 0.1)
        raw = rec.get_traces(0.0, 0.1, unit="int16")
        assert raw.dtype == np.int16
        assert uv.dtype == np.float32
        np.testing.assert_allclose(uv, raw * rec.uv_per_bit, rtol=1e-5)


def test_window_is_requested_in_seconds(recording_dir):
    with RawRecording.open(recording_dir) as rec:
        block = rec.get_traces(1.0, 1.5)
        assert block.shape[0] == pytest.approx(0.5 * FS, abs=1)
        assert rec.duration_s == pytest.approx(4.0, abs=0.01)


def test_a_window_off_the_end_is_an_error_not_an_empty_array(recording_dir):
    with RawRecording.open(recording_dir) as rec:
        with pytest.raises(ValueError):
            rec.get_traces(100.0, 101.0)


# ---------------------------------------------------------------- sorter

def test_spike_times_are_seconds(sorter_output):
    ks = KilosortResults.load(sorter_output)
    t = ks.spike_times_s
    # The whole recording is 4 s. Sample indices would run to ~120,000.
    assert t.max() < 5.0
    assert t.min() >= 0.0


def test_amplitudes_require_a_gain_rather_than_defaulting_to_one(sorter_output):
    from qualitymetrics.ksdata import KsError

    ks = KilosortResults.load(sorter_output)          # no uv_per_bit
    with pytest.raises(KsError):
        _ = ks.unit_amplitude_uv


def test_unit_amplitudes_are_physically_plausible(sorter_output):
    ks = KilosortResults.load(sorter_output, uv_per_bit=UV_PER_BIT)
    amps = np.array(list(ks.unit_amplitude_uv.values()))
    assert amps.size == 5
    # Extracellular spikes are tens to hundreds of microvolts. The Kilosort-2
    # recipe applied to Kilosort 4 output lands in the thousands, which is the
    # error this pins against.
    assert amps.min() > 1.0
    assert amps.max() < 2000.0


def test_template_anchored_amplitudes_agree_with_the_unit_metric(sorter_output):
    """The per-spike scale must reconcile with the per-unit one, by design."""
    ks = KilosortResults.load(sorter_output, uv_per_bit=UV_PER_BIT)
    per_spike = ks.spike_amplitudes_uv("template_anchored")
    for uid, target in ks.unit_amplitude_uv.items():
        mine = per_spike[ks.spike_clusters == uid]
        assert mine.mean() == pytest.approx(target, rel=1e-6)


def test_trough_to_peak_is_measured_and_never_a_bare_zero(sorter_output):
    ks = KilosortResults.load(sorter_output, uv_per_bit=UV_PER_BIT)
    values = np.array(list(ks.trough_to_peak_ms.values()))
    finite = values[np.isfinite(values)]
    assert finite.size > 0
    # Zero would place a unit at the extreme narrow-spiking end of every figure.
    assert (finite > 0).all()
    assert (finite < 3.0).all()


def test_titles_come_from_the_data(sorter_output):
    ks = KilosortResults.load(sorter_output, uv_per_bit=UV_PER_BIT,
                              label="SUBJ 2026-01-01 imec0_shank2")
    title = ks.title("Drift")
    assert "SUBJ" in title and "shank2" in title
    assert "probe00a" not in title


def test_firing_rate_uses_the_duration_it_is_given(sorter_output):
    ks = KilosortResults.load(sorter_output, uv_per_bit=UV_PER_BIT)
    short = ks.firing_rate(10.0)
    long = ks.firing_rate(100.0)
    for uid in short:
        assert short[uid] == pytest.approx(10 * long[uid])


# ---------------------------------------------------------------- metrics

def test_noise_cutoff_pinned_against_upstream():
    """Pin the vendored copy so a divergence from ibllib is visible."""
    rng = np.random.default_rng(3)
    clean = rng.normal(100, 20, 5000)
    passed, cutoff, _low = noise_cutoff(clean)
    assert np.isfinite(cutoff)
    assert passed

    truncated = clean[clean > 95]          # amplitude distribution cut off
    passed_t, cutoff_t, _ = noise_cutoff(truncated)
    assert cutoff_t > cutoff
    assert not passed_t


def test_censored_contamination_sorts_as_worst_rather_than_vanishing():
    """A unit nothing could be confirmed about must not read as unmeasured."""
    assert CONTAMINATION_CENSORED > 35.0


def test_cluster_tsv_is_the_format_phy_reads(tmp_path):
    path = write_cluster_tsv(tmp_path, "AmpuV", {3: 12.345, 1: 7.0}, "{:.1f}")
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    assert path.name == "cluster_AmpuV.tsv"
    assert lines[0] == "cluster_id\tAmpuV"
    # Sorted by cluster id, tab separated, formatted as asked.
    assert lines[1] == "1\t7.0"
    assert lines[2] == "3\t12.3"


def test_cluster_tsv_omits_unmeasurable_units_rather_than_writing_zero(tmp_path):
    """Blank and zero mean different things and must not be conflated."""
    path = write_cluster_tsv(tmp_path, "NoiseCutoff",
                             {0: 1.0, 1: float("nan"), 2: None}, "{:.2f}")
    body = path.read_text(encoding="utf-8").strip().splitlines()[1:]
    assert body == ["0\t1.00"]
