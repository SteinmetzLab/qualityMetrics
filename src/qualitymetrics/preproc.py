"""Filtering applied to raw traces before they are looked at.

Neuropixels 2.0 saves wideband: the AP file carries the LFP and the per-channel
DC offset as well as the spikes. On our data the unfiltered median absolute
voltage is around 1100 uV, so any figure that draws raw traces without a
high-pass is drawing offsets, not spikes. Every viewer in this package filters
before display and says so in the axis label.

Arrays are (n_samples, n_channels) throughout, matching RawRecording.get_traces.
"""
from __future__ import annotations

import numpy as np


def highpass(traces: np.ndarray, fs: float, cutoff_hz: float = 300.0,
             order: int = 3) -> np.ndarray:
    """Zero-phase Butterworth high-pass along time.

    Zero-phase (filtfilt) because a causal filter shifts a spike trough in time,
    and these traces are read next to spike times taken from the sorter.
    """
    from scipy.signal import butter, sosfiltfilt
    sos = butter(order, cutoff_hz / (fs / 2), btype="highpass", output="sos")
    return sosfiltfilt(sos, np.asarray(traces, dtype=np.float64), axis=0)


def lowpass(traces: np.ndarray, fs: float, cutoff_hz: float = 300.0,
            order: int = 3) -> np.ndarray:
    """Zero-phase Butterworth low-pass along time."""
    from scipy.signal import butter, sosfiltfilt
    sos = butter(order, cutoff_hz / (fs / 2), btype="lowpass", output="sos")
    return sosfiltfilt(sos, np.asarray(traces, dtype=np.float64), axis=0)


def bandpass(traces: np.ndarray, fs: float, low_hz: float, high_hz: float,
             order: int = 3) -> np.ndarray:
    """Zero-phase Butterworth band-pass along time."""
    from scipy.signal import butter, sosfiltfilt
    nyq = fs / 2
    sos = butter(order, [low_hz / nyq, min(high_hz / nyq, 0.99)],
                 btype="bandpass", output="sos")
    return sosfiltfilt(sos, np.asarray(traces, dtype=np.float64), axis=0)


def common_median_reference(traces: np.ndarray, groups=None) -> np.ndarray:
    """Subtract the across-channel median at each sample.

    Median rather than mean so a few large spikes or one broken channel do not
    leak into every other channel. Pass groups (one label per channel) to
    reference within shank or within ADC bank instead of globally.
    """
    traces = np.asarray(traces, dtype=np.float64)
    if groups is None:
        return traces - np.median(traces, axis=1, keepdims=True)
    groups = np.asarray(groups)
    out = traces.copy()
    for g in np.unique(groups):
        sel = groups == g
        out[:, sel] -= np.median(traces[:, sel], axis=1, keepdims=True)
    return out


def preprocess_ap(traces: np.ndarray, fs: float, cutoff_hz: float = 300.0,
                  cmr: bool = True, groups=None) -> np.ndarray:
    """The standard look-at-the-spikes chain: high-pass then common reference."""
    out = highpass(traces, fs, cutoff_hz)
    if cmr:
        out = common_median_reference(out, groups=groups)
    return out


def derive_lfp(traces: np.ndarray, fs: float, cutoff_hz: float = 300.0,
               target_fs: float = 2500.0) -> tuple[np.ndarray, float]:
    """Anti-alias low-pass then decimate, returning (traces, new_fs).

    Only needed when no separate LF band was saved. Neuropixels 2.0 has no LF
    stream from the probe, though a pipeline may have derived and archived one
    already, in which case read that instead of calling this.
    """
    from scipy.signal import decimate
    factor = max(1, int(round(fs / target_fs)))
    filtered = lowpass(traces, fs, min(cutoff_hz, 0.4 * fs / factor))
    if factor == 1:
        return filtered, fs
    return decimate(filtered, factor, axis=0, zero_phase=True), fs / factor


def channel_rms(traces: np.ndarray, subtract_median: bool = True
                ) -> np.ndarray:
    """Per-channel RMS in the units of the input.

    The DC offset is removed first, otherwise this measures the offset rather
    than the noise.
    """
    traces = np.asarray(traces, dtype=np.float64)
    if subtract_median:
        traces = traces - np.median(traces, axis=0, keepdims=True)
    return np.sqrt(np.mean(traces ** 2, axis=0))
