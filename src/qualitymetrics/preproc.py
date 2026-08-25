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


def destripe_available() -> bool:
    """Is IBL's destriping importable here?"""
    try:
        import ibldsp.voltage  # noqa: F401
        import neuropixel  # noqa: F401
    except Exception:  # noqa: BLE001 - optional
        return False
    return True


def destripe(traces: np.ndarray, fs: float, geometry=None,
             neuropixel_version: int = 2, k_filter: bool = True) -> np.ndarray:
    """IBL's destriping: ADC phase correction, high-pass, then a spatial filter.

    A stronger clean-up than high-pass plus common median reference, and aimed
    at a different artifact. Common median subtraction removes whatever is
    identical across channels at one instant. Destriping first undoes the ADC
    sampling skew, so channels digitised at different moments line up before
    anything is subtracted, and then applies a spatial filter across the array.
    Line noise and its harmonics arrive with exactly that skew, which is why
    they survive a plain common reference and show up as regular horizontal
    banding down the probe.

    Takes and returns (n_samples, n_channels) in microvolts. IBL works in volts
    with channels first, so the conversion happens here rather than at every
    call site.

    Pass geometry to use the recording's real site positions. Without it, the
    standard header for the probe version is used, whose depth origin is offset
    by 20 um from ours; that offset does not change the spatial filter, but
    relying on it would be an assumption nobody would remember making.
    """
    import ibldsp.voltage as voltage
    from neuropixel import trace_header

    x = np.asarray(traces, dtype=np.float64).T * 1e-6      # volts, channels first
    header = trace_header(version=neuropixel_version, nshank=1)
    if geometry is not None and len(geometry.x) >= x.shape[0]:
        header = dict(header)
        header["x"] = np.asarray(geometry.x[:x.shape[0]], dtype=float)
        header["y"] = np.asarray(geometry.y[:x.shape[0]], dtype=float)
    cleaned = voltage.destripe(x, fs=fs, h=header, k_filter=k_filter)
    return (np.asarray(cleaned).T * 1e6).astype(np.float32)
