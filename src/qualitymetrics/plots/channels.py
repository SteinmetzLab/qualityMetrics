"""Channel-level and band-level checks on the raw signal.

These answer the questions you should ask before sorting: are all the channels
alive, has the file already been filtered, and where along the probe is there
any activity at all.
"""
from __future__ import annotations

import numpy as np

from ..preproc import channel_rms, highpass
from ..style import DEPTH_LABEL, TIME_LABEL, despine, use_lab_style


def channel_health(rec, n_chunks: int = 8, win_s: float = 1.0,
                   highpass_hz: float | None = 300.0,
                   dead_factor: float = 0.1, noisy_factor: float = 3.0,
                   z_threshold: float = 5.0) -> dict:
    """Per-channel RMS across several windows, with bad channels flagged.

    Sampling several short windows rather than one long one keeps the estimate
    from being a statement about one arbitrary second.

    Two detectors, because one is not enough. The multiplicative one (0.1x and
    3x the array median) catches catastrophe: a disconnected site, a shorted
    one. On healthy Neuropixels 2.0 data it never fires at all, because the RMS
    distribution is tight. On our test shank it spans 17 to 27 uV around a
    median of 23.4, so the "noisy" threshold sits at 70 uV and nothing comes
    near it. A detector that cannot fire is not a check.

    So the RMS is also scored as a robust z, (value - median) / (1.4826 * MAD),
    which adapts to however tight the distribution actually is and flags a
    channel that is anomalous for this recording. Both verdicts are returned
    separately, neither is authoritative on its own, and the caller decides.
    """
    rms_per_chunk = []
    for _t0, block in rec.sample_windows(n_chunks, win_s):
        if highpass_hz:
            block = highpass(block, rec.fs, highpass_hz)
        rms_per_chunk.append(channel_rms(block))
    rms = np.median(np.vstack(rms_per_chunk), axis=0)

    median = float(np.median(rms))
    mad = float(np.median(np.abs(rms - median)))
    # 1.4826 makes the MAD a consistent estimator of sigma for Gaussian data.
    scale = 1.4826 * mad
    z = (rms - median) / scale if scale > 0 else np.zeros_like(rms)

    dead = rms < dead_factor * median
    noisy = rms > noisy_factor * median
    outlier = (np.abs(z) > z_threshold) & ~dead & ~noisy
    geom = rec.geometry
    depths = geom.y[:rec.n_data_channels] if geom is not None \
        else np.arange(rec.n_data_channels, dtype=float)

    return {
        "rms_uv": rms,
        "depth_um": depths,
        "median_rms_uv": median,
        "mad_rms_uv": mad,
        "robust_z": z,
        "dead": dead,
        "noisy": noisy,
        "outlier": outlier,
        "n_dead": int(dead.sum()),
        "n_noisy": int(noisy.sum()),
        "n_outlier": int(outlier.sum()),
        "z_threshold": z_threshold,
        "n_chunks": n_chunks,
        "win_s": win_s,
        "highpass_hz": highpass_hz,
    }


def plot_channel_health(rec, result: dict | None = None, title: str | None = None,
                        figsize=(11, 8), **kwargs):
    """Draw the per-channel RMS against depth, flagging bad sites two ways."""
    import matplotlib.pyplot as plt
    use_lab_style()

    r = result if result is not None else channel_health(rec, **kwargs)
    rms, depth = r["rms_uv"], r["depth_um"]
    flagged = r["dead"] | r["noisy"] | r["outlier"]

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(rms[~flagged], depth[~flagged], "o", ms=3.5, color="0.3",
            label=f"Normal ({int((~flagged).sum())})")
    if r["n_dead"]:
        ax.plot(rms[r["dead"]], depth[r["dead"]], "o", ms=6, color="tab:blue",
                label=f"Dead ({r['n_dead']})")
    if r["n_noisy"]:
        ax.plot(rms[r["noisy"]], depth[r["noisy"]], "o", ms=6, color="tab:red",
                label=f"Noisy ({r['n_noisy']})")
    if r["n_outlier"]:
        ax.plot(rms[r["outlier"]], depth[r["outlier"]], "o", ms=6,
                markerfacecolor="none", markeredgecolor="tab:orange", mew=1.5,
                label=f"Outlier, |z| > {r['z_threshold']:.0f} "
                      f"({r['n_outlier']})")
    ax.axvline(r["median_rms_uv"], color="0.55", ls="--", lw=1,
               label=f"Median {r['median_rms_uv']:.1f} µV")
    for sign in (-1, 1):
        edge = r["median_rms_uv"] + sign * r["z_threshold"] * 1.4826 * r["mad_rms_uv"]
        ax.axvline(edge, color="tab:orange", ls=":", lw=1)

    ax.set_xlabel("RMS noise (µV)")
    ax.set_ylabel(DEPTH_LABEL)
    ax.legend(loc="best")
    despine(ax)
    band = (f"high-passed at {r['highpass_hz']:.0f} Hz" if r["highpass_hz"]
            else "wideband")
    ax.set_title(title or f"{rec.path.name}: channel noise "
                          f"({r['n_chunks']} windows of {r['win_s']:.0f} s, "
                          f"{band})")
    fig.tight_layout()
    return fig


def filter_state(rec, n_chunks: int = 4, win_s: float = 2.0,
                 n_channels: int = 32, cutoff_hz: float = 300.0,
                 title: str | None = None, figsize=(10, 6)):
    """Average power spectrum, to reveal whether the file is already filtered.

    A recording that has been high-passed at 300 Hz has almost no power below
    300 Hz, which means it is ready for spike sorting but cannot yield an LFP by
    low-pass filtering. Knowing which of those two you have is the difference
    between deriving a usable LFP band and deriving a plausible-looking one made
    of filter roll-off.
    """
    import matplotlib.pyplot as plt
    from scipy.signal import welch
    use_lab_style()

    chans = np.unique(np.linspace(0, rec.n_data_channels - 1,
                                  n_channels).astype(int))
    spectra = []
    for _t0, block in rec.sample_windows(n_chunks, win_s, channels=chans):
        freqs, pxx = welch(block, fs=rec.fs, nperseg=min(4096, block.shape[0]),
                           axis=0)
        spectra.append(pxx.mean(axis=1))
    pxx = np.mean(np.vstack(spectra), axis=0)

    low = pxx[(freqs > 1) & (freqs < cutoff_hz)].mean()
    high = pxx[(freqs >= cutoff_hz) & (freqs < min(3000, rec.fs / 2))].mean()
    ratio = float(low / high) if high > 0 else float("inf")
    # A wideband file has far more power below 300 Hz than above it. A filtered
    # one does not. The boundary is soft, so the number is reported either way.
    verdict = ("looks already high-passed" if ratio < 2
               else "looks wideband (LFP present)")

    fig, ax = plt.subplots(figsize=figsize)
    ax.loglog(freqs[1:], pxx[1:], lw=1.2, color="0.2")
    ax.axvline(cutoff_hz, color="crimson", ls="--", lw=1,
               label=f"{cutoff_hz:.0f} Hz")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Power spectral density (µV²/Hz)")
    ax.legend(loc="best")
    despine(ax)
    ax.set_title(title or f"{rec.path.name}: low/high power ratio "
                          f"{ratio:.1f}, {verdict}")
    fig.tight_layout()
    return fig, {"freqs": freqs, "psd": pxx, "ratio": ratio, "verdict": verdict}


def band_rms(rec, lf_rec=None, n_time_bins: int = 60, win_s: float = 0.25,
             depth_bin_um: float = 40.0, highpass_hz: float = 300.0,
             title: str | None = None, figsize=(13, 9)):
    """Time by depth map of RMS voltage, for the AP band and optionally the LFP.

    Bins the session in time and along the probe and shows RMS in each block, so
    a period when a region went quiet, or an artifact that swept the array, is
    visible as a stripe. Sampling short windows keeps this affordable: reading
    every sample of a 17 GB shank to make a 120-column image is not a good
    trade.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    panels = [("AP band", rec, highpass_hz)]
    if lf_rec is not None:
        panels.append(("LFP band", lf_rec, None))

    fig, axes = plt.subplots(len(panels), 1, figsize=figsize, sharex=True,
                             squeeze=False)
    axes = axes[:, 0]

    for ax, (name, recording, hp) in zip(axes, panels, strict=False):
        geom = recording.geometry
        n = recording.n_data_channels
        depths = geom.y[:n] if geom is not None else np.arange(n, dtype=float)
        edges = np.arange(depths.min(), depths.max() + depth_bin_um,
                          depth_bin_um)
        bin_of = np.clip(np.digitize(depths, edges) - 1, 0, len(edges) - 2)

        times, image = [], []
        for t0, block in recording.sample_windows(n_time_bins, win_s):
            if hp:
                block = highpass(block, recording.fs, hp)
            rms = channel_rms(block)
            column = np.array([rms[bin_of == b].mean() if np.any(bin_of == b)
                               else np.nan for b in range(len(edges) - 1)])
            times.append(t0)
            image.append(column)
        image = np.array(image).T                    # (depth bins, time bins)

        finite = image[np.isfinite(image)]
        vmax = np.percentile(finite, 99) if finite.size else 1.0
        im = ax.imshow(image, aspect="auto", origin="lower", cmap="viridis",
                       vmin=0, vmax=vmax,
                       extent=[times[0], times[-1] + win_s,
                               edges[0], edges[-1]])
        cbar = fig.colorbar(im, ax=ax, pad=0.02)
        cbar.set_label("RMS voltage (µV)")
        ax.set_ylabel(DEPTH_LABEL)
        ax.set_title(f"{name} ({len(times)} windows of {win_s * 1000:.0f} ms)")
        despine(ax)
    axes[-1].set_xlabel(TIME_LABEL)

    fig.suptitle(title or f"{rec.path.name}: activity over time and depth",
                 fontsize=12)
    fig.tight_layout()
    return fig


def depth_power(rec, t_start_s: float = 100.0, win_s: float = 10.0,
                fmax_hz: float = 300.0, depth_bin_um: float = 40.0,
                normalize: str | None = None, title: str | None = None,
                figsize=(11, 8)):
    """Depth by frequency power map, the standard read on laminar structure.

    Welch PSD per channel, averaged into depth bins. normalize='freq' divides
    each frequency by its across-depth mean, which removes the 1/f slope and
    makes a band that is stronger at one depth stand out; without it the figure
    is dominated by low frequencies everywhere.
    """
    import matplotlib.pyplot as plt
    from scipy.signal import welch
    use_lab_style()

    geom = rec.geometry
    n = rec.n_data_channels
    depths = geom.y[:n] if geom is not None else np.arange(n, dtype=float)

    block = rec.get_traces(t_start_s, t_start_s + win_s)
    freqs, pxx = welch(block, fs=rec.fs,
                       nperseg=min(int(rec.fs), block.shape[0]), axis=0)
    keep = (freqs > 0) & (freqs <= fmax_hz)
    freqs, pxx = freqs[keep], pxx[keep]              # (n_freq, n_channels)

    edges = np.arange(depths.min(), depths.max() + depth_bin_um, depth_bin_um)
    bin_of = np.clip(np.digitize(depths, edges) - 1, 0, len(edges) - 2)
    image = np.array([pxx[:, bin_of == b].mean(axis=1) if np.any(bin_of == b)
                      else np.full(freqs.size, np.nan)
                      for b in range(len(edges) - 1)])   # (depth, freq)

    label = "Power spectral density (dB re 1 µV²/Hz)"
    with np.errstate(divide="ignore", invalid="ignore"):
        shown = 10 * np.log10(image)
    if normalize == "freq":
        shown = shown - np.nanmean(shown, axis=0, keepdims=True)
        label = "Power relative to depth mean (dB)"

    fig, ax = plt.subplots(figsize=figsize)
    finite = shown[np.isfinite(shown)]
    vmin, vmax = np.percentile(finite, [2, 98]) if finite.size else (0, 1)
    im = ax.imshow(shown, aspect="auto", origin="lower", cmap="viridis",
                   vmin=vmin, vmax=vmax,
                   extent=[freqs[0], freqs[-1], edges[0], edges[-1]])
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(label)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel(DEPTH_LABEL)
    despine(ax)
    ax.set_title(title or f"{rec.path.name}: power against depth "
                          f"({win_s:.0f} s from {t_start_s:.0f} s)")
    fig.tight_layout()
    return fig
