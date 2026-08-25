"""Population views along the probe, borrowed from IBL's alignment GUI.

The IBL ephys alignment GUI exists to place a probe against histology, and the
plots it offers are chosen to make anatomical boundaries visible from the
electrophysiology alone. Three of them do work that nothing else in this package
was doing, so they are reproduced here in our conventions:

* a depth-by-depth correlation matrix, which shows regions as blocks along the
  diagonal;
* firing rate as a binned image over time and depth, which the per-spike scatter
  cannot show once the scatter saturates;
* firing rate and amplitude profiles down the probe.

See iblapps/atlaselectrophysiology/plot_data.py for the originals.
"""
from __future__ import annotations

import numpy as np

from ..style import DEPTH_LABEL, RATE_LABEL, TIME_LABEL, despine, use_lab_style


def _binned_rates(ks, depth_bin_um: float, time_bin_s: float,
                  duration_s: float | None = None):
    """Spike counts as (n_depth_bins, n_time_bins), plus the bin edges."""
    z = ks.spike_depths_um
    t = ks.spike_times_s
    duration = duration_s if duration_s is not None else float(t.max())
    depth_edges = np.arange(np.floor(z.min() / depth_bin_um) * depth_bin_um,
                            z.max() + depth_bin_um, depth_bin_um)
    time_edges = np.arange(0, duration + time_bin_s, time_bin_s)
    counts, _, _ = np.histogram2d(z, t, bins=[depth_edges, time_edges])
    return counts, depth_edges, time_edges


def depth_correlation(ks, depth_bin_um: float = 40.0, time_bin_s: float = 0.05,
                      duration_s: float | None = None,
                      min_spikes: int = 50, title: str | None = None,
                      figsize=(9, 8)):
    """Correlation of binned firing rate between every pair of depths.

    Neurons in the same structure share fluctuations, so a brain region shows up
    as a block of high correlation along the diagonal and a boundary shows up as
    the edge of that block. This is the closest thing in the package to a
    functional anatomy read, and it needs no histology.

    Depth bins with almost no spikes are excluded rather than correlated: the
    correlation of two nearly empty bins is noise, and leaving them in scatters
    bright pixels across the figure that look like structure.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    counts, depth_edges, _time_edges = _binned_rates(
        ks, depth_bin_um, time_bin_s, duration_s)
    totals = counts.sum(axis=1)
    keep = totals >= min_spikes
    if keep.sum() < 3:
        raise ValueError("too few occupied depth bins to correlate")

    kept = counts[keep]
    corr = np.corrcoef(kept)
    centres = (depth_edges[:-1] + depth_edges[1:]) / 2
    kept_centres = centres[keep]

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(corr, aspect="auto", origin="lower", cmap="viridis",
                   extent=[kept_centres.min(), kept_centres.max(),
                           kept_centres.min(), kept_centres.max()],
                   vmin=np.percentile(corr, 2), vmax=np.percentile(corr, 98),
                   interpolation="nearest")
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("Correlation of firing rate (r)")
    ax.set_xlabel(DEPTH_LABEL)
    ax.set_ylabel(DEPTH_LABEL)
    despine(ax)
    ax.set_title(title or ks.title(
        f"Firing rate correlation between depths "
        f"({depth_bin_um:.0f} µm bins, {time_bin_s * 1000:.0f} ms windows, "
        f"{int(keep.sum())} of {keep.size} bins occupied)"))
    fig.tight_layout()
    return fig


def firing_rate_image(ks, depth_bin_um: float = 20.0, time_bin_s: float = 1.0,
                      duration_s: float | None = None,
                      title: str | None = None, figsize=(13, 7)):
    """Firing rate over time and depth, as a binned image.

    The per-spike scatter in unit_drift saturates wherever the recording is
    dense, so a change in rate inside a busy band is invisible there. Binning
    into a rate image gives that back, at the cost of no longer showing
    individual units.

    Empty bins are grey rather than the bottom of the colormap, so a stretch of
    shank with nothing in it is distinguishable from one that is merely quiet.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    use_lab_style()

    counts, depth_edges, time_edges = _binned_rates(
        ks, depth_bin_um, time_bin_s, duration_s)
    rate = counts / time_bin_s

    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("0.75")
    shown = np.ma.masked_where(rate <= 0, rate)
    positive = rate[rate > 0]
    norm = (LogNorm(vmin=max(positive.min(), 1e-3),
                    vmax=np.percentile(positive, 99.5))
            if positive.size else None)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(shown, aspect="auto", origin="lower", cmap=cmap, norm=norm,
                   extent=[time_edges[0], time_edges[-1],
                           depth_edges[0], depth_edges[-1]],
                   interpolation="nearest")
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(RATE_LABEL)
    ax.set_xlabel(TIME_LABEL)
    ax.set_ylabel(DEPTH_LABEL)
    despine(ax)
    ax.set_title(title or ks.title(
        f"Firing rate over time and depth "
        f"({depth_bin_um:.0f} µm by {time_bin_s:.0f} s bins)"))
    fig.tight_layout()
    return fig


def depth_profiles(ks, depth_bin_um: float = 20.0,
                   duration_s: float | None = None, min_spikes: int = 50,
                   title: str | None = None, figsize=(11, 8)):
    """Firing rate, spike amplitude and unit count down the probe.

    The marginal summaries of the drift map. Useful precisely because they are
    one-dimensional: a boundary that is hard to see in a dense scatter is often
    obvious as a step in a profile.

    Amplitude is left blank in bins with fewer than min_spikes, following IBL,
    because a mean amplitude from three spikes is not a measurement.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    z = ks.spike_depths_um
    amps = ks.spike_amplitudes_uv()
    duration = duration_s if duration_s is not None else ks.duration_s
    edges = np.arange(np.floor(z.min() / depth_bin_um) * depth_bin_um,
                      z.max() + depth_bin_um, depth_bin_um)
    centres = (edges[:-1] + edges[1:]) / 2
    which = np.clip(np.digitize(z, edges) - 1, 0, len(centres) - 1)

    counts = np.bincount(which, minlength=len(centres)).astype(float)
    rate = counts / duration
    amp_sum = np.bincount(which, weights=amps, minlength=len(centres))
    with np.errstate(invalid="ignore", divide="ignore"):
        amp_mean = amp_sum / counts
    amp_mean[counts < min_spikes] = np.nan

    unit_depth = ks.unit_depth_um
    unit_which = np.clip(
        np.digitize(np.array(list(unit_depth.values())), edges) - 1,
        0, len(centres) - 1)
    n_units = np.bincount(unit_which, minlength=len(centres)).astype(float)

    fig, axes = plt.subplots(1, 3, figsize=figsize, sharey=True)
    axes[0].plot(rate, centres, lw=1.2, color="0.2")
    axes[0].set_xlabel(RATE_LABEL)
    axes[0].set_title("Multi-unit rate")

    axes[1].plot(amp_mean, centres, lw=1.2, color="tab:blue")
    axes[1].set_xlabel("Mean spike amplitude (µV)")
    axes[1].set_title(f"Amplitude (bins with ≥ {min_spikes} spikes)")

    axes[2].barh(centres, n_units, height=depth_bin_um * 0.9, color="0.35")
    axes[2].set_xlabel("Units (count)")
    axes[2].set_title("Sorted units")

    for ax in axes:
        despine(ax)
    axes[0].set_ylabel(DEPTH_LABEL)
    fig.suptitle(title or ks.title(f"Profiles down the probe "
                                   f"({depth_bin_um:.0f} µm bins)"),
                 fontsize=12)
    fig.tight_layout()
    return fig


#: The bands IBL's alignment GUI splits LFP power into. Kept identical so a
#: figure from here and one from that GUI can be compared without arithmetic.
LFP_BANDS = ((0.0, 4.0, "Delta, 0 to 4 Hz"),
             (4.0, 10.0, "Theta, 4 to 10 Hz"),
             (10.0, 30.0, "Alpha/beta, 10 to 30 Hz"),
             (30.0, 80.0, "Gamma, 30 to 80 Hz"),
             (80.0, 200.0, "High, 80 to 200 Hz"))


def lfp_band_profiles(rec, t_start_s: float = 100.0, win_s: float = 10.0,
                      depth_bin_um: float = 40.0, title: str | None = None,
                      figsize=(10, 8)):
    """LFP power in five named bands, as profiles down the probe.

    The depth-by-frequency image shows everything at once and is therefore hard
    to read a boundary off. Collapsing to the conventional bands gives five
    curves whose relative shape is the anatomical signal: hippocampus, for
    instance, is a theta peak that the full spectrogram makes you hunt for.

    Each band is normalised to its own maximum, since absolute power falls
    steeply with frequency and plotting the raw values puts every band except
    delta on the floor.
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

    edges = np.arange(depths.min(), depths.max() + depth_bin_um, depth_bin_um)
    centres = (edges[:-1] + edges[1:]) / 2
    bin_of = np.clip(np.digitize(depths, edges) - 1, 0, len(centres) - 1)

    fig, ax = plt.subplots(figsize=figsize)
    for lo, hi, label in LFP_BANDS:
        if lo >= rec.fs / 2:
            continue
        sel = (freqs >= lo) & (freqs < min(hi, rec.fs / 2))
        if not sel.any():
            continue
        per_channel = pxx[sel].mean(axis=0)
        profile = np.array([per_channel[bin_of == b].mean()
                            if np.any(bin_of == b) else np.nan
                            for b in range(len(centres))])
        peak = np.nanmax(profile)
        if peak and np.isfinite(peak):
            ax.plot(profile / peak, centres, lw=1.4, label=label)

    ax.set_xlabel("Power, normalized to each band's maximum")
    ax.set_ylabel(DEPTH_LABEL)
    ax.legend(loc="best", fontsize=8)
    despine(ax)
    ax.set_title(title or f"{rec.path.name}: LFP band power against depth "
                          f"({win_s:.0f} s from {t_start_s:.0f} s)")
    fig.tight_layout()
    return fig
