"""Drift figures: where each spike was, over the whole session.

These are the figures that decide whether a recording needs motion correction
and whether the sorter split one drifting neuron into several units. Both
questions are about time, so both axes here are in seconds, taken from the
sampling rate. The MATLAB original labelled its x-axis "time (s)" and plotted
sample indices, which put 5.6e7 on an axis whose true extent was 1877.
"""
from __future__ import annotations

import colorsys

import numpy as np

from ..style import DEPTH_LABEL, TIME_LABEL, despine, use_lab_style


def _golden_angle_colors(n: int) -> np.ndarray:
    """n maximally separated hues, so neighbouring units look different.

    The golden angle keeps successive hues far apart for any n, which matters
    because the point of the unit-coloured panel is to notice that a continuous
    depth band changes colour partway through the session.
    """
    phi = 0.6180339887498949
    return np.array([colorsys.hsv_to_rgb((i * phi) % 1.0, 0.85, 0.9)
                     for i in range(n)])


def _subsample(n: int, limit: int | None, seed: int = 0):
    """Indices for at most limit points, or all of them."""
    if limit is None or n <= limit:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n, size=limit, replace=False))


def unit_drift(ks, artifact_z_um: float | None = None,
               max_spikes: int | None = 400_000, title: str | None = None,
               figsize=(13, 7)):
    """Per-spike depth against time, coloured two ways.

    Panel A colours by unit. A depth band that changes colour partway through
    the session is one neuron the sorter re-clustered, which is a drift split
    you can see without any motion estimate at all.

    Panel B colours by spike amplitude. Coherent motion shows up as amplitude
    changing together across a band of depths, which panel A cannot show.

    Two colourings of the same scatter, because certifying drift from one
    colouring alone is how people talk themselves into seeing it.

    artifact_z_um shades everything above a depth, for the band above the brain
    surface. It defaults to no shading: the MATLAB version hardcoded 4850 um,
    which on a 2865 um shank shades nothing and on a shorter one shades
    everything, and it also wrote that number into the annotation text so the
    label disagreed with the line whenever a caller passed something else.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    t = ks.spike_times_s
    z = ks.spike_depths_um
    clusters = ks.spike_clusters
    amps = ks.spike_amplitudes_uv()

    idx = _subsample(t.size, max_spikes)
    shown = f"{idx.size:,} of {t.size:,} spikes"

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=True, sharey=True)

    # -- Panel A: colour = unit, ordered by depth so the palette walks the shank
    uids = ks.unit_ids
    order = np.argsort([ks.unit_depth_um[int(u)] for u in uids])
    colors = _golden_angle_colors(len(uids))
    lookup = {int(uids[o]): colors[i] for i, o in enumerate(order)}
    point_colors = np.array([lookup[int(c)] for c in clusters[idx]])

    axes[0].scatter(t[idx], z[idx], s=1.0, c=point_colors, alpha=0.45,
                    linewidths=0, rasterized=True)
    axes[0].set_title(f"Colored by unit ({len(uids)} units)")

    # -- Panel B: colour = amplitude, log-scaled because the range is wide
    a = amps[idx]
    vmin, vmax = np.percentile(a, [2, 98])
    sc = axes[1].scatter(t[idx], z[idx], s=1.0, c=a, cmap="viridis",
                         vmin=vmin, vmax=vmax, alpha=0.45, linewidths=0,
                         rasterized=True)
    axes[1].set_title("Colored by spike amplitude")
    cbar = fig.colorbar(sc, ax=axes[1], pad=0.02)
    cbar.set_label("Spike amplitude (µV)")
    cbar.solids.set_alpha(1.0)

    for ax in axes:
        ax.set_xlabel(TIME_LABEL)
        despine(ax)
        if artifact_z_um is not None and np.isfinite(artifact_z_um):
            top = max(np.nanmax(z), artifact_z_um)
            if artifact_z_um < top:
                ax.axhspan(artifact_z_um, top, color="0.85", zorder=0)
                ax.axhline(artifact_z_um, color="0.45", lw=0.8, ls="--")
                # Annotation text derived from the argument, never hardcoded.
                ax.text(0.99, artifact_z_um, f" Above {artifact_z_um:.0f} µm",
                        transform=ax.get_yaxis_transform(), ha="right",
                        va="bottom", fontsize=8, color="0.35")
    axes[0].set_ylabel(DEPTH_LABEL)
    axes[0].set_xlim(0, t.max())

    fig.suptitle(title or ks.title(f"Spike depth over time ({shown})"),
                 fontsize=12)
    fig.tight_layout()
    return fig


def drift_map(ks, amp_threshold_uv: float = 0.0, max_spikes: int | None = 600_000,
              title: str | None = None, figsize=(13, 7)):
    """The classic drift map: every spike as a dot, darkness by amplitude.

    Equivalent in intent to plotDriftmap in cortex-lab/spikes, with the
    amplitude scale in real microvolts rather than in the sorter's arbitrary
    units, and with a colorbar that says what the shading means.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    t = ks.spike_times_s
    z = ks.spike_depths_um
    amps = ks.spike_amplitudes_uv()

    keep = amps >= amp_threshold_uv
    t, z, amps = t[keep], z[keep], amps[keep]
    idx = _subsample(t.size, max_spikes)

    fig, ax = plt.subplots(figsize=figsize)
    vmin, vmax = np.percentile(amps, [5, 99])
    sc = ax.scatter(t[idx], z[idx], s=1.5, c=amps[idx], cmap="gray_r",
                    vmin=vmin, vmax=vmax, alpha=0.35, linewidths=0,
                    rasterized=True)
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Spike amplitude (µV)")
    cbar.solids.set_alpha(1.0)

    ax.set_xlabel(TIME_LABEL)
    ax.set_ylabel(DEPTH_LABEL)
    ax.set_xlim(0, ks.spike_times_s.max())
    despine(ax)
    thresh = (f", amplitude ≥ {amp_threshold_uv:.0f} µV"
              if amp_threshold_uv > 0 else "")
    ax.set_title(title or ks.title(
        f"Drift map ({idx.size:,} spikes{thresh})"))
    fig.tight_layout()
    return fig


def amplitude_cdf_over_depth(ks, depth_bin_um: float = 40.0,
                             amp_bin_uv: float = 10.0,
                             max_amp_uv: float | None = None,
                             title: str | None = None, figsize=(12, 6)):
    """Spike-amplitude distribution as a function of depth.

    Left: the rate of spikes at each amplitude and depth, in spikes/s, so the
    colour is a rate rather than a raw count and is comparable between
    recordings of different length.

    Right: the same thing accumulated from the top down, which reads as "how
    many spikes per second above this amplitude live at this depth" and is the
    view that answers whether a stretch of the probe is producing anything worth
    sorting.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    use_lab_style()

    z = ks.spike_depths_um
    amps = ks.spike_amplitudes_uv()
    duration = ks.duration_s

    if max_amp_uv is None:
        max_amp_uv = float(np.percentile(amps, 99.5))
    depth_edges = np.arange(np.floor(z.min() / depth_bin_um) * depth_bin_um,
                            z.max() + depth_bin_um, depth_bin_um)
    amp_edges = np.arange(0, max_amp_uv + amp_bin_uv, amp_bin_uv)

    counts, _, _ = np.histogram2d(z, amps, bins=[depth_edges, amp_edges])
    rate = counts / duration                        # spikes/s per bin
    # Accumulate from the largest amplitude downwards: "at least this big".
    survival = np.cumsum(rate[:, ::-1], axis=1)[:, ::-1]

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    extent = [amp_edges[0], amp_edges[-1], depth_edges[0], depth_edges[-1]]

    # Empty bins get their own colour, off the end of the scale. At high
    # amplitudes a handful of spikes is all there is, and that handful is the
    # interesting part; on a continuous scale anchored at zero it is
    # indistinguishable from a bin where nothing happened. Masking the true
    # zeros to grey makes "rare" and "never" different things to look at.
    empty_color = "0.75"
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(empty_color)

    # Log colour scale for the same reason: one spike and a thousand differ by
    # three orders of magnitude, and a linear scale spends the whole palette on
    # the densest bins.
    for ax, data, name in ((axes[0], rate, "Rate in bin (spikes/s)"),
                           (axes[1], survival,
                            "Rate above amplitude (spikes/s)")):
        shown = np.ma.masked_where(data <= 0, data)
        positive = data[data > 0]
        if positive.size:
            norm = LogNorm(vmin=positive.min(),
                           vmax=max(np.percentile(positive, 99.5),
                                    positive.min() * 10))
        else:
            norm = None
        im = ax.imshow(shown, aspect="auto", origin="lower", extent=extent,
                       cmap=cmap, norm=norm, interpolation="nearest")
        ax.set_xlabel("Spike amplitude (µV)")
        despine(ax)
        cbar = fig.colorbar(im, ax=ax, pad=0.02, extend="min")
        cbar.set_label(name)
        # Name the grey rather than leaving the reader to infer it.
        cbar.ax.plot([0.5], [-0.04], marker="s", markersize=9,
                     color=empty_color, markeredgecolor="0.4",
                     transform=cbar.ax.transAxes, clip_on=False)
        cbar.ax.text(2.2, -0.04, "No spikes", transform=cbar.ax.transAxes,
                     va="center", ha="left", fontsize=8)
    axes[0].set_ylabel(DEPTH_LABEL)
    axes[0].set_title("Amplitude distribution")
    axes[1].set_title("Cumulative from high amplitude")

    fig.suptitle(title or ks.title("Spike amplitude against depth"), fontsize=12)
    fig.tight_layout()
    return fig
