"""Per-unit figures: where the units are, how big, and what shape.

Everything here is derived from the sorter output plus the recording's gain, so
these run on an archived session with no raw data present.
"""
from __future__ import annotations

import numpy as np

from ..style import AMP_LABEL, DEPTH_LABEL, despine, size_legend, use_lab_style


def _spike_count_sizes(counts: np.ndarray, lo: float = 8.0, hi: float = 160.0):
    """Map spike counts to marker areas on a log scale, plus a legend key.

    Returns (sizes, legend_sizes, legend_labels). The legend is the point: the
    MATLAB version scaled markers by spike count and never said so, leaving the
    most visually salient variable in the figure undecodable.
    """
    counts = np.asarray(counts, dtype=float)
    safe = np.clip(counts, 1, None)
    logs = np.log10(safe)
    lo_l, hi_l = logs.min(), logs.max()
    if hi_l <= lo_l:
        sizes = np.full(counts.shape, (lo + hi) / 2)
        return sizes, [(lo + hi) / 2], [f"{int(safe[0]):,}"]
    sizes = lo + (logs - lo_l) / (hi_l - lo_l) * (hi - lo)

    ticks = [10 ** np.ceil(lo_l), 10 ** np.floor((lo_l + hi_l) / 2),
             10 ** np.floor(hi_l)]
    ticks = sorted({int(t) for t in ticks if lo_l <= np.log10(t) <= hi_l})
    key = [lo + (np.log10(t) - lo_l) / (hi_l - lo_l) * (hi - lo) for t in ticks]
    return sizes, key, [f"{t:,}" for t in ticks]


def amp_depth_scatter(ks, duration_s: float | None = None,
                      artifact_z_um: float | None = None,
                      title: str | None = None, figsize=(12, 8)):
    """Unit amplitude against depth, coloured by rate and by waveform width.

    Panel A colours by firing rate, on a log colour axis because rates span
    orders of magnitude. Panel B colours by trough-to-peak duration, the classic
    narrow versus broad spiking split.

    Marker area encodes spike count, and unlike the original there is a legend
    that says so.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    use_lab_style()

    uids = [int(u) for u in ks.unit_ids]
    amp = np.array([ks.unit_amplitude_uv[u] for u in uids])
    depth = np.array([ks.unit_depth_um[u] for u in uids])
    counts = np.array([ks.n_spikes[u] for u in uids], dtype=float)
    rate = np.array([v for v in
                     (ks.firing_rate(duration_s)[u] for u in uids)])
    ttp = np.array([ks.trough_to_peak_ms.get(u, np.nan) for u in uids])

    sizes, key_sizes, key_labels = _spike_count_sizes(counts)

    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)

    # Floor the colour axis at 0.01 spikes/s. Below that a unit fired a handful
    # of times in half an hour, and letting it set vmin spends most of the
    # colour range on units nobody will use while flattening everything else.
    positive = rate[rate > 0]
    norm = LogNorm(vmin=max(positive.min(), 0.01), vmax=rate.max()) \
        if positive.size else None
    sc0 = axes[0].scatter(amp, depth, s=sizes, c=np.clip(rate, 0.01, None),
                          cmap="viridis", norm=norm, alpha=0.8,
                          edgecolors="0.25", linewidths=0.3)
    cb0 = fig.colorbar(sc0, ax=axes[0], pad=0.02)
    cb0.set_label("Firing rate (spikes/s)")
    axes[0].set_title("Colored by firing rate")

    sc1 = axes[1].scatter(amp, depth, s=sizes, c=ttp, cmap="RdYlGn",
                          vmin=0, vmax=np.nanpercentile(ttp, 98) or 1.0,
                          alpha=0.8, edgecolors="0.25", linewidths=0.3)
    cb1 = fig.colorbar(sc1, ax=axes[1], pad=0.02)
    cb1.set_label("Trough to peak (ms)")
    axes[1].set_title("Colored by waveform duration")

    for ax in axes:
        ax.set_xlabel(AMP_LABEL)
        despine(ax)
        if artifact_z_um is not None and np.isfinite(artifact_z_um):
            top = max(depth.max(), artifact_z_um)
            if artifact_z_um < top:
                ax.axhspan(artifact_z_um, top, color="0.88", zorder=0)
                ax.axhline(artifact_z_um, color="0.45", lw=0.8, ls="--")
    axes[0].set_ylabel(DEPTH_LABEL)
    size_legend(axes[0], key_sizes, key_labels, "Spikes", loc="lower right")

    fig.suptitle(title or ks.title(f"Unit amplitude against depth "
                                   f"({len(uids)} units)"), fontsize=12)
    fig.tight_layout()
    return fig


def unit_summary_table(ks, duration_s: float | None = None):
    """The numbers behind amp_depth_scatter, as a list of dicts.

    Useful on its own, and it keeps the plotting function free of anything a
    caller might want without a figure.
    """
    rates = ks.firing_rate(duration_s)
    rows = []
    for u in (int(x) for x in ks.unit_ids):
        rows.append({
            "unit_id": u,
            "n_spikes": ks.n_spikes[u],
            "firing_rate_sps": rates[u],
            "amplitude_uv": ks.unit_amplitude_uv[u],
            "depth_um": ks.unit_depth_um[u],
            "trough_to_peak_ms": ks.trough_to_peak_ms.get(u, float("nan")),
            "ks_label": ks.ks_label.get(u, ""),
        })
    return rows


def templates_grid(ks, unit_ids=None, n_channels: int = 8, n_cols: int = 6,
                   title: str | None = None, scale_uv: float | None = None):
    """A small multiple per unit: the template on its best few channels.

    Channels are chosen around the peak and drawn at their true depth spacing,
    so a unit whose template spreads over 200 um looks different from one
    confined to two rows.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    if unit_ids is None:
        unit_ids = [int(u) for u in ks.unit_ids]
    unit_ids = list(unit_ids)
    if not unit_ids:
        raise ValueError("no units to plot")

    unw = ks.templates_unwhitened
    peak = ks.peak_channel
    pos = ks.channel_positions
    scale = ks.uv_per_bit or 1.0
    n_rows = int(np.ceil(len(unit_ids) / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.0 * n_cols, 2.2 * n_rows),
                             squeeze=False)
    t_ms = np.arange(unw.shape[1]) / ks.fs * 1000.0

    for k, uid in enumerate(unit_ids):
        ax = axes[k // n_cols][k % n_cols]
        if uid >= len(unw):
            ax.axis("off")
            continue
        pk = int(peak[uid])
        # Nearest channels by physical distance, not by index: index order is a
        # wiring detail and jumps across columns on a multi-column probe.
        d = np.hypot(pos[:, 0] - pos[pk, 0], pos[:, 1] - pos[pk, 1])
        chans = np.argsort(d)[:n_channels]
        chans = chans[np.argsort(pos[chans, 1])]

        wave = unw[uid][:, chans] * scale
        step = scale_uv or (np.ptp(wave) * 0.6 or 1.0)
        for j in range(len(chans)):
            ax.plot(t_ms, wave[:, j] + j * step, lw=0.8, color="0.2")
        ax.set_title(f"Unit {uid}", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("left", "bottom"):
            ax.spines[side].set_visible(False)
        despine(ax)

    for k in range(len(unit_ids), n_rows * n_cols):
        axes[k // n_cols][k % n_cols].axis("off")

    # One scale bar for the whole grid beats an unreadable axis on each panel.
    axes[0][0].plot([t_ms[0], t_ms[0] + 1.0], [0, 0], lw=2, color="crimson",
                    clip_on=False)
    axes[0][0].text(t_ms[0], 0, " 1 ms", color="crimson", fontsize=7,
                    va="top")

    fig.suptitle(title or ks.title(f"Templates ({len(unit_ids)} units, "
                                   f"{n_channels} channels each)"), fontsize=12)
    fig.tight_layout()
    return fig


def template_waterfall(ks, unit_id: int, title: str | None = None,
                       figsize=(7, 8)):
    """One unit's template across every channel, as an image over depth.

    Shows how far a unit's footprint spreads, which is the quickest read on
    whether a "unit" is really a shared artifact across the whole shank.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    unw = ks.templates_unwhitened
    if unit_id >= len(unw):
        raise ValueError(f"unit {unit_id} has no template")
    pos = ks.channel_positions
    order = np.argsort(pos[:, 1])
    wave = unw[unit_id][:, order] * (ks.uv_per_bit or 1.0)
    t_ms = np.arange(wave.shape[0]) / ks.fs * 1000.0

    fig, ax = plt.subplots(figsize=figsize)
    lim = np.abs(wave).max() or 1.0
    im = ax.imshow(wave.T, aspect="auto", origin="lower", cmap="RdBu_r",
                   vmin=-lim, vmax=lim,
                   extent=[t_ms[0], t_ms[-1], pos[order, 1].min(),
                           pos[order, 1].max()])
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(AMP_LABEL)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel(DEPTH_LABEL)
    despine(ax)
    amp = ks.unit_amplitude_uv.get(int(unit_id), float("nan"))
    ax.set_title(title or ks.title(
        f"Unit {unit_id} template, {amp:.0f} µV peak to peak"))
    fig.tight_layout()
    return fig
