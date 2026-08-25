"""Figures that read the archived recording.

These are the ones the MATLAB version could not run once a session had been
archived, because it looked for a SpikeInterface scratch file the pipeline
deletes. They read the .cbin instead, which is the copy that is kept.

The traces are wideband as stored, so every viewer here high-passes before
drawing and says so in the title. Drawing NP2.0 traces unfiltered shows the
per-channel DC offset, which on our data is 1100 uV against spikes of 70.
"""
from __future__ import annotations

import numpy as np

from ..preproc import preprocess_ap
from ..style import DEPTH_LABEL, despine, use_lab_style


def _depth_order(rec):
    """Channel indices sorted shallow to deep, and their depths."""
    geom = rec.geometry
    n = rec.n_data_channels
    if geom is None:
        return np.arange(n), np.arange(n, dtype=float)
    y = geom.y[:n]
    order = np.argsort(y)
    return order, y[order]


def _filtered_window(rec, t_start_s, win_ms, highpass_hz, cmr, channels):
    """Filter over the whole shank, then return only the requested channels.

    The common median reference must be taken across the array, not across
    whatever subset is being drawn. Referencing a contiguous 24-channel band
    against its own median subtracts the local population signal, which removes
    exactly the spikes the figure exists to show, and does it invisibly: the
    result still looks like plausible traces.

    Costs nothing to do properly. mtscomp decompresses whole samples across all
    channels regardless, so the wide read is already paid for and only the
    column slice differs.
    """
    t1 = t_start_s + win_ms / 1000.0
    full = rec.get_traces(t_start_s, t1)             # every data channel
    if highpass_hz:
        full = preprocess_ap(full, rec.fs, highpass_hz, cmr=cmr)
    return full[:, np.asarray(channels, dtype=int)]


def wall_heatmap(rec, t_starts_s=None, win_ms: float = 25.0,
                 highpass_hz: float | None = 300.0, cmr: bool = True,
                 clim_uv: float | None = None, clim_pct: float = 99.5,
                 artifact_z_um: float | None = None, title: str | None = None,
                 figsize=(15, 8)):
    """Voltage across the whole shank at three moments in the session.

    Depth on y, time on x, voltage as a diverging colour centred on zero, so a
    spike is a dark spot and a bad channel is a stripe. This is the fastest way
    to see dead channels, a saturating amplifier, or line noise, and it is the
    figure to look at before trusting anything downstream.

    Three windows by default, early, middle and late, because one window says
    nothing about whether a problem was there the whole time. A dead channel is
    dead in all three; an artifact that appears only in the last panel is a
    different problem with a different cause.

    All panels share one colour scale, computed from the three together. Letting
    each autoscale would make them individually prettier and mutually
    incomparable, which defeats the point of showing three.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    order, depths = _depth_order(rec)
    if t_starts_s is None:
        # 10%, 50% and 90% through, so the first window clears any settling
        # transient and the last is not truncated by the end of the file.
        span = max(rec.duration_s - win_ms / 1000.0, 0.0)
        t_starts_s = [span * f for f in (0.10, 0.50, 0.90)]
    elif np.isscalar(t_starts_s):
        t_starts_s = [float(t_starts_s)]
    t_starts_s = [float(t) for t in t_starts_s]

    panels = []
    for t0 in t_starts_s:
        traces = rec.get_traces(t0, t0 + win_ms / 1000.0, channels=order)
        if highpass_hz:
            traces = preprocess_ap(traces, rec.fs, highpass_hz, cmr=cmr)
        panels.append((t0, traces))

    if clim_uv is None:
        pooled = np.concatenate([np.abs(p).ravel() for _t, p in panels])
        clim_uv = float(np.percentile(pooled, clim_pct)) or 1.0

    fig, axes = plt.subplots(1, len(panels), figsize=figsize, sharey=True,
                             squeeze=False)
    axes = axes[0]
    for ax, (t0, traces) in zip(axes, panels, strict=True):
        im = ax.imshow(traces.T, aspect="auto", origin="lower", cmap="RdBu_r",
                       vmin=-clim_uv, vmax=clim_uv,
                       extent=[0, win_ms, depths.min(), depths.max()],
                       interpolation="nearest")
        ax.set_xlabel("Time from window start (ms)")
        despine(ax)
        ax.set_title(f"{t0:.0f} s ({t0 / rec.duration_s * 100:.0f}% in)")
        if artifact_z_um is not None and np.isfinite(artifact_z_um) \
                and artifact_z_um < depths.max():
            ax.axhline(artifact_z_um, color="0.2", lw=1.0, ls="--")
            ax.text(win_ms, artifact_z_um, f" Above {artifact_z_um:.0f} µm ",
                    ha="right", va="bottom", fontsize=8, color="0.2")
    axes[0].set_ylabel(DEPTH_LABEL)

    cbar = fig.colorbar(im, ax=axes, pad=0.02, fraction=0.03)
    cbar.set_label("Voltage (µV)")

    band = f"high-passed at {highpass_hz:.0f} Hz" if highpass_hz else "wideband"
    ref = ", common median referenced" if (highpass_hz and cmr) else ""
    fig.suptitle(title or f"{rec.path.name}: {win_ms:.0f} ms windows "
                          f"({band}{ref}, shared color scale)", fontsize=12)
    return fig


def raw_traces(rec, t_start_s: float = 100.0, win_ms: float = 100.0,
               channels=None, n_channels: int = 40,
               highpass_hz: float | None = 300.0, cmr: bool = True,
               spacing_uv: float | None = None, title: str | None = None,
               figsize=(13, 9)):
    """Stacked traces for a slice of the shank, with a real microvolt scale bar.

    Channels are evenly spread over the shank by default rather than taken from
    one end, so the figure describes the whole probe.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    order, depths = _depth_order(rec)
    if channels is None:
        pick = np.linspace(0, len(order) - 1, min(n_channels, len(order)))
        pick = np.unique(pick.astype(int))
        channels = order[pick]
        chan_depths = depths[pick]
    else:
        channels = np.asarray(channels, dtype=int)
        geom = rec.geometry
        chan_depths = (geom.y[channels] if geom is not None
                       else channels.astype(float))

    traces = _filtered_window(rec, t_start_s, win_ms, highpass_hz, cmr,
                              channels)

    if spacing_uv is None:
        spacing_uv = float(np.percentile(np.abs(traces), 99.5)) * 2.5 or 50.0
    t_ms = np.arange(traces.shape[0]) / rec.fs * 1000.0

    fig, ax = plt.subplots(figsize=figsize)
    for j in range(traces.shape[1]):
        ax.plot(t_ms, traces[:, j] + j * spacing_uv, lw=0.5, color="0.15")

    ax.set_yticks(np.arange(traces.shape[1]) * spacing_uv)
    ax.set_yticklabels([f"{d:.0f}" for d in chan_depths], fontsize=7)
    ax.set_ylabel(DEPTH_LABEL)
    ax.set_xlabel("Time from window start (ms)")
    ax.set_xlim(0, win_ms)
    despine(ax)

    # Scale bar, because a stacked-trace plot has no usable voltage axis.
    bar = _nice_round(spacing_uv / 2.5)
    ax.plot([win_ms * 0.99] * 2, [0, bar], lw=3, color="crimson",
            solid_capstyle="butt", clip_on=False)
    ax.text(win_ms * 0.985, bar / 2, f"{bar:.0f} µV ", ha="right", va="center",
            color="crimson", fontsize=9)

    band = f"high-passed at {highpass_hz:.0f} Hz" if highpass_hz else "wideband"
    ref = ", common median referenced" if (highpass_hz and cmr) else ""
    ax.set_title(title or f"{rec.path.name}: {traces.shape[1]} channels, "
                          f"{win_ms:.0f} ms from {t_start_s:.0f} s "
                          f"({band}{ref})")
    fig.tight_layout()
    return fig


def raw_with_spikes(rec, ks, t_start_s: float = 100.0, win_ms: float = 50.0,
                    n_channels: int = 24, depth_range_um=None,
                    highpass_hz: float | None = 300.0, cmr: bool = True,
                    spacing_uv: float | None = None, title: str | None = None,
                    figsize=(13, 9)):
    """Raw traces with the sorter's spikes marked where it says they are.

    This is the one figure that checks the sorter against the data rather than
    against itself: if the marks do not land on visible deflections, no
    downstream metric means anything. Each unit gets its own colour, and marks
    are drawn at the unit's own depth.

    Channels are a *contiguous* block centred on the busiest part of the shank
    unless depth_range_um says otherwise. Spreading 40 channels evenly over 384
    (the obvious thing, and what the sibling raw_traces does) breaks this figure
    specifically: a spike lands on four or five adjacent sites, so a subsampled
    set shows the mark with no deflection under it, which looks exactly like a
    sorter error and is not one.
    """
    import matplotlib.pyplot as plt

    from .drift import _golden_angle_colors
    use_lab_style()

    order, depths = _depth_order(rec)
    if depth_range_um is not None:
        lo, hi = depth_range_um
        keep = (depths >= lo) & (depths <= hi)
        order, depths = order[keep], depths[keep]
        pick = np.unique(np.linspace(0, len(order) - 1,
                                     min(n_channels, len(order))).astype(int))
    else:
        # Centre on the depth where this sorting actually found spikes, so the
        # window is not spent on a silent stretch of shank.
        busiest = float(np.median(ks.spike_depths_um))
        centre = int(np.argmin(np.abs(depths - busiest)))
        n_take = min(n_channels, len(order))
        start = int(np.clip(centre - n_take // 2, 0, len(order) - n_take))
        pick = np.arange(start, start + n_take)
    channels, chan_depths = order[pick], depths[pick]

    t1 = t_start_s + win_ms / 1000.0
    traces = _filtered_window(rec, t_start_s, win_ms, highpass_hz, cmr,
                              channels)
    if spacing_uv is None:
        # Tight enough that a spike is a visible deflection rather than a
        # wiggle. Traces overlap a little; that is the right trade here.
        spacing_uv = float(np.percentile(np.abs(traces), 99.0)) * 2.2 or 50.0
    t_ms = np.arange(traces.shape[0]) / rec.fs * 1000.0

    fig, ax = plt.subplots(figsize=figsize)
    for j in range(traces.shape[1]):
        ax.plot(t_ms, traces[:, j] + j * spacing_uv, lw=0.5, color="0.55")

    st = ks.spike_times_s
    in_window = (st >= t_start_s) & (st < t1)
    clusters = ks.spike_clusters[in_window]
    times_ms = (st[in_window] - t_start_s) * 1000.0
    unit_depth = ks.unit_depth_um

    uids = np.unique(clusters)
    colors = _golden_angle_colors(max(len(uids), 1))
    lo_d, hi_d = chan_depths.min(), chan_depths.max()
    n_marked = 0
    for i, uid in enumerate(uids):
        d = unit_depth.get(int(uid))
        if d is None or not (lo_d <= d <= hi_d):
            continue
        row = float(np.interp(d, chan_depths,
                              np.arange(len(chan_depths))))
        sel = clusters == uid
        # Filled circles rather than tick marks. A vertical tick sits on top of
        # a trace that is itself a thin vertical squiggle, so it reads as part
        # of the signal; a filled disc with a pale edge separates from it at a
        # glance, which is the whole job of this figure.
        ax.plot(times_ms[sel], np.full(sel.sum(), row * spacing_uv),
                marker="o", linestyle="none", ms=7, mew=0.8,
                markerfacecolor=colors[i], markeredgecolor="white",
                alpha=0.95, label=f"Unit {uid}", zorder=5)
        n_marked += int(sel.sum())

    ax.set_yticks(np.arange(len(chan_depths)) * spacing_uv)
    ax.set_yticklabels([f"{d:.0f}" for d in chan_depths], fontsize=7)
    ax.set_ylabel(DEPTH_LABEL)
    ax.set_xlabel("Time from window start (ms)")
    ax.set_xlim(0, win_ms)
    despine(ax)
    ax.set_title(title or f"{rec.path.name}: sorted spikes on the raw signal "
                          f"({n_marked} spikes from {len(uids)} units, "
                          f"{win_ms:.0f} ms from {t_start_s:.0f} s)")
    fig.tight_layout()
    return fig


def spike_snippets(rec, ks, unit_id: int, n_waveforms: int = 100,
                   win_ms: float = 3.0, n_channels: int = 6,
                   highpass_hz: float | None = 300.0, seed: int = 0):
    """Pull individual spike waveforms for one unit out of the raw file.

    Returns (t_ms, stack, chans), where stack is
    (n_snippets, n_samples, n_channels) in microvolts and chans is ordered
    shallow to deep.

    Separated from the plotting so the standalone figure and the example-neuron
    grid extract snippets exactly the same way.
    """
    st = ks.spike_times_s[ks.spike_clusters == unit_id]
    if st.size == 0:
        raise ValueError(f"unit {unit_id} has no spikes")
    rng = np.random.default_rng(seed)
    take = st if st.size <= n_waveforms else rng.choice(st, n_waveforms, False)
    take = np.sort(take)

    pos = ks.channel_positions
    pk = int(ks.peak_channel[unit_id]) if unit_id < len(ks.peak_channel) else 0
    d = np.hypot(pos[:, 0] - pos[pk, 0], pos[:, 1] - pos[pk, 1])
    chans = np.argsort(d)[:n_channels]
    chans = chans[np.argsort(pos[chans, 1])]

    half = win_ms / 2000.0
    n_samp = int(round(win_ms / 1000.0 * rec.fs))
    snippets = []
    for t in take:
        if t - half < 0 or t + half > rec.duration_s:
            continue
        block = rec.get_traces(t - half, t + half, channels=chans)
        if highpass_hz:
            # No common reference on a snippet: the six channels here are all
            # neighbours of the same spike, so referencing them against each
            # other would subtract the waveform being looked at.
            block = preprocess_ap(block, rec.fs, highpass_hz, cmr=False)
        if block.shape[0] >= n_samp:
            snippets.append(block[:n_samp])
    if not snippets:
        raise ValueError(f"no usable snippets for unit {unit_id}")
    stack = np.stack(snippets)
    t_ms = (np.arange(n_samp) / rec.fs - half) * 1000.0
    return t_ms, stack, chans


def _site_spacing(pos):
    """Column separation and row pitch of a set of sites, in micrometres."""
    xs = np.unique(np.round(pos[:, 0], 3))
    ys = np.unique(np.round(pos[:, 1], 3))
    dx = float(np.min(np.diff(xs))) if xs.size > 1 else float("nan")
    dy = float(np.median(np.diff(ys))) if ys.size > 1 else float("nan")
    if not np.isfinite(dy) or dy <= 0:
        dy = 15.0                       # Neuropixels 2.0 row pitch
    if not np.isfinite(dx) or dx <= 0:
        dx = 2 * dy                     # single-column probe: invent a width
    return dx, dy


def draw_waveforms_on_probe(ax, ks, t_ms, stack, chans, template=None,
                            show_depth_labels=True, x_span_frac: float = 0.9,
                            y_gain_rows: float = 1.5):
    """Draw snippets at the sites they were recorded on, the way phy does.

    Each channel's waveform is drawn centred on that channel's real position on
    the shank, so a two-column probe produces two columns of waveforms and the
    spatial spread of a unit is visible directly. A single vertical stack, which
    is what this used to do, throws away the across-shank coordinate and makes a
    unit centred on the left column look identical to one centred on the right.

    Pass template as (t_ms, wave) in microvolts to overlay it. That overlay is
    the check that matters: the template and the mean of the raw snippets are
    computed from different things by different code, on the same channels and
    the same time origin, so if they disagree something is genuinely wrong
    rather than merely conventional.
    """
    pos = np.asarray(ks.channel_positions)[np.asarray(chans, dtype=int)]
    dx, dy = _site_spacing(pos)

    span_ms = float(t_ms[-1] - t_ms[0]) or 1.0
    x_scale = (x_span_frac * dx) / span_ms            # micrometres per ms
    peak_uv = float(np.percentile(np.abs(stack), 99.8)) or 1.0
    if template is not None:
        peak_uv = max(peak_uv, float(np.abs(template[1]).max()) or 1.0)
    y_scale = (y_gain_rows * dy) / peak_uv            # micrometres per microvolt

    for j in range(len(chans)):
        px, py = pos[j, 0], pos[j, 1]
        x = px + t_ms * x_scale
        ax.plot(x, py + stack[:, :, j].T * y_scale, lw=0.3, color="0.7",
                alpha=0.10, zorder=1)
        ax.plot(x, py + stack[:, :, j].mean(0) * y_scale, lw=1.5,
                color="crimson", zorder=3)
        if template is not None:
            tt, tw = template
            ax.plot(px + tt * x_scale, py + tw[:, j] * y_scale, lw=1.3,
                    color="tab:blue", alpha=0.9, zorder=2)

    ax.set_xticks(np.unique(np.round(pos[:, 0], 3)))
    ax.set_xlabel("Position across shank (µm)")
    ax.set_yticks(np.unique(np.round(pos[:, 1], 3)))
    if not show_depth_labels:
        ax.set_yticklabels([])
    ax.tick_params(labelsize=7)
    despine(ax)

    # Time is folded into the x coordinate and voltage into y, so neither axis
    # can be read as a duration or a voltage. A scale bar for each is the only
    # honest way to put units on this view.
    bar_uv = _nice_round(peak_uv / 2)
    bar_ms = _nice_round(span_ms / 3)
    x0 = pos[:, 0].min() - 0.75 * dx
    y0 = pos[:, 1].min() - 1.2 * dy
    ax.plot([x0, x0], [y0, y0 + bar_uv * y_scale], lw=2.5, color="0.15",
            solid_capstyle="butt", clip_on=False)
    ax.text(x0, y0 + bar_uv * y_scale / 2, f" {bar_uv:.0f} µV", fontsize=7,
            va="center", ha="left")
    ax.plot([x0, x0 + bar_ms * x_scale], [y0, y0], lw=2.5, color="0.15",
            solid_capstyle="butt", clip_on=False)
    ax.text(x0 + bar_ms * x_scale / 2, y0, f"{bar_ms:g} ms ", fontsize=7,
            va="top", ha="center")

    ax.set_xlim(x0 - 0.35 * dx, pos[:, 0].max() + 0.75 * dx)
    ax.set_ylim(y0 - 0.9 * dy, pos[:, 1].max() + 1.9 * dy)
    return y_scale


def sample_waveforms(rec, ks, unit_id: int, n_waveforms: int = 100,
                     win_ms: float = 3.0, n_channels: int = 12,
                     highpass_hz: float | None = 300.0, show_template: bool = True,
                     title: str | None = None, figsize=(8, 9)):
    """Individual spike waveforms for one unit, drawn on the probe layout.

    The template is an average and hides everything an average hides. Drawing
    the individual snippets shows amplitude spread, and shows immediately when a
    "unit" is really two.

    The sorter's template is overlaid in blue by default. Both are in
    microvolts on the same channels and the same time origin, so agreement is a
    real check on the whole chain: gain, channel mapping and spike timing.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    t_ms, stack, chans = spike_snippets(rec, ks, unit_id, n_waveforms, win_ms,
                                        n_channels, highpass_hz)
    n_total = int((ks.spike_clusters == unit_id).sum())
    template = None
    if show_template:
        try:
            template = ks.template_uv(unit_id, chans)
        except Exception:  # noqa: BLE001 - the snippets are worth having alone
            template = None

    fig, ax = plt.subplots(figsize=figsize)
    draw_waveforms_on_probe(ax, ks, t_ms, stack, chans, template=template)
    ax.set_ylabel(DEPTH_LABEL)
    legend = "mean in red, template in blue" if template is not None \
        else "mean in red"
    ax.set_title(title or ks.title(
        f"Unit {unit_id}: {len(stack)} spikes of {n_total:,} ({legend})"))
    fig.tight_layout()
    return fig


def _nice_round(value: float) -> float:
    """Round to something a person would put on a scale bar."""
    if value <= 0:
        return 1.0
    exp = np.floor(np.log10(value))
    for m in (1, 2, 5, 10):
        candidate = m * 10 ** exp
        if candidate >= value * 0.8:
            return float(candidate)
    return float(10 ** (exp + 1))
