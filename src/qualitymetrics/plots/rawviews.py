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


def wall_heatmap(rec, t_start_s: float = 100.0, win_ms: float = 25.0,
                 highpass_hz: float | None = 300.0, cmr: bool = True,
                 clim_uv: float | None = None, clim_pct: float = 99.5,
                 artifact_z_um: float | None = None, title: str | None = None,
                 figsize=(13, 8)):
    """Voltage across the whole shank for a short window: the "wall" view.

    Depth on y, time on x, voltage as a diverging colour centred on zero, so a
    spike is a dark spot and a bad channel is a stripe. This is the fastest way
    to see dead channels, a saturating amplifier, or line noise, and it is the
    figure to look at before trusting anything downstream.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    order, depths = _depth_order(rec)
    t1 = t_start_s + win_ms / 1000.0
    traces = rec.get_traces(t_start_s, t1, channels=order)
    if highpass_hz:
        traces = preprocess_ap(traces, rec.fs, highpass_hz, cmr=cmr)

    if clim_uv is None:
        clim_uv = float(np.percentile(np.abs(traces), clim_pct)) or 1.0

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(traces.T, aspect="auto", origin="lower", cmap="RdBu_r",
                   vmin=-clim_uv, vmax=clim_uv,
                   extent=[0, win_ms, depths.min(), depths.max()],
                   interpolation="nearest")
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("Voltage (µV)")
    ax.set_xlabel("Time from window start (ms)")
    ax.set_ylabel(DEPTH_LABEL)
    despine(ax)

    if artifact_z_um is not None and np.isfinite(artifact_z_um) \
            and artifact_z_um < depths.max():
        ax.axhline(artifact_z_um, color="0.2", lw=1.0, ls="--")
        ax.text(win_ms, artifact_z_um, f" Above {artifact_z_um:.0f} µm ",
                ha="right", va="bottom", fontsize=8, color="0.2")

    band = f"high-passed at {highpass_hz:.0f} Hz" if highpass_hz else "wideband"
    ref = ", common median referenced" if (highpass_hz and cmr) else ""
    ax.set_title(title or f"{rec.path.name}: {win_ms:.0f} ms from "
                          f"{t_start_s:.0f} s ({band}{ref})")
    fig.tight_layout()
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
        ax.plot(times_ms[sel], np.full(sel.sum(), row * spacing_uv),
                marker="|", linestyle="none", ms=9, mew=1.4,
                color=colors[i], label=f"Unit {uid}")
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


def sample_waveforms(rec, ks, unit_id: int, n_waveforms: int = 100,
                     win_ms: float = 3.0, n_channels: int = 6,
                     highpass_hz: float | None = 300.0,
                     title: str | None = None, figsize=(9, 7)):
    """Individual spike waveforms for one unit, pulled from the raw file.

    The template is an average and hides everything an average hides. Drawing
    the individual snippets shows amplitude spread, and shows immediately when a
    "unit" is really two.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    st = ks.spike_times_s[ks.spike_clusters == unit_id]
    if st.size == 0:
        raise ValueError(f"unit {unit_id} has no spikes")
    rng = np.random.default_rng(0)
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
            block = preprocess_ap(block, rec.fs, highpass_hz, cmr=False)
        if block.shape[0] >= n_samp:
            snippets.append(block[:n_samp])
    if not snippets:
        raise ValueError(f"no usable snippets for unit {unit_id}")
    stack = np.stack(snippets)                       # (n, samples, channels)

    t_ms = (np.arange(n_samp) / rec.fs - half) * 1000.0
    step = float(np.percentile(np.abs(stack), 99.8)) * 2.2 or 50.0

    fig, ax = plt.subplots(figsize=figsize)
    for j in range(len(chans)):
        ax.plot(t_ms, stack[:, :, j].T + j * step, lw=0.35, color="0.65",
                alpha=0.5)
        ax.plot(t_ms, stack[:, :, j].mean(0) + j * step, lw=1.6, color="crimson")
    ax.set_yticks(np.arange(len(chans)) * step)
    ax.set_yticklabels([f"{pos[c, 1]:.0f}" for c in chans], fontsize=8)
    ax.set_ylabel(DEPTH_LABEL)
    ax.set_xlabel("Time from spike (ms)")
    despine(ax)

    bar = _nice_round(step / 2.2)
    ax.plot([t_ms[-1]] * 2, [0, bar], lw=3, color="black",
            solid_capstyle="butt", clip_on=False)
    ax.text(t_ms[-1], bar / 2, f" {bar:.0f} µV", va="center", fontsize=9)

    ax.set_title(title or ks.title(
        f"Unit {unit_id}: {len(snippets)} spikes of {st.size:,} "
        f"(mean in red)"))
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
