"""Diagnostics: figures for working out why a number came out the way it did.

Two here. A wideband version of the depth-power map, which shows the noise the
300 Hz version cuts off. And a per-unit view of the noise cutoff metric, drawn
after the one in IBL's data portal
(int-brain-lab/website, plots/static_plots.py, DataLoader.plot_noise_cutoff)
with the intermediate quantities of the calculation added, because on our data
the metric returns values in the thousands and the point is to see why.
"""
from __future__ import annotations

import numpy as np

from ..style import DEPTH_LABEL, TIME_LABEL, despine, use_lab_style

#: The metric's own constants, repeated here so the drawing and the calculation
#: cannot drift apart. They match qualitymetrics.metrics.
N_BINS = 100
QUANTILE_LENGTH = 0.25
PERCENT_THRESHOLD = 0.10
NC_THRESHOLD = 5.0


def depth_power_wide(rec, t_start_s: float = 100.0, win_s: float = 4.0,
                     depth_bin_um: float = 40.0, fmin_hz: float = 1.0,
                     normalize: str | None = None, title: str | None = None,
                     figsize=(12, 8)):
    """Power against depth across the whole band, on a logarithmic frequency axis.

    The companion depth_power stops at 300 Hz, which is where the interesting
    anatomy is but also where the interesting problems are not. Line noise and
    its harmonics, switching supplies and amplifier trouble all live above it,
    and on a linear axis the decades below 100 Hz are squeezed into nothing.

    Running to Nyquist on a log axis puts both in one figure: laminar structure
    on the left, noise on the right, and any narrow artifact between them
    visible as a vertical stripe crossing every depth.

    Drawn with pcolormesh rather than imshow because a log axis needs real bin
    edges. imshow would lay a log-spaced image on a linear grid and mislabel
    every frequency on it.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from scipy.signal import welch
    use_lab_style()

    geom = rec.geometry
    n = rec.n_data_channels
    depths = geom.y[:n] if geom is not None else np.arange(n, dtype=float)

    block = rec.get_traces(t_start_s, t_start_s + win_s)
    freqs, pxx = welch(block, fs=rec.fs,
                       nperseg=min(int(rec.fs), block.shape[0]), axis=0)
    keep = freqs >= fmin_hz
    freqs, pxx = freqs[keep], pxx[keep]

    edges = np.arange(depths.min(), depths.max() + depth_bin_um, depth_bin_um)
    bin_of = np.clip(np.digitize(depths, edges) - 1, 0, len(edges) - 2)
    image = np.array([pxx[:, bin_of == b].mean(axis=1) if np.any(bin_of == b)
                      else np.full(freqs.size, np.nan)
                      for b in range(len(edges) - 1)])       # (depth, freq)

    label = "Power spectral density (dB re 1 µV²/Hz)"
    with np.errstate(divide="ignore", invalid="ignore"):
        shown = 10 * np.log10(image)
    if normalize == "freq":
        shown = shown - np.nanmean(shown, axis=0, keepdims=True)
        label = "Power relative to depth mean (dB)"

    # pcolormesh wants one more edge than value on each axis. Geometric midpoints
    # are the right interior edges for a log axis.
    f_edges = np.concatenate([[freqs[0]], np.sqrt(freqs[:-1] * freqs[1:]),
                              [freqs[-1]]])

    fig, ax = plt.subplots(figsize=figsize)
    finite = shown[np.isfinite(shown)]
    vmin, vmax = np.percentile(finite, [2, 98]) if finite.size else (0.0, 1.0)
    mesh = ax.pcolormesh(f_edges, edges, np.ma.masked_invalid(shown),
                         cmap="viridis", norm=Normalize(vmin=vmin, vmax=vmax),
                         shading="auto")
    ax.set_xscale("log")
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(label)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel(DEPTH_LABEL)
    ax.set_xlim(freqs[0], freqs[-1])
    despine(ax)
    ax.set_title(title or f"{rec.path.name}: power against depth to "
                          f"{freqs[-1] / 1000:.0f} kHz "
                          f"({win_s:.0f} s from {t_start_s:.0f} s)")
    fig.tight_layout()
    return fig


def _noise_cutoff_parts(amps):
    """Every intermediate the noise cutoff formula uses, for drawing.

    Recomputed here exactly as in metrics.noise_cutoff rather than returned from
    it, because the metric is vendored verbatim from ibllib and must stay that
    way. If the two ever disagree the figure is wrong, so the caller checks the
    reported value against the metric.
    """
    amps = np.asarray(amps, dtype=float)
    edges = np.linspace(0, np.max(amps), N_BINS)
    n, _ = np.histogram(amps, bins=edges)
    centres = (edges[:-1] + edges[1:]) / 2

    idx_peak = int(np.argmax(n))
    length_top_half = len(np.where(n[idx_peak:-1] > 0)[0])
    start = int(np.ceil(2 * QUANTILE_LENGTH * length_top_half + idx_peak))
    high_idx = np.arange(start, len(n))
    used = high_idx[n[high_idx] >= 1] if high_idx.size else high_idx

    parts = {
        "counts": n, "centres": centres, "edges": edges,
        "idx_peak": idx_peak, "peak_height": int(np.max(n)) if n.size else 0,
        "high_start": start, "high_idx": high_idx, "used_idx": used,
        "mean_high": float(n[used].mean()) if used.size else float("nan"),
        "std_high": float(n[used].std()) if used.size else float("nan"),
    }
    nonzero = n[n != 0]
    parts["first_low"] = float(nonzero[1]) if nonzero.size > 1 else float("nan")
    if np.isfinite(parts["std_high"]) and parts["std_high"] > 0:
        parts["cutoff"] = ((parts["first_low"] - parts["mean_high"])
                           / parts["std_high"])
    else:
        parts["cutoff"] = float("nan")
    return parts


def noise_cutoff_diagnostic(ks, unit_ids=None, n_units: int = 5,
                            time_bins: int = 120, amp_bins: int = 100,
                            title: str | None = None, figsize=None):
    """Show how the noise cutoff metric arrives at its number, per unit.

    Two panels per unit, following IBL's data portal. Left: spike amplitude
    against time in the session as a density, so a unit that fades or steps is
    obvious. Right: the amplitude histogram the metric actually works on, drawn
    horizontally against the same axis, with the pieces of the formula marked.

    Marked on the histogram:

    * the peak bin, and the 10% of peak line that is the metric's second
      criterion;
    * the high-amplitude quantile the metric takes as its reference, shaded,
      with its mean and one standard deviation;
    * the "first low bin", the second occupied bin from zero, which is the value
      the metric compares against that reference.

    Those last two are what make the number legible. The metric is
    (first low bin - mean of the reference) / std of the reference, so a unit
    whose distribution is strongly peaked near its low end and has a thin
    high-amplitude tail produces an enormous value without anything being wrong
    with the histogram, which is exactly the regime our data sits in.

    Density is drawn with a 2D histogram rather than mpl-scatter-density, which
    is what IBL uses, to avoid the extra dependency.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    use_lab_style()

    if unit_ids is None:
        unit_ids = units_spanning_noise_cutoff(ks, n_units)
    unit_ids = [int(u) for u in unit_ids]
    if not unit_ids:
        raise ValueError("no units to show")

    stored = ks.metrics_by_unit()
    amps_all = ks.spike_amplitudes_uv()
    times_all = ks.spike_times_s
    n = len(unit_ids)
    figsize = figsize or (11, 3.1 * n)
    fig, axes = plt.subplots(n, 2, figsize=figsize, squeeze=False,
                             gridspec_kw={"width_ratios": [3, 1.5],
                                          "wspace": 0.05})

    for row, uid in enumerate(unit_ids):
        sel = ks.spike_clusters == uid
        amps, times = amps_all[sel], times_all[sel]
        ax_t, ax_h = axes[row]
        if amps.size < 3:
            ax_t.text(0.5, 0.5, f"Unit {uid}: too few spikes",
                      transform=ax_t.transAxes, ha="center", va="center")
            despine(ax_t)
            despine(ax_h)
            continue

        parts = _noise_cutoff_parts(amps)
        top = float(np.max(amps))

        # -- amplitude against time, as a density
        counts, xe, ye = np.histogram2d(
            times, amps, bins=[time_bins, amp_bins],
            range=[[0, max(times.max(), 1.0)], [0, top]])
        shown = np.ma.masked_where(counts.T <= 0, counts.T)
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad("white")
        ax_t.pcolormesh(xe, ye, shown, cmap=cmap,
                        norm=LogNorm(vmin=1, vmax=max(counts.max(), 2)),
                        shading="auto")
        ax_t.set_ylim(0, top)
        ax_t.set_ylabel("Spike amplitude (µV)")
        despine(ax_t)

        # -- the histogram the metric works on
        ax_h.barh(parts["centres"], parts["counts"],
                  height=np.diff(parts["edges"]).mean(), color="0.35",
                  linewidth=0)
        ax_h.axvline(PERCENT_THRESHOLD * parts["peak_height"], color="tab:blue",
                     lw=1.2, ls="--",
                     label=f"{PERCENT_THRESHOLD:.0%} of peak bin")
        if parts["high_idx"].size:
            lo = parts["centres"][parts["high_idx"][0]]
            ax_h.axhspan(lo, top, color="tab:orange", alpha=0.13, zorder=0)
            ax_h.axvline(parts["mean_high"], color="tab:orange", lw=1.2,
                         label="Mean of reference tail")
            if np.isfinite(parts["std_high"]):
                ax_h.axvspan(parts["mean_high"] - parts["std_high"],
                             parts["mean_high"] + parts["std_high"],
                             color="tab:orange", alpha=0.22, zorder=0)
        if np.isfinite(parts["first_low"]):
            ax_h.axvline(parts["first_low"], color="crimson", lw=1.2,
                         label="First low bin")
        ax_h.set_ylim(0, top)
        ax_h.set_yticklabels([])
        ax_h.set_xlabel("Spikes in bin (count)")
        despine(ax_h)
        if row == 0:
            ax_h.legend(loc="upper right", fontsize=7)

        reported = (stored.get(uid, {}) or {}).get("noise_cutoff")
        drawn = parts["cutoff"]
        # If the drawing and the pipeline disagree, say so on the figure rather
        # than letting a reader trust a picture of a different calculation.
        agree = (reported is None or not np.isfinite(reported)
                 or not np.isfinite(drawn)
                 or abs(reported - drawn) <= 0.02 * max(abs(reported), 1.0))
        pct = (100.0 * parts["first_low"] / parts["peak_height"]
               if parts["peak_height"] else float("nan"))
        verdict = "pass" if (np.isfinite(drawn) and drawn <= NC_THRESHOLD) \
            else "fail"
        colour = "tab:green" if verdict == "pass" else "crimson"
        note = "" if agree else f"  (pipeline recorded {reported:.2f})"
        ax_t.set_title(
            f"Unit {uid}: noise cutoff {drawn:,.1f} ({verdict} at "
            f"{NC_THRESHOLD:.0f}), first low bin {pct:.0f}% of peak, "
            f"{amps.size:,} spikes{note}", fontsize=9, color=colour)

    axes[-1][0].set_xlabel(TIME_LABEL)
    fig.suptitle(title or ks.title("Noise cutoff, and how it is computed"),
                 fontsize=12)
    # subplots_adjust rather than tight_layout: the two columns are deliberately
    # butted together by the gridspec, and tight_layout would pull them apart
    # again and warn that it cannot honour the width ratios.
    fig.subplots_adjust(left=0.09, right=0.98, top=1 - 0.5 / max(n, 1),
                        bottom=0.06, hspace=0.55, wspace=0.05)
    return fig


def units_spanning_noise_cutoff(ks, n_units: int = 5, min_spikes: int = 100):
    """Units evenly spread across the recorded range of noise cutoff values.

    Picked from what the pipeline recorded, so the examples represent the
    distribution actually in the data rather than the ones that happen to be
    interesting.
    """
    stored = ks.metrics_by_unit()
    have = [(u, v.get("noise_cutoff")) for u, v in stored.items()
            if v.get("noise_cutoff") is not None
            and np.isfinite(v["noise_cutoff"])
            and ks.n_spikes.get(u, 0) >= min_spikes]
    if not have:
        counts = ks.n_spikes
        return [u for u, _ in sorted(counts.items(),
                                     key=lambda kv: -kv[1])[:n_units]]
    have.sort(key=lambda uv: uv[1])
    idx = np.linspace(0, len(have) - 1, min(n_units, len(have)))
    return [have[int(round(i))][0] for i in idx]
