"""Saturation: how many samples sat at the ADC's rails, channel by channel.

Samples at the rail carry no signal, and a spike there is lost whatever the
pipeline does. On FD_013 2026-07-29 an optogenetic laser drove the ADC to its
rail for up to ~1 ms on 9-22% of pulse-channels, and those samples were what
no artifact subtraction could recover. Other causes (licking, grounding, a
channel whose DC offset sits near the rail) look the same here, so the figure
is made for every session, light or not.

The counts are made by SortingManager's worker while it reads every raw sample
for the LFP (sortingmanager/saturation.py) and saved as saturation.npz in the
sort folder: per channel, samples at the negative and positive rail and the
number of separate events; and per second, saturated samples per channel.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..style import TIME_LABEL, color_legend, use_lab_style

SATURATION_FILE = "saturation.npz"
NEG_COLOR, POS_COLOR, ALL_COLOR = "#1f77b4", "#d62728", "#222222"


def load_counts(sorting_dir) -> dict | None:
    """The counts in a sort folder, or None if it has none (sorts made before
    the worker counted saturation, 2026-10)."""
    path = Path(sorting_dir) / SATURATION_FILE
    if not path.exists():
        return None
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def depth_order(counts: dict):
    """Channel indices from the deepest site to the shallowest, or None when the
    counts carry no positions. On NP2 the channel order is not the depth order
    (on a Quad shank it runs 0, 2880, 15, 2895, ... um), so plotting by channel
    number would turn one band of saturation into a comb."""
    if "z_um" not in counts:
        return None
    z = np.asarray(counts["z_um"], float)
    x = np.asarray(counts.get("x_um", np.zeros_like(z)), float)
    if z.size != np.asarray(counts["negative"]).size:
        return None
    return np.lexsort((x, z))


def saturation_report(counts: dict, title: str | None = None, figsize=(12, 8)):
    """Three panels: saturated samples per channel (log scale, zero drawn at the
    bottom), in depth order when the counts carry site positions; then, over the
    session, saturated samples per second summed over channels, and the number
    of channels saturating in each second."""
    import matplotlib.pyplot as plt

    use_lab_style()
    neg = np.asarray(counts["negative"], np.int64)
    pos = np.asarray(counts["positive"], np.int64)
    per_bin = np.asarray(counts["per_bin"], np.int64).reshape(-1, neg.size)
    bin_s = float(counts["bin_s"])
    n_ch = neg.size
    total = neg + pos
    fs = float(counts["fs"])
    n_samples = int(counts["n_samples"])
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize,
                                        gridspec_kw={"height_ratios": [1.4, 1, 1]})
    order = depth_order(counts)
    pos_x = np.arange(n_ch)
    shown = order if order is not None else pos_x
    ax1.plot(pos_x, total[shown], color=ALL_COLOR, lw=0.6, alpha=0.5)
    ax1.plot(pos_x, neg[shown], "o", ms=3, color=NEG_COLOR, label="At the negative rail")
    ax1.plot(pos_x, pos[shown], "o", ms=3, color=POS_COLOR, label="At the positive rail")
    # Symmetric log, linear below 1, so channels that never saturated sit at zero
    # instead of vanishing off a log axis.
    ax1.set_yscale("symlog", linthresh=1, linscale=0.5)
    ax1.set_ylim(-0.3, max(10, total.max(initial=0) * 3))
    ax1.set_xlim(-1, n_ch)
    ax1.set_ylabel("Saturated samples (count)")
    if order is not None:
        # Ticks name the channels at their place in depth order, and the depth
        # itself goes underneath, so a band at one depth reads as a band.
        z = np.asarray(counts["z_um"], float)[order]
        ticks = np.unique(np.linspace(0, n_ch - 1, min(n_ch, 9)).round().astype(int))
        ax1.set_xticks(ticks)
        ax1.set_xticklabels([f"{order[t]}\n{z[t]:.0f}" for t in ticks])
        ax1.set_xlabel("Channel, in depth order (channel number above, depth in µm below)")
    else:
        ax1.set_xlabel("Channel")
    seconds = total / fs
    pinned = np.flatnonzero(total >= 0.9 * n_samples) if n_samples else np.zeros(0, int)
    note = (f"{int((total > 0).sum())} of {n_ch} channels saturated at some point; "
            f"{seconds.sum():.3f} channel-seconds in all ({total.sum() / max(1, n_samples * n_ch) * 100:.4f}% of samples)")
    if pinned.size:
        note += f"; pinned at the rail throughout (dead): {pinned.tolist()}"
    if "events" in counts:
        note += f"; {int(np.asarray(counts['events']).sum())} separate events"
    if not total.any():
        note = "No sample reached the ADC's rails on any channel"
    if "minimum" in counts and "maximum" in counts:
        # The threshold assumes where the rails are; the extremes show whether
        # the data ever got there, so a zero here cannot hide a lower rail.
        note += (f"\nThreshold +/-{float(counts['threshold_bits']):.0f} bits; most extreme raw "
                 f"values {int(np.min(counts['minimum']))} and {int(np.max(counts['maximum']))} bits")
    ax1.set_title(note, fontsize=9)
    if total.any():
        color_legend(ax1, loc="upper right", fontsize=8)
    # Over time: leave out channels pinned throughout, which would swamp it.
    live = np.setdiff1d(np.arange(n_ch), pinned)
    t = (np.arange(per_bin.shape[0]) + 0.5) * bin_s
    per_s = per_bin[:, live].sum(axis=1)
    n_sat = (per_bin[:, live] > 0).sum(axis=1)
    left_out = " (channels pinned throughout left out)" if pinned.size else ""
    ax2.plot(t, per_s, color=ALL_COLOR, lw=0.7)
    ax2.set_yscale("symlog", linthresh=1, linscale=0.5)
    ax2.set_ylim(-0.3, max(10, per_s.max(initial=0) * 3))
    ax2.set_ylabel(f"Saturated samples\nper {bin_s:g} s (count)")
    ax2.set_title("Saturated samples over the session, all channels" + left_out, fontsize=9)
    ax3.plot(t, n_sat, color=ALL_COLOR, lw=0.7)
    ax3.set_ylim(-0.5, max(5, n_sat.max(initial=0) * 1.1))
    ax3.set_ylabel(f"Channels saturating\nper {bin_s:g} s (count)")
    ax3.set_title("Channels at the rail in each second" + left_out, fontsize=9)
    for ax in (ax2, ax3):
        ax.set_xlabel(TIME_LABEL)
        if t.size:
            ax.set_xlim(0, t[-1] + bin_s / 2)
    if title:
        fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    return fig
