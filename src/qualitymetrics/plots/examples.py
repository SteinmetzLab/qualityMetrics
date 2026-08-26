"""One figure that says everything about three representative units.

A QC report full of population summaries never answers the question a person
actually has, which is "what does a unit from this recording look like". Three
columns, one per unit, three rows: the raw snippets, the template across the
whole shank, and the autocorrelogram. The metrics are printed on the figure, so
the numbers and the thing they describe are in the same place.

The three units are chosen at percentiles of the amplitude distribution rather
than by rank. The single largest unit in a recording is usually atypical, which
is why the 95th percentile is a better "big unit" than the maximum.
"""
from __future__ import annotations

import numpy as np

from ..style import AMP_LABEL, DEPTH_LABEL, despine, use_lab_style
from .rawviews import draw_waveforms_on_probe, spike_snippets


def draw_template_waterfall(ax, ks, unit_id: int, depth_pad_um: float = 250.0):
    """One unit's template over depth, as an image, into an existing axes.

    Cropped to a window around the unit rather than the whole shank: a template
    occupies maybe 200 um of a 2865 um probe, and drawing all of it makes the
    part that matters two pixels tall.
    """
    unw = ks.templates_unwhitened
    if unit_id >= len(unw):
        raise ValueError(f"unit {unit_id} has no template")
    pos = ks.channel_positions
    order = np.argsort(pos[:, 1])
    wave = unw[unit_id][:, order] * (ks.uv_per_bit or 1.0)
    depths = pos[order, 1]

    centre = ks.unit_depth_um.get(int(unit_id), float(np.median(depths)))
    keep = np.abs(depths - centre) <= depth_pad_um
    if keep.sum() < 4:
        keep = np.ones_like(depths, dtype=bool)
    wave, depths = wave[:, keep], depths[keep]

    t_ms = np.arange(wave.shape[0]) / ks.fs * 1000.0
    lim = float(np.abs(wave).max()) or 1.0
    im = ax.imshow(wave.T, aspect="auto", origin="lower", cmap="RdBu_r",
                   vmin=-lim, vmax=lim,
                   extent=[t_ms[0], t_ms[-1], depths.min(), depths.max()],
                   interpolation="nearest")
    ax.set_xlabel("Time (ms)")
    despine(ax)
    return im, lim


def draw_acg(ax, ks, unit_id: int, window_ms: float = 50.0,
             bin_ms: float = 0.1, log_x: bool = True,
             rp_min_val_ms: float | None = None):
    """Autocorrelogram into an existing axes, on a logarithmic lag axis.

    Log-spaced lags after the manner of myACG in cortex-lab/spikes. A linear
    axis out to 50 ms puts the entire refractory question into the first two
    pixels, and whether a dip is 0.25 ms wide or 2 ms wide is exactly what
    decides if it is a neuron, a deduplication window, or structured noise.
    Only the positive half is drawn, since the function is symmetric by
    construction and mirroring it wastes half the axis.

    rp_min_val_ms, from slidingRP, marks the refractory duration at which that
    unit's contamination estimate was minimised. That is the width the metric
    itself settled on, which beats any fixed band drawn by eye.

    Returns nothing to be interpreted as a metric. An earlier version returned
    a ratio of counts inside a 2 ms window to counts in the flanks and printed
    it on the figure, which was a worse version of what slidingRP computes
    properly, and it invited exactly the wrong reading.
    """
    lags, counts = ks.autocorrelogram(unit_id, window_ms, bin_ms)
    if lags.size == 0:
        ax.text(0.5, 0.5, "Too few spikes", transform=ax.transAxes,
                ha="center", va="center", fontsize=9, color="0.4")
        despine(ax)
        return

    positive = lags > 0
    lag, count = lags[positive], counts[positive]
    ax.step(lag, count, where="mid", color="0.25", lw=1.0)
    ax.fill_between(lag, count, step="mid", color="0.35", alpha=0.5)

    if rp_min_val_ms is not None and np.isfinite(rp_min_val_ms)             and rp_min_val_ms > 0:
        ax.axvline(rp_min_val_ms, color="crimson", lw=1.2, ls="--")
        ax.text(rp_min_val_ms, ax.get_ylim()[1], f" {rp_min_val_ms:.2f} ms",
                color="crimson", fontsize=7, va="top", ha="left")

    if log_x:
        ax.set_xscale("log")
        ax.set_xlim(max(bin_ms / 2, 0.05), window_ms)
    else:
        ax.set_xlim(0, window_ms)
    ax.set_xlabel("Lag (ms)")
    despine(ax)


def example_neurons(rec, ks, unit_ids=None, percentiles=(95, 75, 50),
                    n_waveforms: int = 150, win_ms: float = 3.0,
                    n_channels: int = 12, acg_window_ms: float = 50.0,
                    passing_only: bool = True, title: str | None = None,
                    figsize=(15, 14)):
    """Three representative units, three views each.

    Rows are raw snippets on the probe layout, template over depth, and
    autocorrelogram. Columns are units taken at the given percentiles of the
    amplitude distribution, largest on the left, restricted to units that passed
    quality control. A rejected unit can be anything at all, so it says nothing
    about what the recording will actually yield.

    Each column is annotated with that unit's metrics, and the top panel carries
    the sorter's template in blue over the mean snippet in red.
    """
    import matplotlib.pyplot as plt
    use_lab_style()

    if unit_ids is None:
        unit_ids = ks.units_at_amplitude_percentiles(
            percentiles, passing_only=passing_only)
        labels = [f"{p}th percentile of QC-passing amplitude"
                  if passing_only else f"{p}th percentile amplitude"
                  for p in percentiles]
    else:
        unit_ids = [int(u) for u in unit_ids]
        labels = [""] * len(unit_ids)
    if not unit_ids:
        raise ValueError("no units to show")

    n = len(unit_ids)
    fig, axes = plt.subplots(3, n, figsize=figsize, squeeze=False)

    for col, uid in enumerate(unit_ids):
        # -- row 0: raw snippets on the probe layout, template overlaid
        ax = axes[0][col]
        try:
            t_ms, stack, chans = spike_snippets(
                rec, ks, uid, n_waveforms=n_waveforms, win_ms=win_ms,
                n_channels=n_channels)
            try:
                template = ks.template_uv(uid, chans)
            except Exception:  # noqa: BLE001 - snippets still worth drawing
                template = None
            draw_waveforms_on_probe(ax, ks, t_ms, stack, chans,
                                    template=template,
                                    show_depth_labels=(col == 0))
            shown = f"{len(stack)} of {ks.n_spikes.get(uid, 0):,} spikes"
        except Exception as exc:  # noqa: BLE001 - one bad unit must not stop it
            ax.text(0.5, 0.5, f"No snippets\n({type(exc).__name__})",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="0.4")
            despine(ax)
            shown = "no snippets"
        header = f"Unit {uid}"
        if labels[col]:
            header += f"\n{labels[col]}"
        ax.set_title(f"{header}\n{shown}: mean red, template blue", fontsize=9)

        summary = ks.qc_summary(uid)
        if summary:
            ax.text(0.02, 0.98, "\n".join(summary), transform=ax.transAxes,
                    va="top", ha="left", fontsize=7.5, color="0.15",
                    bbox={"facecolor": "white", "alpha": 0.78,
                          "edgecolor": "0.75", "boxstyle": "round,pad=0.35"})

        # -- row 1: the template across depth
        ax = axes[1][col]
        try:
            im, lim = draw_template_waterfall(ax, ks, uid)
            cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.05)
            cbar.set_label(AMP_LABEL, fontsize=8)
            cbar.ax.tick_params(labelsize=7)
            # Largest absolute excursion, which is not the peak-to-peak figure
            # printed in the metrics box. Saying which avoids the two numbers
            # on the same column looking like a contradiction.
            ax.set_title(f"Template, largest excursion {lim:.0f} µV",
                         fontsize=9)
        except Exception as exc:  # noqa: BLE001
            ax.text(0.5, 0.5, f"No template\n({type(exc).__name__})",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="0.4")
            despine(ax)

        # -- row 2: the autocorrelogram, annotated from slidingRP rather than
        # from anything measured off the picture.
        ax = axes[2][col]
        rp = (ks.metrics_by_unit().get(uid, {}) or {})
        rp_min = rp.get("sliding_rp_rp_min_val_ms")
        draw_acg(ax, ks, uid, window_ms=acg_window_ms,
                 rp_min_val_ms=rp_min)
        cont = rp.get("sliding_rp_contamination")
        if cont is None:
            note = ""
        elif not np.isfinite(cont):
            note = ", contamination above the 35% cap"
        else:
            note = f", contamination {cont:.1f}%"
        ax.set_title(f"Autocorrelogram{note}", fontsize=9)

    axes[0][0].set_ylabel(DEPTH_LABEL)
    axes[1][0].set_ylabel(DEPTH_LABEL)
    axes[2][0].set_ylabel("Spike pairs (count)")

    fig.suptitle(title or ks.title(
        f"Example units at the {', '.join(str(p) for p in percentiles)} "
        f"percentiles of amplitude"), fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig
