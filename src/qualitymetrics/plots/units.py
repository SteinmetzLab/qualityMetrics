"""Per-unit figures: where the units are, how big, and what shape.

Everything here is derived from the sorter output plus the recording's gain, so
these run on an archived session with no raw data present.
"""
from __future__ import annotations

import numpy as np

from ..ksdata import unit_summary_table  # noqa: F401  Moved to ``ksdata``; re-exported.
from ..style import AMP_LABEL, DEPTH_LABEL, despine, size_legend, use_lab_style

#: Marker areas for a unit that passed quality control and one that did not.
#: Binary, following the MATLAB original, which used 50 and 10. Size encoding a
#: continuous quantity here would be wasted: colour already carries firing rate
#: and waveform duration, and what a reader needs from the third channel is the
#: one bit that says whether to trust the unit at all.
QC_PASS_SIZE = 55.0
QC_FAIL_SIZE = 11.0


def amp_depth_scatter(ks, duration_s: float | None = None,
                      artifact_z_um: float | None = None,
                      title: str | None = None, figsize=(12, 8)):
    """Unit amplitude against depth, coloured by rate and by waveform width.

    Panel A colours by firing rate, on a log colour axis because rates span
    orders of magnitude. Panel B colours by trough-to-peak duration, the classic
    narrow versus broad spiking split.

    Marker size is binary: large for a unit that passed quality control, small
    for one that did not. The legend says which criterion was used, since it is
    our own gate when quality_metrics.json is present and Kilosort's good/mua
    label otherwise, and those are not the same question.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    use_lab_style()

    uids = [int(u) for u in ks.unit_ids]
    amp = np.array([ks.unit_amplitude_uv[u] for u in uids], dtype=float)
    # A unit with no measurable amplitude cannot sit on the log axis: drawn at
    # 1 uV, below any real unit, rather than silently dropped.
    amp = np.where(np.isfinite(amp) & (amp > 0), amp, 1.0)
    depth = np.array([ks.unit_depth_um[u] for u in uids])
    rate = np.array([v for v in
                     (ks.firing_rate(duration_s)[u] for u in uids)])
    ttp = np.array([ks.trough_to_peak_ms.get(u, np.nan) for u in uids])

    passed = np.array([bool(ks.qc_pass.get(u, False)) for u in uids])
    sizes = np.where(passed, QC_PASS_SIZE, QC_FAIL_SIZE)
    n_pass = int(passed.sum())

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
        # Log, because amplitudes span 20 to 1,000 uV and most units sat in
        # the lowest tenth of a linear axis (Nick, 2026-11-15).
        ax.set_xscale("log")
        despine(ax)
        if artifact_z_um is not None and np.isfinite(artifact_z_um):
            top = max(depth.max(), artifact_z_um)
            if artifact_z_um < top:
                ax.axhspan(artifact_z_um, top, color="0.88", zorder=0)
                ax.axhline(artifact_z_um, color="0.45", lw=0.8, ls="--")
    axes[0].set_ylabel(DEPTH_LABEL)
    size_legend(axes[0], [QC_PASS_SIZE, QC_FAIL_SIZE],
                [f"Pass ({n_pass})", f"Fail ({len(uids) - n_pass})"],
                f"QC, by {ks.qc_source}", loc="lower right")

    fig.suptitle(title or ks.title(f"Unit amplitude against depth "
                                   f"({len(uids)} units, {n_pass} passing QC)"),
                 fontsize=12)
    fig.tight_layout()
    return fig


#: A site counts toward a unit's footprint when its template is larger than
#: this, peak to peak (Nick, 2026-11-15).
FOOTPRINT_UV = 30.0
#: The template grid's picks, two of each, in the order drawn (Nick,
#: 2026-11-15), then RANDOM_PICKS more at random.
DIVERSE_PICKS = (
    ("Highest rate", "rate", True),
    ("Largest amplitude", "amplitude", True),
    ("Largest footprint", "footprint", True),
    ("Smallest footprint", "footprint", False),
    ("Broadest", "width", True),
    ("Narrowest", "width", False),
)
PICKS_EACH = 2
RANDOM_PICKS = 6


def template_footprints(ks, threshold_uv: float = FOOTPRINT_UV) -> dict[int, int]:
    """Number of sites on which each unit's template exceeds ``threshold_uv``
    peak to peak. Needs microvolts per bit."""
    if not ks.uv_per_bit:
        raise ValueError("footprints in microvolts need uv_per_bit")
    unw = ks.templates_unwhitened
    out = {}
    for uid in ks.unit_ids:
        uid = int(uid)
        if uid < len(unw):
            ptp = np.ptp(unw[uid], axis=0) * ks.uv_per_bit
            out[uid] = int((ptp > threshold_uv).sum())
    return out


def diverse_template_units(ks, duration_s: float | None = None,
                           seed: int = 0) -> list[tuple[int, str]]:
    """``(unit, why)`` for the template grid: two units at each extreme of
    :data:`DIVERSE_PICKS`, then :data:`RANDOM_PICKS` others at random.

    From the units that passed QC, where any did: among every unit, the
    smallest footprints and narrowest waveforms are mostly noise. No unit is
    picked twice; one already taken goes to the next in line. Ties go to the
    unit with more spikes. The random picks are seeded, so a redraw shows
    the same units.
    """
    have_template = {int(u) for u in ks.unit_ids if int(u) < len(ks.templates_unwhitened)}
    passing = {int(u) for u, ok in ks.qc_pass.items() if ok} & have_template
    pool = sorted(passing or have_template)
    if not pool:
        raise ValueError("no units with templates")
    counts = ks.n_spikes
    rate = ks.firing_rate(duration_s)
    measures = {
        "rate": {u: rate.get(u, np.nan) for u in pool},
        "amplitude": {u: ks.unit_amplitude_uv.get(u, np.nan) for u in pool},
        "footprint": template_footprints(ks),
        "width": ks.trough_to_peak_ms,
    }
    chosen: list[tuple[int, str]] = []
    taken: set[int] = set()
    for why, key, largest in DIVERSE_PICKS:
        values = measures[key]
        ranked = [u for u in pool if u not in taken
                  and np.isfinite(values.get(u, np.nan))]
        # Largest or smallest first; ties to the unit with more spikes.
        ranked.sort(key=lambda u: ((-values[u] if largest else values[u]),
                                   -counts.get(u, 0), u))
        for u in ranked[:PICKS_EACH]:
            chosen.append((u, f"{why}: {_value_text(key, values[u])}"))
            taken.add(u)
    rest = [u for u in pool if u not in taken]
    rng = np.random.default_rng(seed)
    for u in sorted(rng.choice(rest, size=min(RANDOM_PICKS, len(rest)),
                               replace=False).tolist()) if rest else []:
        chosen.append((int(u), "Random"))
    return chosen


def _value_text(key: str, value: float) -> str:
    return {"rate": f"{value:.1f} spikes/s",
            "amplitude": f"{value:.0f} µV",
            "footprint": (f"{value:.0f} site{'' if value == 1 else 's'} over "
                          f"{FOOTPRINT_UV:.0f} µV"),
            "width": f"{value:.2f} ms trough to peak"}[key]


def templates_grid(ks, unit_ids=None, n_channels: int = 32, n_cols: int = 6,
                   title: str | None = None, scale_uv: float | None = None,
                   reasons: dict | None = None):
    """A small multiple per unit: the template on its nearest sites, each
    waveform drawn where its site sits on the probe.

    ``reasons`` maps a unit to why it is shown (:func:`diverse_template_units`),
    written above its panel.

    The ``n_channels`` sites nearest the peak site by physical distance (32 is
    sixteen rows of two columns on NP2, 240 um of shank), each waveform at its
    site's position across and along the shank, the way phy draws them. A unit
    whose template spreads over 200 um then looks different from one confined
    to two rows, and one centered on the left column from one on the right.

    Each panel has its own vertical scale, the unit's largest excursion filling
    about three rows, so a 40 uV unit is as legible as a 600 uV one; its
    peak-to-peak amplitude is in its title. ``scale_uv`` fixes one scale for
    every panel instead (microvolts per three rows).
    """
    import matplotlib.pyplot as plt

    from .rawviews import _site_spacing
    use_lab_style()

    if unit_ids is None:
        unit_ids = [int(u) for u in ks.unit_ids]
    unit_ids = list(unit_ids)
    if not unit_ids:
        raise ValueError("no units to plot")

    unw = ks.templates_unwhitened
    peak = ks.peak_channel
    pos = np.asarray(ks.channel_positions, dtype=float)
    scale = ks.uv_per_bit or 1.0
    n_rows = int(np.ceil(len(unit_ids) / n_cols))
    t_ms = np.arange(unw.shape[1]) / ks.fs * 1000.0
    span_ms = float(t_ms[-1] - t_ms[0]) or 1.0
    dx, dy = _site_spacing(pos)
    x_scale = 0.85 * dx / span_ms                     # micrometres per ms
    n_columns = max(1, np.unique(np.round(pos[:, 0], 3)).size)
    site_rows = int(np.ceil(min(n_channels, len(pos)) / n_columns))

    fig, axes = plt.subplots(
        n_rows, n_cols, squeeze=False,
        figsize=(2.0 * n_cols, 0.8 + max(2.2, 0.21 * site_rows) * n_rows))

    for k, uid in enumerate(unit_ids):
        ax = axes[k // n_cols][k % n_cols]
        if uid >= len(unw):
            ax.axis("off")
            continue
        pk = int(peak[uid])
        # Nearest sites by physical distance, not by index: index order is a
        # wiring detail and jumps across columns on a multi-column probe.
        d = np.hypot(pos[:, 0] - pos[pk, 0], pos[:, 1] - pos[pk, 1])
        chans = np.argsort(d, kind="stable")[:n_channels]
        wave = unw[uid][:, chans] * scale
        biggest = scale_uv or (float(np.abs(wave).max()) or 1.0)
        y_scale = 3.0 * dy / biggest                  # micrometres per microvolt
        for j, c in enumerate(chans):
            ax.plot(pos[c, 0] + (t_ms - t_ms[0]) * x_scale,
                    pos[c, 1] + wave[:, j] * y_scale, lw=0.7,
                    color="crimson" if c == pk else "0.2")
        ax.set_xlim(pos[chans, 0].min() - 0.15 * dx,
                    pos[chans, 0].max() + 1.0 * dx)
        ax.set_ylim(pos[chans, 1].min() - 3.2 * dy, pos[chans, 1].max() + 3.2 * dy)
        # The amplitude only where microvolts are known: without a gain the
        # grid is still worth drawing, in arbitrary units.
        if ks.uv_per_bit:
            amp = ks.unit_amplitude_uv.get(int(uid), float("nan"))
            name = f"Unit {uid}, {amp:.0f} µV"
        else:
            name = f"Unit {uid}"
        why = (reasons or {}).get(int(uid))
        ax.set_title(f"{why}\n{name}" if why else name, fontsize=8)
        ax.set_xticks([])
        ax.tick_params(labelsize=6)
        ax.spines["bottom"].set_visible(False)
        despine(ax)
        if k % n_cols == 0:
            ax.set_ylabel(DEPTH_LABEL, fontsize=7)

    for k in range(len(unit_ids), n_rows * n_cols):
        axes[k // n_cols][k % n_cols].axis("off")

    # Time is folded into x, so the x axis is not a duration: one time bar.
    ax = axes[0][0]
    x_lo, _x_hi = ax.get_xlim()
    y_lo = ax.get_ylim()[0]
    ax.plot([x_lo, x_lo + 1.0 * x_scale], [y_lo, y_lo], lw=2, color="0.15",
            clip_on=False, solid_capstyle="butt")
    ax.text(x_lo + 0.5 * x_scale, y_lo, "1 ms", fontsize=7, va="top",
            ha="center", clip_on=False)

    scaling = (f"{scale_uv:g} µV per three rows" if scale_uv
               else "each scaled to its own peak, peak site in red")
    which = (f"{len(unit_ids)} units chosen for variety" if reasons
             else f"{len(unit_ids)} units")
    fig.suptitle(title or ks.title(
        f"Templates of {which}, {n_channels} nearest sites at their "
        f"positions ({scaling})"), fontsize=11)
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
