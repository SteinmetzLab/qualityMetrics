"""Session noise figures: how each processing stage changes what the probe sees.

Ported from the lab's quad-base noise testing (MATLAB, ``quad_noiseTesting``),
which compared grounding and reference configurations recorded one after
another. Here the comparison is across the shanks of one session instead, so
the figures answer: how noisy is each shank, in which frequency bands, at which
depths, and how much of it is common to the shank (removed by a global median
reference) or common to channels digitised at the same instant (removed by the
"demux" reference, one median per group of simultaneously sampled channels).

The four stages, in the order of the original's ``rmsCalcNPQuad``:

1. **Raw**: each channel's median removed.
2. **Band-pass**: 500 Hz to 10 kHz, 3rd-order Butterworth. Zero-phase here
   (the original filtered causally), with a margin read on each side and
   dropped, so the one-second windows have no edge transient.
3. **+ Global CAR**: the shank's across-channel median removed sample by
   sample, from the raw data, *then* band-passed.
4. **+ Demux CAR**: after the global CAR, within each group of channels sampled
   at the same instant, that group's median removed, then band-passed. Groups
   come from the NP2 multiplexing, the same rule the sorting pipeline's ADC
   phase correction uses (see :func:`adc_groups`), rather than a hard-coded
   channel list. Within a shank only, as in the original.

Referencing happens before the filter, not after, deliberately. A median
across channels is not a linear operation, so taking it after band-passing
puts power back outside the band: on FD_008 it raised the spectrum below
500 Hz from the filter's floor to about -35 dB. (The narrow lines at 3214.4
and 4790.0 Hz on FD_008 are not that either: they are in the raw data, on
every channel, about 31 dB over the floor, and no referencing removes them.
Fixed in frequency across a month and seen only on one rig, so interference
from equipment rather than anything in the processing.)

A shank is its own column or panel everywhere. Stages are overlaid within a
shank; shanks are never overlaid on one another, because their channels, depths
and referencing are not comparable one to one.

Data: one second each from early, middle and late in the recording for the
spectra, snippets and per-channel RMS, and twelve one-second windows spread
across it for the RMS over time. Read from the archived ``.cbin``.
"""
from __future__ import annotations

import numpy as np

from ..style import DEPTH_LABEL, color_legend, despine, use_lab_style

STAGES = ("Raw", "Band-pass 0.5-10 kHz", "+ Global CAR", "+ Demux CAR")
STAGE_COLORS = ("#8c8c8c", "#1f77b4", "#2ca02c", "#d62728")
BAND_HZ = (500.0, 10_000.0)
#: Where in the recording the three detailed windows sit (early, middle, late).
WINDOW_FRACTIONS = (0.1, 0.5, 0.9)
WINDOW_S = 1.0
#: Read on each side of a window and dropped, so the zero-phase filter's edges
#: fall outside what is measured.
MARGIN_S = 0.1
#: One-second windows spread across the recording for the RMS over time.
N_TIME_WINDOWS = 12
#: Channels read in blocks of 32, two at a time: within a block a pair shares
#: an instant, successive pairs step by one sixteenth of a sample (NP2).
ADC_BLOCK = 32
#: Log-spaced frequency bins for the depth-power images, as the original's
#: max pooling into 300 bins from 2 Hz.
FREQ_BINS = 300
FREQ_MIN_HZ = 2.0
#: Welch segment length: 1.8 Hz resolution at 30 kHz.
WELCH_SAMPLES = 16384
#: Color limits of the depth-power images and snippets, as in the original.
POWER_LIMITS_DB = (-25.0, 40.0)
SNIPPET_LIMITS_UV = 50.0
SNIPPET_SAMPLES = 1000


def adc_groups(n_channels: int) -> np.ndarray:
    """Group id per channel: channels with the same id are sampled together.

    NP2 multiplexes 24 ADCs over a shank's 384 channels, each ADC reading 16 in
    turn, so 24 channels share each of 16 instants. The index within a block
    of 32, halved, is which instant: the same rule as SortingManager's
    inter_sample_shift and ibl-neuropixel's trace_header(version=2), and the
    same grouping as the original MATLAB's hard-coded ADC table.
    """
    return (np.arange(n_channels) % ADC_BLOCK) // 2


def _group_median(x: np.ndarray) -> np.ndarray:
    """Each channel's ADC-group median, sample by sample, shaped like ``x``."""
    groups = adc_groups(x.shape[1])
    counts = np.bincount(groups)
    if counts.min() == counts.max():
        # Every group the same size (a full shank): one median over a
        # (samples, groups, members) view, rather than a median per group.
        order = np.argsort(groups, kind="stable")
        medians = np.median(x[:, order].reshape(x.shape[0], counts.size, -1), axis=2)
        return medians[:, groups]
    out = np.empty_like(x)
    for group in np.unique(groups):
        members = groups == group
        out[:, members] = np.median(x[:, members], axis=1, keepdims=True)
    return out


def _referenced(block_uv: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Raw (channel medians removed), + global CAR, + demux CAR; unfiltered."""
    x = np.asarray(block_uv, dtype=np.float32)
    # The per-channel offset from every tenth sample: an offset estimate, and
    # a tenth of the cost of the one step that dominated this function.
    raw = x - np.median(x[::10], axis=0, keepdims=True)
    global_car = raw - np.median(raw, axis=1, keepdims=True)
    demux = global_car - _group_median(global_car)
    return raw, global_car, demux


def _band_pass(x: np.ndarray, fs: float) -> np.ndarray:
    from scipy import signal

    sos = signal.butter(3, BAND_HZ, btype="bandpass", fs=fs, output="sos")
    return signal.sosfiltfilt(sos, x, axis=0).astype(np.float32)


def processing_stages(block_uv: np.ndarray, fs: float, *, margin: int = 0) -> list:
    """The four stages of one block, (samples, channels) microvolts.

    ``margin`` samples at each end are used by the filter and then dropped.
    """
    raw, global_car, demux = _referenced(block_uv)
    keep = slice(margin, raw.shape[0] - margin if margin else None)
    return [x[keep] for x in (raw, _band_pass(raw, fs), _band_pass(global_car, fs),
                              _band_pass(demux, fs))]


def demux_band(block_uv: np.ndarray, fs: float, *, margin: int = 0) -> np.ndarray:
    """Only the last stage, for the RMS over time, which needs nothing else."""
    demux = _referenced(block_uv)[2]
    keep = slice(margin, demux.shape[0] - margin if margin else None)
    return _band_pass(demux, fs)[keep]


def log_bins(freqs: np.ndarray, n_bins: int = FREQ_BINS,
             f_min: float = FREQ_MIN_HZ) -> tuple[np.ndarray, np.ndarray]:
    edges = np.logspace(np.log10(f_min), np.log10(freqs.max()), n_bins + 1)
    return edges, np.sqrt(edges[:-1] * edges[1:])


def pool_max(power: np.ndarray, freqs: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Max over the frequencies in each log bin; a bin with none repeats the last.

    Max rather than mean, as in the original, so a narrow line (60 Hz and its
    harmonics) keeps its height when many frequencies share a bin. The lowest
    bins are narrower than the frequency resolution; those before the first
    frequency take the first populated bin's value, where the original drew
    zeros, so the image has no blank strip at its low end.
    """
    out = np.empty(power.shape[:-1] + (edges.size - 1,), dtype=np.float32)
    last = None
    first = None
    for i in range(edges.size - 1):
        inside = (freqs >= edges[i]) & (freqs < edges[i + 1])
        if inside.any():
            last = power[..., inside].max(axis=-1)
            if first is None:
                first = i
                out[..., :i] = last[..., None]
        if last is not None:
            out[..., i] = last
    if first is None:
        out[...] = np.nan
    return out


def _windows(total: int, fs: float, fractions, width: int, margin: int) -> list[int]:
    starts = []
    for fraction in fractions:
        start = int(fraction * total) - width // 2
        starts.append(int(np.clip(start, margin, total - width - margin)))
    return starts


def measure_shank(shank, *, reader_factory=None) -> dict:
    """Everything the noise figures need for one shank, already reduced.

    ``shank`` is a session.SessionShank. Returns per-stage pooled power,
    median spectra, per-channel RMS, a snippet, and the RMS over time.
    ``reader_factory(shank)`` returns an object with ``shape`` and row slicing
    (an mtscomp Reader); tests pass their own.
    """
    from scipy import signal

    from .session import _ap_paths

    if reader_factory is None:
        import mtscomp

        def reader_factory(s):
            cbin, ch = _ap_paths(s)
            r = mtscomp.Reader()
            r.open(str(cbin), str(ch))
            return r

    reader = reader_factory(shank)
    try:
        total, n_stored = reader.shape
        positions = np.load(shank.ks.path / "channel_positions.npy")
        n_ap = positions.shape[0]
        if n_stored < n_ap:
            raise ValueError(f"{shank.label}: the file stores {n_stored} channels "
                             f"but the geometry describes {n_ap}")
        fs = float(shank.ks.fs)
        width = int(round(WINDOW_S * fs))
        margin = int(round(MARGIN_S * fs))

        def read(start: int) -> np.ndarray:
            # The trailing stored channel is the SY sync word, never a signal.
            rows = np.asarray(reader[start - margin:start + width + margin],
                              dtype=np.float32)[:, :n_ap]
            return rows * shank.ks.uv_per_bit

        detail = [processing_stages(read(s), fs, margin=margin)
                  for s in _windows(total, fs, WINDOW_FRACTIONS, width, margin)]
        freqs = None
        pooled, median_psd, rms = [], [], []
        for stage in range(len(STAGES)):
            psd = []
            for window in detail:
                freqs, p = signal.welch(window[stage], fs=fs, axis=0,
                                        nperseg=min(WELCH_SAMPLES, width))
                psd.append(p)
            psd = np.mean(psd, axis=0).T                      # (channels, freqs)
            keep = freqs >= FREQ_MIN_HZ
            db = 10 * np.log10(np.maximum(psd[:, keep], 1e-12))
            edges, centers = log_bins(freqs[keep])
            pooled.append(pool_max(db, freqs[keep], edges))
            median_psd.append(np.median(db, axis=0))
            joined = np.concatenate([w[stage] for w in detail], axis=0)
            rms.append(np.sqrt(np.mean(joined.astype(np.float64) ** 2, axis=0)))
        middle = detail[len(detail) // 2]
        snippet_start = max(0, middle[0].shape[0] // 2 - SNIPPET_SAMPLES // 2)
        snippets = [s[snippet_start:snippet_start + SNIPPET_SAMPLES] for s in middle]

        times, over_time = [], []
        fractions = np.linspace(0.05, 0.95, N_TIME_WINDOWS)
        for start in _windows(total, fs, fractions, width, margin):
            demux = demux_band(read(start), fs, margin=margin)
            over_time.append(float(np.median(np.sqrt(np.mean(
                demux.astype(np.float64) ** 2, axis=0)))))
            times.append((start + width / 2) / fs)
    finally:
        close = getattr(reader, "close", None)
        if close:
            close()

    return {
        "label": shank.label, "short": f"imec{shank.probe} sh{shank.shank}",
        "probe": shank.probe, "shank": shank.shank, "fs": fs,
        "depth_um": positions[:, 1].astype(float),
        "freq_centers": centers, "freqs": freqs[freqs >= FREQ_MIN_HZ],
        "pooled_db": pooled, "median_psd_db": median_psd, "rms_uv": rms,
        "snippets": snippets,
        "window_times_s": [(s + width / 2) / fs for s in
                           _windows(total, fs, WINDOW_FRACTIONS, width, margin)],
        "times_s": np.array(times), "rms_over_time_uv": np.array(over_time),
    }


# ------------------------------------------------------------------ figures
def by_depth(values: np.ndarray, depth: np.ndarray, *, axis: int = 0):
    """Average the channels at each depth, for display; returns (depths, values).

    NP2 has two sites at every depth. Drawn as rows of an image they make
    rows of zero height, which showed as thin dark lines across every panel;
    as a line they zigzag between the two columns. The pair also shares an ADC
    instant, so averaging it keeps the 32-channel structure the demux
    reference is about.
    """
    depths, inverse = np.unique(np.round(depth, 3), return_inverse=True)
    moved = np.moveaxis(np.asarray(values, dtype=np.float64), axis, 0)
    sums = np.zeros((depths.size,) + moved.shape[1:])
    np.add.at(sums, inverse, moved)
    counts = np.bincount(inverse, minlength=depths.size).astype(float)
    means = sums / counts.reshape((-1,) + (1,) * (moved.ndim - 1))
    return depths, np.moveaxis(means, 0, axis)


def _grid(n_rows, n_cols, *, panel_w=1.55, panel_h=2.2, **kw):
    import matplotlib.pyplot as plt

    use_lab_style()
    return plt.subplots(n_rows, n_cols, figsize=(0.9 + panel_w * n_cols,
                                                 0.9 + panel_h * n_rows),
                        squeeze=False, **kw)


def _caption(fig, title: str, subtitle: str) -> float:
    """Title and subtitle placed in inches from the top, whatever the height.

    Returns the fraction of the figure's height they use, for tight_layout. A
    fraction fixed for a tall grid put the subtitle through the title of the
    short summary figure.
    """
    height = fig.get_figheight()
    fig.suptitle(title, x=0.01, y=1 - 0.12 / height, ha="left", va="top",
                 fontsize=11, fontweight="bold")
    if subtitle:
        fig.text(0.01, 1 - 0.42 / height, subtitle, ha="left", va="top",
                 fontsize=8, color="#555555", wrap=True)
    return 1 - 0.95 / height


def _label_grid(ax, row, col, m, *, xlabel):
    despine(ax)
    ax.tick_params(labelsize=7)
    if row == 0:
        ax.set_title(m["short"], fontsize=9)
    if col == 0:
        ax.set_ylabel(f"{STAGES[row]}\n{DEPTH_LABEL}", fontsize=8)
    else:
        # Every shank has the same geometry, so one column of depth labels
        # says it; repeated, they crowded into the neighbouring panel.
        ax.tick_params(labelleft=False)
    if row == len(STAGES) - 1:
        ax.set_xlabel(xlabel, fontsize=8)
    else:
        ax.tick_params(labelbottom=False)


def depth_power_grid(measured: list[dict], *, title="", subtitle=""):
    """Rows: processing stage. Columns: shank. Depth by log frequency."""
    fig, axes = _grid(len(STAGES), len(measured), sharey=True)
    image = None
    for col, m in enumerate(measured):
        f = m["freq_centers"]
        for row in range(len(STAGES)):
            ax = axes[row, col]
            depths, power = by_depth(m["pooled_db"][row], m["depth_um"])
            image = ax.pcolormesh(f, depths, power, shading="nearest",
                                  cmap="viridis", vmin=POWER_LIMITS_DB[0],
                                  vmax=POWER_LIMITS_DB[1], rasterized=True)
            ax.set_xscale("log")
            ax.set_xlim(FREQ_MIN_HZ, f.max())
            _label_grid(ax, row, col, m, xlabel="Frequency (Hz)")
    bar = fig.colorbar(image, ax=axes, shrink=0.6, pad=0.01)
    bar.set_label("Power (dB re 1 µV²/Hz)", fontsize=8)
    _caption(fig, title, subtitle)
    return fig


def snippet_grid(measured: list[dict], *, title="", subtitle=""):
    """The same grid, as 33 ms of voltage from the middle of the recording.

    The raw row has its own color scale, as in the original: it carries the
    LFP, hundreds of microvolts, and on the filtered rows' +-50 uV it was
    solid color.
    """
    fig, axes = _grid(len(STAGES), len(measured), sharey=True)
    raw_limit = float(np.percentile(np.abs(np.concatenate(
        [m["snippets"][0].ravel() for m in measured])), 99))
    images = [None, None]
    for col, m in enumerate(measured):
        ms = np.arange(m["snippets"][0].shape[0]) / m["fs"] * 1000
        for row in range(len(STAGES)):
            ax = axes[row, col]
            limit = raw_limit if row == 0 else SNIPPET_LIMITS_UV
            depths, volts = by_depth(m["snippets"][row], m["depth_um"], axis=1)
            image = ax.pcolormesh(ms, depths, volts.T, shading="nearest",
                                  cmap="RdBu_r", vmin=-limit, vmax=limit,
                                  rasterized=True)
            images[0 if row == 0 else 1] = image
            _label_grid(ax, row, col, m, xlabel="Time (ms)")
    raw_bar = fig.colorbar(images[0], ax=axes[0, :], shrink=0.9, pad=0.01)
    raw_bar.set_label("Raw (µV)", fontsize=8)
    bar = fig.colorbar(images[1], ax=axes[1:, :], shrink=0.6, pad=0.01)
    bar.set_label("Filtered (µV)", fontsize=8)
    _caption(fig, title, subtitle)
    return fig


def _per_shank_axes(measured, **kw):
    probes = sorted({m["probe"] for m in measured})
    shanks = sorted({m["shank"] for m in measured})
    fig, axes = _grid(len(probes), len(shanks), **kw)
    where = {(m["probe"], m["shank"]): axes[probes.index(m["probe"]),
                                            shanks.index(m["shank"])]
             for m in measured}
    used = set(id(a) for a in where.values())
    for ax in axes.ravel():
        if id(ax) not in used:
            ax.set_visible(False)
    return fig, axes, where


def rms_by_channel(measured: list[dict], *, title="", subtitle=""):
    """One panel per shank, depth by RMS, the four stages overlaid.

    RMS on a log axis: raw (with the LFP) is ten times the filtered stages,
    and on a linear axis it pressed them against the left edge.
    """
    fig, axes, where = _per_shank_axes(measured, panel_w=2.3, panel_h=2.8,
                                       sharex=True)
    for m in measured:
        ax = where[(m["probe"], m["shank"])]
        for stage, color, values in zip(STAGES, STAGE_COLORS, m["rms_uv"]):
            depths, rms = by_depth(values, m["depth_um"])
            # A dead channel's zero cannot sit on a log axis; 1 uV is below
            # anything a working channel shows.
            ax.plot(np.maximum(rms, 1.0), depths, color=color, lw=0.8, label=stage)
        ax.set_xscale("log")
        ax.set_title(m["short"], fontsize=9)
        ax.set_xlabel("RMS (µV)", fontsize=8)
        ax.set_ylabel(DEPTH_LABEL, fontsize=8)
        ax.tick_params(labelsize=7)
        despine(ax)
    first = where[(measured[0]["probe"], measured[0]["shank"])]
    color_legend(first, fontsize=7, loc="best")
    top = _caption(fig, title, subtitle)
    fig.tight_layout(rect=(0, 0, 1, top))
    return fig


def median_spectra(measured: list[dict], *, title="", subtitle=""):
    """One panel per shank: the median spectrum across channels, per stage."""
    fig, axes, where = _per_shank_axes(measured, panel_w=2.3, panel_h=2.1,
                                       sharey=True)
    for m in measured:
        ax = where[(m["probe"], m["shank"])]
        for stage, color, values in zip(STAGES, STAGE_COLORS, m["median_psd_db"]):
            ax.semilogx(m["freqs"], values, color=color, lw=0.8, label=stage)
        ax.set_xlim(FREQ_MIN_HZ, m["freqs"].max())
        # The band-pass takes the stopband to -120 dB; drawn to there, every
        # difference that matters shared the top fifth of the panel.
        ax.set_ylim(-60, 45)
        ax.set_title(m["short"], fontsize=9)
        ax.set_xlabel("Frequency (Hz)", fontsize=8)
        ax.set_ylabel("Power (dB re 1 µV²/Hz)", fontsize=8)
        ax.tick_params(labelsize=7)
        despine(ax)
    first = where[(measured[0]["probe"], measured[0]["shank"])]
    color_legend(first, fontsize=7, loc="best")
    top = _caption(fig, title, subtitle)
    fig.tight_layout(rect=(0, 0, 1, top))
    return fig


def rms_over_time(measured: list[dict], *, title="", subtitle=""):
    """One panel per shank: median RMS across channels, after demux CAR, over time."""
    fig, axes, where = _per_shank_axes(measured, panel_w=2.3, panel_h=2.0,
                                       sharey=True)
    for m in measured:
        ax = where[(m["probe"], m["shank"])]
        ax.plot(m["times_s"], m["rms_over_time_uv"], "o-", color=STAGE_COLORS[-1],
                ms=3, lw=0.8)
        ax.set_title(m["short"], fontsize=9)
        ax.set_xlabel("Time in recording (s)", fontsize=8)
        ax.set_ylabel("Median RMS (µV)", fontsize=8)
        ax.tick_params(labelsize=7)
        despine(ax)
    top = _caption(fig, title, subtitle)
    fig.tight_layout(rect=(0, 0, 1, top))
    return fig


def bootstrap_ci(values: np.ndarray, *, n: int = 10_000, seed: int = 0):
    """95% CI of the median, by resampling the windows (as bootci in the original)."""
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    medians = np.median(rng.choice(values, size=(n, values.size)), axis=1)
    return np.percentile(medians, [2.5, 97.5])


def rms_summary(measured: list[dict], *, title="", subtitle=""):
    """Each shank's median RMS over time, with a bootstrap 95% CI."""
    import matplotlib.pyplot as plt

    use_lab_style()
    fig, ax = plt.subplots(figsize=(max(4.0, 0.45 * len(measured) + 1.5), 3.2))
    for i, m in enumerate(measured):
        values = m["rms_over_time_uv"]
        low, high = bootstrap_ci(values)
        ax.plot([i, i], [low, high], color="#444444", lw=1.5)
        ax.plot(i, np.median(values), "o", color=STAGE_COLORS[-1], ms=5)
    ax.set_xticks(range(len(measured)))
    ax.set_xticklabels([m["short"] for m in measured], rotation=60, fontsize=7,
                       ha="right")
    ax.set_xlabel("Shank")
    ax.set_ylabel("Median RMS after demux CAR (µV)")
    despine(ax)
    top = _caption(fig, title, subtitle)
    fig.tight_layout(rect=(0, 0, 1, top))
    return fig
