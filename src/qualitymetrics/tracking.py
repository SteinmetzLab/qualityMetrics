"""Descriptive firing-continuity QC, never a reason to relabel a unit.

Ten-second rates are split at every eligible bin edge. A strong, sustained
rate change is a tracking-loss *candidate*, not evidence that drift caused it.
Only sorted spikes and small sidecars are read; raw traces and ops pickles are
not needed. The caller contains failures so a picture cannot cost a sort.
"""
from __future__ import annotations

import ast
import csv
import json
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Settings:
    bin_s: float = 10.0
    minimum_side_s: float = 600.0
    minimum_side_exposure_s: float = 300.0
    minimum_bin_exposure_s: float = 5.0
    minimum_active_rate: float = 1.0  # spikes/s
    quiet_active_ratio: float = 0.1
    active_bin_fraction: float = 0.2
    persistence_fraction: float = 0.8
    light_padding_s: float = 0.3


SETTINGS = Settings()
CHUNK_SPIKES = 250_000
MAX_BIN_CELLS = 20_000_000
FILES = ("tracking_continuity.png", "tracking_continuity_units.csv",
         "tracking_continuity_bins.npz", "tracking_continuity.json")


def merge_intervals(intervals, duration_s):
    values = np.asarray(intervals, dtype=float).reshape(-1, 2)
    if not np.isfinite(values).all() or np.any(values[:, 1] < values[:, 0]):
        raise ValueError("invalid exclusion intervals")
    merged = []
    for a, b in values[np.argsort(values[:, 0])]:
        a, b = max(0.0, a), min(duration_s, b)
        if b <= a:
            continue
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(b, merged[-1][1])
        else:
            merged.append([a, b])
    return np.asarray(merged, dtype=float).reshape(-1, 2)


def exposure(edges, intervals):
    seconds = np.diff(edges).astype(float)
    for a, b in intervals:
        seconds -= np.maximum(0.0, np.minimum(edges[1:], b) - np.maximum(edges[:-1], a))
    return np.maximum(seconds, 0.0)


def outside_intervals(times, intervals):
    if not len(intervals):
        return np.ones(len(times), dtype=bool)
    j = np.searchsorted(intervals[:, 0], times, side="right") - 1
    return (j < 0) | (times >= intervals[np.maximum(j, 0), 1])


def screen(counts, seconds, edges, settings=SETTINGS):
    """Best qualifying split, otherwise best descriptive split (no p-value).

    Poisson gain is a ranking score, not a test: splits were searched, spikes
    need not be Poisson, and real biology can also produce sustained changes.
    """
    usable = seconds >= settings.minimum_bin_exposure_s
    e = seconds * usable
    c = counts * usable
    cs, es = np.cumsum(c), np.cumsum(e)
    rates = np.divide(c, e, out=np.zeros_like(e), where=e > 0)
    rows = []
    for k in range(1, len(c)):
        split = float(edges[k])
        if min(split, edges[-1] - split) < settings.minimum_side_s:
            continue
        e1, e2 = es[k - 1], es[-1] - es[k - 1]
        if min(e1, e2) < settings.minimum_side_exposure_s:
            continue
        n1, n2 = cs[k - 1], cs[-1] - cs[k - 1]
        r1, r2 = n1 / e1, n2 / e2
        high, low = max(r1, r2), min(r1, r2)
        if high <= 0:
            continue
        left = np.arange(len(c)) < k
        hi = (left if r1 >= r2 else ~left) & usable
        lo = (~left if r1 >= r2 else left) & usable
        high_fraction = np.mean(rates[hi] >= settings.active_bin_fraction * high)
        low_fraction = np.mean(rates[lo] < settings.active_bin_fraction * high)
        overall = (n1 + n2) / (e1 + e2)
        gain = 2 * sum(n * np.log(r / overall) if n else 0.0
                       for n, r in ((n1, r1), (n2, r2)))
        flagged = (high >= settings.minimum_active_rate
                   and low / high <= settings.quiet_active_ratio
                   and high_fraction >= settings.persistence_fraction
                   and low_fraction >= settings.persistence_fraction)
        rows.append(dict(assessed=True, flagged=bool(flagged), split_s=split,
                         before_spikes_per_s=float(r1), after_spikes_per_s=float(r2),
                         low_high_ratio=float(low / high),
                         high_bin_fraction=float(high_fraction),
                         low_bin_fraction=float(low_fraction), poisson_gain=float(gain),
                         change="disappears" if r1 > r2 else "appears"))
    if not rows:
        return dict(assessed=False, flagged=False, split_s=None,
                    before_spikes_per_s=None, after_spikes_per_s=None,
                    low_high_ratio=None, high_bin_fraction=None, low_bin_fraction=None,
                    poisson_gain=None, change="unassessed")
    return max([r for r in rows if r["flagged"]] or rows,
               key=lambda r: r["poisson_gain"])


def _sample_rate(sorter):
    # params.py is user-editable; read the literal, never execute it.
    for node in ast.parse((sorter / "params.py").read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "sample_rate" for t in node.targets):
            return float(ast.literal_eval(node.value))
    raise ValueError("sample_rate is missing from params.py")


def _exclusions(sort_dir, duration, settings):
    intervals, sources = [], []
    path = sort_dir / "light_events.npz"
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            on, off = data["onset_s"].ravel(), data["offset_s"].ravel()
            if on.shape != off.shape or not np.isfinite(on).all() or not np.isfinite(off).all() or np.any(off < on):
                raise ValueError("invalid light event pairs")
            intervals.extend(zip(on - settings.light_padding_s, off + settings.light_padding_s))
        sources.append(path.name)
    path = sort_dir / "mask_intervals.csv"
    if path.exists():
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if not {"start_s", "stop_s"}.issubset(reader.fieldnames or []):
                raise ValueError("mask_intervals.csv needs start_s and stop_s")
            intervals.extend((float(r["start_s"]), float(r["stop_s"])) for r in reader)
        sources.append(path.name)
    return merge_intervals(intervals, duration), sources


def _labels(sorter):
    path = sorter / "cluster_KSLabel.tsv"
    if not path.exists():
        return {}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return {int(r["cluster_id"]): r["KSLabel"] for r in csv.DictReader(stream, delimiter="\t")}


def bin_spikes(sorter, fs, edges, intervals):
    spikes = np.load(sorter / "spike_times.npy", mmap_mode="r", allow_pickle=False).reshape(-1)
    clusters = np.load(sorter / "spike_clusters.npy", mmap_mode="r", allow_pickle=False).reshape(-1)
    if spikes.size != clusters.size or spikes.dtype.kind not in "iu" or clusters.dtype.kind not in "iu":
        raise ValueError("spike samples and cluster IDs must be paired integer arrays")
    # Do not allocate one array indexed by the largest (possibly sparse) ID, or
    # copy/sort all spikes. Both passes use bounded chunks and a small ID set.
    unit_ids = set()
    for a in range(0, len(clusters), CHUNK_SPIKES):
        unit_ids.update(map(int, np.unique(clusters[a:a + CHUNK_SPIKES])))
    ids = np.array(sorted(unit_ids), dtype=clusters.dtype)
    if len(ids) * (len(edges) - 1) > MAX_BIN_CELLS:
        raise ValueError("continuity matrix exceeds diagnostic memory limit")
    counts = np.zeros((len(ids), len(edges) - 1), dtype=np.int64)
    for a in range(0, len(spikes), CHUNK_SPIKES):
        times = spikes[a:a + CHUNK_SPIKES].astype(float) / fs
        if np.any(times < 0) or np.any(times > edges[-1]):
            raise ValueError("spikes outside provenance duration; not truncating or inferring duration")
        row = np.searchsorted(ids, clusters[a:a + CHUNK_SPIKES])
        # Match numpy.histogram: the final bin includes the exact right edge.
        # Kilosort exports can contain a spike at the acquisition endpoint.
        col = np.minimum(np.searchsorted(edges, times, side="right") - 1, len(edges) - 2)
        # Treat an exported endpoint spike as part of the last sample for
        # masking too; clipping an end-reaching mask must not expose it.
        mask_times = np.minimum(times, np.nextafter(edges[-1], -np.inf))
        keep = outside_intervals(mask_times, intervals)
        np.add.at(counts, (row[keep], col[keep]), 1)
    return ids, counts


def normalized_rates(counts, seconds):
    usable = seconds >= SETTINGS.minimum_bin_exposure_s
    rates = np.divide(counts, seconds, out=np.full(counts.shape, np.nan), where=usable)
    scale = np.zeros(len(counts))
    if np.any(usable):
        scale = np.percentile(rates[:, usable], 95, axis=1)
        # Sustained activity can occupy under 5% of a long recording. A zero
        # percentile must not erase that activity from a tracking-loss map.
        scale = np.where(scale > 0, scale, np.max(rates[:, usable], axis=1))
    norm = np.divide(rates, scale[:, None], out=np.zeros_like(rates), where=scale[:, None] > 0)
    norm[:, ~usable] = np.nan
    return norm


def write(sort_dir: Path, out_dir: Path, *, label="", duration_s=None) -> dict:
    """Stage a complete generation before publishing; remove stale outputs on failure."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    error_path = out_dir / "tracking_continuity_error.json"
    try:
        with tempfile.TemporaryDirectory(prefix=".tracking-", dir=out_dir) as temp:
            stage = Path(temp)
            report = _write(sort_dir, stage, label=label, duration_s=duration_s)
            for name in FILES:
                (stage / name).replace(out_dir / name)
        error_path.unlink(missing_ok=True)
        return report
    except Exception as exc:
        # A failed replot must not present an earlier map as the current one.
        failure = {"status": "failed", "error": f"{type(exc).__name__}: {exc}",
                   "stale_files": []}
        error_path.write_text(json.dumps(failure), encoding="utf-8")
        for name in FILES:
            try:
                (out_dir / name).unlink(missing_ok=True)
            except OSError:
                failure["stale_files"].append(name)
        error_path.write_text(json.dumps(failure), encoding="utf-8")
        raise


def _write(sort_dir: Path, out_dir: Path, *, label="", duration_s=None) -> dict:
    """Build one shank's figure/tables; report.build_report catches failures.

    Duration comes from the acquisition provenance, never the final spike.
    Short/empty/all-masked recordings still get a figure and explicit status.
    An optional duration is useful for standalone, externally sorted data.
    """
    sort_dir, out_dir = Path(sort_dir), Path(out_dir)
    from .report import find_sorter_output
    sorter = find_sorter_output(sort_dir)
    if duration_s is None:
        duration_s = json.loads((sort_dir / "provenance.json").read_text(encoding="utf-8"))["duration_s"]
    duration_s, fs = float(duration_s), _sample_rate(sorter)
    if not np.isfinite([duration_s, fs]).all() or min(duration_s, fs) <= 0:
        raise ValueError("duration and sample rate must be finite and positive")
    if duration_s / SETTINGS.bin_s > MAX_BIN_CELLS:
        raise ValueError("recording duration exceeds diagnostic memory limit")
    edges = np.r_[np.arange(0, duration_s, SETTINGS.bin_s), duration_s]
    intervals, exclusion_sources = _exclusions(sort_dir, duration_s, SETTINGS)
    seconds = exposure(edges, intervals)
    ids, counts = bin_spikes(sorter, fs, edges, intervals)
    labels = _labels(sorter)
    rows = [dict(cluster_id=int(cid), ks_label=labels.get(int(cid), "unknown"),
                 **screen(count, seconds, edges)) for cid, count in zip(ids, counts)]
    order = sorted(range(len(rows)), key=lambda i: (
        not rows[i]["flagged"], rows[i]["change"],
        rows[i]["split_s"] if rows[i]["split_s"] is not None else float("inf"), rows[i]["cluster_id"]))
    ordered = [dict(plot_row=k, **rows[i]) for k, i in enumerate(order)]
    norm = normalized_rates(counts, seconds)
    flagged = sum(r["flagged"] for r in rows)
    assessed = sum(r["assessed"] for r in rows)
    summary = dict(label=label, duration_s=duration_s, sample_rate=fs, settings=asdict(SETTINGS),
                   n_units=len(ids), n_candidates=flagged, n_assessed=assessed,
                   n_ks_good=sum(r["ks_label"] == "good" for r in rows),
                   n_candidates_ks_good=sum(r["flagged"] and r["ks_label"] == "good" for r in rows),
                   valid_exposure_s=float(seconds.sum()), exclusion_sources=exclusion_sources,
                   exclusions_s=intervals.tolist(),
                   rate_basis="saved light/mask intervals excluded" if exclusion_sources else "full recording; no exclusion sidecars available",
                   row_order="candidates first, appears/disappears, split time, cluster ID",
                   normalization="unit 95th percentile of usable-bin rates; maximum if percentile is zero; clipped at 1",
                   interpretation="descriptive rate-change candidates, not confirmed tracking failures; no labels changed",
                   source_files={name: {"size_bytes": (sorter / name).stat().st_size,
                                        "mtime_ns": (sorter / name).stat().st_mtime_ns}
                                 for name in ("spike_times.npy", "spike_clusters.npy", "params.py")})
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / FILES[1]).open("w", encoding="utf-8", newline="") as stream:
        columns = ["plot_row", "cluster_id", "ks_label", *screen(np.zeros(len(seconds)), seconds, edges)]
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(ordered)
    np.savez_compressed(out_dir / FILES[2], cluster_ids=ids, counts=counts,
                        edges_s=edges, exposure_s=seconds, plot_order=np.asarray(order, dtype=int))
    (out_dir / FILES[3]).write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    # FigureCanvasAgg avoids global backend changes in a worker with Qt open.
    from matplotlib import rc_context
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    with rc_context({"font.family": "Arial", "font.size": 10}):
        fig = Figure(figsize=(10, 8), layout="constrained")
        FigureCanvasAgg(fig)
        grid = fig.add_gridspec(2, 1, height_ratios=[1, .07])
        ax = fig.add_subplot(grid[0])
        if len(ids):
            cmap = __import__("matplotlib").colormaps["magma"].with_extremes(bad="0.75")
            im = ax.pcolormesh(edges / 60, np.arange(len(ids) + 1),
                               np.clip(norm[order], 0, 1), shading="flat",
                               cmap=cmap, vmin=0, vmax=1, rasterized=True)
            ax.set_ylim(len(ids), 0)
            if 0 < flagged < len(ids):
                ax.axhline(flagged, color="cyan", linewidth=0.8)
            fig.colorbar(im, ax=ax, label="Rate / unit 95th percentile (max if zero; clipped at 1)")
        else:
            ax.set_xlim(0, duration_s / 60)
            ax.text(.5, .5, "No sorted spikes", transform=ax.transAxes, ha="center")
        ax.set(xlabel="Recording time (min)", ylabel="Unit row (not cluster ID)",
               title=f"{label + ': ' if label else ''}firing continuity\n{flagged}/{len(ids)} candidates; {assessed} assessed")
        ax.spines[["top", "right"]].set_visible(False)
        footer = fig.add_subplot(grid[1])
        footer.set_axis_off()
        footer.text(.5, .5, "Candidates first, then direction and split time; gray = insufficient exposure\n"
                    + summary["rate_basis"] + "; candidates are not confirmed tracking failures",
                    transform=footer.transAxes, ha="center", va="center", fontsize=8)
        fig.savefig(out_dir / FILES[0], dpi=180, bbox_inches="tight")
    return {"made": list(FILES), "n_candidates": flagged, "n_units": len(ids)}
