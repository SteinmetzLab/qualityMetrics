"""Load a Kilosort 4 output directory into arrays with real units.

Two things this fixes relative to the MATLAB version.

Time is seconds. spike_times.npy holds sample indices, and plotting those on an
axis labelled "time (s)" reads 5.6e7 where the answer is 1877. Nothing in this
module exposes a sample index.

Amplitude is microvolts. Getting there needs the gain from the recording's
.meta, which the sorter output does not carry, so amplitude accessors take a
uv_per_bit argument and there is no default. See spike_amplitudes_uv for why
the obvious Kilosort-2-era formula is wrong for Kilosort 4.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np


class KsError(RuntimeError):
    """Raised when a sorter output directory is unusable."""


def _read_tsv(path: Path) -> dict[int, str]:
    """Read a phy cluster_*.tsv into {cluster_id: value}."""
    out: dict[int, str] = {}
    if not path.exists():
        return out
    lines = path.read_text(encoding="utf-8").splitlines()
    for line in lines[1:]:                      # [0] is the header
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0].strip():
            try:
                out[int(parts[0])] = parts[1].strip()
            except ValueError:
                continue
    return out


@dataclass
class KilosortResults:
    """One sorted shank, with units attached.

    Construct with KilosortResults.load(). Arrays are read lazily, so opening a
    directory to ask one question does not pay for the 100 MB of PC features.
    """

    path: Path
    fs: float
    uv_per_bit: float | None = None
    #: Free-form provenance for figure titles: subject, date, probe, shank.
    label: str = ""
    _cache: dict = field(default_factory=dict, repr=False)

    # -- construction ------------------------------------------------------

    @classmethod
    def load(cls, path: str | Path, uv_per_bit: float | None = None,
             label: str = "") -> KilosortResults:
        """Open a sorter_output directory (or its parent, which is also fine)."""
        path = Path(path)
        if not (path / "spike_times.npy").exists():
            nested = path / "sorter_output"
            if (nested / "spike_times.npy").exists():
                path = nested
            else:
                raise KsError(f"{path} holds no spike_times.npy")

        fs = 30000.0
        params = path / "params.py"
        if params.exists():
            for line in params.read_text(encoding="utf-8").splitlines():
                if line.startswith("sample_rate"):
                    fs = float(line.split("=", 1)[1])
        return cls(path=path, fs=fs, uv_per_bit=uv_per_bit, label=label)

    def _npy(self, name: str) -> np.ndarray:
        if name not in self._cache:
            p = self.path / f"{name}.npy"
            if not p.exists():
                raise KsError(f"{self.path.name} has no {name}.npy")
            self._cache[name] = np.load(p)
        return self._cache[name]

    # -- per-spike ---------------------------------------------------------

    @cached_property
    def spike_times_s(self) -> np.ndarray:
        """Spike times in seconds. Never sample indices."""
        return self._npy("spike_times").ravel().astype(np.float64) / self.fs

    @cached_property
    def spike_clusters(self) -> np.ndarray:
        return self._npy("spike_clusters").ravel().astype(np.int64)

    @cached_property
    def spike_templates(self) -> np.ndarray:
        return self._npy("spike_templates").ravel().astype(np.int64)

    @cached_property
    def spike_depths_um(self) -> np.ndarray:
        """Per-spike depth along the shank, from the sorter's localization."""
        return self._npy("spike_positions")[:, 1].astype(np.float64)

    @cached_property
    def spike_x_um(self) -> np.ndarray:
        return self._npy("spike_positions")[:, 0].astype(np.float64)

    @cached_property
    def spike_amplitudes_raw(self) -> np.ndarray:
        """Kilosort's per-spike template scaling, in its own arbitrary units."""
        return self._npy("amplitudes").ravel().astype(np.float64)

    def spike_amplitudes_uv(self, method: str = "template_anchored"
                            ) -> np.ndarray:
        """Per-spike amplitude in microvolts.

        The Kilosort-2-era recipe, amplitudes.npy times the unwhitened template
        peak-to-peak, is wrong for Kilosort 4: its templates are not unit-norm
        (median L2 norm ~11 on our data), so the product double-counts and lands
        near 1300 uV where the true median is tens of uV.

        Two defensible options are offered instead:

        'template_anchored' (default) rescales each unit's spikes so their mean
        equals that unit's unwhitened template peak-to-peak in uV. Within-unit
        variation over time, which is the whole point of a drift map, is
        preserved exactly; the absolute scale is anchored to the same quantity
        the per-unit amplitude metric reports, so the two agree by construction.

        'kilosort' is amplitudes.npy times uv_per_bit, with no template
        involved. Simpler and closer to the sorter's own bookkeeping, but its
        scale is a peak-like quantity rather than peak-to-peak, so it runs about
        1.4x below the per-unit amplitudes on our data.

        Neither is more "correct" in the abstract; they answer slightly
        different questions. The default is the one that makes the figures in
        this package internally consistent.
        """
        if self.uv_per_bit is None:
            raise KsError(
                "amplitudes in microvolts need uv_per_bit, which comes from the "
                "recording's .meta and is not stored in the sorter output. "
                "Pass it to KilosortResults.load().")
        raw = self.spike_amplitudes_raw
        if method == "kilosort":
            return raw * self.uv_per_bit
        if method != "template_anchored":
            raise ValueError(f"unknown method {method!r}")

        out = np.empty_like(raw)
        unit_amp = self.unit_amplitude_uv
        clusters = self.spike_clusters
        for uid in np.unique(clusters):
            sel = clusters == uid
            mean_raw = raw[sel].mean()
            target = unit_amp.get(int(uid))
            if target is None or mean_raw == 0:
                out[sel] = raw[sel] * self.uv_per_bit
            else:
                out[sel] = raw[sel] * (target / mean_raw)
        return out

    # -- templates ---------------------------------------------------------

    @cached_property
    def templates(self) -> np.ndarray:
        """Whitened templates, (n_templates, n_samples, n_channels)."""
        return self._npy("templates")

    @cached_property
    def templates_unwhitened(self) -> np.ndarray:
        """Templates back in data space, which is where amplitudes live."""
        winv = self._npy("whitening_mat_inv")
        return np.einsum("ijk,kl->ijl", self.templates, winv)

    @cached_property
    def channel_positions(self) -> np.ndarray:
        """(n_channels, 2) site coordinates in micrometres, x then depth."""
        return self._npy("channel_positions")

    @cached_property
    def peak_channel(self) -> np.ndarray:
        """Index of each template's largest-amplitude channel."""
        unw = self.templates_unwhitened
        return np.ptp(unw, axis=1).argmax(axis=1)

    # -- per-unit ----------------------------------------------------------

    @cached_property
    def unit_ids(self) -> np.ndarray:
        return np.unique(self.spike_clusters)

    @cached_property
    def n_spikes(self) -> dict[int, int]:
        uid, counts = np.unique(self.spike_clusters, return_counts=True)
        return {int(u): int(c) for u, c in zip(uid, counts, strict=False)}

    @cached_property
    def duration_s(self) -> float:
        """Recording duration, from the last spike unless told otherwise.

        This underestimates by however long the recording ran silent at the end.
        Pass the true duration from the recording when firing rates matter.
        """
        return float(self.spike_times_s.max())

    def firing_rate(self, duration_s: float | None = None) -> dict[int, float]:
        """Mean rate per unit, in spikes per second."""
        dur = duration_s if duration_s is not None else self.duration_s
        return {u: n / dur for u, n in self.n_spikes.items()}

    @cached_property
    def unit_amplitude_uv(self) -> dict[int, float]:
        """Peak-to-peak of each unit's unwhitened template, in microvolts.

        Keyed by cluster id. Clusters that Kilosort merged carry the template of
        their dominant contributor, which is the usual convention and is what
        the per-spike anchoring above assumes.
        """
        if self.uv_per_bit is None:
            raise KsError("unit amplitudes in microvolts need uv_per_bit")
        ptp = np.ptp(self.templates_unwhitened, axis=1).max(axis=1)
        out: dict[int, float] = {}
        for uid in self.unit_ids:
            # A cluster id indexes the template table directly in Kilosort 4
            # output; fall back to the modal template of its spikes if not.
            if uid < len(ptp):
                out[int(uid)] = float(ptp[uid] * self.uv_per_bit)
            else:
                tmpl = self.spike_templates[self.spike_clusters == uid]
                if tmpl.size:
                    modal = np.bincount(tmpl).argmax()
                    out[int(uid)] = float(ptp[modal] * self.uv_per_bit)
        return out

    @cached_property
    def unit_depth_um(self) -> dict[int, float]:
        """Median depth of each unit's spikes."""
        out = {}
        for uid in self.unit_ids:
            out[int(uid)] = float(
                np.median(self.spike_depths_um[self.spike_clusters == uid]))
        return out

    @cached_property
    def trough_to_peak_ms(self) -> dict[int, float]:
        """Trough-to-peak duration of each unit's template, in milliseconds.

        The classic narrow/broad spiking split. Measured on the peak channel of
        the unwhitened template, trough first then the following maximum, so a
        template whose peak precedes its trough does not report a negative
        duration.
        """
        unw = self.templates_unwhitened
        peak = self.peak_channel
        out = {}
        for uid in self.unit_ids:
            idx = uid if uid < len(unw) else None
            if idx is None:
                continue
            wave = unw[idx, :, peak[idx]]
            trough = int(np.argmin(wave))
            after = wave[trough:]
            # A trough in the last few samples leaves nothing to find a peak in,
            # and argmax on that stub returns 0. Reporting 0.0 ms would put the
            # unit at the extreme narrow-spiking end of every figure, which is a
            # confident wrong answer where the honest one is "not measurable".
            samples = int(np.argmax(after)) if after.size >= 3 else 0
            out[int(uid)] = (float("nan") if samples == 0
                             else float(samples / self.fs * 1000.0))
        return out

    # -- labels and metrics ------------------------------------------------

    @cached_property
    def ks_label(self) -> dict[int, str]:
        """Kilosort's own good/mua call."""
        return _read_tsv(self.path / "cluster_KSLabel.tsv")

    @cached_property
    def group(self) -> dict[int, str]:
        """Curator's label from phy, when a human has been through it."""
        return _read_tsv(self.path / "cluster_group.tsv")

    @cached_property
    def quality_metrics(self) -> dict:
        """Our pipeline's quality_metrics.json, if it sits nearby.

        Looked for beside the sorter output and two levels up, which covers the
        SortingManager layout (shank dir holds the json, kilosort4/sorter_output
        holds the arrays).
        """
        for cand in (self.path / "quality_metrics.json",
                     self.path.parent / "quality_metrics.json",
                     self.path.parent.parent / "quality_metrics.json"):
            if cand.exists():
                return json.loads(cand.read_text(encoding="utf-8"))
        return {}

    def metrics_by_unit(self) -> dict[int, dict]:
        """quality_metrics.json reshaped to {unit_id: {...}}."""
        qm = self.quality_metrics
        return {int(u["unit_id"]): u for u in qm.get("units", [])}

    # -- provenance --------------------------------------------------------

    def title(self, suffix: str = "") -> str:
        """A figure title derived from the data, never hardcoded.

        The MATLAB figures all said "probe00a full KS4" regardless of what they
        were plotting, which silently mislabels every figure made from anything
        else.
        """
        base = self.label or self.path.parent.parent.name or self.path.name
        return f"{base} - {suffix}" if suffix else base


def label_from_path(path: str | Path) -> str:
    """Guess subject / date / probe / shank from a SortingManager path.

    Best-effort: used only for figure titles, and an odd layout degrades to the
    directory name rather than to a wrong label.
    """
    parts = [p for p in Path(path).parts if p not in ("sorter_output",)]
    shank = next((p for p in reversed(parts) if "shank" in p.lower()), None)
    date = next((p for p in parts if len(p) == 10 and p.count("-") == 2), None)
    subject = None
    if date and date in parts:
        i = parts.index(date)
        if i > 0:
            subject = parts[i - 1]
    bits = [b for b in (subject, date, shank) if b]
    return " ".join(bits) if bits else Path(path).name
