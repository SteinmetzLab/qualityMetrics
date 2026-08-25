"""Random access to an archived recording, compressed or not.

The MATLAB version of this repository could only read a SpikeInterface scratch
file (traces_cached_seg0.raw) that a sorting pipeline deletes when it finishes,
so half its figures were unrunnable on any recording that had already been
archived. Reading the archived .cbin instead removes that dependency: the
compressed file is the copy that is kept forever, and mtscomp decompresses on
the fly, so a 17 GB shank costs a tenth of a second to open and a fraction of a
second per window.

Everything here returns microvolts. Raw integers are meaningless across probes,
and not applying the gain is why the MATLAB amplitude axes read in the tens of
thousands where the real values were tens.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


class RawError(RuntimeError):
    """Raised when a recording cannot be opened or described."""


# --------------------------------------------------------------------------
# SpikeGLX metadata
# --------------------------------------------------------------------------

def parse_meta(path: str | Path) -> dict[str, str]:
    """Read a SpikeGLX .meta into a plain dict of strings."""
    path = Path(path)
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def paren_entries(value: str) -> list[str]:
    """Split a SpikeGLX (a)(b)(c) list into its parts."""
    return re.findall(r"\(([^)]*)\)", value)


@dataclass
class Geometry:
    """Site coordinates aligned 1:1 with the data channels."""

    x: np.ndarray          # micrometres, across the shank
    y: np.ndarray          # micrometres, depth along the shank
    shank: np.ndarray      # shank index per channel
    used: np.ndarray       # bool, site connected

    @property
    def n_channels(self) -> int:
        return int(self.x.size)

    def pitch(self) -> float:
        """Vertical spacing between adjacent rows, in micrometres."""
        uniq = np.unique(np.round(self.y, 3))
        return float(np.median(np.diff(uniq))) if uniq.size > 1 else float("nan")


def geometry_from_meta(meta: dict[str, str]) -> Geometry:
    """Geometry from ~snsGeomMap, whose entries are shank:x:z:used."""
    if "~snsGeomMap" not in meta:
        raise RawError(
            "meta has no ~snsGeomMap, so channel positions are unknown. Pass "
            "geometry explicitly, or take it from the sorter's "
            "channel_positions.npy."
        )
    entries = paren_entries(meta["~snsGeomMap"])[1:]  # [0] is the header
    rows = [e.split(":") for e in entries]
    return Geometry(
        x=np.array([float(r[1]) for r in rows]),
        y=np.array([float(r[2]) for r in rows]),
        shank=np.array([int(r[0]) for r in rows]),
        used=np.array([bool(int(r[3])) for r in rows]),
    )


def uv_per_bit(meta: dict[str, str], band: str = "ap") -> float:
    """Microvolts per integer step, from the range, the ADC depth and the gain.

    Every amplitude in every figure depends on this one number, so this raises
    rather than quietly falling back to 1.0. A silent fallback is exactly how
    the MATLAB figures came to label raw integers as microvolts.
    """
    try:
        ai_range = float(meta["imAiRangeMax"])
        max_int = int(meta["imMaxInt"])
        # The band's own gain key first, then the AP gain. An LF band derived
        # from wideband data (Neuropixels 2.0 has no LF stream from the probe)
        # is written with only imChan0apGain, because the AP gain is the one
        # that actually applied when the samples were digitised. A reader that
        # insists on imChan0lfGain cannot open those files at all.
        keys = (["imChan0apGain"] if band == "ap"
                else ["imChan0lfGain", "imChan0apGain"])
        gain = next((float(meta[k]) for k in keys if k in meta), None)
        if gain is None:
            raise KeyError(" or ".join(keys))
    except (KeyError, ValueError) as exc:
        raise RawError(f"cannot derive uV/bit from this meta: {exc}") from exc
    if gain == 0 or max_int == 0:
        raise RawError("meta reports a zero gain or zero max integer")
    return ai_range / max_int / gain * 1e6


# --------------------------------------------------------------------------
# the recording
# --------------------------------------------------------------------------

@dataclass
class RawRecording:
    """A view of one shank's traces that never loads the whole file.

    Open with RawRecording.open(). Use as a context manager, or call close();
    mtscomp holds a file handle and a decompression cache.
    """

    path: Path
    fs: float
    n_samples: int
    n_channels: int
    uv_per_bit: float
    geometry: Geometry | None
    meta: dict[str, str]
    band: str
    _reader: object = None
    _memmap: np.ndarray | None = None

    # -- construction ------------------------------------------------------

    @classmethod
    def open(cls, path: str | Path, meta_path: str | Path | None = None,
             band: str | None = None) -> RawRecording:
        """Open a .cbin (compressed) or .bin/.raw (flat) recording.

        The sibling .meta is found automatically when it sits next to the data
        under the usual SpikeGLX name.
        """
        path = Path(path)
        if not path.exists():
            raise RawError(f"no such recording: {path}")
        if band is None:
            band = "lf" if ".lf." in path.name else "ap"

        meta_path = Path(meta_path) if meta_path else path.with_suffix(".meta")
        if not meta_path.exists():
            raise RawError(
                f"no meta beside {path.name} (looked for {meta_path.name})")
        meta = parse_meta(meta_path)

        n_channels = int(meta["nSavedChans"])
        fs = float(meta["imSampRate"])
        scale = uv_per_bit(meta, band)
        try:
            geom = geometry_from_meta(meta)
        except RawError:
            geom = None

        reader = memmap = None
        if path.suffix == ".cbin":
            from mtscomp import decompress
            ch = path.with_suffix(".ch")
            if not ch.exists():
                raise RawError(f"{path.name} has no .ch index beside it; a "
                               ".cbin cannot be read without one")
            reader = decompress(path, ch)
            n_samples = int(reader.shape[0])
            if int(reader.shape[1]) != n_channels:
                raise RawError(
                    f"meta says {n_channels} channels, the compressed file has "
                    f"{reader.shape[1]}. One of them is not this recording.")
        else:
            itemsize = np.dtype(np.int16).itemsize
            n_samples = path.stat().st_size // (itemsize * n_channels)
            memmap = np.memmap(path, dtype=np.int16, mode="r",
                               shape=(n_samples, n_channels))

        # Geometry describes the data channels; a sync channel has no site.
        if geom is not None and geom.n_channels > n_channels:
            geom = Geometry(geom.x[:n_channels], geom.y[:n_channels],
                            geom.shank[:n_channels], geom.used[:n_channels])

        return cls(path=path, fs=fs, n_samples=n_samples, n_channels=n_channels,
                   uv_per_bit=scale, geometry=geom, meta=meta, band=band,
                   _reader=reader, _memmap=memmap)

    # -- basic facts -------------------------------------------------------

    @property
    def duration_s(self) -> float:
        return self.n_samples / self.fs

    @property
    def n_data_channels(self) -> int:
        """Channels carrying voltage, i.e. excluding the sync channel."""
        counts = self.meta.get("snsApLfSy")
        if counts:
            ap, lf, _sy = (int(v) for v in counts.split(","))
            return ap or lf or self.n_channels
        return self.n_channels

    def describe(self) -> str:
        g = self.geometry
        span = f"{g.y.min():.0f}-{g.y.max():.0f} um" if g else "unknown"
        return (f"{self.path.name}: {self.n_data_channels} channels, "
                f"{self.duration_s / 60:.1f} min at {self.fs / 1000:.1f} kHz, "
                f"{self.uv_per_bit:.4f} uV/bit, depth span {span}")

    # -- data --------------------------------------------------------------

    def get_traces(self, t0: float, t1: float, channels=None,
                   unit: str = "uV") -> np.ndarray:
        """Pull the half-open window [t0, t1) seconds as (samples, channels).

        Returns microvolts by default. unit='int16' returns the stored integers,
        which is what a checksum or a saturation test wants.
        """
        if t1 <= t0:
            raise ValueError(f"empty window: t0={t0} t1={t1}")
        i0 = max(0, int(round(t0 * self.fs)))
        i1 = min(self.n_samples, int(round(t1 * self.fs)))
        if i1 <= i0:
            raise ValueError(f"window {t0}-{t1} s lies outside this "
                             f"{self.duration_s:.1f} s recording")
        if channels is None:
            channels = np.arange(self.n_data_channels)
        channels = np.asarray(channels, dtype=int)

        source = self._reader if self._reader is not None else self._memmap
        block = np.asarray(source[i0:i1, :])[:, channels]
        if unit == "int16":
            return block
        if unit != "uV":
            raise ValueError(f"unit must be 'uV' or 'int16', got {unit!r}")
        return block.astype(np.float32) * self.uv_per_bit

    def sample_windows(self, n: int, win_s: float, unit: str = "uV",
                       channels=None, skip_edges_s: float = 1.0):
        """Yield (t0, traces) for n windows spread across the recording.

        Sampling several short windows rather than one is how every per-channel
        statistic in this package avoids being a statement about one arbitrary
        second of the session.
        """
        usable = self.duration_s - 2 * skip_edges_s - win_s
        if usable <= 0:
            starts = [0.0]
        else:
            starts = list(np.linspace(skip_edges_s, skip_edges_s + usable, n))
        for t0 in starts:
            yield t0, self.get_traces(t0, t0 + win_s, channels=channels,
                                      unit=unit)

    # -- lifecycle ---------------------------------------------------------

    def close(self) -> None:
        if self._reader is not None:
            try:
                self._reader.close()
            except Exception:  # noqa: BLE001 - closing must not raise
                pass
            self._reader = None
        self._memmap = None

    def __enter__(self) -> RawRecording:
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def find_band(session_dir: str | Path, probe: str, shank: int,
              band: str = "ap") -> Path:
    """Locate the archived recording for one shank of one probe.

    Written against the SortingManager archive layout, but falls back to a glob
    so a differently named archive still resolves.
    """
    session_dir = Path(session_dir)
    roots = [session_dir / "raw_ephys_data", session_dir]
    patterns = [f"*{probe}*sh{shank}.{band}.cbin",
                f"*{probe}*sh{shank}.{band}.bin",
                f"*sh{shank}.{band}.cbin", f"*sh{shank}.{band}.bin"]
    for root in roots:
        if not root.exists():
            continue
        for pattern in patterns:
            hits = sorted(root.glob(f"**/{pattern}"))
            if hits:
                return hits[0]
    raise RawError(
        f"no {band} recording for {probe} shank {shank} under {session_dir}")
