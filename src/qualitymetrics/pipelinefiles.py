"""Files SortingManager's pipeline writes beside a sort, read without matplotlib.

``destriped_windows.npz`` (``ks4driver.save_destriped_windows``) holds three
windows of exactly what Kilosort received: destriped, with light artifacts
removed when the session had them, before motion interpolation. Each window is
a core second plus a margin on either side, kept for the figures' filters.

``light_events.npz`` (``sortingmanager/lightevents.py``) holds the laser pulse
edges, in seconds on the probe's own clock, when the session had light
stimulation. Whether the sort removed their artifacts is taken from the sort's
``provenance.json`` (``destripe.light_artifacts.applied``), not inferred from
the file being there.

Both the noise figures and ``wall_heatmap`` read these, through this module, so
the one place that turns a recording sample into a row of the saved traces is
:meth:`DestripedWindows.samples_uv`.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

#: Written beside the sort by SortingManager's ks4driver.save_destriped_windows.
PIPELINE_FILE = "destriped_windows.npz"
#: Written beside the sort by SortingManager when the session had light pulses.
LIGHT_EVENTS_FILE = "light_events.npz"


@dataclass(frozen=True)
class DestripedWindows:
    """The pipeline's kept windows of one shank.

    ``starts[i]`` is the recording sample where window i's core second begins.
    The saved traces begin one ``margin`` earlier: row 0 of ``traces[i]`` is
    recording sample ``starts[i] - margin``. Reading row ``s - starts[i]``
    for sample ``s`` draws the data ``margin`` samples (100 ms) late, a mistake
    that has been made once already; use :meth:`samples_uv`.
    """

    traces: np.ndarray            # int16 counts, (windows, samples, channels)
    uv_per_count: float
    starts: np.ndarray            # sample of each core second's first sample
    margin: int                   # samples saved before and after each core
    fs: float
    #: Recording channel index of each column, from ``channel_ids``. Some
    #: channels (above the brain edge, dead or listed) are absent. None when
    #: the file has no ``channel_ids`` or names channels another way (a
    #: probe-wide sort names them ``AP<n>``, the acquired channel; see
    #: :attr:`channel_names`).
    channels: np.ndarray | None
    fractions: np.ndarray | None = None
    adc_group_median: bool | None = None
    description: str = ""
    path: Path | None = None
    #: ``channel_ids`` as written, or None when the file has none.
    channel_names: tuple[str, ...] | None = None
    #: From the sort's provenance.json (destripe.light_artifacts.applied):
    #: True or False when it says, None when there is no record.
    light_artifacts_removed: bool | None = None

    @property
    def n_windows(self) -> int:
        return int(self.traces.shape[0])

    @property
    def width(self) -> int:
        """Samples in each window's core, without the margins."""
        return int(self.traces.shape[1]) - 2 * self.margin

    def first_sample(self, window: int) -> int:
        """The recording sample held in row 0 of ``traces[window]``."""
        return int(self.starts[window]) - self.margin

    def core(self, window: int) -> tuple[int, int]:
        """The half-open span of recording samples in a window's core second."""
        start = int(self.starts[window])
        return start, start + self.width

    def window_covering(self, i0: int, i1: int) -> int | None:
        """The window whose core second holds samples [i0, i1), or None.

        The core rather than the whole saved span, deliberately: the margin
        exists for filters and is not promised to be anything else.
        """
        for window in range(self.n_windows):
            lo, hi = self.core(window)
            if lo <= i0 and i1 <= hi:
                return window
        return None

    def samples_uv(self, i0: int, i1: int) -> np.ndarray | None:
        """Recording samples [i0, i1) in microvolts, (samples, columns), or None
        when no window's core covers them. Columns are in file order; map them
        to recording channels through :attr:`channels`."""
        if i1 <= i0:
            raise ValueError(f"empty sample span {i0}-{i1}")
        window = self.window_covering(i0, i1)
        if window is None:
            return None
        # The margin correction: rows start one margin before the core.
        row = i0 - self.first_sample(window)
        block = self.traces[window, row:row + (i1 - i0)]
        return block.astype(np.float32) * np.float32(self.uv_per_count)


def load_destriped_windows(where: str | Path) -> DestripedWindows | None:
    """The kept windows in a sort folder (or at a file path), or None if absent.

    Raises ValueError for a file whose parts disagree with one another, rather
    than returning something that would be drawn misaligned.
    """
    path = Path(where)
    if path.is_dir():
        path = path / PIPELINE_FILE
    if not path.exists():
        return None
    with np.load(path) as kept:
        files = set(kept.files)
        traces = kept["traces"]
        # Files from before uv_per_count was written were already microvolts.
        scale = float(kept["uv_per_count"]) if "uv_per_count" in files else 1.0
        starts = np.asarray(kept["starts"], dtype=np.int64).ravel()
        margin, fs = int(kept["margin"]), float(kept["fs"])
        ids = kept["channel_ids"] if "channel_ids" in files else None
        fractions = (np.asarray(kept["fractions"], dtype=float)
                     if "fractions" in files else None)
        group_median = (bool(kept["adc_group_median"])
                        if "adc_group_median" in files else None)
        description = str(kept["description"]) if "description" in files else ""
    if traces.ndim != 3:
        raise ValueError(f"{path.name}: traces are {traces.ndim}-D, expected "
                         "(windows, samples, channels)")
    if starts.size != traces.shape[0]:
        raise ValueError(f"{path.name}: {starts.size} starts for "
                         f"{traces.shape[0]} windows")
    if margin < 0 or traces.shape[1] <= 2 * margin:
        raise ValueError(f"{path.name}: windows of {traces.shape[1]} samples "
                         f"cannot hold a margin of {margin} on each side")
    channels = names = None
    if ids is not None:
        names = tuple(str(c) for c in np.ravel(ids))
        if len(names) != traces.shape[2]:
            raise ValueError(f"{path.name}: {len(names)} channel_ids for "
                             f"{traces.shape[2]} columns")
        if len(set(names)) != len(names):
            raise ValueError(f"{path.name}: a channel is listed twice")
        # A single-shank sort names each channel by its column in the shank's
        # file ("0" to "383"), so the name is the recording's index. Any other
        # naming (a probe-wide sort's "AP<n>") is kept as names, not refused:
        # the noise figures need only the count, and wall_heatmap maps "AP<n>"
        # through the recording's own meta.
        if all(n.isdigit() for n in names):
            channels = np.array([int(n) for n in names], dtype=np.int64)
    return DestripedWindows(traces=traces, uv_per_count=scale, starts=starts,
                            margin=margin, fs=fs, channels=channels,
                            fractions=fractions, adc_group_median=group_median,
                            description=description, path=path,
                            channel_names=names,
                            light_artifacts_removed=_light_removed(path.parent))


def _light_removed(sort_dir: Path) -> bool | None:
    """Whether the sort removed light artifacts, by its own provenance record.

    The record rather than the presence of light_events.npz, which says only
    that there were pulses. None when the sort recorded nothing either way.
    """
    import json

    try:
        provenance = json.loads((sort_dir / "provenance.json")
                                .read_text(encoding="utf-8"))
        applied = provenance["destripe"]["light_artifacts"]["applied"]
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return bool(applied)


@dataclass(frozen=True)
class LightEvents:
    """Laser pulse edges, in seconds on the recording's clock."""

    onset_s: np.ndarray
    offset_s: np.ndarray

    def edges_between(self, t0: float, t1: float) -> list[tuple[str, float]]:
        """("onset" or "offset", time) for every edge in [t0, t1), in time order."""
        edges = [(kind, float(t))
                 for kind, times in (("onset", self.onset_s),
                                     ("offset", self.offset_s))
                 for t in times if t0 <= t < t1]
        return sorted(edges, key=lambda e: e[1])


def load_light_events(where: str | Path) -> LightEvents | None:
    """The light pulse edges in a sort folder (or at a file path), or None."""
    path = Path(where)
    if path.is_dir():
        path = path / LIGHT_EVENTS_FILE
    if not path.exists():
        return None
    with np.load(path) as z:
        return LightEvents(onset_s=np.asarray(z["onset_s"], dtype=float).ravel(),
                           offset_s=np.asarray(z["offset_s"], dtype=float).ravel())
