"""Figures that describe a whole session rather than one shank.

Every other module here takes one :class:`~qualitymetrics.ksdata.KilosortResults`
and draws what that shank shows. Two questions are not answerable that way:

* **How noisy is this recording?** Per-site AP band RMS for every shank at once,
  drawn as distributions. A shank with a dozen bad channels and a shank that is
  uniformly noisy have the same median and look nothing alike, so a bar of
  medians hides the two things a noise measurement is read for.
* **What did the whole probe see?** Every unit on one time axis, hue by probe
  and lightness by shank, so drift, a dead shank or a session-wide event is one
  glance rather than sixteen.

The unit figures come in two forms, every sorted unit and only those that passed
quality control, and the second needs nothing computed: ``quality_metrics.json``
is written beside each sort and already carries the per-unit verdicts. Nothing
here recomputes a metric. A figure that disagreed with the pipeline's own
numbers would be worse than no figure.
"""
from __future__ import annotations

import glob
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..ksdata import KilosortResults
from ..style import despine, save, use_lab_style

#: The directory that one sort writes. Its name is ``<run>_imec<P>_shank<S>``
#: with a run such as ``p1_g0``. Sorts made before 2026-10-28 are named
#: ``imec<P>_shank<S>``. Two computers that record one session each number
#: their probes from imec0. The run tells their sorts apart.
SHANK_DIR_RE = re.compile(
    r"^(?:(?P<run>.+_g\d+)_)?imec(?P<probe>\d+)_shank(?P<shank>\d+)$")
#: Matches the run and trigger in a sort's source file name such as
#: ``p1_g0_t0.imec0.sh2.ap.bin``.
_SOURCE_RUN_RE = re.compile(r"^(?P<run>.+_g\d+)_(?P<trigger>t(?:\d+|cat))\.imec\d+")


def file_run(run: str | None) -> str | None:
    """Returns the run as its file names spell it.

    SortingManager's queue names the run of a second trigger file ``p4_t1_g0``.
    Its file names spell the run ``p4_g0``. This function returns ``p4_g0``
    for that run.
    """
    return re.sub(r"_t\d+(?=_g\d+$)", "", run) if run else run


def sort_run(folder_run: str | None, recording: dict) -> str | None:
    """Returns the run of one sort.

    The source file name in ``provenance.json`` names the run without doubt.
    A folder name can be read two ways. ``sess_t3_g0`` is either a run named
    ``sess_t3`` or the third trigger file of run ``sess``. The folder name is
    used only if the provenance names no source.
    """
    named = _SOURCE_RUN_RE.match(_source_name(recording))
    return named.group("run") if named else file_run(folder_run)


def sort_trigger(folder_run: str | None, recording: dict) -> str | None:
    """Returns the trigger file of one sort, such as ``t1``, or None.

    One run can hold several trigger files. Each has its own raw files. The
    source file name in ``provenance.json`` names the trigger. A folder such
    as ``p4_t1_g0`` names it too. A folder such as ``p4_g0`` does not: its
    trigger may be ``t0`` or ``tcat``. Without a source, a folder such as
    ``sess_t3_g0`` reads as trigger ``t3`` of run ``sess``, as in
    ``sort_run``.
    """
    named = _SOURCE_RUN_RE.match(_source_name(recording))
    if named:
        return named.group("trigger")
    in_folder = re.search(r"_(t\d+)_g\d+$", folder_run or "")
    return in_folder.group(1) if in_folder else None


def _source_name(recording: dict) -> str:
    return Path(str(recording.get("source") or "").replace("\\", "/")).name

#: Hue identifies the probe, lightness identifies the shank.
#:
#: Generated in OKLCH and checked all-pairs for separation under protanopia and
#: deuteranopia (Machado, Oliveira and Fernandes 2009, severity 1.0), and within
#: each hue for monotone lightness so shank order is visible in the color.
#: Sixteen colors can never be pairwise distinct, so color is never the only
#: cue here: every figure also labels each probe.
#:
#: Probes 4 to 7 use four more hues on the same lightness ramp. The hues are
#: 110, 200, 290 and 15 degrees. A session on two computers can hold eight
#: probes. No set of four more hues keeps all 32 colors apart under
#: protanopia and deuteranopia. So these four hues are not checked. The labels
#: carry the probe identity. They are drawn in ink at a readable size.
PALETTE: dict[tuple[int, int], str] = {
    (0, 0): "#006eae", (0, 1): "#0090e2", (0, 2): "#50b3ff", (0, 3): "#a2d4ff",
    (1, 0): "#9f5000", (1, 1): "#ce6a00", (1, 2): "#f28e42", (1, 3): "#ffbe91",
    (2, 0): "#007d51", (2, 1): "#00a36a", (2, 2): "#38c789", (2, 3): "#8ae4b5",
    (3, 0): "#9b418b", (3, 1): "#c25faf", (3, 2): "#e484d0", (3, 3): "#fbb2ea",
    (4, 0): "#6d6e01", (4, 1): "#8f9001", (4, 2): "#b1b22b", (4, 3): "#d1d487",
    (5, 0): "#01787d", (5, 1): "#059ca3", (5, 2): "#0ac2ca", (5, 3): "#77e1e7",
    (6, 0): "#6959ae", (6, 1): "#8a76df", (6, 2): "#aa9aff", (6, 3): "#ccc6ff",
    (7, 0): "#a74350", (7, 1): "#d65b6b", (7, 2): "#fa808c", (7, 3): "#ffb8bc",
}
#: Grey, for a layout beyond the 32 colors the palette names. A wrong hue would
#: read as a probe identity that does not exist.
UNKNOWN_COLOR = "#8a8a8a"
INK = "#222222"
MUTED = "#666666"
GRID = "#dddddd"

#: The AP band: above the LFP, below where a spike has no energy left.
SPIKE_BAND_HZ = 300.0
#: Four one-second windows per shank measures a noise floor perfectly well and
#: keeps a 16-shank session to about a minute of reading over a network share.
NOISE_WINDOWS = 4
NOISE_WINDOW_S = 1.0
#: Skipped at each end: the first moments of a recording carry settling
#: transients that are not the noise floor.
NOISE_EDGE_SKIP_S = 20.0
#: Point size of the label on each shank's block in the raster. It is large
#: enough to read without color cues. A color-blind reader may not tell
#: probes 4 to 7 from probes 0 to 3 by color.
SHANK_LABEL_PT = 10.0
#: Height of the raster axes in inches, near enough, after tight_layout on the
#: 8 inch figure below. Only used to convert a label height into rows.
RASTER_AXES_IN = 6.2
#: Width of both session figures: 18 cm, a standard page figure width.
FIGURE_WIDTH_IN = 18 / 2.54
#: Raster height in inches for a session of more than 16 shanks. Each block
#: holds a label. 32 blocks do not fit in 8 inches at SHANK_LABEL_PT.
TALL_RASTER_IN = 10.5


def min_block_rows(max_rows: int, *, label_pt: float = SHANK_LABEL_PT,
                   axes_in: float = RASTER_AXES_IN, spacing: float = 1.6) -> int:
    """Rows a shank's block needs before its label stops hitting its neighbour.

    Computed rather than picked, because the two quantities are in different
    units and guessing converts between them badly. A block's height is a share
    of ``max_rows`` spread over the axes; a label's height is a share of an
    inch. Fixed guesses of 8 and then 20 rows both looked reasonable and both
    collided, the second one by about a pixel.

    This bites only when quality control leaves a shank with almost no units,
    which is exactly when the figure most needs to say which shank that is.
    """
    import math

    return int(math.ceil(max_rows * (label_pt / 72.0) / axes_in * spacing))


def color_for(probe: int, shank: int) -> str:
    return PALETTE.get((probe, shank), UNKNOWN_COLOR)


# --------------------------------------------------------------------------
# Finding the sorts
# --------------------------------------------------------------------------
@dataclass
class SessionShank:
    """One sorted shank of a session, and which shank of which probe it is."""

    probe: int
    shank: int
    directory: Path
    ks: KilosortResults
    duration_s: float
    #: The run such as ``p1_g0``. It comes from the folder name or from
    #: ``provenance.json``. None if neither names it.
    run: str | None = None
    #: The trigger file such as ``t1``. None if neither names it.
    trigger: str | None = None
    #: The probe's position in the session for color and order. It is the
    #: probe number in a session with one run. Each run numbers its probes
    #: from 0. So with several runs it is the index of the (run, probe) pair.
    #: See find_session_shanks.
    slot: int | None = None
    #: True if the session has several runs. Labels then include the run.
    several_runs: bool = False

    @property
    def label(self) -> str:
        return self.directory.name

    @property
    def recording(self) -> str | None:
        """The run, with the trigger file if it is a later one: ``p4_g0``, or
        ``p4_g0 t1``. Each trigger file is its own recording."""
        if self.run and self.trigger not in (None, "t0", "tcat"):
            return f"{self.run} {self.trigger}"
        return self.run

    @property
    def run_token(self) -> str:
        """The recording without its gate: ``p1`` for ``p1_g0``, ``p4 t1``
        for trigger ``t1`` of ``p4_g0``. ``?`` if the run is not known."""
        return re.sub(r"_g\d+(?= |$)", "", self.recording or "") or "?"

    @property
    def probe_name(self) -> str:
        """The probe name such as ``imec0``. It is ``p1 imec0`` if the
        session has several runs."""
        return (f"{self.run_token} imec{self.probe}" if self.several_runs
                else f"imec{self.probe}")

    @property
    def name(self) -> str:
        """The probe and shank such as ``imec0 sh2`` or ``p1 imec0 sh2``."""
        return f"{self.probe_name} sh{self.shank}"

    @property
    def short_label(self) -> str:
        """A label that fits on a crowded axis: ``0.2`` for imec0 shank 2. It
        is ``p1:0.2`` if the session has several runs."""
        prefix = f"{self.run_token}:" if self.several_runs else ""
        return f"{prefix}{self.probe}.{self.shank}"

    @property
    def color(self) -> str:
        return color_for(self.color_probe, self.shank)

    @property
    def color_probe(self) -> int:
        return self.probe if self.slot is None else self.slot

    def quality_metrics(self) -> dict:
        """The per-unit verdicts the pipeline recorded, or ``{}`` if absent.

        Read rather than recomputed. The sliding refractory test, the noise
        cutoff and the amplitude are all decided at sort time and written here;
        deciding them again in a plotting library invites two answers to one
        question.
        """
        path = self.directory / "quality_metrics.json"
        if not path.exists():
            return {}
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return {}

    def passing_units(self, criterion: str = "pass_rescued") -> set[int] | None:
        """Unit ids that passed, or None when this sort recorded no verdicts.

        None and an empty set mean different things: nothing was recorded, as
        against nothing passed. A figure must not present the first as the
        second.
        """
        metrics = self.quality_metrics()
        units = metrics.get("units")
        if not units:
            return None
        if not any(criterion in unit for unit in units):
            return None
        return {int(unit["unit_id"]) for unit in units if unit.get(criterion)}


@dataclass
class SessionSorts:
    """What was found under a session, including what was not usable."""

    session_dir: Path
    shanks: list[SessionShank] = field(default_factory=list)
    #: Directories that look sorted but cannot be read, with the reason.
    incomplete: list[str] = field(default_factory=list)

    @property
    def probes(self) -> list[int]:
        """The color slot of each probe in the session. Two runs' imec0 count
        as two probes."""
        return sorted({s.color_probe for s in self.shanks})

    @property
    def describe(self) -> str:
        return (f"{len(self.shanks)} shanks across "
                f"{len(self.probes)} probe{'s' if len(self.probes) != 1 else ''}")


def find_session_shanks(session_dir: str | Path) -> SessionSorts:
    """Every sorted shank under a session, ordered by probe then shank.

    ``session_dir`` may be the session folder or its ``sorting`` directory.

    Unusable directories are collected rather than raised on, so the figures
    are still drawn for the shanks that are there. They are reported so the
    caller can say so *on the figure*: a session quietly becoming a subset of
    itself is the failure worth avoiding, and a subtitle naming what is missing
    avoids it without throwing the picture away.
    """
    session_dir = Path(session_dir)
    sorting_dir = (session_dir if session_dir.name == "sorting"
                   else session_dir / "sorting")
    found = SessionSorts(session_dir=session_dir)
    if not sorting_dir.is_dir():
        found.incomplete.append(f"no sorting directory under {session_dir}")
        return found

    for entry in sorted(sorting_dir.iterdir()):
        match = SHANK_DIR_RE.match(entry.name)
        if not (match and entry.is_dir()):
            continue
        provenance_path = entry / "provenance.json"
        if not provenance_path.exists():
            found.incomplete.append(f"{entry.name} (no provenance.json)")
            continue
        try:
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            found.incomplete.append(f"{entry.name} (provenance unreadable: {exc})")
            continue
        recording = provenance.get("recording", {})
        gain = provenance.get("gain_to_uv", recording.get("uv_per_bit"))
        if gain is None:
            # Every voltage in the noise figure would otherwise be ADC counts
            # wearing a microvolt label.
            found.incomplete.append(f"{entry.name} (no gain_to_uv or uv_per_bit)")
            continue
        label = entry.name
        run = sort_run(match.group("run"), recording)
        trigger = sort_trigger(match.group("run"), recording)
        try:
            # Imported here rather than at module scope: report.py imports the
            # plots package, so a top-level import would close a cycle.
            from ..report import find_sorter_output

            # A SortingManager sort nests its output under kilosort4/, one
            # level deeper than KilosortResults.load looks by itself.
            ks = KilosortResults.load(find_sorter_output(entry),
                                      uv_per_bit=float(gain), label=label)
        except Exception as exc:  # noqa: BLE001 - reported, never raised
            found.incomplete.append(f"{entry.name} ({type(exc).__name__}: {exc})")
            continue
        found.shanks.append(SessionShank(
            probe=int(match.group("probe")),
            shank=int(match.group("shank")),
            directory=entry,
            ks=ks,
            duration_s=float(provenance.get("duration_s", 0.0)),
            run=run,
            trigger=trigger,
        ))
    # A sort whose provenance names no run belongs to the session's only
    # named run. The session has several runs only if two or more are named.
    # Two trigger files of one run are two recordings, labeled apart.
    runs = {s.recording for s in found.shanks if s.run is not None}
    if len(runs) == 1:
        for shank in found.shanks:
            shank.run = shank.run or next(s.run for s in found.shanks if s.run)
    elif len(runs) > 1:
        pairs = sorted({(s.recording or "", s.probe) for s in found.shanks})
        for shank in found.shanks:
            shank.slot = pairs.index((shank.recording or "", shank.probe))
            shank.several_runs = True
    found.shanks.sort(key=lambda s: (s.color_probe, s.shank))
    return found


# --------------------------------------------------------------------------
# Noise
# --------------------------------------------------------------------------
def _ap_paths(shank: SessionShank) -> tuple[Path, Path]:
    """The archived ``.cbin`` and ``.ch`` this shank was sorted from."""
    raw_root = shank.directory.parent.parent / "raw_ephys_data"
    # Match on the run and trigger if they are known. Two runs of one session
    # can each have an imec0, and one run can have several trigger files. The
    # file names include both even if the folder name spells them
    # differently. The folder p4_t1_g0_imec0 of a second trigger file holds
    # p4_g0_t1 files.
    files = (f"{glob.escape(shank.run)}_{glob.escape(shank.trigger)}."
             if shank.run and shank.trigger
             else f"{glob.escape(shank.run)}_*" if shank.run else "*")
    pattern = f"*imec{shank.probe}/{files}imec{shank.probe}.sh{shank.shank}.ap.cbin"
    matches = sorted(raw_root.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"no AP file matching {pattern} under {raw_root}")
    if len(matches) > 1:
        raise ValueError(f"{len(matches)} AP files match {pattern}; refusing to guess")
    return matches[0], matches[0].with_suffix(".ch")


def channel_rms(shank: SessionShank, *, n_windows: int = NOISE_WINDOWS,
                window_s: float = NOISE_WINDOW_S) -> dict:
    """Per-site AP band RMS for one shank, in microvolts.

    Two numbers per site, because they answer different questions. The
    unreferenced RMS is the conventional figure quoted for Neuropixels and
    includes whatever the whole shank picks up in common. The residual, after
    subtracting the across-site median sample by sample, is what sets the spike
    detection floor. A grounding problem moves the first and leaves the second
    alone, which is how it is told apart from a genuinely noisy probe.

    Median across windows, per site, so one transient in one window cannot set
    a site's level.
    """
    import mtscomp
    from scipy import signal

    cbin, ch = _ap_paths(shank)
    reader = mtscomp.Reader()
    reader.open(str(cbin), str(ch))
    try:
        total, n_stored = reader.shape
        positions = np.load(shank.ks.path / "channel_positions.npy")
        n_ap = positions.shape[0]
        if n_stored < n_ap:
            raise ValueError(
                f"{shank.label}: the AP file stores {n_stored} channels but the "
                f"geometry describes {n_ap}. Refusing to guess which is which.")

        fs = float(shank.ks.fs)
        per_window = int(round(window_s * fs))
        skip = int(round(NOISE_EDGE_SKIP_S * fs))
        start, stop = skip, total - skip - per_window
        if stop <= start:                       # a recording shorter than 40 s
            start, stop = 0, max(total - per_window, 1)
        starts = np.linspace(start, stop, n_windows).astype(np.int64)

        sos = signal.butter(3, SPIKE_BAND_HZ, "hp", fs=fs, output="sos")
        common_rms, residual_rms, raw_rms = [], [], []
        for begin in starts:
            # The trailing stored channel is the SY sync word, a bit field
            # rather than a signal, so it is never part of a noise estimate.
            block = np.asarray(reader[begin:begin + per_window],
                               dtype=np.float32)[:, :n_ap] * shank.ks.uv_per_bit
            block = signal.sosfiltfilt(sos, block, axis=0)
            common = np.median(block, axis=1)
            common_rms.append(float(np.sqrt(np.mean(common ** 2))))
            residual_rms.append(np.std(block - common[:, None], axis=0))
            raw_rms.append(np.std(block, axis=0))
    finally:
        reader.close()

    residual = np.median(np.stack(residual_rms), axis=0)
    raw = np.median(np.stack(raw_rms), axis=0)
    return {
        "label": shank.label,
        "probe": shank.probe,
        "probe_name": shank.probe_name,
        "run": shank.run_token if shank.several_runs else None,
        "slot": shank.color_probe,
        "short": shank.short_label,
        "shank": shank.shank,
        "n_sites": int(n_ap),
        "common_mode_uv": float(np.median(common_rms)),
        "residual_uv": float(np.median(residual)),
        "unreferenced_uv": float(np.median(raw)),
        "per_channel": {
            "depth_um": positions[:, 1].astype(float),
            "unreferenced_uv": raw,
            "residual_uv": residual,
        },
    }


def rms_across_shanks(measurements: list[dict], *, title: str = "",
                      subtitle: str = ""):
    """Per-site AP band RMS for every shank, as distributions.

    One violin per shank rather than a bar of medians: how tightly the sites
    agree, and whether a handful are far out, are both invisible in a median
    and are the two things this figure is read for.
    """
    import matplotlib.pyplot as plt

    use_lab_style()
    rows = [
        ("unreferenced_uv", "Unreferenced"),
        ("residual_uv", "Common-median referenced"),
    ]
    # Each probe has one slot per shank. A shank without a measurement leaves
    # its slot empty. The shanks are not labeled, so a slot's place says which
    # shank it is.
    probes = _probe_groups(measurements)
    per_probe = max(m["shank"] for m in measurements) + 1
    # A gap between probes, in slots, keeps neighboring probe names apart.
    stride = per_probe + 0.7
    positions = np.zeros(len(measurements))
    for place, (start, end, _name) in enumerate(probes):
        for index in range(start, end):
            positions[index] = place * stride + measurements[index]["shank"]
    n_slots = len(probes) * stride - 0.7
    # Only each probe is labeled, so 32 violins fit in the standard width.
    figure, axes = plt.subplots(len(rows), 1, figsize=(FIGURE_WIDTH_IN, 5.5),
                                sharex=True, sharey=True, squeeze=False)
    # One y-scale for both panels, so the drop from referencing is the drop
    # between them. The axis starts at 0, so the drop looks as large as it is.
    # A few far-out sites would leave the panel mostly empty. So the axis ends
    # at the tallest violin once each violin's own far-out sites are set
    # aside. Sites above that are left out, and a count over their violin
    # says how many. The headroom holds the counts.
    panels = [[np.asarray(m["per_channel"][measure], dtype=float)
               for m in measurements] for measure, _r in rows]
    fenced = [_without_far_out(group) for groups in panels for group in groups]
    # A value of 0 at every site would give an axis of zero height.
    high = max(max((float(np.max(f)) for f in fenced if f.size), default=0.0), 1.0)
    shown_panels = [[group[np.isfinite(group) & (group <= high)] for group in groups]
                    for groups in panels]
    hidden = [[int(np.sum(group > high)) for group in groups] for groups in panels]
    lines = [_count_lines(positions, counts, n_slots) for counts in hidden]
    n_lines = max((max(rows_, default=-1) + 1 for rows_ in lines), default=0)
    several_runs = any(m.get("run") for m in measurements)
    headroom = 0.05 + 0.08 * n_lines
    bottom, top, step = _round_axis(0.0, high * (1 + headroom), max_intervals=5)
    span = top - bottom
    count_line = 0.065 * span

    for row, (_measure, reference) in enumerate(rows):
        axis = axes[row][0]
        groups = panels[row]
        # Each shape is drawn from the sites it shows. Far-out sites in the
        # density estimate would flatten the shape into a block.
        shown = shown_panels[row]
        drawn = [index for index, group in enumerate(shown) if group.size]
        parts = axis.violinplot([shown[index] for index in drawn],
                                positions=positions[drawn],
                                widths=0.85, showextrema=False, showmedians=True)
        for body, index in zip(parts["bodies"], drawn, strict=True):
            measurement = measurements[index]
            body.set_facecolor(color_for(measurement.get("slot", measurement["probe"]),
                                         measurement["shank"]))
            body.set_alpha(0.85)
            body.set_edgecolor("none")
        parts["cmedians"].set_color(INK)
        parts["cmedians"].set_linewidth(1.2)
        for position, count, line in zip(positions, hidden[row], lines[row],
                                         strict=True):
            if count:
                unit = "site" if count == 1 else "sites"
                # The math-text arrow sits on the baseline with the digits.
                # Arial's own arrow drops below it.
                axis.text(position, high + 0.01 * span + line * count_line,
                          rf"$\uparrow${count} {unit}", ha="center", va="bottom",
                          fontsize=9, color=INK)

        overall = float(np.median(np.concatenate(groups)))
        axis.axhline(overall, color=INK, linewidth=0.9, linestyle="--", alpha=0.6)
        axis.text(n_slots - 0.4, overall, f" {overall:.1f}", va="center",
                  ha="left", fontsize=8, color=INK, clip_on=False)
        axis.set_title(reference, fontsize=11)
        axis.set_ylabel("AP band RMS (µV)", fontsize=11)
        axis.tick_params(labelsize=10)
        axis.grid(axis="y", color=GRID, linewidth=0.6)
        axis.set_axisbelow(True)
        axis.set_yticks(np.arange(bottom, top + step / 2, step))
        axis.set_ylim(bottom, top)
        despine(axis)

    # The probe names sit under the lower panel. Both panels share its slots.
    # A bracket under each probe's slots carries its name. It runs from the
    # center of the first shank's slot to the center of the last. With several
    # runs, a second bracket under each run's probes carries the run name once.
    # A long run name then never crowds a probe name.
    lower = axes[-1][0]
    lower.set_xticks([])
    lower.set_xlim(-0.6, n_slots - 0.4)
    for place, (start, _end, name) in enumerate(probes):
        if several_runs:
            name = f"imec{measurements[start]['probe']}"
        _bracket(lower, place * stride, place * stride + per_probe - 1, -0.04, name)
    if several_runs:
        runs = [measurements[start].get("run") or "?" for start, _e, _n in probes]
        first = 0
        for last in range(1, len(runs) + 1):
            if last < len(runs) and runs[last] == runs[first]:
                continue
            _bracket(lower, first * stride, (last - 1) * stride + per_probe - 1,
                     -0.19, runs[first])
            first = last
    _caption(figure, title, subtitle)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    return figure


#: A site is far out above its shank's third quartile plus this many
#: interquartile ranges: Tukey's far-out fence.
FAR_OUT_IQR = 3.0


def _without_far_out(group: np.ndarray) -> np.ndarray:
    """Returns the finite values of ``group`` at or below its far-out fence."""
    finite = group[np.isfinite(group)]
    if not finite.size:
        return finite
    q1, q3 = np.percentile(finite, [25, 75])
    return finite[finite <= q3 + FAR_OUT_IQR * (q3 - q1)]


def _count_lines(positions: np.ndarray, counts: list[int],
                 n_slots: float) -> list[int]:
    """Returns the line each hidden-site count goes on. Counts that would
    overlap go on separate lines.

    A label's width is estimated from its length: a character of 9-point Arial
    is about 5 points wide. The axes are about 85% of the figure width.
    """
    slot_pt = FIGURE_WIDTH_IN * 72 * 0.85 / (n_slots + 1)
    lines, ends = [], []
    for position, count in zip(positions, counts, strict=True):
        if not count:
            lines.append(-1)
            continue
        half = 5.0 * len(f"x{count} sites") / 2 / slot_pt
        line = next((index for index, end in enumerate(ends)
                     if end < position - half), len(ends))
        if line == len(ends):
            ends.append(0.0)
        ends[line] = position + half
        lines.append(line)
    return lines


def _bracket(axis, left: float, right: float, level: float, name: str) -> None:
    """Draws a bracket under the axes from ``left`` to ``right``, with ``name``
    under it. ``level`` is a fraction of the axes height, below 0."""
    along = axis.get_xaxis_transform()
    axis.plot([left, left, right, right],
              [level + 0.03, level, level, level + 0.03],
              transform=along, color=INK, linewidth=1.0, clip_on=False)
    # 11 points: eight names fit across 18 cm without touching.
    axis.text((left + right) / 2, level - 0.015, name, transform=along,
              ha="center", va="top", fontsize=11, color=INK)


def _round_axis(low: float, high: float, *,
                max_intervals: int) -> tuple[float, float, float]:
    """Returns ``(bottom, top, step)`` for an axis that covers ``low`` to
    ``high``. The step is the smallest of 1, 2, 2.5 or 5 times a power of 10
    that needs at most ``max_intervals`` intervals. Both ends fall on a tick.
    """
    scale = 10.0 ** math.floor(math.log10(max(high - low, 1e-9) / max_intervals))
    for factor in (1, 2, 2.5, 5, 10, 20):
        step = factor * scale
        bottom = step * math.floor(low / step)
        top = step * math.ceil(high / step)
        if round((top - bottom) / step) <= max_intervals:
            return bottom, top, step
    return bottom, top, step


def _probe_groups(measurements: list[dict]) -> list[tuple[int, int, str]]:
    """Returns ``(start, end, probe name)`` for each run of one probe's
    measurements. Measurements arrive grouped by probe."""
    slots = [m.get("slot", m["probe"]) for m in measurements]
    found, start = [], 0
    for end in range(1, len(slots) + 1):
        if end < len(slots) and slots[end] == slots[start]:
            continue
        first = measurements[start]
        found.append((start, end, first.get("probe_name") or f"imec{first['probe']}"))
        start = end
    return found


# --------------------------------------------------------------------------
# The raster
# --------------------------------------------------------------------------
def unit_rates(shank: SessionShank, edges: np.ndarray,
               keep: set[int] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Binned spike counts per unit, ordered by depth from the shank tip.

    Depth order puts neighbouring rows near each other in the brain rather than
    near each other in cluster numbering. ``keep`` selects unit ids.
    """
    times = shank.ks.spike_times_s
    clusters = shank.ks.spike_clusters
    if times.size == 0:
        return np.zeros((0, len(edges) - 1), dtype=np.float32), np.zeros(0)

    unit_ids = np.unique(clusters)
    if keep is not None:
        unit_ids = np.array([u for u in unit_ids if int(u) in keep], dtype=unit_ids.dtype)
    if unit_ids.size == 0:
        return np.zeros((0, len(edges) - 1), dtype=np.float32), np.zeros(0)

    try:
        spike_depths = shank.ks.spike_depths_um
    except Exception:  # noqa: BLE001 - depth is an ordering nicety, not the data
        spike_depths = None

    counts = np.zeros((unit_ids.size, len(edges) - 1), dtype=np.float32)
    depths = np.zeros(unit_ids.size, dtype=float)
    for row, unit in enumerate(unit_ids):
        mine = clusters == unit
        counts[row], _ = np.histogram(times[mine], bins=edges)
        if spike_depths is not None and mine.any():
            depths[row] = float(np.median(spike_depths[mine]))

    if spike_depths is not None:
        order = np.argsort(depths, kind="stable")
        counts, depths = counts[order], depths[order]
    return counts, depths


def _fit_rows(counts: np.ndarray, target: int) -> np.ndarray:
    """Make a block exactly ``target`` rows tall, both directions.

    Shrinking averages neighbouring units: handing ``imshow`` more rows than it
    has pixels lets it drop rows silently, which on a five thousand unit
    session means units missing from the picture with nothing saying so.

    Growing repeats rows, which is what gives a shank with one surviving unit
    enough height to carry its own label. A row is therefore not a unit in
    either direction, which is why the axis says "units, by depth" rather
    than labelling row counts.
    """
    n = counts.shape[0]
    if n == target or n == 0:
        return counts
    if n > target:
        edges = np.linspace(0, n, target + 1).astype(int)
        return np.stack([counts[a:b].mean(axis=0) if b > a else counts[a]
                         for a, b in zip(edges[:-1], edges[1:], strict=False)])
    return counts[np.linspace(0, n - 1, target).round().astype(int)]


def session_raster(shanks: list[SessionShank], *, duration_s: float | None = None,
                   criterion: str | None = None, bin_s: float = 1.0,
                   max_rows: int = 1400, title: str = "", subtitle: str = ""):
    """Every unit of every shank on one time axis.

    ``criterion`` names a column of ``quality_metrics.json`` (``pass_strict``
    or ``pass_rescued``) to filter on. None draws every sorted unit, which is
    the right input to "is this recording worth curating".

    Each unit is a row shaded by its rate normalised to its own maximum, so a
    quiet unit is as visible as a loud one and the figure shows *when* things
    fired rather than how hard.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgb

    use_lab_style()
    duration_s = duration_s or max((s.duration_s for s in shanks), default=0.0)
    if duration_s <= 0:
        raise ValueError("the session has no recorded duration to plot against")
    edges = np.linspace(0.0, duration_s, max(int(round(duration_s / bin_s)), 1) + 1)

    blocks, unverdicted = [], []
    for shank in shanks:
        keep = None
        if criterion is not None:
            keep = shank.passing_units(criterion)
            if keep is None:
                unverdicted.append(shank.label)
                continue
        counts, _ = unit_rates(shank, edges, keep=keep)
        if counts.shape[0]:
            blocks.append((shank, counts))
    if not blocks:
        raise ValueError("no units to draw for this session")

    total_units = sum(counts.shape[0] for _, counts in blocks)
    height_in = 8.0 if len(blocks) <= 16 else TALL_RASTER_IN
    floor_rows = min_block_rows(max_rows,
                                axes_in=RASTER_AXES_IN + height_in - 8.0)
    image_rows, boundaries, row_colors = [], [], []
    for shank, counts in blocks:
        # Proportional to unit count, so a shank that produced few units looks
        # like it, but floored so every shank still has room for its label.
        share = max(int(round(max_rows * counts.shape[0] / total_units)),
                    floor_rows)
        reduced = _fit_rows(counts, share)
        peak = reduced.max(axis=1, keepdims=True)
        peak[peak == 0] = 1.0
        image_rows.append(reduced / peak)
        row_colors.extend([shank.color] * reduced.shape[0])
        boundaries.append((shank, reduced.shape[0], counts.shape[0]))

    intensity = np.vstack(image_rows)
    rgb = np.ones((*intensity.shape, 3), dtype=float)
    for row in range(intensity.shape[0]):
        color = np.asarray(to_rgb(row_colors[row]))
        # Composite onto white: a silent bin is the page, a busy one is the
        # shank's color at full strength.
        rgb[row] = 1.0 - intensity[row][:, None] * (1.0 - color)[None, :]

    figure, axis = plt.subplots(figsize=(FIGURE_WIDTH_IN, height_in))
    axis.imshow(rgb, aspect="auto", interpolation="nearest", origin="upper",
                extent=(0.0, duration_s / 60.0, intensity.shape[0], 0))

    position = 0
    for index, (shank, height, _n) in enumerate(boundaries):
        # On the block itself, so no reader has to match a color.
        axis.text(0.004, position + height / 2, shank.name,
                  transform=axis.get_yaxis_transform(), ha="left",
                  va="center", fontsize=SHANK_LABEL_PT, color=INK,
                  bbox={"boxstyle": "round,pad=0.15", "facecolor": "white",
                        "edgecolor": "none", "alpha": 0.85})
        position += height
        if position < intensity.shape[0]:
            following = (boundaries[index + 1][0]
                         if index + 1 < len(boundaries) else None)
            # A heavier rule between probes, where the shade ramp restarts.
            if following is not None and following.color_probe != shank.color_probe:
                axis.axhline(position, color=INK, linewidth=1.0, alpha=0.5)
            else:
                axis.axhline(position, color="white", linewidth=0.8)

    axis.tick_params(labelsize=10)
    axis.set_xlabel("Time in the recording (minutes)", fontsize=11)
    axis.set_ylabel(f"{total_units} units, by depth", fontsize=11)
    axis.set_yticks([])
    axis.spines["left"].set_visible(False)
    despine(axis)

    if unverdicted:
        subtitle += (f"   No recorded verdicts, so left out: "
                     f"{', '.join(unverdicted)}")
    _caption(figure, title, subtitle)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    return figure


def _caption(figure, title: str, subtitle: str) -> None:
    if title:
        figure.text(0.006, 0.992, title, fontsize=13, fontweight="bold", va="top")
    if subtitle:
        figure.text(0.006, 0.960, subtitle, fontsize=9, color=MUTED, va="top")


def noise_table(measurements: list[dict], path: str | Path) -> str:
    """The medians behind the violins, so two sessions compare as numbers."""
    lines = ["label,probe,shank,n_sites,unreferenced_uv,residual_uv,common_mode_uv"]
    for m in measurements:
        lines.append(f"{m['label']},{m['probe']},{m['shank']},{m['n_sites']},"
                     f"{m['unreferenced_uv']:.3f},{m['residual_uv']:.3f},"
                     f"{m['common_mode_uv']:.3f}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


__all__ = ["SessionShank", "SessionSorts", "find_session_shanks", "channel_rms",
           "rms_across_shanks", "session_raster", "unit_rates", "noise_table",
           "color_for", "PALETTE", "save"]
