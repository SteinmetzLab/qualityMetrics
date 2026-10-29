"""Figures that describe a whole session rather than one shank.

Most of what can go wrong here is discovery and partial sessions, so that is
what most of these are about. The drawing itself is checked for "it produced a
figure and did not raise", because a test cannot tell whether a picture is any
good. A person has to look, and doing that is what caught the three defects
this module had before it worked: an axis label painted through the shank
labels, a label floor that was never applied, and a floor whose units were
wrong.
"""
from __future__ import annotations

import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from qualitymetrics.plots import session as S  # noqa: E402
from qualitymetrics.report import build_session_report  # noqa: E402

FS = 30000.0
N_SAMPLES = 60_000


def write_shank(session_dir, probe: int, shank: int, *, n_units: int = 6,
                gain: float = 3.0, provenance: bool = True,
                sorter_output: bool = True, verdicts: bool = True,
                passing: int = 2, run: str | None = None,
                source_run: str | None = None, trigger: str = "t0") -> None:
    """Writes one shank's sort output in SortingManager's layout.

    ``run`` such as ``p1_g0`` names the folder as SortingManager has since
    2026-10-28. Without it the folder has the older name. ``source_run`` is
    the run that ``provenance.json`` records. Every sort records one.
    ``trigger`` is the trigger file of the recorded source.
    """
    name = f"{run}_imec{probe}_shank{shank}" if run else f"imec{probe}_shank{shank}"
    directory = session_dir / "sorting" / name
    directory.mkdir(parents=True, exist_ok=True)
    if provenance:
        recording = {"sample_rate": FS}
        if source_run or run:
            recording["source"] = (f"D:\\scratch\\{source_run or run}_{trigger}."
                                   f"imec{probe}.sh{shank}.ap.bin")
        (directory / "provenance.json").write_text(json.dumps({
            "gain_to_uv": gain, "duration_s": N_SAMPLES / FS,
            "recording": recording,
        }), encoding="utf-8")
    if verdicts:
        (directory / "quality_metrics.json").write_text(json.dumps({
            "n_units_total": n_units,
            "criteria": {"min_amplitude_uv": 50.0},
            "units": [{"unit_id": u, "pass_strict": u < passing,
                       "pass_rescued": u < passing} for u in range(n_units)],
        }), encoding="utf-8")
    if not sorter_output:
        return

    out = directory / "kilosort4" / "sorter_output"
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(probe * 10 + shank)
    per_unit = 50
    times = np.sort(rng.integers(0, N_SAMPLES, n_units * per_unit))
    clusters = np.repeat(np.arange(n_units), per_unit)
    np.save(out / "spike_times.npy", times.astype(np.int64))
    np.save(out / "spike_clusters.npy", clusters.astype(np.int64))
    np.save(out / "spike_positions.npy", np.column_stack([
        np.full(clusters.size, 40.0),
        np.repeat(np.linspace(0, 700, n_units), per_unit)]).astype(np.float32))
    np.save(out / "channel_positions.npy", np.column_stack([
        np.zeros(8), np.linspace(0, 700, 8)]))
    (out / "params.py").write_text(f"sample_rate = {FS}\n", encoding="utf-8")


@pytest.fixture
def session(tmp_path):
    """Two probes of two shanks each, all sorted."""
    session_dir = tmp_path / "SUBJ_001" / "2026-08-19" / "001"
    for probe in (0, 1):
        for shank in (0, 1):
            write_shank(session_dir, probe, shank)
    return session_dir


# --------------------------------------------------------------------------
# Finding the sorts
# --------------------------------------------------------------------------
def test_every_sorted_shank_is_found_in_probe_then_shank_order(session):
    found = S.find_session_shanks(session)

    assert [s.label for s in found.shanks] == [
        "imec0_shank0", "imec0_shank1", "imec1_shank0", "imec1_shank1"]
    assert found.probes == [0, 1]
    assert not found.incomplete
    assert found.shanks[0].ks.uv_per_bit == 3.0
    assert found.shanks[0].ks.fs == FS


def test_the_sorting_directory_itself_is_also_accepted(session):
    assert len(S.find_session_shanks(session / "sorting").shanks) == 4


def test_a_shank_that_cannot_be_read_is_reported_rather_than_dropped(session):
    """A session quietly becoming a subset of itself is the failure to avoid:
    it would be left out of the rasters and the unit count with nothing
    anywhere saying so."""
    write_shank(session, 2, 0, sorter_output=False)

    found = S.find_session_shanks(session)

    assert len(found.shanks) == 4
    assert any("imec2_shank0" in line for line in found.incomplete)


def test_a_sort_with_no_gain_is_refused_rather_than_guessed(session):
    """Without microvolts per bit every number on the noise figure would be
    ADC counts wearing a microvolt label."""
    write_shank(session, 2, 0)
    (session / "sorting" / "imec2_shank0" / "provenance.json").write_text(
        json.dumps({"duration_s": 1.0, "recording": {"sample_rate": FS}}),
        encoding="utf-8")

    found = S.find_session_shanks(session)

    assert all(s.label != "imec2_shank0" for s in found.shanks)
    assert any("gain_to_uv" in line for line in found.incomplete)


def test_a_session_with_no_sorting_directory_is_not_an_exception(tmp_path):
    found = S.find_session_shanks(tmp_path / "nothing")

    assert found.shanks == []
    assert found.incomplete


# --------------------------------------------------------------------------
# The verdicts, which are read and never recomputed
# --------------------------------------------------------------------------
def test_passing_units_come_from_the_file_the_pipeline_wrote(session):
    shank = S.find_session_shanks(session).shanks[0]

    assert shank.passing_units("pass_rescued") == {0, 1}
    assert shank.passing_units("pass_strict") == {0, 1}


def test_a_sort_with_no_verdicts_is_distinguished_from_one_where_none_passed(
        tmp_path):
    """None and an empty set mean different things, and a figure that showed
    the first as the second would claim the sorter found nothing good."""
    session_dir = tmp_path / "SUBJ" / "2026-08-19" / "001"
    write_shank(session_dir, 0, 0, verdicts=False)
    write_shank(session_dir, 0, 1, passing=0)
    shanks = S.find_session_shanks(session_dir).shanks

    assert shanks[0].passing_units("pass_rescued") is None
    assert shanks[1].passing_units("pass_rescued") == set()


def test_a_criterion_the_file_predates_reads_as_no_verdict(session):
    """quality_metrics.json grew fields over time. Asking for one an older file
    never had must not silently select nothing."""
    shank = S.find_session_shanks(session).shanks[0]

    assert shank.passing_units("pass_a_criterion_from_the_future") is None


# --------------------------------------------------------------------------
# Rows, which are not units in either direction
# --------------------------------------------------------------------------
def test_more_units_than_rows_are_averaged_rather_than_dropped(session):
    """imshow drops rows it has no pixels for, silently."""
    counts = np.arange(100 * 4, dtype=np.float32).reshape(100, 4)

    fitted = S._fit_rows(counts, 10)

    assert fitted.shape == (10, 4)
    assert fitted.sum() == pytest.approx(counts.mean(axis=0).sum() * 10)


def test_fewer_units_than_the_floor_are_stretched_to_it(session):
    """The bug this replaced: the fitting only ever shrank, so a shank with one
    surviving unit kept one row and its label landed on its neighbour's."""
    counts = np.ones((1, 4), dtype=np.float32)

    assert S._fit_rows(counts, 30).shape == (30, 4)


def test_the_label_floor_is_computed_from_the_label_height(session):
    """Fixed guesses of 8 and then 20 rows both looked reasonable and both
    collided, because a block is measured in rows and a label in points."""
    floor = S.min_block_rows(1400)

    assert floor == S.min_block_rows(1400, label_pt=S.SHANK_LABEL_PT)
    assert S.min_block_rows(1400, label_pt=15.0) > floor, "bigger type, more room"
    # Within one row: min_block_rows rounds up.
    assert abs(S.min_block_rows(2800) - 2 * floor) <= 1, "more rows, proportionally more"


def test_binning_puts_every_spike_in_a_bin(session):
    shank = S.find_session_shanks(session).shanks[0]
    edges = np.linspace(0.0, N_SAMPLES / FS, 11)

    counts, depths = S.unit_rates(shank, edges)

    assert counts.shape == (6, 10)
    assert counts.sum() == 6 * 50
    assert list(depths) == sorted(depths), "rows are ordered by depth"


def test_a_filter_selects_only_the_units_asked_for(session):
    shank = S.find_session_shanks(session).shanks[0]
    edges = np.linspace(0.0, N_SAMPLES / FS, 11)

    counts, _ = S.unit_rates(shank, edges, keep={0, 1})

    assert counts.shape[0] == 2


# --------------------------------------------------------------------------
# The figures
# --------------------------------------------------------------------------
def test_the_raster_draws_every_unit_by_default(session):
    shanks = S.find_session_shanks(session).shanks

    figure = S.session_raster(shanks)

    assert "24 units" in figure.axes[0].get_ylabel()


def test_the_filtered_raster_draws_only_the_passing_units(session):
    shanks = S.find_session_shanks(session).shanks

    figure = S.session_raster(shanks, criterion="pass_rescued")

    assert "8 units" in figure.axes[0].get_ylabel()


def test_every_shank_gets_a_legible_label_however_few_units_it_kept(session):
    """The defect a person found by looking, now measured instead.

    A shank left with one unit by quality control still has to say which shank
    it is, and the label is the only cue a colorblind reader has.
    """
    write_shank(session, 1, 2, n_units=40, passing=40)   # one big block
    write_shank(session, 1, 3, n_units=6, passing=1)     # and one tiny one
    shanks = S.find_session_shanks(session).shanks

    figure = S.session_raster(shanks, criterion="pass_rescued")
    figure.canvas.draw()
    axis = figure.axes[0]
    positions = sorted(text.get_transform().transform(text.get_position())[1]
                       for text in axis.texts)
    gaps = np.diff(positions)
    label_px = S.SHANK_LABEL_PT / 72.0 * figure.dpi

    assert len(positions) == len(shanks), "one label per shank"
    assert gaps.min() > label_px, (
        f"labels {gaps.min():.1f} px apart but {label_px:.1f} px tall")


def test_a_session_with_no_units_says_so_rather_than_drawing_nothing(tmp_path):
    session_dir = tmp_path / "SUBJ" / "2026-08-19" / "001"
    write_shank(session_dir, 0, 0, passing=0)
    shanks = S.find_session_shanks(session_dir).shanks

    with pytest.raises(ValueError, match="no units to draw"):
        S.session_raster(shanks, criterion="pass_rescued")


# --------------------------------------------------------------------------
# The report
# --------------------------------------------------------------------------
def test_the_report_makes_both_rasters(session, tmp_path):
    result = build_session_report(session, tmp_path / "out", noise=False)

    assert "session_raster_all_units" in result.made
    assert "session_raster_passing" in result.made
    assert result.skipped["session_rms_across_shanks"] == "not requested"


def test_the_noise_figure_is_skipped_with_a_reason_when_the_files_are_absent(
        session, tmp_path):
    """No archived .cbin here, which is the ordinary case for a sort somebody
    copied off the server without the raw data."""
    result = build_session_report(session, tmp_path / "out", noise=True)

    assert "session_raster_all_units" in result.made
    assert "session_rms_across_shanks" in result.skipped
    assert "imec0_shank0" in result.skipped["session_rms_across_shanks"]


def test_a_session_with_nothing_sorted_reports_instead_of_raising(tmp_path):
    result = build_session_report(tmp_path / "empty", tmp_path / "out")

    assert not result.made
    assert len(result.skipped) == 3


def test_every_quad_base_probe_and_shank_has_its_own_color():
    seen = {S.color_for(p, s) for p in range(4) for s in range(4)}

    assert len(seen) == 16


def test_a_layout_the_palette_does_not_name_still_gets_a_color():
    """Grey on purpose: a wrong hue would read as a probe identity that does
    not exist."""
    assert S.color_for(9, 9) == S.UNKNOWN_COLOR


# --------------------------------------------------------------------------
# Two runs in one session: two computers each number their probes from 0
# --------------------------------------------------------------------------
@pytest.fixture
def two_runs(tmp_path):
    """XL_009 on 2026-10-28. Run p1 was sorted before folder names included
    the run. Run p0 was sorted after. Both have an imec0."""
    session_dir = tmp_path / "XL_009" / "2026-10-28" / "1"
    for shank in (0, 1):
        write_shank(session_dir, 0, shank, source_run="p1_g0")
        write_shank(session_dir, 0, shank, run="p0_g0")
    return session_dir


def test_two_runs_imec0_are_two_probes(two_runs):
    found = S.find_session_shanks(two_runs)

    assert {(s.run, s.probe, s.shank) for s in found.shanks} == {
        ("p0_g0", 0, 0), ("p0_g0", 0, 1), ("p1_g0", 0, 0), ("p1_g0", 0, 1)}
    assert len(found.probes) == 2
    assert len({s.color for s in found.shanks}) == 4
    assert sorted(s.short_label for s in found.shanks) == [
        "p0:0.0", "p0:0.1", "p1:0.0", "p1:0.1"]
    assert {s.label for s in found.shanks} == {
        "imec0_shank0", "imec0_shank1", "p0_g0_imec0_shank0", "p0_g0_imec0_shank1"}


def test_one_run_named_either_way_reads_as_before(tmp_path):
    for run in (None, "p1_g0"):
        session_dir = tmp_path / str(run) / "1"
        write_shank(session_dir, 0, 0, run=run, source_run="p1_g0")
        write_shank(session_dir, 1, 0, run=run, source_run="p1_g0")

        found = S.find_session_shanks(session_dir)

        assert [s.short_label for s in found.shanks] == ["0.0", "1.0"]
        assert [s.color for s in found.shanks] == [S.color_for(0, 0), S.color_for(1, 0)]


def test_a_shanks_raw_file_is_its_own_runs(two_runs):
    for run in ("p0_g0", "p1_g0"):
        folder = two_runs / "raw_ephys_data" / f"{run}_imec0"
        folder.mkdir(parents=True)
        (folder / f"{run}_t0.imec0.sh0.ap.cbin").write_bytes(b"")
    found = S.find_session_shanks(two_runs)

    for shank in (s for s in found.shanks if s.shank == 0):
        cbin, _ch = S._ap_paths(shank)
        assert cbin.name.startswith(shank.run)


def test_find_band_takes_a_run_and_refuses_to_guess_between_two(two_runs):
    from qualitymetrics.raw import RawError, find_band

    for run in ("p0_g0", "p1_g0"):
        folder = two_runs / "raw_ephys_data" / f"{run}_imec0"
        folder.mkdir(parents=True)
        (folder / f"{run}_t0.imec0.sh0.ap.cbin").write_bytes(b"")
    assert find_band(two_runs, "imec0", 0, run="p1_g0").name.startswith("p1_g0")
    assert find_band(two_runs, "imec0", 0, run="p0_g0").name.startswith("p0_g0")
    with pytest.raises(RawError, match="Pass the run"):
        find_band(two_runs, "imec0", 0)


def test_find_band_without_a_run_finds_a_single_runs_file(tmp_path):
    from qualitymetrics.raw import find_band

    folder = tmp_path / "raw_ephys_data" / "p1_g0_imec2"
    folder.mkdir(parents=True)
    (folder / "p1_g0_t0.imec2.sh3.ap.cbin").write_bytes(b"")
    assert find_band(tmp_path, "imec2", 3).name == "p1_g0_t0.imec2.sh3.ap.cbin"


def test_a_sort_naming_no_run_beside_named_ones_is_the_same_run(tmp_path):
    session_dir = tmp_path / "1"
    write_shank(session_dir, 0, 0)                 # Older sort: no source.
    write_shank(session_dir, 0, 1, source_run="p1_g0")
    found = S.find_session_shanks(session_dir)
    assert {s.run for s in found.shanks} == {"p1_g0"}
    assert len(found.probes) == 1
    assert [s.short_label for s in found.shanks] == ["0.0", "0.1"]


def test_find_recording_reads_both_names_and_stays_in_its_session(tmp_path):
    from qualitymetrics.report import find_recording

    session_dir = tmp_path / "XL_009" / "2026-10-28" / "1"
    write_shank(session_dir, 1, 3, run="p1_g0")
    write_shank(session_dir, 0, 0, source_run="p1_g0")
    for probe, shank in ((1, 3), (0, 0)):
        folder = session_dir / "raw_ephys_data" / f"p1_g0_imec{probe}"
        folder.mkdir(parents=True, exist_ok=True)
        (folder / f"p1_g0_t0.imec{probe}.sh{shank}.lf.cbin").write_bytes(b"")
    found = find_recording(session_dir / "sorting" / "p1_g0_imec1_shank3", "lf")
    assert found.name == "p1_g0_t0.imec1.sh3.lf.cbin"
    found = find_recording(session_dir / "sorting" / "imec0_shank0", "lf")
    assert found.name == "p1_g0_t0.imec0.sh0.lf.cbin"
    # Another session of the subject has the band that this session lacks.
    # find_recording must not return it.
    other = tmp_path / "XL_009" / "2026-10-27" / "1" / "raw_ephys_data" / "p1_g0_imec1"
    other.mkdir(parents=True)
    (other / "p1_g0_t0.imec1.sh3.ap.cbin").write_bytes(b"")
    assert find_recording(session_dir / "sorting" / "p1_g0_imec1_shank3", "ap") is None


def test_a_second_trigger_files_folder_finds_its_files(tmp_path):
    from qualitymetrics.raw import find_band

    session_dir = tmp_path / "1"
    write_shank(session_dir, 0, 0, run="p4_t1_g0", source_run="p4_g0", trigger="t1")
    folder = session_dir / "raw_ephys_data" / "p4_t1_g0_imec0"
    folder.mkdir(parents=True)
    (folder / "p4_g0_t1.imec0.sh0.ap.cbin").write_bytes(b"")
    [shank] = S.find_session_shanks(session_dir).shanks
    assert shank.run == "p4_g0"
    assert S._ap_paths(shank)[0].name == "p4_g0_t1.imec0.sh0.ap.cbin"
    odd = session_dir / "raw_ephys_data" / "a[1]_g0_imec0"
    odd.mkdir()
    (odd / "a[1]_g0_t0.imec0.sh2.ap.cbin").write_bytes(b"")
    assert find_band(session_dir, "imec0", 2, run="a[1]_g0").name.startswith("a[1]")


def test_a_run_named_like_a_trigger_file_finds_its_own_files(tmp_path):
    from qualitymetrics.report import find_recording

    session_dir = tmp_path / "1"
    write_shank(session_dir, 0, 0, run="sess_t3_g0")
    folder = session_dir / "raw_ephys_data" / "sess_t3_g0_imec0"
    folder.mkdir(parents=True)
    for band in ("ap", "lf"):
        (folder / f"sess_t3_g0_t0.imec0.sh0.{band}.cbin").write_bytes(b"")
    [shank] = S.find_session_shanks(session_dir).shanks
    assert shank.run == "sess_t3_g0"
    assert S._ap_paths(shank)[0].name == "sess_t3_g0_t0.imec0.sh0.ap.cbin"
    found = find_recording(session_dir / "sorting" / "sess_t3_g0_imec0_shank0", "lf")
    assert found.name == "sess_t3_g0_t0.imec0.sh0.lf.cbin"


def test_eight_probes_get_eight_hues():
    """Two computers can record eight probes in one session."""
    colors = {S.color_for(probe, shank) for probe in range(8) for shank in range(4)}
    assert len(colors) == 32
    assert S.UNKNOWN_COLOR not in colors


def test_a_missing_shank_leaves_its_violin_slot_empty():
    """The shanks are not labeled, so a violin's place says which shank it is."""
    rng = np.random.default_rng(0)
    measurements = [{
        "probe": probe, "probe_name": f"imec{probe}", "slot": probe, "shank": shank,
        "per_channel": {"unreferenced_uv": rng.normal(12, 1, 20),
                        "residual_uv": rng.normal(7, 1, 20)},
    } for probe, shank in ((0, 0), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (1, 3))]

    figure = S.rms_across_shanks(measurements)

    from matplotlib.collections import PolyCollection

    bodies = [c for c in figure.axes[0].collections if isinstance(c, PolyCollection)]
    centers = [round(float(np.mean(c.get_paths()[0].vertices[:, 0])), 1) for c in bodies]
    stride = 4 + 0.7
    assert centers == [0.0, 2.0, 3.0, stride, stride + 1, stride + 2, stride + 3]


def test_find_band_does_not_take_another_probes_file(tmp_path):
    from qualitymetrics.raw import RawError, find_band

    folder = tmp_path / "raw_ephys_data" / "p1_g0_imec1"
    folder.mkdir(parents=True)
    (folder / "p1_g0_t0.imec1.sh0.ap.cbin").write_bytes(b"")
    with pytest.raises(RawError, match="no ap recording"):
        find_band(tmp_path, "imec0", 0, run="p1_g0")


def test_runs_with_a_shared_prefix_get_different_labels(tmp_path):
    session_dir = tmp_path / "1"
    write_shank(session_dir, 0, 0, run="m_grating_g0")
    write_shank(session_dir, 0, 0, run="m_gfp_g0")
    names = {s.name for s in S.find_session_shanks(session_dir).shanks}
    assert names == {"m_grating imec0 sh0", "m_gfp imec0 sh0"}


def test_violins_of_identical_values_still_draw():
    measurements = [{"probe": 0, "shank": shank, "per_channel": {
        "unreferenced_uv": np.full(20, 10.0), "residual_uv": np.full(20, 10.0)}}
        for shank in range(4)]
    import warnings

    with warnings.catch_warnings():
        # Matplotlib only warns about an axis of zero height.
        warnings.simplefilter("error")
        bottom, top = S.rms_across_shanks(measurements).axes[0].get_ylim()
    assert top > bottom


def test_each_run_is_named_once_over_its_probes():
    rng = np.random.default_rng(0)
    measurements = [{
        "probe": probe, "probe_name": f"{run} imec{probe}", "run": run,
        "slot": slot, "shank": shank,
        "per_channel": {"unreferenced_uv": rng.normal(12, 1, 20),
                        "residual_uv": rng.normal(7, 1, 20)},
    } for slot, (run, probe) in enumerate((("JRS_0053_rec", 0), ("JRS_0053_rec", 1),
                                           ("p1", 0)))
        for shank in range(4)]

    # The names sit under the lower panel.
    texts = [t.get_text() for t in S.rms_across_shanks(measurements).axes[-1].texts]

    assert texts.count("JRS_0053_rec") == 1
    assert texts.count("p1") == 1
    assert texts.count("imec0") == 2


def test_sites_above_the_cut_are_counted_over_their_violin():
    rng = np.random.default_rng(0)
    measurements = [{"probe": 0, "shank": shank, "per_channel": {
        "unreferenced_uv": np.r_[rng.normal(12, 1, 2000), [400.0] * (shank == 2) * 3],
        "residual_uv": rng.normal(7, 1, 2000)}} for shank in range(4)]

    texts = [t.get_text() for t in S.rms_across_shanks(measurements).axes[0].texts]

    assert r"$\uparrow$3 sites" in texts


def test_a_sorts_raw_file_is_its_own_trigger_files(tmp_path):
    from qualitymetrics.report import find_recording

    session_dir = tmp_path / "1"
    write_shank(session_dir, 0, 0, run="p4_t1_g0", source_run="p4_g0", trigger="t1")
    shank_dir = session_dir / "sorting" / "p4_t1_g0_imec0_shank0"
    for folder, trigger in (("p4_g0_imec0", "t0"), ("p4_t1_g0_imec0", "t1")):
        raw = session_dir / "raw_ephys_data" / folder
        raw.mkdir(parents=True)
        for band in ("ap", "lf"):
            (raw / f"p4_g0_{trigger}.imec0.sh0.{band}.cbin").write_bytes(b"")

    [shank] = S.find_session_shanks(session_dir).shanks

    assert S._ap_paths(shank)[0].name == "p4_g0_t1.imec0.sh0.ap.cbin"
    assert find_recording(shank_dir, "lf").name == "p4_g0_t1.imec0.sh0.lf.cbin"


def test_find_band_sees_two_runs_stored_in_different_formats(two_runs):
    from qualitymetrics.raw import RawError, find_band

    for run, suffix in (("p0_g0", "cbin"), ("p1_g0", "bin")):
        folder = two_runs / "raw_ephys_data" / f"{run}_imec0"
        folder.mkdir(parents=True)
        (folder / f"{run}_t0.imec0.sh0.ap.{suffix}").write_bytes(b"")
    with pytest.raises(RawError, match="Pass the run"):
        find_band(two_runs, "imec0", 0)
    assert find_band(two_runs, "imec0", 0, run="p1_g0").suffix == ".bin"


def test_two_trigger_files_of_one_run_are_labeled_apart(tmp_path):
    session_dir = tmp_path / "1"
    write_shank(session_dir, 0, 0, run="p4_g0")
    write_shank(session_dir, 0, 0, run="p4_t1_g0", source_run="p4_g0", trigger="t1")

    shanks = S.find_session_shanks(session_dir).shanks

    assert {s.name for s in shanks} == {"p4 imec0 sh0", "p4 t1 imec0 sh0"}
    assert len({s.color for s in shanks}) == 2
