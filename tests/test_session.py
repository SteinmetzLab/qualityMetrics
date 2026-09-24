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
                passing: int = 2) -> None:
    """One shank's sort output, shaped the way SortingManager leaves it."""
    directory = session_dir / "sorting" / f"imec{probe}_shank{shank}"
    directory.mkdir(parents=True, exist_ok=True)
    if provenance:
        (directory / "provenance.json").write_text(json.dumps({
            "gain_to_uv": gain, "duration_s": N_SAMPLES / FS,
            "recording": {"sample_rate": FS},
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
    assert S.min_block_rows(2800) == 2 * floor, "more rows, proportionally more"


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
