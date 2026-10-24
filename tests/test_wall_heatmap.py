"""wall_heatmap's second row: the same moments as Kilosort received them.

The raw row alone misled once (a laser artifact in it was read as a residual of
the artifact removal), and the first attempt at a second row drew the sorted
data 100 ms late, because the kept windows begin one margin before their
``starts``. These tests pin both: the second row holds the recording's own
samples at the raw row's times, channels sit at their own depths whatever the
file leaves out, and the figure says plainly when it has no second row.
"""
from __future__ import annotations

import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from qualitymetrics import RawRecording, plots  # noqa: E402
from qualitymetrics import report as qm_report  # noqa: E402
from qualitymetrics.pipelinefiles import (  # noqa: E402
    LIGHT_EVENTS_FILE,
    PIPELINE_FILE,
    DestripedWindows,
    load_destriped_windows,
    load_light_events,
)
from qualitymetrics.plots.rawviews import (  # noqa: E402
    RAW_ROW_LABEL,
    SORTED_ROW_LABEL,
    _depth_order,
)

FS = 30000.0
N_CH = 16
DURATION_S = 4.0
N_SAMPLES = int(FS * DURATION_S)
UV_PER_BIT = 0.62 / 2048 / 100 * 1e6       # what the meta below implies
UV_PER_COUNT = 0.1                          # SortingManager's DESTRIPED_UV_PER_COUNT
WIDTH = int(FS)                             # one second, as SortingManager keeps
MARGIN = int(0.1 * FS)
FRACTIONS = (0.1, 0.5, 0.9)
WIN_MS = 25.0
SPIKE_AT = 300                              # samples into each wall window
SPIKE_CHANNELS = (3, 4, 5)
#: Depth row of each channel, deliberately not in channel order, the way an
#: NP2 Quad shank runs 0, 2880, 15, 2895 um: a figure that confused the column
#: of a file with a depth would put things visibly in the wrong place.
DEPTH_ROW = np.random.default_rng(7).permutation(N_CH)
#: The acquired channel (SpikeGLX's AP<n>) of each column, as on a Quad Base
#: shank file: neither the column number nor in its order.
ACQUIRED = 1535 - 5 * np.arange(N_CH)


@pytest.fixture(autouse=True)
def _close_figures():
    import matplotlib.pyplot as plt
    yield
    plt.close("all")


def _wall_starts() -> list[int]:
    """The first sample of each default wall_heatmap window."""
    span = DURATION_S - WIN_MS / 1000.0
    return [int(round(span * f * FS)) for f in FRACTIONS]


@pytest.fixture
def session(tmp_path):
    """A flat recording with a big spike planted in each wall window, and a
    sort folder beside it. Returns (recording path, sort folder, data)."""
    rng = np.random.default_rng(0)
    data = rng.normal(0, 8, size=(N_SAMPLES, N_CH)).astype(np.int16)
    # A sharp trough with its minimum exactly SPIKE_AT samples into each window.
    kernel = np.array([-100, -200, -400, -200, -100], dtype=np.int16)
    for start in _wall_starts():
        s = start + SPIKE_AT
        data[s - 2:s + 3, list(SPIKE_CHANNELS)] += kernel[:, None]
    path = tmp_path / "test_g0_t0.imec0.sh0.ap.bin"
    data.tofile(path)
    geom = "".join(f"(0:{32 * (i % 2)}:{DEPTH_ROW[i] * 15}:1)" for i in range(N_CH))
    path.with_suffix(".meta").write_text("\n".join([
        f"nSavedChans={N_CH}", f"imSampRate={FS}",
        f"fileSizeBytes={N_SAMPLES * N_CH * 2}", "imAiRangeMax=0.62",
        "imMaxInt=2048", "imChan0apGain=100", f"snsApLfSy={N_CH},0,0",
        f"~snsGeomMap=(NP2021,4,250,70){geom}",
        f"~snsChanMap=({N_CH},0,1)" + "".join(
            f"(AP{h};{i}:{i})" for i, h in enumerate(ACQUIRED))]) + "\n",
        encoding="utf-8")
    sort_dir = tmp_path / "imec0_shank0"
    sort_dir.mkdir()
    return path, sort_dir, data


def _keep(sort_dir, data, *, channels=None, columns_shuffled=False,
          constant_per_channel=False, acquired_names=False, **extra):
    """Write destriped_windows.npz the way ks4driver.save_destriped_windows does:
    starts centered on the fractions, traces from one margin before them."""
    channels = np.arange(N_CH) if channels is None else np.asarray(channels)
    if columns_shuffled:
        channels = np.random.default_rng(3).permutation(channels)
    starts, traces = [], []
    for f in FRACTIONS:
        start = int(np.clip(int(f * N_SAMPLES) - WIDTH // 2, MARGIN,
                            N_SAMPLES - WIDTH - MARGIN))
        starts.append(start)
        uv = data[start - MARGIN:start + WIDTH + MARGIN][:, channels] * UV_PER_BIT
        if constant_per_channel:
            uv = np.broadcast_to(10.0 * (channels + 1), uv.shape)
        traces.append(uv)
    counts = np.round(np.stack(traces) / UV_PER_COUNT).astype(np.int16)
    np.savez_compressed(
        sort_dir / PIPELINE_FILE, traces=counts, uv_per_count=UV_PER_COUNT,
        starts=np.asarray(starts), margin=MARGIN, fs=FS,
        fractions=np.asarray(FRACTIONS),
        channel_ids=np.asarray([f"AP{ACQUIRED[c]}" if acquired_names else str(c)
                                for c in channels]),
        adc_group_median=True,
        description=np.asarray("pipeline destripe (test), before motion "
                               "interpolation; traces times uv_per_count is "
                               "microvolts"), **extra)
    return load_destriped_windows(sort_dir)


def _images(fig, row_count=2):
    """The heatmap axes as rows, leaving out the colorbar."""
    heat = [ax for ax in fig.axes if ax.get_label() != "<colorbar>"]
    n_cols = len(heat) // row_count
    return [heat[r * n_cols:(r + 1) * n_cols] for r in range(row_count)]


def _texts(fig):
    return [t.get_text() for t in fig.texts] + [
        t.get_text() for ax in fig.axes for t in ax.texts] + [
        fig._suptitle.get_text() if fig._suptitle else ""]


# ----------------------------------------------------------------- the indexing

def test_saved_rows_begin_one_margin_before_the_core_second():
    """Row r of a window holds sample starts - margin + r, not starts + r."""
    start, margin, width = 500, 100, 1000
    rows = np.arange(start - margin, start + width + margin, dtype=np.int16)
    kept = DestripedWindows(traces=rows[None, :, None], uv_per_count=1.0,
                            starts=np.array([start]), margin=margin, fs=FS,
                            channels=np.array([0]))
    got = kept.samples_uv(700, 750)[:, 0]
    np.testing.assert_array_equal(got, np.arange(700, 750))
    # Only the core second is promised; the margins are for filters.
    assert kept.samples_uv(start - 10, start + 10) is None
    assert kept.samples_uv(start + width - 10, start + width + 10) is None
    assert kept.samples_uv(start, start + width) is not None


def test_the_as_sorted_row_is_the_recordings_own_samples_at_the_raw_rows_times(session):
    """Sample for sample: the test that a margin mistake cannot pass."""
    path, sort_dir, _data = session
    with RawRecording.open(path) as rec:
        kept = _keep(sort_dir, rec.get_traces(0, DURATION_S, unit="int16"))
        fig = plots.wall_heatmap(rec, win_ms=WIN_MS, cleaned=kept)
        order, _depths = _depth_order(rec)
        raw_axes, sorted_axes = _images(fig)
        span = rec.duration_s - WIN_MS / 1000.0
        for f, ax in zip(FRACTIONS, sorted_axes, strict=True):
            t0 = span * f
            want = rec.get_traces(t0, t0 + WIN_MS / 1000.0, channels=order)
            got = np.ma.filled(ax.images[0].get_array(), np.nan).T
            assert got.shape == want.shape
            # Within the int16 rounding of the saved copy.
            np.testing.assert_allclose(got, want, atol=UV_PER_COUNT / 2 + 1e-4)


def test_a_planted_spike_is_at_the_same_time_in_both_rows(session):
    path, sort_dir, data = session
    with RawRecording.open(path) as rec:
        kept = _keep(sort_dir, data)
        fig = plots.wall_heatmap(rec, win_ms=WIN_MS, cleaned=kept)
        order, _depths = _depth_order(rec)
        row = int(np.flatnonzero(order == SPIKE_CHANNELS[1])[0])
        for raw_ax, sorted_ax in zip(*_images(fig), strict=True):
            raw = raw_ax.images[0].get_array()[row]
            clean = np.ma.filled(sorted_ax.images[0].get_array(), np.nan)[row]
            assert clean.min() < -500, "the spike is missing from the sorted row"
            assert abs(int(np.argmin(raw)) - int(np.argmin(clean))) <= 1
            # And where it was planted, on the time axis the figure labels.
            x0, x1 = sorted_ax.images[0].get_extent()[:2]
            t_ms = x0 + (np.argmin(clean) + 0.5) * (x1 - x0) / clean.size
            assert t_ms == pytest.approx(SPIKE_AT / FS * 1000, abs=1e-6)


@pytest.mark.parametrize("acquired_names", [False, True],
                         ids=["named by column", "named AP<n> (probe-wide sort)"])
def test_dropped_channels_land_on_the_same_rows_as_in_the_raw_row(session,
                                                                  acquired_names):
    """The file leaves channels out and lists the rest in any order; each must
    still be drawn on the image row its channel has in the raw panel, and the
    missing ones show as missing (NaN, drawn gray). A probe-wide sort names its
    channels by acquired channel, which the recording's meta maps to columns."""
    path, sort_dir, data = session
    dropped = {int(np.argmin(DEPTH_ROW)), 6, 11}      # the shallowest, two inside
    channels = [c for c in range(N_CH) if c not in dropped]
    kept = _keep(sort_dir, data, channels=channels, columns_shuffled=True,
                 constant_per_channel=True, acquired_names=acquired_names)
    assert (kept.channels is None) == acquired_names
    names = [f"AP{ACQUIRED[c]}" if acquired_names else str(c) for c in channels]
    assert list(kept.channel_names) != names, "the test needs shuffled columns"
    with RawRecording.open(path) as rec:
        fig = plots.wall_heatmap(rec, win_ms=WIN_MS, cleaned=kept)
        order, depths = _depth_order(rec)
        raw_axes, sorted_axes = _images(fig)
        for raw_ax, ax in zip(raw_axes, sorted_axes, strict=True):
            image = ax.images[0].get_array()
            assert ax.images[0].get_extent() == raw_ax.images[0].get_extent()
            for k, ch in enumerate(order):
                values = np.ma.filled(image[k], np.nan)
                if ch in dropped:
                    assert np.isnan(values).all(), f"channel {ch} was not sorted"
                else:
                    np.testing.assert_allclose(values, 10.0 * (ch + 1), atol=0.06)
        # The gray rows must not poison the shared scale: finite, and the
        # percentile of every drawn value from both rows.
        pooled = np.concatenate([
            np.abs(np.ma.compressed(ax.images[0].get_array()))
            for row in (raw_axes, sorted_axes) for ax in row])
        lo, hi = sorted_axes[0].images[0].get_clim()
        assert np.isfinite(hi) and lo == -hi
        assert hi == pytest.approx(np.percentile(pooled, 99.5), rel=1e-5)
        # The depth each image row is drawn at is the depth of its channel.
        assert np.all(np.diff(depths) >= 0)
        np.testing.assert_allclose(depths, np.sort(DEPTH_ROW * 15.0))


# --------------------------------------------------------------- when it is absent

def test_without_the_kept_windows_the_figure_is_the_raw_row_alone(session):
    path, _sort_dir, _data = session
    with RawRecording.open(path) as rec:
        fig = plots.wall_heatmap(rec, win_ms=WIN_MS)
        heat = [ax for ax in fig.axes if ax.get_label() != "<colorbar>"]
        assert len(heat) == 3 and all(len(ax.images) == 1 for ax in heat)
        assert not any(SORTED_ROW_LABEL in t for t in _texts(fig))


def test_kept_windows_that_miss_every_panel_leave_the_raw_row_and_say_why(session):
    path, sort_dir, data = session
    with RawRecording.open(path) as rec:
        kept = _keep(sort_dir, data)
        fig = plots.wall_heatmap(rec, t_starts_s=[0.0, 3.9], win_ms=WIN_MS,
                                 cleaned=kept)
        heat = [ax for ax in fig.axes if ax.get_label() != "<colorbar>"]
        assert len(heat) == 2
        assert any("No as-sorted row" in t for t in _texts(fig))


def test_a_window_the_sort_did_not_keep_is_blank_and_says_so(session):
    path, sort_dir, data = session
    with RawRecording.open(path) as rec:
        kept = _keep(sort_dir, data)
        covered = (rec.duration_s - WIN_MS / 1000.0) * 0.5
        fig = plots.wall_heatmap(rec, t_starts_s=[covered, 0.0], win_ms=WIN_MS,
                                 cleaned=kept)
        raw_axes, sorted_axes = _images(fig)
        assert len(sorted_axes[0].images) == 1
        assert len(sorted_axes[1].images) == 0
        assert any("Not among" in t.get_text() for t in sorted_axes[1].texts)


def test_a_file_without_channel_ids_is_not_drawn_by_column_position(session):
    path, sort_dir, data = session
    kept = _keep(sort_dir, data)
    anonymous = DestripedWindows(traces=kept.traces, uv_per_count=kept.uv_per_count,
                                 starts=kept.starts, margin=kept.margin,
                                 fs=kept.fs, channels=None)
    with RawRecording.open(path) as rec:
        fig = plots.wall_heatmap(rec, win_ms=WIN_MS, cleaned=anonymous)
        assert len([a for a in fig.axes if a.get_label() != "<colorbar>"]) == 3
        assert any("channel_ids" in t for t in _texts(fig))


# ------------------------------------------------------------------- the figure

def test_all_six_panels_share_one_color_scale_and_every_axis_is_labeled(session):
    path, sort_dir, data = session
    with RawRecording.open(path) as rec:
        fig = plots.wall_heatmap(rec, win_ms=WIN_MS, cleaned=_keep(sort_dir, data))
        rows = _images(fig)
        images = [im for row in rows for ax in row for im in ax.images]
        assert len(images) == 6
        assert len({im.get_clim() for im in images}) == 1
        for row in rows:
            assert row[0].get_ylabel() == "Depth (µm)"
        assert all(ax.get_xlabel() == "Time from window start (ms)"
                   for ax in rows[-1])
        texts = _texts(fig)
        assert RAW_ROW_LABEL in texts and SORTED_ROW_LABEL in texts
        cbar = next(ax for ax in fig.axes if ax.get_label() == "<colorbar>")
        assert cbar.get_ylabel() == "Voltage (µV)"
        for row in rows:
            for ax in row:
                assert not ax.spines["top"].get_visible()
                assert not ax.spines["right"].get_visible()


def test_laser_edges_inside_a_window_are_marked_in_both_rows(session):
    path, sort_dir, data = session
    _keep(sort_dir, data)
    t_mid = (DURATION_S - WIN_MS / 1000.0) * 0.5
    i0 = int(round(t_mid * FS))
    onset, offset = (i0 + 150) / FS, (i0 + 450) / FS      # 5 and 15 ms in
    np.savez(sort_dir / LIGHT_EVENTS_FILE,
             onset_s=np.array([0.05, onset, 3.0]),
             offset_s=np.array([0.06, offset, 3.01]))
    with RawRecording.open(path) as rec:
        fig = plots.wall_heatmap(rec, win_ms=WIN_MS,
                                 cleaned=load_destriped_windows(sort_dir),
                                 light_events=load_light_events(sort_dir))
        rows = _images(fig)
        for col, expected in ((0, []), (1, [5.0, 15.0]), (2, [])):
            for row in rows:
                ax = row[col]
                xs = sorted(float(line.get_xdata()[0]) for line in ax.lines
                            if line.get_linestyle() == ":")
                np.testing.assert_allclose(xs, expected, atol=1e-6)
        labels = [t.get_text() for t in rows[0][1].texts]
        assert " Laser onset" in labels and " Laser offset" in labels
        assert any("laser onsets and offsets" in t for t in _texts(fig))


# ---------------------------------------------------------------- the files

def test_the_report_reads_both_files_from_the_sort_folder(session):
    path, sort_dir, data = session
    _keep(sort_dir, data)
    (sort_dir / "provenance.json").write_text(json.dumps(
        {"destripe": {"light_artifacts": {"applied": True}}}), encoding="utf-8")
    with RawRecording.open(path) as rec:
        fig = qm_report._wall_heatmap(rec, sort_dir, None)
        assert len([a for a in fig.axes if a.get_label() != "<colorbar>"]) == 6
        assert any("light artifacts removed" in t for t in _texts(fig))


def test_the_report_says_when_the_sort_kept_no_windows(session):
    path, sort_dir, _data = session
    with RawRecording.open(path) as rec:
        fig = qm_report._wall_heatmap(rec, sort_dir, None)
        assert len([a for a in fig.axes if a.get_label() != "<colorbar>"]) == 3
        assert any(f"no {PIPELINE_FILE}" in t for t in _texts(fig))


def test_a_file_whose_parts_disagree_is_refused(tmp_path):
    np.savez(tmp_path / PIPELINE_FILE, traces=np.zeros((3, 100, 4), np.int16),
             uv_per_count=0.1, starts=np.array([20, 40, 60]), margin=10, fs=FS,
             channel_ids=np.array(["0", "1", "2"]))
    with pytest.raises(ValueError, match="channel_ids"):
        load_destriped_windows(tmp_path)


def test_light_removal_is_read_from_the_record_not_guessed(session):
    _path, sort_dir, data = session
    assert _keep(sort_dir, data).light_artifacts_removed is None
    (sort_dir / "provenance.json").write_text(json.dumps(
        {"destripe": {"light_artifacts": {"applied": False}}}), encoding="utf-8")
    assert load_destriped_windows(sort_dir).light_artifacts_removed is False


def test_a_probe_wide_sorts_ap_names_load_and_still_feed_the_noise_figures(session):
    """Names other than column numbers are kept, not refused: the noise figures
    use only the count, and refusing dropped such shanks from them entirely."""
    from types import SimpleNamespace

    from qualitymetrics.plots import noise

    _path, sort_dir, data = session
    kept = _keep(sort_dir, data, acquired_names=True)
    assert kept.channels is None and kept.channel_names[0] == f"AP{ACQUIRED[0]}"
    _starts, blocks, note = noise.pipeline_windows(
        SimpleNamespace(directory=sort_dir), n_channels=N_CH, fs=FS,
        width=WIDTH, margin=MARGIN)
    assert note is None and len(blocks) == 3


def test_ap_names_the_recording_does_not_have_are_refused_not_guessed(session):
    path, sort_dir, data = session
    kept = _keep(sort_dir, data, acquired_names=True)
    stranger = DestripedWindows(
        traces=kept.traces, uv_per_count=kept.uv_per_count, starts=kept.starts,
        margin=kept.margin, fs=kept.fs, channels=None,
        channel_names=tuple(f"AP{9000 + i}" for i in range(N_CH)))
    with RawRecording.open(path) as rec:
        fig = plots.wall_heatmap(rec, win_ms=WIN_MS, cleaned=stranger)
        assert len([a for a in fig.axes if a.get_label() != "<colorbar>"]) == 3
        assert any("not in this recording" in t for t in _texts(fig))


def test_a_light_session_whose_sort_did_not_remove_them_says_so(session):
    path, sort_dir, data = session
    _keep(sort_dir, data)
    (sort_dir / "provenance.json").write_text(json.dumps(
        {"destripe": {"light_artifacts": {"applied": False}}}), encoding="utf-8")
    with RawRecording.open(path) as rec:
        kept = load_destriped_windows(sort_dir)
        quiet = plots.wall_heatmap(rec, win_ms=WIN_MS, cleaned=kept)
        assert not any("light artifacts" in t for t in _texts(quiet))
        np.savez(sort_dir / LIGHT_EVENTS_FILE, onset_s=np.array([1.0]),
                 offset_s=np.array([1.05]))
        loud = plots.wall_heatmap(rec, win_ms=WIN_MS, cleaned=kept,
                                  light_events=load_light_events(sort_dir))
        assert any("light artifacts NOT removed" in t for t in _texts(loud))
