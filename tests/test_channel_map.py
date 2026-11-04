"""A sort that left channels out: Kilosort's channel j is the recording's
column ``channel_map[j]``, and every figure that reads the raw file for a
Kilosort channel must read that column.

The recording here has one large spike, on one column only (``SPIKE_COLUMN``).
The sort left column ``LEFT_OUT`` out, so Kilosort's channel ``SPIKE_COLUMN - 1``
is that column. Reading Kilosort's channel number as the column, as everything
before 2.3.0 did, reads the flat column beside it.

Also here: the uV/bit a caller gives for a meta that states none, and site
positions from ``~snsShankMap``, both for SpikeGLX 20201024's NP2010 metas.
"""
from __future__ import annotations

from types import SimpleNamespace

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import qualitymetrics  # noqa: E402
from qualitymetrics import (  # noqa: E402
    KilosortResults,
    KsError,
    RawError,
    RawRecording,
)
from qualitymetrics.plots import noise as N  # noqa: E402
from qualitymetrics.plots import rawviews  # noqa: E402
from qualitymetrics.plots import session as S  # noqa: E402
from qualitymetrics.raw import (  # noqa: E402
    geometry_from_meta,
    parse_meta,
    resolve_uv_per_bit,
)

FS = 30000.0
N_FILE = 16                 # data columns in the recording
N_SAMPLES = int(FS * 3)
LEFT_OUT = 3                # the column the sort left out
SPIKE_COLUMN = 9            # the only column with the spike
SPIKE_COUNTS = -400         # its size, in counts
SPIKE_TIMES = np.arange(3000, N_SAMPLES - 3000, 6000)
UV = 0.62 / 2048 / 100 * 1e6
#: The sort's channels: every column but LEFT_OUT, in order.
CHANNEL_MAP = np.array([c for c in range(N_FILE) if c != LEFT_OUT])
#: Kilosort's channel of SPIKE_COLUMN.
SPIKE_CHANNEL = int(np.flatnonzero(CHANNEL_MAP == SPIKE_COLUMN)[0])


def _site(column):
    """(x, depth) of a column, as the meta below states it."""
    return (27.0 if column % 2 else 59.0, (column // 2) * 15.0)


def _meta(*, gain=True, geometry="geom", part="NP2021", n=N_FILE, extra=()):
    lines = [f"nSavedChans={n}", f"imSampRate={FS}",
             f"fileSizeBytes={N_SAMPLES * n * 2}", "imAiRangeMax=0.62",
             "imMaxInt=2048", f"snsApLfSy={n},0,0", f"imDatPrb_pn={part}"]
    if gain:
        lines.append("imChan0apGain=100")
    if geometry == "geom":
        lines.append("~snsGeomMap=(NP2021,4,250,70)" + "".join(
            f"(0:{_site(c)[0]:g}:{_site(c)[1]:g}:1)" for c in range(n)))
    elif geometry == "shank":
        lines.append(f"~snsShankMap=(1,2,{n // 2})" + "".join(
            f"(0:{c % 2}:{c // 2}:1)" for c in range(n)))
    lines.extend(extra)
    return "\n".join(lines) + "\n"


def _recording(tmp_path, **meta):
    rng = np.random.default_rng(0)
    data = rng.normal(0, 4, (N_SAMPLES, N_FILE)).astype(np.int16)
    for t in SPIKE_TIMES:
        data[t:t + 12, SPIKE_COLUMN] += SPIKE_COUNTS
    raw = tmp_path / "raw_ephys_data" / "p0_g0_imec0"
    raw.mkdir(parents=True)
    path = raw / "p0_g0_t0.imec0.sh0.ap.bin"
    data.tofile(path)
    path.with_suffix(".meta").write_text(_meta(**meta), encoding="utf-8")
    return path


def _sort(tmp_path, *, channel_map=CHANNEL_MAP, positions=None, write_map=True):
    """A one-unit Kilosort output whose template peaks on SPIKE_CHANNEL."""
    out = tmp_path / "sorting" / "imec0_shank0" / "sorter_output"
    out.mkdir(parents=True)
    n = len(channel_map)
    np.save(out / "spike_times.npy", SPIKE_TIMES.astype(np.int64))
    np.save(out / "spike_clusters.npy", np.zeros(SPIKE_TIMES.size, np.int32))
    np.save(out / "spike_templates.npy", np.zeros(SPIKE_TIMES.size, np.int64))
    np.save(out / "amplitudes.npy", np.full(SPIKE_TIMES.size, 20.0, np.float32))
    depth = _site(SPIKE_COLUMN)[1]
    np.save(out / "spike_positions.npy",
            np.column_stack([np.full(SPIKE_TIMES.size, 40.0),
                             np.full(SPIKE_TIMES.size, depth)]).astype(np.float32))
    t = np.arange(61)
    templates = np.zeros((1, 61, n), np.float32)
    templates[0, :, SPIKE_CHANNEL] = -10 * np.exp(-((t - 20) ** 2) / 8.0)
    np.save(out / "templates.npy", templates)
    np.save(out / "whitening_mat_inv.npy", np.eye(n, dtype=np.float32))
    if positions is None:
        positions = np.array([_site(c) for c in channel_map])
    np.save(out / "channel_positions.npy", np.asarray(positions, dtype=float))
    if write_map:
        np.save(out / "channel_map.npy", np.asarray(channel_map, dtype=np.int32))
    (out / "params.py").write_text(f"sample_rate = {FS}\n", encoding="utf-8")
    (out / "cluster_KSLabel.tsv").write_text("cluster_id\tKSLabel\n0\tgood\n",
                                             encoding="utf-8")
    return out


# ------------------------------------------------------------- the map itself
def test_the_installed_package_says_it_reads_the_channel_map():
    assert {"ks_channel_map", "uv_per_bit_override"} <= qualitymetrics.FEATURES


def test_the_channel_map_is_read_and_identity_without_one(tmp_path):
    ks = KilosortResults.load(_sort(tmp_path / "a"))
    assert ks.channel_map.tolist() == CHANNEL_MAP.tolist() and ks.channels_remapped
    assert ks.recording_columns([SPIKE_CHANNEL]).tolist() == [SPIKE_COLUMN]
    plain = KilosortResults.load(_sort(tmp_path / "b", write_map=False))
    assert plain.channel_map.tolist() == list(range(N_FILE - 1))
    assert not plain.channels_remapped


@pytest.mark.parametrize("bad", [CHANNEL_MAP[:-1], np.r_[CHANNEL_MAP[:-1], 0],
                                 np.r_[-1, CHANNEL_MAP[1:]]])
def test_a_map_that_cannot_be_one_is_refused(tmp_path, bad):
    out = _sort(tmp_path)
    np.save(out / "channel_map.npy", bad)
    with pytest.raises(KsError):
        _ = KilosortResults.load(out).channel_map


# ------------------------------------------------------- the figures' reading
def test_snippets_are_read_from_the_column_the_map_names(tmp_path):
    rec = RawRecording.open(_recording(tmp_path))
    ks = KilosortResults.load(_sort(tmp_path), uv_per_bit=rec.uv_per_bit)
    with rec:
        _t, stack, chans = rawviews.spike_snippets(rec, ks, 0, n_channels=4,
                                                   highpass_hz=None)
    mean = stack.mean(axis=0)
    peak = int(np.flatnonzero(chans == SPIKE_CHANNEL)[0])
    assert mean[:, peak].min() < 0.8 * SPIKE_COUNTS * UV
    others = [i for i in range(len(chans)) if i != peak]
    assert np.abs(mean[:, others]).max() < 0.1 * abs(SPIKE_COUNTS) * UV


def test_a_map_whose_depths_are_not_the_files_is_refused(tmp_path):
    # Positions of the first 15 columns, as if nothing were left out: the map
    # and the positions then name different sites from column 4 on.
    rec = RawRecording.open(_recording(tmp_path))
    ks = KilosortResults.load(_sort(
        tmp_path, positions=[_site(c) for c in range(N_FILE - 1)]))
    with rec, pytest.raises(ValueError, match="do not fit"):
        rawviews.spike_snippets(rec, ks, 0)


def test_a_map_naming_a_column_the_file_lacks_is_refused(tmp_path):
    rec = RawRecording.open(_recording(tmp_path))
    channel_map = np.r_[CHANNEL_MAP[:-1], N_FILE]
    ks = KilosortResults.load(_sort(tmp_path, channel_map=channel_map,
                                    positions=[_site(c) for c in CHANNEL_MAP]))
    with rec, pytest.raises(ValueError, match="names column 16"):
        rawviews.sorted_columns(rec, ks)


@pytest.fixture
def undrawn(monkeypatch):
    """build_report deciding what to draw, without drawing it.

    A full report draws about twenty figures, 45 s on this data, nearly all
    of it matplotlib; these tests are about which figures it attempts and the
    reasons for the rest. test_an_np2010_shank_draws_its_unit_figures draws.
    """
    from qualitymetrics import report

    def attempt(result, name, fn):
        result.made[name] = "not drawn"

    monkeypatch.setattr(report, "_attempt", attempt)
    return report.build_report


def _report(tmp_path, build, raw, **kw):
    return build(tmp_path / "sorting" / "imec0_shank0", tmp_path / "qc",
                 raw_path=raw, write_phy=False, example_units=1, **kw)


def test_the_report_attempts_the_unit_figures_of_a_remapped_sort(tmp_path, undrawn):
    raw = _recording(tmp_path)
    _sort(tmp_path)
    result = _report(tmp_path, undrawn, raw)
    for name in ("raw_with_spikes", "example_neurons", "raw_traces"):
        assert name in result.made, result.skipped.get(name)


def test_the_report_skips_the_unit_figures_when_the_map_does_not_fit(tmp_path,
                                                                     undrawn):
    raw = _recording(tmp_path)
    _sort(tmp_path, positions=[_site(c) for c in range(N_FILE - 1)])
    result = _report(tmp_path, undrawn, raw)
    for name in ("raw_with_spikes", "example_neurons"):
        assert "do not fit" in result.skipped[name]
    assert "raw_traces" in result.made, "reads no Kilosort channel"


# ----------------------------------------------------- the session noise rows
def test_adc_groups_follow_the_columns():
    assert N.adc_groups(3, [0, 2, 35]).tolist() == [0, 1, 1]
    with pytest.raises(ValueError):
        N.adc_groups(2, [0, 1, 2])


def test_spikes_are_not_placed_on_a_recording_without_positions(tmp_path):
    # Before, its channel numbers stood in for depths in micrometres, and
    # the figure showed traces with no unit on them.
    rec = RawRecording.open(_recording(tmp_path, geometry=None))
    ks = KilosortResults.load(_sort(tmp_path), uv_per_bit=UV)
    with rec, pytest.raises(ValueError, match="no site positions"):
        rawviews.raw_with_spikes(rec, ks, t_start_s=0.05, win_ms=20)


def test_the_noise_rows_are_the_sorted_channels_columns(tmp_path):
    """measure_shank reads channel_map's columns, not the first n.

    About 6 s: every stage of three 1 s windows, with Welch at the figures'
    resolution, which is what this function is.
    """
    out = _sort(tmp_path)
    ks = SimpleNamespace(path=out, fs=FS, uv_per_bit=1.0,
                         channel_map=CHANNEL_MAP)
    shank = SimpleNamespace(probe=0, shank=0, label="imec0_shank0", ks=ks,
                            directory=out.parent)
    n_samples = int(FS * 1.5)
    t = np.arange(n_samples) / FS
    data = np.random.default_rng(2).normal(0, 2, (n_samples, N_FILE + 1))
    data[:, SPIKE_COLUMN] += 300 * np.sin(2 * np.pi * 2000 * t)

    class Reader:
        shape = (n_samples, N_FILE + 1)

        def __getitem__(self, rows):
            return data[rows].astype(np.int16)

    saved = N.N_TIME_WINDOWS
    N.N_TIME_WINDOWS = 2
    try:
        m = N.measure_shank(shank, reader_factory=lambda _s: Reader())
    finally:
        N.N_TIME_WINDOWS = saved
    band = m["rms_uv"][1]                    # band-passed, no reference
    assert int(np.argmax(band)) == SPIKE_CHANNEL
    assert band[SPIKE_CHANNEL] > 150 and np.delete(band, SPIKE_CHANNEL).max() < 10
    assert m["depth_um"][SPIKE_CHANNEL] == _site(SPIKE_COLUMN)[1]


def test_the_session_rms_reads_the_sorted_channels_columns(tmp_path, monkeypatch):
    mtscomp = pytest.importorskip("mtscomp")
    out = _sort(tmp_path)
    n_samples = int(FS * 2)
    t = np.arange(n_samples) / FS
    data = np.random.default_rng(3).normal(0, 2, (n_samples, N_FILE + 1))
    data[:, SPIKE_COLUMN] += 300 * np.sin(2 * np.pi * 2000 * t)
    flat = tmp_path / "x.bin"
    data.astype(np.int16).tofile(flat)
    cbin, ch = tmp_path / "x.cbin", tmp_path / "x.ch"
    mtscomp.compress(flat, cbin, ch, sample_rate=FS, n_channels=N_FILE + 1,
                     dtype=np.int16)
    monkeypatch.setattr(S, "_ap_paths", lambda _s: (cbin, ch))
    ks = KilosortResults.load(out, uv_per_bit=1.0)
    shank = SimpleNamespace(probe=0, shank=0, label="imec0_shank0", ks=ks,
                            directory=out.parent, probe_name="imec0",
                            run_token="", several_runs=False, color_probe=0,
                            short_label="imec0 sh0")
    m = S.channel_rms(shank, n_windows=2, window_s=0.5)
    raw = m["per_channel"]["unreferenced_uv"]
    assert int(np.argmax(raw)) == SPIKE_CHANNEL
    assert m["per_channel"]["depth_um"][SPIKE_CHANNEL] == _site(SPIKE_COLUMN)[1]


def test_kept_windows_naming_other_channels_are_refused(tmp_path):
    out = _sort(tmp_path)
    shank = SimpleNamespace(directory=out.parent)
    counts = np.zeros((3, 100, N_FILE - 1), np.int16)
    counts[:, :, 0] = 7                          # column 0's samples
    kept = dict(traces=counts, uv_per_count=1.0, starts=np.array([200, 400, 600]),
                margin=10, fs=FS)
    np.savez(out.parent / N.PIPELINE_FILE, **kept,
             channel_ids=np.array([str(c) for c in CHANNEL_MAP[::-1]]))
    starts, blocks, why = N.pipeline_windows(shank, n_channels=N_FILE - 1, fs=FS,
                                             width=80, margin=10,
                                             columns=CHANNEL_MAP)
    assert why is None
    # Put in the sort's order: column 0 was the file's last channel.
    assert blocks[0][0, -1] == 7 and blocks[0][0, 0] == 0
    np.savez(out.parent / N.PIPELINE_FILE, **kept,
             channel_ids=np.array([str(c) for c in range(N_FILE - 1)]))
    assert N.pipeline_windows(shank, n_channels=N_FILE - 1, fs=FS, width=80,
                              margin=10, columns=CHANNEL_MAP)[2] == \
        "kept other channels than the sort's"


# --------------------------------------------------- the caller's uV per bit
def test_a_meta_with_no_gain_opens_with_the_callers(tmp_path):
    path = _recording(tmp_path, gain=False)
    with pytest.raises(RawError):
        RawRecording.open(path)
    with RawRecording.open(path, uv_per_bit=0.7629) as rec:
        assert rec.uv_per_bit == 0.7629


def test_a_callers_gain_that_disagrees_with_the_meta_is_refused(tmp_path):
    path = _recording(tmp_path)
    with RawRecording.open(path, uv_per_bit=UV * (1 + 1e-6)) as rec:
        assert rec.uv_per_bit == pytest.approx(UV)
    with pytest.raises(RawError, match="disagrees"):
        RawRecording.open(path, uv_per_bit=0.7629)
    with pytest.raises(RawError, match="positive"):
        resolve_uv_per_bit({}, "ap", 0.0)


def test_the_report_of_a_meta_with_no_gain(tmp_path, undrawn):
    """Without a gain the report raised before 2.3.0; it now says why, and
    with the caller's attempts everything."""
    raw = _recording(tmp_path, gain=False)
    _sort(tmp_path)
    without = _report(tmp_path, undrawn, raw)
    assert "imChan0apGain" in without.skipped["raw_traces"]
    assert "imChan0apGain" in without.skipped["unit_drift"]
    given = _report(tmp_path, undrawn, raw, uv_per_bit=UV)
    for name in ("raw_traces", "example_neurons", "unit_drift", "raw_with_spikes"):
        assert name in given.made, given.skipped.get(name)


def test_the_report_without_a_raw_file_takes_the_callers_gain(tmp_path, undrawn):
    _sort(tmp_path)
    # An empty session root, so find_recording stops here rather than
    # searching every other test's folder (48 s, and it could find one).
    (tmp_path / "raw_ephys_data").mkdir()
    result = undrawn(tmp_path / "sorting" / "imec0_shank0", tmp_path / "qc",
                     write_phy=False, example_units=1, uv_per_bit=UV)
    assert "unit_drift" in result.made
    assert "no archived recording" in result.skipped["raw_traces"]


def test_the_report_refuses_a_gain_the_meta_contradicts(tmp_path, undrawn):
    raw = _recording(tmp_path)
    _sort(tmp_path)
    result = _report(tmp_path, undrawn, raw, uv_per_bit=2 * UV)
    for name in ("raw_traces", "example_neurons", "unit_drift", "depth_power_wide"):
        assert "disagrees" in result.skipped[name]


# ----------------------------------------------------------- ~snsShankMap
def test_shank_map_positions_for_an_np2_part(tmp_path):
    meta = parse_meta(_recording(tmp_path, geometry="shank", part="NP2010")
                      .with_suffix(".meta"))
    geom = geometry_from_meta(meta)
    assert geom.y.tolist() == [(c // 2) * 15.0 for c in range(N_FILE)]
    assert geom.x.tolist() == [32.0 * (c % 2) for c in range(N_FILE)]


def test_shank_map_positions_take_the_metas_pitch_and_never_guess(tmp_path):
    stated = _meta(geometry="shank", part="NP1000", extra=("imZPitch=20",
                                                           "imXPitch=16"))
    path = tmp_path / "a.meta"
    path.write_text(stated, encoding="utf-8")
    assert geometry_from_meta(parse_meta(path)).y[2] == 20.0
    path.write_text(_meta(geometry="shank", part="NP1000"), encoding="utf-8")
    with pytest.raises(RawError, match="pitch"):
        geometry_from_meta(parse_meta(path))


def test_an_np2010_shank_draws_its_unit_over_its_own_column(tmp_path, undrawn):
    """No gain and no ~snsGeomMap, as SpikeGLX 20201024 writes, and a channel
    left out, as the NP2010's reference channel is: every raw figure is
    attempted, and the unit's snippets are its own column's."""
    raw = _recording(tmp_path, gain=False, geometry="shank", part="NP2010")
    out = _sort(tmp_path, positions=[(32.0 * (c % 2), (c // 2) * 15.0)
                                     for c in CHANNEL_MAP])
    result = _report(tmp_path, undrawn, raw, uv_per_bit=UV)
    for name in ("raw_traces", "raw_with_spikes", "example_neurons",
                 "wall_heatmap", "channel_health", "unit_drift"):
        assert name in result.made, result.skipped.get(name)
    ks = KilosortResults.load(out, uv_per_bit=UV)
    with RawRecording.open(raw, uv_per_bit=UV) as rec:
        _t, stack, chans = rawviews.spike_snippets(rec, ks, 0, n_channels=4,
                                                   highpass_hz=None)
        # Drawn once for real (2 s): the figure itself runs on such a shank.
        rawviews.raw_with_spikes(rec, ks, t_start_s=0.05, win_ms=20)
    peak = int(np.flatnonzero(chans == SPIKE_CHANNEL)[0])
    assert stack.mean(axis=0)[:, peak].min() < 0.8 * SPIKE_COUNTS * UV
