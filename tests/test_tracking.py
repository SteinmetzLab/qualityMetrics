"""Descriptive continuity flags must never change sorting or censor its duration."""
import csv
import json

import numpy as np
import pytest

from qualitymetrics import tracking as tp


def test_sustained_step_in_either_direction_has_the_right_split():
    edges = np.arange(0, 2401, 10, dtype=float)
    for reverse in [False, True]:
        counts = np.r_[np.zeros(120), np.full(120, 50)]
        if reverse:
            counts = counts[::-1]
        row = tp.screen(counts, np.diff(edges), edges)
        assert row['flagged'] and row['assessed']
        assert row['split_s'] == 1200
        assert row['change'] == ('disappears' if reverse else 'appears')


def test_constant_bursty_short_or_unexposed_units_are_not_candidates():
    edges = np.arange(0, 2401, 10, dtype=float)
    seconds = np.diff(edges)
    assert not tp.screen(np.full(240, 50), seconds, edges)['flagged']
    bursty = np.zeros(240)
    bursty[120::10] = 500
    assert not tp.screen(bursty, seconds, edges)['flagged']
    assert not tp.screen(np.zeros(240), seconds, edges)['assessed']
    assert not tp.screen(np.ones(240), np.zeros(240), edges)['assessed']
    assert not tp.screen(np.ones(60), seconds[:60], edges[:61])['assessed']


def test_best_qualifying_split_matches_a_brute_force_reference():
    """Independent slice-based reference catches prefix-sum off-by-one errors."""
    rng = np.random.default_rng(41)
    edges = np.arange(0, 2401, 10, dtype=float)
    seconds = rng.uniform(0, 10, 240)
    counts = rng.poisson(seconds * np.r_[np.full(120, .01), np.full(120, 6)])
    usable = seconds >= 5
    candidates = []
    for k in range(60, 181):
        a = counts[:k][usable[:k]]; b = counts[k:][usable[k:]]
        ea = seconds[:k][usable[:k]]; eb = seconds[k:][usable[k:]]
        if min(ea.sum(), eb.sum()) < 300:
            continue
        ra, rb = a.sum()/ea.sum(), b.sum()/eb.sum()
        hi, lo = max(ra, rb), min(ra, rb)
        if not hi:
            continue
        hr, lr = (a/ea, b/eb) if ra >= rb else (b/eb, a/ea)
        flag = hi >= 1 and lo/hi <= .1 and np.mean(hr >= .2*hi) >= .8 and np.mean(lr < .2*hi) >= .8
        rate = (a.sum()+b.sum())/(ea.sum()+eb.sum())
        score = 2*sum(n*np.log(r/rate) if n else 0 for n,r in [(a.sum(),ra),(b.sum(),rb)])
        candidates.append((flag, score, edges[k]))
    qualifying = [x for x in candidates if x[0]]
    expected = max(qualifying or candidates, key=lambda x:x[1])
    got = tp.screen(counts, seconds, edges)
    assert got['flagged'] == expected[0]
    assert got['split_s'] == expected[2]
    assert got['poisson_gain'] == pytest.approx(expected[1])


def test_overlapping_exclusions_are_counted_once_with_half_open_boundaries():
    intervals = tp.merge_intervals([[-2,2],[1,5],[5,6],[8,20]],10)
    np.testing.assert_equal(intervals,[[0,6],[8,10]])
    np.testing.assert_equal(tp.exposure(np.array([0,5,10]),intervals),[0,2])
    np.testing.assert_equal(tp.outside_intervals(np.array([0,5,6,7,8,10]),intervals),[False,False,True,True,False,True])


def _sort(tmp_path, duration=2400, times=None, clusters=None):
    root = tmp_path/'sort'; sorter = root/'kilosort4/sorter_output'
    sorter.mkdir(parents=True)
    (root/'provenance.json').write_text(json.dumps(dict(duration_s=duration)))
    (sorter/'params.py').write_text('sample_rate = 1000.0\n')
    if times is None:
        times = np.arange(1200000,2400000,200,dtype=np.int64)
        clusters = np.full(len(times),1000000,dtype=np.int64)
    np.save(sorter/'spike_times.npy',np.asarray(times,dtype=np.int64))
    np.save(sorter/'spike_clusters.npy',np.asarray(clusters,dtype=np.int64))
    (sorter/'cluster_KSLabel.tsv').write_text('cluster_id\tKSLabel\n1000000\tgood\n')
    (sorter/'cluster_group.tsv').write_text('cluster_id\tgroup\n1000000\tgood\n')
    return root,sorter


def test_sparse_ids_and_chunk_boundaries_preserve_exact_counts(tmp_path,monkeypatch):
    root,sorter = _sort(tmp_path,duration=30,times=[0,9999,10000,19999,20000,29999],
                        clusters=[2,1000000,2,2,1000000,2])
    monkeypatch.setattr(tp,'CHUNK_SPIKES',2)
    ids,counts = tp.bin_spikes(sorter,1000,np.array([0,10,20,30]),np.array([[10,20]]))
    np.testing.assert_equal(ids,[2,1000000])
    np.testing.assert_equal(counts,[[1,0,1],[1,0,1]])


def test_late_silence_is_not_truncated_to_the_final_spike(tmp_path):
    root,sorter = _sort(tmp_path,times=np.arange(0,1200000,200),clusters=np.full(6000,1000000))
    before = {p.name:p.read_bytes() for p in sorter.iterdir()}
    report = tp.write(root,tmp_path/'qc',label='shank 0')
    assert report['n_candidates']==1
    saved = np.load(tmp_path/'qc/tracking_continuity_bins.npz')
    assert saved['edges_s'][-1]==2400
    assert saved['counts'][0,120:].sum()==0
    assert before=={p.name:p.read_bytes() for p in sorter.iterdir()}
    with (tmp_path/'qc/tracking_continuity_units.csv').open() as stream:
        row = next(csv.DictReader(stream))
    assert row['cluster_id']=='1000000' and row['change']=='disappears'
    assert set(report['made'])==set(tp.FILES)


def test_mask_and_light_union_removes_spikes_and_exposure(tmp_path):
    root,sorter = _sort(tmp_path,duration=30,times=[1000,10000,15500,22000],clusters=[7]*4)
    np.savez(root/'light_events.npz',onset_s=[10.],offset_s=[15.])
    (root/'mask_intervals.csv').write_text('start_s,stop_s\n14,16\n')
    tp.write(root,tmp_path/'qc')
    data = np.load(tmp_path/'qc/tracking_continuity_bins.npz')
    np.testing.assert_equal(data['counts'],[[1,0,1]])
    assert data['exposure_s'].sum()==pytest.approx(23.7)
    report=json.loads((tmp_path/'qc/tracking_continuity.json').read_text())
    assert set(report['exclusion_sources'])=={'light_events.npz','mask_intervals.csv'}


@pytest.mark.parametrize('kind',['empty','masked','short'])
def test_unassessable_recordings_have_explicit_status_and_a_plot(tmp_path,kind):
    root,_ = _sort(tmp_path,duration=30,times=[] if kind=='empty' else [1000],clusters=[] if kind=='empty' else [7])
    if kind=='masked':
        (root/'mask_intervals.csv').write_text('start_s,stop_s\n0,30\n')
    tp.write(root,tmp_path/'qc')
    report=json.loads((tmp_path/'qc/tracking_continuity.json').read_text())
    assert report['n_assessed']==report['n_candidates']==0
    assert (tmp_path/'qc/tracking_continuity.png').stat().st_size>1000


def test_bad_metadata_does_not_silently_change_recording_length(tmp_path):
    root,_ = _sort(tmp_path,duration=10,times=[10001],clusters=[0])
    with pytest.raises(ValueError,match='outside provenance'):
        tp.write(root,tmp_path/'qc')


def test_missing_duration_is_not_inferred_from_spikes(tmp_path):
    root,_ = _sort(tmp_path)
    (root/'provenance.json').unlink()
    with pytest.raises(FileNotFoundError):
        tp.write(root,tmp_path/'qc')


def test_params_file_is_not_executed(tmp_path):
    root,sorter = _sort(tmp_path)
    (sorter/'params.py').write_text("raise RuntimeError('must not execute')\nsample_rate = 1000.0\n")
    assert tp._sample_rate(sorter)==1000


def test_exact_endpoint_matches_numpy_histogram(tmp_path):
    _, sorter = _sort(tmp_path, duration=10, times=[0, 10000], clusters=[0, 0])
    _, counts = tp.bin_spikes(sorter, 1000, np.array([0., 10.]), np.empty((0, 2)))
    np.testing.assert_equal(counts[0], np.histogram([0., 10.], [0., 10.])[0])


def test_endpoint_spike_is_masked_with_the_final_sample(tmp_path):
    _, sorter = _sort(tmp_path, duration=30, times=[29000, 30000], clusters=[7, 7])
    intervals = tp.merge_intervals([[26, 40]], 30)
    _, counts = tp.bin_spikes(sorter, 1000, np.array([0.,10.,20.,30.]), intervals)
    assert counts.sum() == 0


def test_sorter_output_directly_under_sort_is_supported(tmp_path):
    root, sorter = _sort(tmp_path)
    sorter.rename(root/'sorter_output')
    assert tp.write(root, tmp_path/'qc')['n_units'] == 1


def test_locked_old_figure_does_not_prevent_failure_manifest(tmp_path, monkeypatch):
    from pathlib import Path
    from matplotlib.figure import Figure
    root, _ = _sort(tmp_path)
    out = tmp_path/'qc'
    tp.write(root, out)
    unlink = Path.unlink
    def locked(self, *args, **kwargs):
        if self == out/'tracking_continuity.png':
            raise PermissionError('figure is locked')
        return unlink(self, *args, **kwargs)
    def fail(*args, **kwargs):
        raise OSError('redraw failed')
    monkeypatch.setattr(Path, 'unlink', locked)
    monkeypatch.setattr(Figure, 'savefig', fail)
    with pytest.raises(OSError, match='redraw failed'):
        tp.write(root, out)
    error = json.loads((out/'tracking_continuity_error.json').read_text())
    assert 'tracking_continuity.png' in error['stale_files']
    assert not (out/'tracking_continuity.json').exists()


def test_short_active_period_remains_visible_even_when_95th_percentile_is_zero():
    counts = np.r_[np.full(60, 50), np.zeros(1380)][None, :]
    seconds = np.full(1440, 10.)
    assert tp.screen(counts[0], seconds, np.arange(1441)*10.)['flagged']
    norm = tp.normalized_rates(counts, seconds)
    np.testing.assert_equal(norm[0, :60], 1.)
    np.testing.assert_equal(norm[0, 60:], 0.)


def test_failed_redraw_cannot_leave_an_old_plot_beside_new_tables(tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    root, _ = _sort(tmp_path)
    out = tmp_path/'qc'
    tp.write(root, out)
    (root/'mask_intervals.csv').write_text('start_s,stop_s\n0,2400\n')
    def fail(*args, **kwargs):
        raise OSError('cannot save figure')
    monkeypatch.setattr(Figure, 'savefig', fail)
    with pytest.raises(OSError, match='cannot save'):
        tp.write(root, out)
    assert not any((out/name).exists() for name in tp.FILES)
    assert 'cannot save figure' in (out/'tracking_continuity_error.json').read_text()


def test_partial_final_bin_uses_its_true_time_edges(tmp_path, monkeypatch):
    from matplotlib.axes import Axes
    original = Axes.pcolormesh
    captured = []
    def capture(self, *args, **kwargs):
        if len(args) == 3 and np.ndim(args[0]) == 1:
            captured.append(np.asarray(args[0]))
        return original(self, *args, **kwargs)
    monkeypatch.setattr(Axes, 'pcolormesh', capture)
    root, _ = _sort(tmp_path, duration=25, times=[1000, 21000], clusters=[7, 7])
    tp.write(root, tmp_path/'qc')
    assert any(np.array_equal(x, np.array([0,10,20,25])/60) for x in captured)


def test_report_includes_continuity_without_raw_and_survives_its_failure(sorter_output, tmp_path, monkeypatch):
    from qualitymetrics import report
    monkeypatch.setattr(report, 'find_recording', lambda *a: None)
    monkeypatch.setattr(report, '_attempt', lambda *a: None)
    (sorter_output/'provenance.json').write_text('{"duration_s": 4}')
    result = report.build_report(sorter_output, tmp_path/'qc', write_phy=False)
    assert 'tracking_continuity' in result.made
    assert 'tracking_continuity_bins.npz' in result.made
    def fail(*a, **k):
        raise RuntimeError('failed continuity diagnostic')
    monkeypatch.setattr(tp, 'write', fail)
    result = report.build_report(sorter_output, tmp_path/'qc2', write_phy=False)
    assert 'failed continuity diagnostic' in result.skipped['tracking_continuity']
    assert 'raw_traces' in result.skipped
