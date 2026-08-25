"""A synthetic sorter output and recording, so the suite needs no server.

Small enough to build in a fraction of a second, real enough that the loaders,
the unit conversions and the figure code all run against it.
"""
from __future__ import annotations

import numpy as np
import pytest

FS = 30000.0
N_CHANNELS = 16
N_SAMPLES = int(FS * 4)          # 4 seconds
N_TEMPLATES = 5
UV_PER_BIT = 3.0273


def _meta_text(n_channels: int, n_samples: int, fs: float, band: str) -> str:
    """A SpikeGLX meta with the fields this package actually reads."""
    geom = "".join(f"(0:{27 if i % 2 else 59}:{(i // 2) * 15}:1)"
                   for i in range(n_channels))
    lines = [
        f"nSavedChans={n_channels}",
        f"imSampRate={fs}",
        f"fileSizeBytes={n_samples * n_channels * 2}",
        "imAiRangeMax=0.62",
        "imMaxInt=2048",
        "imChan0apGain=100",
        f"snsApLfSy={n_channels},0,0" if band == "ap" else f"snsApLfSy=0,{n_channels},0",
        f"~snsGeomMap=(NP2021,4,250,70){geom}",
    ]
    return "\n".join(lines) + "\n"


@pytest.fixture
def recording_dir(tmp_path):
    """A flat .bin recording plus its meta. Noise with a few big deflections."""
    rng = np.random.default_rng(0)
    data = (rng.normal(0, 8, size=(N_SAMPLES, N_CHANNELS))).astype(np.int16)
    # Plant spikes on channel 8 so a trace figure has something to show.
    for t in range(1000, N_SAMPLES - 1000, 3000):
        data[t:t + 20, 7:10] -= 60
    path = tmp_path / "test_g0_t0.imec0.sh0.ap.bin"
    data.tofile(path)
    path.with_suffix(".meta").write_text(
        _meta_text(N_CHANNELS, N_SAMPLES, FS, "ap"), encoding="utf-8")
    return path


@pytest.fixture
def sorter_output(tmp_path):
    """A Kilosort-4-shaped output directory."""
    rng = np.random.default_rng(1)
    out = tmp_path / "sorter_output"
    out.mkdir()

    n_spikes = 4000
    times = np.sort(rng.integers(100, N_SAMPLES - 100, n_spikes)).astype(np.int64)
    clusters = rng.integers(0, N_TEMPLATES, n_spikes).astype(np.int32)
    np.save(out / "spike_times.npy", times)
    np.save(out / "spike_clusters.npy", clusters)
    np.save(out / "spike_templates.npy", clusters.astype(np.int64))
    np.save(out / "amplitudes.npy",
            rng.gamma(4, 4, n_spikes).astype(np.float32) + 5)

    depths = np.array([(c + 1) * 20.0 for c in clusters])
    positions = np.column_stack([np.full(n_spikes, 40.0), depths])
    np.save(out / "spike_positions.npy", positions.astype(np.float32))

    # A biphasic template: trough then peak, so trough-to-peak is measurable.
    t = np.arange(61)
    wave = -np.exp(-((t - 20) ** 2) / 8.0) + 0.4 * np.exp(-((t - 32) ** 2) / 20.0)
    templates = np.zeros((N_TEMPLATES, 61, N_CHANNELS), dtype=np.float32)
    for k in range(N_TEMPLATES):
        for offset, weight in ((0, 1.0), (1, 0.5), (-1, 0.5)):
            ch = np.clip(2 * k + offset, 0, N_CHANNELS - 1)
            templates[k, :, ch] += wave * (10 + 4 * k) * weight
    np.save(out / "templates.npy", templates)
    np.save(out / "whitening_mat_inv.npy",
            np.eye(N_CHANNELS, dtype=np.float32))

    pos = np.column_stack([
        np.array([27.0 if i % 2 else 59.0 for i in range(N_CHANNELS)]),
        np.array([(i // 2) * 15.0 for i in range(N_CHANNELS)]),
    ])
    np.save(out / "channel_positions.npy", pos)
    np.save(out / "channel_map.npy", np.arange(N_CHANNELS, dtype=np.int32))

    (out / "params.py").write_text(
        f"n_channels_dat = {N_CHANNELS}\nsample_rate = {FS}\n"
        "dtype = 'int16'\nhp_filtered = False\noffset = 0\n", encoding="utf-8")
    (out / "cluster_KSLabel.tsv").write_text(
        "cluster_id\tKSLabel\n" +
        "".join(f"{k}\tgood\n" for k in range(N_TEMPLATES)), encoding="utf-8")
    return out
