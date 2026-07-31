# Session 2026-07-16 — raw-activity + spike-overlay figure (`plotRawWithSpikes`)

Built a MATLAB tool that recreates the IBL "good raw data" QC figure for the
`probe00a` sort: a grayscale channels×time image of baselined AP activity with
Kilosort spikes overlaid (green = spike from a `good` unit, red = `mua`). Then
added a manual mask to exclude the top artifact band.

- **Deliverable:** `C:\Users\stein\alex\autoqc\qc\plotRawWithSpikes.m`
  (started life as `sandbox\analysis\plot_raw_with_spikes.m`; user moved it into
  the new `autoqc` project tree and renamed it to the repo's camelCase convention).
- **Figure output (function's own `print`):** `C:\Users\stein\alex\autoqc\figures\rawWithSpikes.png`
- **Data source:** the finished sort `C:\spikesort_scratch\probe00a_full_ks4\`.
- **Env:** MATLAB **R2025b** (`C:\Program Files\MATLAB\R2025b\bin\matlab.exe`).
  Batch-run headless: `matlab -batch "run('<script.m>')"`. No toolboxes required.

---

## What the figure shows / how it's built

For a chosen time window (default 40 ms at T=1800 s):
1. Read that window from `recording.dat` (never the whole 143 GB file — `fseek`
   to the sample offset, then `fread [nChan × nSamp]`).
2. Convert to µV (× gain 3.02734375) and **baseline each channel** (subtract the
   per-channel median over the window) → quiescence = mid-gray, negative
   deflection = black, positive = white. Displayed via `imagesc` + `colormap(gray)`
   with symmetric `clim = ±clim_uV`.
3. Sort channel rows by **depth** (the file's row order is interleaved, NOT
   depth-sorted — z goes 0, 2880, 15, 2895, …).
4. Overlay every KS4 spike in the window at its **unit's peak channel** (from
   template peak-to-peak argmax), colored green/red by `cluster_group.tsv` label.
   mua plotted first, good on top so good markers win visually.

### Key data facts this relies on (verified this session)
- `recording.dat`: **int16, 383 ch, sample-major** (all 383 channels contiguous
  per time sample), offset 0, fs 30000, 186,752,122 samples (~6225 s). Already
  hp-filtered + CMR'd + drift-corrected (it's phy's display source).
  MATLAB `fread(fid,[nChan,nS],'*int16')` fills column-major = one sample per
  column = correct for this layout.
- `channel_map.npy` is **identity 0–382**; `channel_positions.npy` x∈{0,32} µm,
  z 0–5730 µm, one contact per depth after masking → clean monotonic y-axis once
  sorted.
- Per-cluster **peak channel** = argmax over channels of `ptp(templates[u])`;
  `templates.npy` row index == cluster_id for raw KS4 output.
- Labels (`cluster_group.tsv`): **160 good, 365 mua** (of 525 clusters).
- Self-contained `.npy` reader (`readNPY`, also now `autoqc\io\readNPY.m`):
  parses the npy header, handles C-order + int/uint/float/bool, reshapes
  reversed-dims + permute to get numpy shape into MATLAB column-major.

---

## The top-band artifact mask (main addition this session)

Per the 2026-07-15 revision, **z > 4850 µm (59 sites) is artifact/noise**, but the
existing 383-ch sort still includes it (the deferred re-sort hasn't run). Added a
display/overlay mask mirroring the LFP pipeline's `load_lf_masked` convention:

- **New param `z_max` (default 4850; `Inf` disables).**
- Channels filtered to `z <= z_max` *before* the depth-sort
  (`keepChan → idxKeep → sort(z(idxKeep))`), so only real rows are drawn.
- Spikes whose peak channel is in the masked band are dropped
  (`keepSpk = keepChan(spk_chan)`).
- Report line prints the effect, e.g.:
  `324/383 channels kept (z<=4850), 17 spikes dropped in artifact band`.
- Also removed the old `rank` roundtrip — spike depth is now just `z(spk_chan)`
  (it was computing the same value the long way).

**Verified on the T=1800 s / 40 ms window:** 181 → **164 spikes** (61 good, 103
mua), **324/383 channels kept**, 17 artifact-band spikes dropped. The striped
artifact band is visibly gone and the y-axis tops out at the last real channel.

> **Important scope note:** this is a *display/overlay* mask only. The underlying
> sort is unchanged, so a unit whose peak sits just below 4850 µm still carries
> template footprint energy from the artifact band. The real fix remains the
> deferred re-sort on the 324-ch set (see `SKILL_probe00a.md` "Deferred items").

---

## User edits folded in (their versions kept)

- `clim_uV` default **100** (was 30); `markersize` default **10** (was 30) — reads
  much closer to the reference figure (spike footprints visible under markers).
- Left a **commented-out dF/F baseline block** with a note that they'd prefer a
  proper *reference window* normalization rather than in-window baseline
  ("not task specific so I guess it's ok" for now). **Open question, not resolved.**
- Added a `print(saveDir,'-dpng')` writing to `autoqc\figures\rawWithSpikes.png`.

---

## Open / follow-up items

- **Reference-window normalization** for the baseline (user's flagged preference)
  vs. the current per-window median subtraction — not yet decided/implemented.
- **`print` filename is fixed** (`rawWithSpikes.png`) → every call overwrites the
  previous figure. Fold `t_start_s` / `z_max` into the filename if batching windows.
- **Per-spike depth option:** markers currently sit at the unit's static peak
  channel; `spike_positions.npy` (per-spike x, depth) would scatter them vertically
  more like the IBL figure but is noisier. Left as peak-channel for robustness.
- The **deferred 324-ch re-sort** is the actual data fix behind the display mask.

## `autoqc` project layout (as of this session)
User has reorganized into a MATLAB project: `io/` (readNPY, loadRawRecording,
getTraces, recordingGeometry), `preproc/` (preprocessAP, preprocessLFP), `qc/`
(plotRawWithSpikes, checkChannels, checkFilterState), `spikes/` (cortex-lab
spikes toolbox), `viz/`, `figures/`.
