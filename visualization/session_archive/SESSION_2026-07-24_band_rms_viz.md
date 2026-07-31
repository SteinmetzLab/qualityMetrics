# Session 2026-07-24 — time×depth RMS-activity map for AP + LFP bands (`plotBandRMS`)

Built a QC figure that shows **RMS voltage binned by time and by depth along the probe**,
for the AP and LFP bands, as a viridis heatmap (**purple = low activity, yellow = high**,
per the user's explicit colour request). Fills the empty "RMS over time" slot in
`autoqc_working.m`.

- **Deliverable:** `C:\Users\stein\alex\autoqc\qc\plotBandRMS.m`
- **Figure output:** `C:\Users\stein\alex\autoqc\figures\bandRMS\` (PNG when `SaveToggle`).
  Smoke-test render this session: `bandRMS_both_t1000_w10.png`.
- **Data sources:** `C:\spikesort_scratch\probe00a_pre` (AP, 30 kHz, 383 ch) and
  `C:\spikesort_scratch\probe00a_lf_pre` (LFP, 2500 Hz, 324 ch) — both already CMR'd.
- **Env:** MATLAB **R2025b**, no toolboxes required. Uses autoqc infra
  (`loadRawRecording`, `getTraces`, `recordingGeometry`, `viz\viridisColors.m`).

---

## What the figure shows

`'Band'` selects the layout: **`'both'`** (default) stacks **AP (top) over LFP (bottom)**
in a `tiledlayout(2,1)` sharing the time x-axis; `'ap'` / `'lfp'` draw one panel.

- **y = depth along probe (µm)**, `YDir` normal (deep at bottom, up = toward tip surface).
- **x = time (s)**.
- **colour = RMS voltage (µV)**, **viridis** (`viridisColors(256)`), each band with its
  **own** colour limits (AP and LFP live on very different scales).
- Empty depth bins render **transparent** (`AlphaData`) over a gray axis.

### Binning + RMS definition
- **Depth bins:** channel sites grouped into `depthBin_um`-wide bins (default 40 µm) via
  `discretize`; bin centre is the y-coordinate.
- **Time bins:** the window is split into time bins; RMS is accumulated **one time bin at
  a time** (read → reduce → discard) so long spans stay memory-bounded even at 30 kHz.
- **RMS:** within each (depth-bin, time-bin) block, per-channel DC is removed (mean over
  that bin), then `rms = sqrt(mean(x.^2))` is taken over **all samples of all channels
  pooled in the block**. Implemented as `accumarray` sum-of-squares per depth bin divided
  by (nCh·nSamp), then sqrt — exact, not an average-of-per-channel-RMS.
- **`z_max = 4850` mask** (the probe00a convention) drops the top artifact band before
  binning → AP contributes 324/383 ch. Display-time exclusion only (same scope caveat as
  the rest of the toolset — the underlying 383-ch sort still includes that band).

### What the smoke test showed (T = 1000 s, physiologically sensible)
- **LFP RMS hotspot at ~2400–3050 µm** — coincides with the broadband power hotspot the
  Morlet (`plotDepthPower`) and Welch (`plotLFPWelch`) tools independently flagged at the
  same depth → good cross-validation.
- **AP RMS highest ~500–2000 µm.**
- Value ranges: AP ≈ 9–42 µV RMS, LFP ≈ 16–235 µV RMS.

---

## Interface (current state, incl. user edits folded in)

**Signature changed by the user to take loaded recordings positionally**, so the driver
can memory-map once and reuse:

```matlab
info = plotBandRMS(apRec, lfpRec, Name,Value, ...)
```

- `apRec`, `lfpRec` = structs from `loadRawRecording` (positional).
- `computeBandRMS(rec, o)` calls `recordingGeometry(rec)` on the rec struct directly.
- **Whole-recording scan:** time axis now starts at `t0 = 0` and uses a **fixed
  `nBin = 100`** time bins spanning the full recording duration
  (`timebinSize_s = (nSamples/fs)/100`), rather than a user-specified `t_start_s`/`win_s`
  window. (The original draft took `apDir`/`lfpDir` + `t_start_s`/`win_s`; those params
  were removed in favour of this.)

**Name/Value options that remain:** `'Band'`, `'timeBin_s'`, `'depthBin_um'`, `'z_max'`,
`'CAR'`, `'climPct'`, `'clim'` (scalar `[lo hi]` for all, or `2x2 [apLo apHi; lfpLo lfpHi]`
for `'both'`), `'SaveToggle'`.

**Output `info`:** one field per band drawn, each a struct with `.depth` (bin centres),
`.t` (bin centres, s), `.rms` ([nDepth×nTime], NaN for empty bins), `.fs`, `.nKept`,
`.nChan`.

Self-contained helpers in the file: `computeBandRMS`, `resolveClim`, `prctile_` (no
Statistics Toolbox).

---

## Open / follow-up items
- **Save block references removed params.** The `SaveToggle` filename still uses
  `o.t_start_s` / `o.win_s`, which no longer exist after the switch to the whole-recording
  `nBin=100` scan → `print` path will error if `SaveToggle=true`. Needs the filename
  updated (e.g. `bandRMS_%s_%dbins.png` with `bands`, `nBin`) before saving works again.
  The smoke-test PNG above was rendered from the earlier windowed draft.
- **`nBin` is hardcoded to 100** and `timeBin_s` is currently unused by the time-axis
  math (bin width is derived from `nBin`, not from `timeBin_s`). Decide whether the knob
  should be bin **count** or bin **width** and wire one through consistently.
- **`clim` on the whole-recording scan** uses robust `[2 98]` percentiles across all bins;
  fine, but a shared/fixed scale may be wanted if comparing sessions.
- Artifact band is **display-only** — real fix is the deferred 324-ch re-sort
  (`SKILL_probe00a.md` "Deferred items").
- Not yet wired into `autoqc_working.m` line ~59 ("RMS over time") — offered.
