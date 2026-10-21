# qualityMetrics

Quality metrics and QC figures for Neuropixels recordings, in Python.

Reads a Kilosort 4 output directory and the archived `.cbin` beside it, and
produces the figure set you look at before trusting a sorting: where the units
are, how they drift, whether the channels are alive, and whether the sorter's
spikes land on real deflections in the raw signal.

## Install

```bash
pip install -e ".[compressed,refractory,dev]"
```

`compressed` pulls mtscomp, which is what reads an archived `.cbin`. `refractory`
pulls slidingRP for the contamination metric. Without either, the package still
works and the figures that needed them report why they were skipped.

## Use

Whole figure set for one shank:

```bash
qm report /path/to/session/sorting/imec0_shank0 --out ./figures
```

It finds the archived recording from the session layout, takes the gain from its
`.meta`, and writes every figure it can build. Anything it cannot build is listed
with the reason, so a missing figure is visible rather than silently absent.

Just the phy curation columns:

```bash
qm phy /path/to/session/sorting/imec0_shank0
```

From Python:

```python
from qualitymetrics import KilosortResults, RawRecording, plots

ks = KilosortResults.load(sorter_output, uv_per_bit=3.0273, label="KM_077 sh0")
rec = RawRecording.open(cbin_path)

plots.unit_drift(ks)
plots.raw_with_spikes(rec, ks, t_start_s=600)
```

## The phy columns

`qm report` and `qm phy` write three `cluster_*.tsv` files into the sorter
output. phy discovers them by filename and shows each as a sortable column:

| Column | Meaning |
| --- | --- |
| `AmpuV` | Peak-to-peak amplitude of the unwhitened template, in microvolts |
| `SlidingCont` | Lowest refractory contamination confirmable at 90% confidence, as a percentage |
| `NoiseCutoff` | How far the amplitude distribution is truncated at the detection floor, in standard deviations |

Two things worth knowing about `SlidingCont`. It is a number rather than a
pass/fail, which is what makes it useful for sorting a cluster table. And
slidingRP caps its search at 35%: a unit for which nothing better than 35% could
be confirmed comes back as NaN, which is a verdict, not a missing measurement.
Those units are written as **100**, a value the estimator can never produce, so
they sort to the worst end instead of vanishing into blanks. On our test shank
that covers 172 of 209 units, every one of which failed, so the distinction is
not academic.

`NoiseCutoff` blanks are genuinely missing: a unit with one spike has no
amplitude distribution to truncate.

## The figures

Sorter output only, so these run on an archived session with no raw data:

| Figure | Question |
| --- | --- |
| `unit_drift` | Did anything move, and did the sorter split a drifting neuron? |
| `drift_map` | The classic drift map, amplitude as darkness |
| `amp_depth_scatter` | Where are the units, how big, how fast, how narrow? |
| `amplitude_cdf_over_depth` | Is any stretch of shank producing spikes worth sorting? |
| `templates_grid` | What do the waveforms look like? |
| `template_waterfall` | How far does one unit's footprint spread? |

Needing the archived recording:

| Figure | Question |
| --- | --- |
| `wall_heatmap` | The whole shank for a few milliseconds: dead channels, artifacts |
| `raw_traces` | What does the signal actually look like? |
| `raw_with_spikes` | Do the sorter's spikes land on real deflections? |
| `sample_waveforms` | Individual snippets, not just the averaged template |
| `channel_health` | Is any channel dead or anomalous? |
| `filter_state` | Is this file wideband, or already high-passed? |
| `band_rms` | Activity over time and depth |
| `depth_power` | Laminar structure in the power spectrum |

## Two things to know about the numbers

**Amplitudes.** The Kilosort-2-era recipe for per-spike amplitude, `amplitudes.npy`
times the unwhitened template peak-to-peak, is wrong for Kilosort 4, whose
templates are not unit-norm. On our data it lands near 1300 uV where the true
median is about 70. `spike_amplitudes_uv` defaults to `template_anchored`, which
rescales each unit's spikes so their mean equals that unit's template
peak-to-peak in microvolts: within-unit variation over time is preserved exactly,
and the absolute scale agrees with the per-unit metric by construction. Pass
`method="kilosort"` for the sorter's own scale instead.

**Channel health.** Two detectors, because a multiplicative threshold (0.1x and
3x the array median) never fires on healthy Neuropixels 2.0 data: the RMS
distribution is tight enough that "noisy" sits far above anything real. A robust
z score against the MAD is reported alongside it, which adapts to however tight
the distribution actually is. Neither is authoritative; the caller decides.

## Conventions

Arial, no top or right spines, editable text in vector output, sentence-case
labels with units on every axis, and firing rates in spikes/s rather than Hz.
Set by `qualitymetrics.style.use_lab_style`, which every figure calls.

## Relation to the MATLAB version

This replaces a MATLAB implementation, preserved on the
[`version1_deprecated`](../../tree/version1_deprecated) branch. If you were using
it, take what you need from that branch and move over rather than maintaining
both.

The rewrite was not cosmetic. The MATLAB figures had, among other things, a time
axis labelled seconds that plotted sample indices; an amplitude axis labelled
microvolts that plotted raw integers, off by roughly a factor of 74 because no
gain was applied; probe names and an artifact-band annotation hardcoded into
every title regardless of what was being plotted; and a hard dependency on a
SpikeInterface scratch file that the sorting pipeline deletes, which made eight
of its twelve figures unrunnable on any archived session.

The choice of figures was good, and is kept. `unit_drift`, whose left panel
reveals a drift split as a depth band that changes colour partway through the
session, is the best idea in the original and is reproduced faithfully.

## Tests

```bash
pytest
```

The suite builds a synthetic sorter output and recording, so it needs no data
server. Most tests pin a specific defect from the MATLAB version: that an axis
labelled seconds carries seconds, that amplitudes are physically plausible, that
no title is hardcoded, that a censored metric does not silently vanish.

## Credits

`noise_cutoff` is vendored from
[ibllib](https://github.com/int-brain-lab/ibllib) (International Brain
Laboratory, MIT licence), unchanged, because importing it would pull ONE, Alyx
and pandas in for thirty lines of numpy.
# Firing continuity (2.2)

The per-shank report now adds `tracking_continuity.png`: all sorted clusters,
candidates first, then appearance/disappearance and selected split time.
Rows are not cluster IDs. Rates are divided by each unit's 95th percentile
(maximum if that is zero), clipped at one; unavailable bins are gray.

This is a descriptive tracking-loss screen, not a significance test and not
a unit-label rule. Ten-second bins are tested at every eligible split, with
at least 600 elapsed seconds and 300 valid seconds per side. Candidates need
an active-side mean of at least 1 spikes/s, quiet/active mean at most 0.1,
and at least 80% of bins on each side respectively above/below 20% of the
active mean. The best qualifying split maximizes two-rate Poisson gain.
Biological rate changes can also qualify. Short or unassessable records are
marked explicitly rather than counted as stable.

Acquisition duration must be provided in `provenance.json` (`duration_s`),
never inferred from the last spike. `qualitymetrics.tracking.write` also
accepts an explicit `duration_s` for externally sorted data. Optional
`light_events.npz` (`onset_s`, `offset_s`, padded by 0.3 s) and
`mask_intervals.csv` (`start_s`, `stop_s`) beside the sort are excluded from
both spike counts and exposure. No sidecars means full-recording rates;
the figure states this. Bins with less than 5 valid seconds are missing.

The same folder receives `tracking_continuity_units.csv` (plot row, cluster
ID, KS label, split, flags and scores), `tracking_continuity_bins.npz`
(cluster IDs, counts, exposure, bin edges and plot order), and
`tracking_continuity.json` (parameters, exclusions and input fingerprints).
The final bin includes the exact acquisition endpoint, as in NumPy's
histogram convention. No source arrays or unit labels are changed.
An endpoint spike follows the final sample's mask, including masks reaching
the end of acquisition.
Outputs are staged together. A failed redraw removes earlier continuity
outputs and records `tracking_continuity_error.json`, so stale results are
not presented as current; locked files that cannot be removed are named in
that failure manifest. Other QC figures continue.
