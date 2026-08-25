"""Single-unit quality metrics, and getting them in front of a curator.

Three metrics are computed here:

* amplitude, in microvolts, from the unwhitened template;
* sliding refractory period contamination, via the slidingRP package;
* noise cutoff, vendored from the International Brain Laboratory.

and one export: phy reads any cluster_<name>.tsv in a sorter output directory
and shows it as a sortable column, which is the cheapest way to put a metric in
front of the person doing the curation.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

#: A unit fails noise cutoff if its amplitude distribution is truncated at the
#: detection threshold by more than this many standard deviations.
NOISE_CUTOFF_THRESHOLD = 5.0
NOISE_CUTOFF_PERCENT = 0.10

#: Sliding refractory period defaults, from the slidingRP paper.
RP_CONF_THRESH = 90.0
RP_CONT_THRESH = 10.0
#: Stand-in written to the phy column when slidingRP could not confirm any
#: contamination level within its 35% cap. Deliberately outside the estimator's
#: range so it cannot be mistaken for an estimate, and high so that sorting the
#: column descending puts the uncertifiable units where a curator will see them.
CONTAMINATION_CENSORED = 100.0


# --------------------------------------------------------------------------
# noise cutoff
# --------------------------------------------------------------------------

def noise_cutoff(amps, quantile_length: float = 0.25, n_bins: int = 100,
                 nc_threshold: float = NOISE_CUTOFF_THRESHOLD,
                 percent_threshold: float = NOISE_CUTOFF_PERCENT):
    """Is a unit's amplitude distribution cut off at the detection floor?

    A unit whose low-amplitude spikes fall below the sorter's threshold has a
    truncated amplitude histogram, which means spikes are missing and the unit
    is not what it appears to be. This compares the lowest occupied histogram
    bin against the spread of the bins in the upper part of the distribution, in
    standard deviations, without assuming the distribution is Gaussian.

    Returns (passed, cutoff, first_low_quantile). A unit fails when the cutoff
    exceeds nc_threshold and the lowest bin holds more than percent_threshold of
    the peak bin, so a distribution that merely tapers is not punished for it.

    Vendored verbatim from ibllib.brainbox.metrics.single_units.noise_cutoff
    (International Brain Laboratory, MIT licence), behaviour unchanged.
    """
    amps = np.asarray(amps)
    cutoff = np.float64(np.nan)
    first_low_quantile = np.float64(np.nan)
    fail_criteria = np.ones(1).astype(bool)[0]

    if amps.size > 1:
        bins_list = np.linspace(0, np.max(amps), n_bins)
        n, bins = np.histogram(amps, bins=bins_list)
        idx_peak = np.argmax(n)
        length_top_half = len(np.where(n[idx_peak:-1] > 0)[0])
        high_quantile = 2 * quantile_length
        high_quantile_start_ind = int(
            np.ceil(high_quantile * length_top_half + idx_peak))
        indices_bins_high_quantile = np.arange(high_quantile_start_ind, len(n))
        idx_use = np.where(n[indices_bins_high_quantile] >= 1)[0]

        if len(n[indices_bins_high_quantile]) > 0:
            mean_high_quantile = np.mean(n[indices_bins_high_quantile][idx_use])
            std_high_quantile = np.std(n[indices_bins_high_quantile][idx_use])
            if std_high_quantile > 0:
                first_low_quantile = n[(n != 0)][1]
                cutoff = ((first_low_quantile - mean_high_quantile)
                          / std_high_quantile)
                peak_bin_height = np.max(n)
                percent_of_peak = percent_threshold * peak_bin_height
                fail_criteria = ((cutoff > nc_threshold)
                                 & (first_low_quantile > percent_of_peak))

    nc_pass = ~fail_criteria
    return bool(nc_pass), float(cutoff), float(first_low_quantile)


def noise_cutoff_per_unit(spike_clusters, spike_amplitudes_uv,
                          unit_ids=None) -> dict[int, dict]:
    """Run noise_cutoff for every unit. Amplitudes in microvolts."""
    spike_clusters = np.asarray(spike_clusters)
    amps = np.asarray(spike_amplitudes_uv, dtype=float)
    if unit_ids is None:
        unit_ids = np.unique(spike_clusters)
    out = {}
    for unit in unit_ids:
        mine = amps[spike_clusters == unit]
        passed, cutoff, low = noise_cutoff(mine)
        out[int(unit)] = {"pass": passed, "cutoff": cutoff,
                          "first_low_bin": low, "n_spikes": int(mine.size)}
    return out


# --------------------------------------------------------------------------
# sliding refractory period
# --------------------------------------------------------------------------

def sliding_rp(spike_times_s, spike_clusters,
               conf_thresh: float = RP_CONF_THRESH,
               cont_thresh: float = RP_CONT_THRESH) -> tuple[dict, str]:
    """Sliding refractory period contamination, per unit.

    Returns ({unit_id: {...}}, status). The interesting field is
    min_contamination: the lowest contamination level that can be ruled out at
    the given confidence, as a percentage. It is a continuous number, which
    makes it far more useful in a curation table than the pass/fail bit that
    most pipelines keep, because a curator can sort by it.

    Degrades to an empty dict with a status string rather than raising, so a
    missing optional dependency never costs a finished analysis.
    """
    try:
        from slidingRP.metrics import slidingRP_all
    except Exception as exc:  # noqa: BLE001 - optional dependency
        return {}, f"slidingRP unavailable ({type(exc).__name__}): {exc}"
    try:
        result = slidingRP_all(np.asarray(spike_times_s),
                               np.asarray(spike_clusters),
                               conf_thresh=conf_thresh,
                               cont_thresh=cont_thresh)
    except Exception as exc:  # noqa: BLE001 - report, never fail the caller
        return {}, f"slidingRP failed: {type(exc).__name__}: {exc}"

    out = {}
    ids = np.asarray(result["cidx"])
    for i, uid in enumerate(ids):
        raw = float(result["min_contamination"][i])
        # slidingRP returns NaN when the minimum confirmable contamination
        # exceeds its internal 35% cap (compute_min_contamination's max_contam,
        # which the public API does not expose). That NaN is a verdict, not a
        # missing measurement: it means nothing better than 35% could be ruled
        # out, and on our test shank it covers 172 of 209 units, every one of
        # which failed. Treating it as missing would hide the worst units from
        # anyone sorting by this column, which is precisely backwards.
        censored = not np.isfinite(raw)
        out[int(uid)] = {
            "min_contamination": raw,
            "min_contamination_censored": censored,
            "min_contamination_for_sorting": CONTAMINATION_CENSORED if censored
            else raw,
            "max_confidence": float(result["max_confidence"][i]),
            "rp_min_val_ms": float(result["rp_min_val"][i]) * 1000.0,
            "pass": bool(result["value"][i] == 1),
        }
    return out, "ok"


# --------------------------------------------------------------------------
# phy export
# --------------------------------------------------------------------------

def write_cluster_tsv(sorter_output: str | Path, name: str,
                      values: dict[int, object], fmt: str = "{}") -> Path:
    """Write one cluster_<name>.tsv that phy will show as a column.

    phy discovers these by filename, uses the second column header as the field
    name, and offers it for sorting and filtering in the cluster view. Missing
    units are simply omitted, which phy renders as blank rather than as zero:
    that distinction matters, because a unit with no measurement is not a unit
    that scored nothing.
    """
    sorter_output = Path(sorter_output)
    path = sorter_output / f"cluster_{name}.tsv"
    lines = [f"cluster_id\t{name}"]
    for uid in sorted(values):
        value = values[uid]
        if value is None:
            continue
        if isinstance(value, float) and not np.isfinite(value):
            continue
        lines.append(f"{uid}\t{fmt.format(value)}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_phy_metrics(sorter_output: str | Path, *,
                      amplitude_uv: dict[int, float] | None = None,
                      sliding_cont: dict[int, float] | None = None,
                      noise_cutoff_value: dict[int, float] | None = None
                      ) -> list[Path]:
    """Write the lab's three curation columns: AmpuV, SlidingCont, NoiseCutoff.

    Names are deliberately short: phy shows the column header verbatim in a
    narrow cluster table, and a long name is truncated to uselessness.
    """
    written = []
    if amplitude_uv:
        written.append(write_cluster_tsv(sorter_output, "AmpuV",
                                         amplitude_uv, "{:.1f}"))
    if sliding_cont:
        written.append(write_cluster_tsv(sorter_output, "SlidingCont",
                                         sliding_cont, "{:.2f}"))
    if noise_cutoff_value:
        written.append(write_cluster_tsv(sorter_output, "NoiseCutoff",
                                         noise_cutoff_value, "{:.2f}"))
    return written


def phy_metrics_from_results(ks, write: bool = True) -> dict[str, dict]:
    """Compute the three curation columns from a KilosortResults and write them.

    Convenience for retrofitting sorts that finished before the pipeline learned
    to write these files.
    """
    amps = ks.unit_amplitude_uv
    rp, rp_status = sliding_rp(ks.spike_times_s, ks.spike_clusters)
    cont = {u: v["min_contamination_for_sorting"] for u, v in rp.items()}

    stored = ks.metrics_by_unit()
    if stored and all("noise_cutoff" in v for v in stored.values()):
        nc = {u: v["noise_cutoff"] for u, v in stored.items()}
    else:
        computed = noise_cutoff_per_unit(
            ks.spike_clusters, ks.spike_amplitudes_uv())
        nc = {u: v["cutoff"] for u, v in computed.items()}

    if write:
        write_phy_metrics(ks.path, amplitude_uv=amps, sliding_cont=cont,
                          noise_cutoff_value=nc)
    return {"AmpuV": amps, "SlidingCont": cont, "NoiseCutoff": nc,
            "sliding_rp_status": rp_status}
