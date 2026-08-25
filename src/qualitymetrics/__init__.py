"""Quality metrics and QC figures for Neuropixels recordings.

Two entry points cover most use:

    from qualitymetrics import KilosortResults, RawRecording
    ks = KilosortResults.load(sorter_output, uv_per_bit=3.0273)
    rec = RawRecording.open(cbin_path)

and then anything in qualitymetrics.plots. For a whole session at once, see
qualitymetrics.report.build_report, or the command line:

    qm report <sorting_dir> --out <figures_dir>
"""
from .ksdata import KilosortResults, KsError, label_from_path
from .metrics import (
                      noise_cutoff,
                      noise_cutoff_per_unit,
                      sliding_rp,
                      write_cluster_tsv,
                      write_phy_metrics,
)
from .raw import Geometry, RawError, RawRecording, find_band, parse_meta
from .style import use_lab_style

__version__ = "2.0.0"

__all__ = [
    "Geometry", "KilosortResults", "KsError", "RawError", "RawRecording",
    "find_band", "label_from_path", "noise_cutoff", "noise_cutoff_per_unit",
    "parse_meta", "sliding_rp", "use_lab_style", "write_cluster_tsv",
    "write_phy_metrics",
]
