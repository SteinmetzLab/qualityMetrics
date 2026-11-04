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

__version__ = "2.3.0"

#: What this version can do that a caller may need to know before relying on
#: it, checked by name, since a machine can lag behind main:
#:
#: ``"ks_channel_map"``: every figure that reads the raw file for a Kilosort
#: channel reads the column ``channel_map.npy`` names
#: (``KilosortResults.channel_map``), and refuses a map that does not fit the
#: recording. Before 2.3.0 they took Kilosort's channel number for the column,
#: wrong whenever channels were left out of the sort.
#:
#: ``"uv_per_bit_override"``: ``build_report``, ``RawRecording.open`` and
#: ``qm report`` take ``uv_per_bit`` for a meta that states no gain, and
#: refuse one that disagrees with a meta that does.
#:
#: ``"shank_map_geometry"``: site positions from ``~snsShankMap`` where a meta
#: has no ``~snsGeomMap`` (SpikeGLX 20201024), for the Neuropixels 2.0 parts
#: or a meta stating its pitch.
FEATURES = frozenset({"ks_channel_map", "uv_per_bit_override", "shank_map_geometry"})

__all__ = [
    "FEATURES", "Geometry", "KilosortResults", "KsError", "RawError", "RawRecording",
    "find_band", "label_from_path", "noise_cutoff", "noise_cutoff_per_unit",
    "parse_meta", "sliding_rp", "use_lab_style", "write_cluster_tsv",
    "write_phy_metrics",
]


def __getattr__(name):
    # Imported on first access, not at the top: ``style`` imports matplotlib,
    # and the calculations must work without it (``tests/test_compute_api.py``).
    if name == "use_lab_style":
        from .style import use_lab_style

        return use_lab_style
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
