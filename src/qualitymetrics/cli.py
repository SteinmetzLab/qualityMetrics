"""Command line: qm report, qm phy, qm info."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _add_report(sub):
    p = sub.add_parser("report", help="build the QC figure set for one shank")
    p.add_argument("sorting_dir", help="shank directory, or a sorter_output")
    p.add_argument("--out", required=True, help="directory for the figures")
    p.add_argument("--raw", default=None,
                   help="archived .cbin/.bin (found automatically if omitted)")
    p.add_argument("--lf", default=None, help="archived LF band recording")
    p.add_argument("--artifact-z", type=float, default=None,
                   help="shade everything above this depth in micrometres")
    p.add_argument("--t-start", type=float, default=None,
                   help="window start for the raw views, in seconds")
    p.add_argument("--no-phy", action="store_true",
                   help="do not write the cluster_*.tsv curation columns")


def _add_phy(sub):
    p = sub.add_parser("phy", help="write the cluster_*.tsv columns only")
    p.add_argument("sorting_dir")
    p.add_argument("--uv-per-bit", type=float, default=None,
                   help="override the gain instead of reading it from the meta")


def _add_info(sub):
    p = sub.add_parser("info", help="describe a recording or a sorter output")
    p.add_argument("path")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="qm", description="Neuropixels quality metrics and QC figures")
    sub = parser.add_subparsers(dest="command", required=True)
    _add_report(sub)
    _add_phy(sub)
    _add_info(sub)
    args = parser.parse_args(argv)

    if args.command == "report":
        from .report import build_report
        result = build_report(
            args.sorting_dir, args.out, raw_path=args.raw, lf_path=args.lf,
            artifact_z_um=args.artifact_z, t_start_s=args.t_start,
            write_phy=not args.no_phy)
        print(result.summary())
        return 0

    if args.command == "phy":
        from .ksdata import KilosortResults, label_from_path
        from .metrics import phy_metrics_from_results
        from .raw import parse_meta, uv_per_bit
        from .report import find_recording, find_sorter_output

        sorter_output = find_sorter_output(args.sorting_dir)
        scale = args.uv_per_bit
        if scale is None:
            raw = find_recording(args.sorting_dir, "ap")
            if raw is None:
                print("cannot find the recording, so the gain is unknown; "
                      "pass --uv-per-bit", file=sys.stderr)
                return 2
            scale = uv_per_bit(parse_meta(Path(raw).with_suffix(".meta")), "ap")
        ks = KilosortResults.load(sorter_output, uv_per_bit=scale,
                                  label=label_from_path(args.sorting_dir))
        out = phy_metrics_from_results(ks, write=True)
        print(f"wrote cluster_AmpuV.tsv, cluster_SlidingCont.tsv and "
              f"cluster_NoiseCutoff.tsv for {len(out['AmpuV'])} units "
              f"into {sorter_output}")
        if out["sliding_rp_status"] != "ok":
            print(f"note: {out['sliding_rp_status']}", file=sys.stderr)
        return 0

    if args.command == "info":
        path = Path(args.path)
        if path.suffix in (".cbin", ".bin", ".raw"):
            from .raw import RawRecording
            with RawRecording.open(path) as rec:
                print(rec.describe())
        else:
            from .ksdata import KilosortResults
            from .report import find_sorter_output
            ks = KilosortResults.load(find_sorter_output(path))
            print(f"{ks.path}: {len(ks.unit_ids)} units, "
                  f"{ks.spike_times_s.size:,} spikes, "
                  f"{ks.duration_s / 60:.1f} min at {ks.fs / 1000:.0f} kHz")
        return 0

    parser.error(f"unknown command {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
