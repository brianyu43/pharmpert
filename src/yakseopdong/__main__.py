"""Command-line interface for reproducible project gates."""

from __future__ import annotations

import argparse
import json

from yakseopdong.checks import print_smoke_report, repository_root


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="yakseopdong")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("smoke", help="validate the Stage 0 repository and environment")
    probe = subparsers.add_parser("data-probe", help="load and summarize the McFarland dataset")
    probe.add_argument(
        "--download",
        action="store_true",
        help="explicitly allow the pertpy loader to download/cache the dataset",
    )
    subparsers.add_parser("core-audit", help="audit authoritative experiment 3 metadata")
    subparsers.add_parser(
        "build-pseudobulk",
        help="build experiment 3 pseudobulk matrices and QC evidence",
    )
    subparsers.add_parser("metadata-audit", help="join author cell-line annotations")
    subparsers.add_parser("landscape", help="build exploratory control/response PCA")
    subparsers.add_parser("build-splits", help="freeze nested cell-line splits")
    subparsers.add_parser("run-baselines", help="run nested-CV B0-B4 baselines")
    subparsers.add_parser("run-cclr", help="run nested-CV CCLR main model")
    subparsers.add_parser("validate-cclr", help="independently validate W6 outputs")
    subparsers.add_parser("run-ablation", help="run frozen W7 diagnostic ablations")
    subparsers.add_parser("validate-ablation", help="independently validate W7 outputs")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "smoke":
        return print_smoke_report()
    if args.command == "data-probe":
        from yakseopdong.data_probe import run_probe

        report = run_probe(repository_root(), download=args.download)
        print(json.dumps({"ok": True, "shape": report["shape"]}, ensure_ascii=False))
        return 0
    if args.command == "core-audit":
        from yakseopdong.core_audit import run_core_audit

        report = run_core_audit(repository_root())
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 0
    if args.command == "build-pseudobulk":
        from yakseopdong.pseudobulk import run_pseudobulk

        report = run_pseudobulk(repository_root())
        print(
            json.dumps(
                {
                    "ok": True,
                    "primary_lines": report["primary_lines"],
                    "sensitivity_lines": report["sensitivity_lines"],
                    "gene_count": report["gene_count"],
                },
                ensure_ascii=False,
            )
        )
        return 0
    if args.command == "metadata-audit":
        from yakseopdong.metadata import run_metadata_audit

        print(json.dumps(run_metadata_audit(repository_root()), indent=2))
        return 0
    if args.command == "landscape":
        from yakseopdong.landscape import run_landscape

        print(json.dumps(run_landscape(repository_root()), indent=2))
        return 0
    if args.command == "build-splits":
        from yakseopdong.splits import run_splits

        print(json.dumps(run_splits(repository_root()), indent=2))
        return 0
    if args.command == "run-baselines":
        from yakseopdong.benchmark import run_baselines

        report = run_baselines(repository_root())
        print(json.dumps(report, indent=2))
        return 0
    if args.command == "run-cclr":
        from yakseopdong.cclr import run_cclr

        report = run_cclr(repository_root())
        print(json.dumps(report, indent=2))
        return 0
    if args.command == "validate-cclr":
        from yakseopdong.cclr_validation import validate_cclr

        report = validate_cclr(repository_root())
        print(json.dumps(report, indent=2))
        return 0
    if args.command == "run-ablation":
        from yakseopdong.ablation import run_ablation

        report = run_ablation(repository_root())
        print(json.dumps(report, indent=2))
        return 0
    if args.command == "validate-ablation":
        from yakseopdong.ablation_validation import validate_ablation

        report = validate_ablation(repository_root())
        print(json.dumps(report, indent=2))
        return 0
    raise RuntimeError(f"unknown command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
