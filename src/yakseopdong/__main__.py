"""Command-line interface for reproducible project gates."""

from __future__ import annotations

import argparse
import json

from yakseopdong.checks import print_smoke_report, repository_root
from yakseopdong.core_audit import run_core_audit
from yakseopdong.data_probe import run_probe


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
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "smoke":
        return print_smoke_report()
    if args.command == "data-probe":
        report = run_probe(repository_root(), download=args.download)
        print(json.dumps({"ok": True, "shape": report["shape"]}, ensure_ascii=False))
        return 0
    if args.command == "core-audit":
        report = run_core_audit(repository_root())
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 0
    raise RuntimeError(f"unknown command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
