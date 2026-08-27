#!/usr/bin/env python3
"""Freeze the external MSigDB Hallmark sets used by the W7 ablation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import yaml

ENDPOINT = "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp"
SOURCE_SETS = {
    "HALLMARK_KRAS_SIGNALING_UP": "MAPK_KRAS_signaling",
    "HALLMARK_KRAS_SIGNALING_DN": "MAPK_KRAS_signaling",
    "HALLMARK_E2F_TARGETS": "E2F_targets",
    "HALLMARK_G2M_CHECKPOINT": "G2M_checkpoint",
    "HALLMARK_APOPTOSIS": "apoptosis",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": (
        "epithelial_mesenchymal_transition"
    ),
}
MANUAL_IMMEDIATE_EARLY = [
    "EGR1",
    "ETV4",
    "ETV5",
    "DUSP4",
    "DUSP5",
    "DUSP6",
    "SPRY2",
    "SPRY4",
]


def fetch_set(name: str) -> dict[str, object]:
    """Download and validate one official MSigDB JSON gene set."""
    url = f"{ENDPOINT}?{urlencode({'fileType': 'json', 'geneSetName': name})}"
    with urlopen(url, timeout=60) as response:  # noqa: S310 - fixed official host
        payload = json.load(response)
    record = payload.get(name)
    if not isinstance(record, dict):
        raise ValueError(f"MSigDB response is missing {name}")
    symbols = sorted({str(value) for value in record.get("geneSymbols", []) if value})
    if len(symbols) < 8:
        raise ValueError(f"MSigDB response for {name} has too few genes")
    return {
        "systematic_name": str(record.get("systematicName", "")),
        "pmid": str(record.get("pmid", "")),
        "gene_count": len(symbols),
        "genes": symbols,
    }


def build_config() -> dict[str, object]:
    """Build the deterministic W7 pathway configuration."""
    downloaded = {name: fetch_set(name) for name in SOURCE_SETS}
    collections: dict[str, dict[str, object]] = {}
    for source_name, logical_name in SOURCE_SETS.items():
        collection = collections.setdefault(logical_name, {"source_sets": [], "genes": []})
        collection["source_sets"].append(source_name)
        collection["genes"].extend(downloaded[source_name]["genes"])
    for collection in collections.values():
        collection["source_sets"] = sorted(collection["source_sets"])
        collection["genes"] = sorted(set(collection["genes"]))
        collection["gene_count"] = len(collection["genes"])
    collections["immediate_early_response"] = {
        "source_sets": ["manual_predefined_markers"],
        "gene_count": len(MANUAL_IMMEDIATE_EARLY),
        "genes": MANUAL_IMMEDIATE_EARLY,
    }
    return {
        "status": "frozen",
        "frozen_for": "W7_ablation_before_outer_test_evaluation",
        "source": {
            "name": "Molecular Signatures Database Hallmark collection",
            "version": "2026.1.Hs",
            "release_date": "2026-01",
            "url": "https://www.gsea-msigdb.org/gsea/msigdb/",
            "download_endpoint": ENDPOINT,
            "license": "CC BY 4.0",
            "license_url": "https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp",
        },
        "source_set_metadata": downloaded,
        "collections": collections,
        "input_panel_collections": sorted(collections),
        "input_panel_rule": (
            "union_then_dataset_intersection_then_outer_training_only_expression_filter"
        ),
        "target_gene_rule": "all_32738_genes_unchanged",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("config/genesets.yaml"),
        help="path for the deterministic YAML snapshot",
    )
    args = parser.parse_args()
    config = build_config()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        yaml.safe_dump(config, sort_keys=False, allow_unicode=True, width=100),
        encoding="utf-8",
    )
    print(f"wrote {args.output} with {len(config['source_set_metadata'])} MSigDB sets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
