"""Audit the authoritative MIX-seq experiment 3 metadata archives."""

from __future__ import annotations

import hashlib
import json
from io import TextIOWrapper
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

ARCHIVES = {
    "dmso_6h": "DMSO_6hr_expt3.zip",
    "dmso_24h": "DMSO_24hr_expt3.zip",
    "trametinib_24h": "Trametinib_24hr_expt3.zip",
}


def _member_name(archive: ZipFile, suffix: str) -> str:
    matches = [name for name in archive.namelist() if name.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"expected one {suffix} member, found {matches}")
    return matches[0]


def load_archive_metadata(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    """Read classifications/barcodes directly from a source zip and validate alignment."""
    with ZipFile(path) as archive:
        classifications_name = _member_name(archive, "/classifications.csv")
        barcodes_name = _member_name(archive, "/barcodes.tsv")
        genes_name = _member_name(archive, "/genes.tsv")

        with archive.open(classifications_name) as handle:
            classifications = pd.read_csv(handle)
        with archive.open(barcodes_name) as handle:
            barcodes = pd.read_csv(TextIOWrapper(handle, encoding="utf-8"), header=None)[0]
        with archive.open(genes_name) as handle:
            genes_bytes = handle.read()

    required = {"barcode", "singlet_ID", "DepMap_ID", "cell_quality"}
    missing = required - set(classifications.columns)
    if missing:
        raise ValueError(f"missing classification columns: {sorted(missing)}")
    if len(barcodes) != len(classifications):
        raise ValueError("barcode and classification row counts differ")
    if barcodes.astype(str).tolist() != classifications["barcode"].astype(str).tolist():
        raise ValueError("barcode order differs between barcodes.tsv and classifications.csv")

    summary = {
        "archive": path.name,
        "cells_total": int(len(classifications)),
        "cells_normal": int((classifications["cell_quality"] == "normal").sum()),
        "quality_counts": {
            str(key): int(value)
            for key, value in classifications["cell_quality"].value_counts().items()
        },
        "genes_tsv_sha256": hashlib.sha256(genes_bytes).hexdigest(),
    }
    return classifications, summary


def build_cell_count_matrix(
    dmso_6h: pd.DataFrame, dmso_24h: pd.DataFrame, trametinib_24h: pd.DataFrame
) -> pd.DataFrame:
    """Build per-cell-line normal-cell counts for strict and pooled-control cohorts."""
    inputs = {
        "dmso_6h_normal": dmso_6h,
        "dmso_24h_normal": dmso_24h,
        "trametinib_24h_normal": trametinib_24h,
    }
    count_frames: list[pd.DataFrame] = []
    identity_frames: list[pd.DataFrame] = []

    for count_name, frame in inputs.items():
        normal = frame.loc[frame["cell_quality"] == "normal", ["singlet_ID", "DepMap_ID"]]
        count_frames.append(normal["singlet_ID"].value_counts().rename(count_name).to_frame())
        identity_frames.append(normal.drop_duplicates())

    identities = pd.concat(identity_frames, ignore_index=True).drop_duplicates()
    conflicts = identities.groupby("singlet_ID")["DepMap_ID"].nunique()
    if (conflicts > 1).any():
        raise ValueError(f"cell-line to DepMap conflicts: {conflicts[conflicts > 1].to_dict()}")
    depmap = identities.drop_duplicates("singlet_ID").set_index("singlet_ID")["DepMap_ID"]

    counts = pd.concat(count_frames, axis=1).fillna(0).astype(int)
    counts.insert(0, "depmap_id", depmap.reindex(counts.index))
    counts["dmso_pooled_normal"] = counts["dmso_6h_normal"] + counts["dmso_24h_normal"]
    counts["primary_strict_eligible"] = (
        (counts["dmso_24h_normal"] >= 20) & (counts["trametinib_24h_normal"] >= 20)
    )
    counts["reproduction_pooled_eligible"] = (
        (counts["dmso_pooled_normal"] >= 20) & (counts["trametinib_24h_normal"] >= 20)
    )
    return counts.rename_axis("cell_line").reset_index().sort_values("cell_line")


def run_core_audit(root: Path) -> dict[str, object]:
    """Run the core metadata audit and save its tracked count matrix."""
    raw_dir = root / "data" / "raw"
    frames: dict[str, pd.DataFrame] = {}
    archive_summaries: dict[str, dict[str, object]] = {}
    for key, filename in ARCHIVES.items():
        frames[key], archive_summaries[key] = load_archive_metadata(raw_dir / filename)

    gene_hashes = {summary["genes_tsv_sha256"] for summary in archive_summaries.values()}
    if len(gene_hashes) != 1:
        raise ValueError("genes.tsv differs across the three experiment 3 archives")

    counts = build_cell_count_matrix(
        frames["dmso_6h"], frames["dmso_24h"], frames["trametinib_24h"]
    )
    counts.to_csv(root / "cell_count_matrix.csv", index=False)

    report = {
        "archives": archive_summaries,
        "cell_lines_in_normal_union": int(len(counts)),
        "primary_strict_eligible_lines": int(counts["primary_strict_eligible"].sum()),
        "reproduction_pooled_eligible_lines": int(
            counts["reproduction_pooled_eligible"].sum()
        ),
        "genes_tsv_sha256": next(iter(gene_hashes)),
    }
    output = root / "results" / "logs" / "core_audit_summary.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report
