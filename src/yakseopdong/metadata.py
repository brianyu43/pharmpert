"""Authoritative cell-line annotations for the experiment 3 cohort."""

from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import pandas as pd
import rdata

from yakseopdong.pseudobulk import read_vector_parquet

RDS_FILE = "all_CL_features.rds"
RDS_OBJECT = "Trametinib_24hr_expt3"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_author_annotations(root: Path) -> pd.DataFrame:
    """Join the frozen 94-line cohort to the authors' sensitivity/omics RDS."""
    response_meta, _ = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet",
        "delta_log1p_cpm",
    )
    cohort = response_meta[["cell_line", "depmap_id"]].copy()
    if len(cohort) != 94 or not cohort["cell_line"].is_unique:
        raise ValueError("expected the frozen 94-line response cohort")

    rds_path = root / "data" / "raw" / RDS_FILE
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Missing constructor for R class")
        objects = rdata.read_rds(rds_path)

    features = objects[RDS_OBJECT].copy()
    metadata = objects["metadata"].copy()
    features = features.rename(
        columns={
            "DEPMAP_ID": "depmap_id",
            "CCLE_ID": "author_ccle_id",
            "AUC_avg": "trametinib_auc",
            "sens": "trametinib_sensitivity",
            "PRISM_AUC": "prism_auc",
            "GDSC_AUC": "gdsc_auc",
        }
    )
    metadata = metadata.rename(
        columns={
            "DEPMAP_ID": "depmap_id",
            "Disease": "lineage",
            "Subtype": "lineage_subtype",
        }
    )

    feature_columns = [
        "depmap_id",
        "author_ccle_id",
        "trametinib_auc",
        "trametinib_sensitivity",
        "prism_auc",
        "gdsc_auc",
        "BRAF_MUT",
        "KRAS_MUT",
        "HRAS_MUT",
        "NRAS_MUT",
        "in_pool",
    ]
    annotation_columns = ["depmap_id", "lineage", "lineage_subtype"]
    features = features[feature_columns].drop_duplicates("depmap_id")
    metadata = metadata[annotation_columns].drop_duplicates("depmap_id")
    annotations = cohort.merge(features, on="depmap_id", how="left", validate="one_to_one")
    annotations = annotations.merge(metadata, on="depmap_id", how="left", validate="one_to_one")

    required = [
        "author_ccle_id",
        "trametinib_auc",
        "trametinib_sensitivity",
        "lineage",
        "BRAF_MUT",
        "KRAS_MUT",
        "HRAS_MUT",
        "NRAS_MUT",
    ]
    missing = annotations[required].isna().sum()
    if missing.any():
        raise ValueError(f"missing required author annotations: {missing[missing > 0].to_dict()}")
    if not annotations["in_pool"].fillna(False).all():
        raise ValueError("a strict-cohort cell line is not marked as experiment 3 in-pool")
    if annotations["lineage"].nunique() != 21:
        raise ValueError("expected 21 author-defined disease lineages")

    for source, target in [
        ("BRAF_MUT", "braf_mut"),
        ("KRAS_MUT", "kras_mut"),
        ("HRAS_MUT", "hras_mut"),
        ("NRAS_MUT", "nras_mut"),
    ]:
        annotations[target] = annotations[source].gt(0)
    annotations["lineage_subtype"] = annotations["lineage_subtype"].fillna("unspecified")
    annotations["annotation_source"] = "Figshare all_CL_features.rds v3"
    annotations["sensitivity_use"] = "interpretation_only_not_model_input"

    return annotations[
        [
            "cell_line",
            "depmap_id",
            "author_ccle_id",
            "lineage",
            "lineage_subtype",
            "trametinib_auc",
            "trametinib_sensitivity",
            "prism_auc",
            "gdsc_auc",
            "braf_mut",
            "kras_mut",
            "hras_mut",
            "nras_mut",
            "annotation_source",
            "sensitivity_use",
        ]
    ].sort_values("cell_line", ignore_index=True)


def run_metadata_audit(root: Path) -> dict[str, object]:
    """Write the joined annotations and a compact provenance/coverage audit."""
    annotations = read_author_annotations(root)
    output = root / "cell_line_annotations.csv"
    annotations.to_csv(output, index=False)
    report = {
        "source": "https://figshare.com/articles/dataset/MIX-seq_data/10298696",
        "source_file": RDS_FILE,
        "source_sha256": _sha256(root / "data" / "raw" / RDS_FILE),
        "rds_object": RDS_OBJECT,
        "cell_lines": int(len(annotations)),
        "depmap_ids": int(annotations["depmap_id"].nunique()),
        "lineages": int(annotations["lineage"].nunique()),
        "sensitivity_complete": bool(annotations["trametinib_sensitivity"].notna().all()),
        "mutation_complete": bool(
            annotations[["braf_mut", "kras_mut", "hras_mut", "nras_mut"]]
            .notna()
            .all()
            .all()
        ),
        "model_inputs": ["control_24h_log1p_cpm", "training_lineage_for_B2_only"],
        "interpretation_only": [
            "trametinib_sensitivity",
            "trametinib_auc",
            "braf_mut",
            "kras_mut",
            "hras_mut",
            "nras_mut",
        ],
    }
    log_path = root / "results" / "logs" / "metadata_audit.json"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
