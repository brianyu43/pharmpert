"""Independent W8 validation from frozen artifacts and source metadata."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from yakseopdong.landscape import _aligned_rows
from yakseopdong.metrics import rmse
from yakseopdong.pseudobulk import _sha256, read_vector_parquet
from yakseopdong.temporal import TIMEPOINTS, parse_timecourse_tag


def validate_temporal(root: Path) -> dict[str, object]:
    """Recompute high-impact W8 counts, deltas, metrics, and leakage invariants."""
    source = root / "data" / "mcfarland_2020.h5ad"
    adata = ad.read_h5ad(source, backed="r")
    try:
        parsed = adata.obs["hash_tag"].astype(object).map(parse_timecourse_tag)
        normal_assigned = int((adata.obs["cell_quality"].eq("normal") & parsed.notna()).sum())
    finally:
        adata.file.close()

    counts = pd.read_csv(root / "results" / "tables" / "timecourse_cell_counts.csv")
    cohort = pd.read_csv(root / "results" / "tables" / "timecourse_cohort.csv")
    pseudobulk_meta, pseudobulk = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_timecourse.parquet", "log1p_cpm"
    )
    response_meta, response = read_vector_parquet(
        root / "data" / "processed" / "response_timecourse.parquet", "delta_log1p_cpm"
    )
    if normal_assigned != 13_713 or int(counts["n_cells"].sum()) != normal_assigned:
        raise ValueError("raw assigned-normal cell count differs from the W8 contract")
    if pseudobulk.shape != (264, 32_738) or response.shape != (120, 32_738):
        raise ValueError("W8 processed matrix shapes differ from the frozen contract")
    if len(cohort) != 24 or int(cohort["eligible_t10"].sum()) != 22:
        raise ValueError("W8 24-line coverage or 22-line primary cohort differs")
    if not all(response_meta.groupby("time_hours").size().reindex(TIMEPOINTS).eq(24)):
        raise ValueError("each time point must contain all 24 response rows")

    primary_lines = sorted(cohort.loc[cohort["eligible_t10"], "cell_line"].astype(str))
    max_delta_error = 0.0
    for hour in TIMEPOINTS:
        control = _aligned_rows(
            pseudobulk_meta,
            pseudobulk,
            primary_lines,
            condition="control",
            time_hours=hour,
        )
        treated = _aligned_rows(
            pseudobulk_meta,
            pseudobulk,
            primary_lines,
            condition="trametinib",
            time_hours=hour,
        )
        observed = _aligned_rows(response_meta, response, primary_lines, time_hours=hour)
        max_delta_error = max(
            max_delta_error,
            float(np.max(np.abs(observed - (treated - control)))),
        )
    if max_delta_error > 2e-6:
        raise ValueError(f"time-matched response reconstruction error is {max_delta_error}")

    train_response_meta, _ = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet", "delta_log1p_cpm"
    )
    train_depmap = set(train_response_meta["depmap_id"].astype(str))
    prediction_meta, predictions = read_vector_parquet(
        root / "results" / "predictions" / "temporal_transfer_predictions.parquet",
        "predicted_delta_log1p_cpm",
    )
    if not prediction_meta["external_test"].all():
        raise ValueError("temporal transfer includes a row not marked external")
    if prediction_meta["depmap_id"].astype(str).isin(train_depmap).any():
        raise ValueError("temporal transfer leaks a training-cohort DepMap ID")
    if prediction_meta["temporal_response_used_for_fit"].any():
        raise ValueError("temporal response was marked as used for fitting")
    external_lines = sorted(prediction_meta["cell_line"].astype(str).unique())
    if len(prediction_meta) != len(external_lines) * len(TIMEPOINTS) * 3:
        raise ValueError("temporal prediction grid is incomplete")

    metrics = pd.read_csv(root / "results" / "tables" / "temporal_transfer_metrics_by_line.csv")
    metric_keys = ["cell_line", "time_hours", "model"]
    metric_lookup = metrics.set_index(metric_keys)["rmse_delta"]
    recomputed_errors: list[float] = []
    for index, row in prediction_meta.iterrows():
        observed = _aligned_rows(
            response_meta,
            response,
            [str(row.cell_line)],
            time_hours=int(row.time_hours),
        )[0]
        recomputed = rmse(observed, predictions[index])
        reported = float(metric_lookup.loc[(row.cell_line, row.time_hours, row.model)])
        recomputed_errors.append(abs(recomputed - reported))
    max_metric_error = float(max(recomputed_errors))
    if max_metric_error > 1e-9:
        raise ValueError(f"reported temporal RMSE differs by {max_metric_error}")

    manifest = pd.read_csv(root / "processed_manifest.csv").set_index("path")
    artifact_hashes: dict[str, str] = {}
    for relative in (
        "data/processed/pseudobulk_timecourse.parquet",
        "data/processed/response_timecourse.parquet",
    ):
        actual = _sha256(root / relative)
        if str(manifest.loc[relative, "sha256"]) != actual:
            raise ValueError(f"processed manifest hash differs for {relative}")
        artifact_hashes[relative] = actual

    required_figures = [
        root / "results" / "figures" / "temporal_component_trajectories.png",
        root / "results" / "figures" / "temporal_heterogeneity_transfer.png",
    ]
    if any(not path.exists() or path.stat().st_size < 20_000 for path in required_figures):
        raise ValueError("a required W8 figure is absent or unexpectedly small")

    report = {
        "stage": "W8_temporal",
        "status": "passed",
        "raw_normal_assigned_cells": normal_assigned,
        "pseudobulk_shape": list(pseudobulk.shape),
        "response_shape": list(response.shape),
        "all_lines": len(cohort),
        "primary_lines_t10": int(cohort["eligible_t10"].sum()),
        "external_transfer_lines": len(external_lines),
        "max_response_reconstruction_error": max_delta_error,
        "max_metric_recomputation_error": max_metric_error,
        "temporal_response_used_for_fit": False,
        "training_depmap_overlap": 0,
        "artifact_hashes": artifact_hashes,
    }
    output = root / "results" / "logs" / "temporal_validation.json"
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
