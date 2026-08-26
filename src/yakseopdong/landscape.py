"""Exploratory control/response landscape and perturbation-direction diagnostics."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from yakseopdong.metadata import run_metadata_audit
from yakseopdong.metrics import pearson_or_nan, spearman_or_nan
from yakseopdong.plots import write_landscape_figures
from yakseopdong.pseudobulk import read_vector_parquet


def _aligned_rows(
    metadata: pd.DataFrame,
    matrix: np.ndarray,
    cell_lines: list[str],
    **filters: object,
) -> np.ndarray:
    mask = np.ones(len(metadata), dtype=bool)
    for column, value in filters.items():
        mask &= metadata[column].eq(value).to_numpy()
    selected = metadata.loc[mask, ["cell_line"]].copy()
    selected["row_index"] = np.flatnonzero(mask)
    if selected["cell_line"].duplicated().any():
        raise ValueError(f"duplicate rows for filters {filters}")
    positions = selected.set_index("cell_line")["row_index"]
    missing = sorted(set(cell_lines) - set(positions.index))
    if missing:
        raise ValueError(f"missing aligned rows for filters {filters}: {missing}")
    return np.asarray(matrix)[positions.loc[cell_lines].to_numpy()]


def _descriptive_pca(
    matrix: np.ndarray, seed: int, max_genes: int = 5_000
) -> tuple[np.ndarray, PCA, int]:
    variances = np.asarray(matrix, dtype=np.float64).var(axis=0)
    eligible = np.flatnonzero(variances > 1e-10)
    order = np.argsort(-variances[eligible], kind="stable")
    selected = eligible[order[:max_genes]]
    pca = PCA(n_components=2, svd_solver="randomized", random_state=seed)
    scores = pca.fit_transform(np.asarray(matrix)[:, selected])
    return scores, pca, int(len(selected))


def run_landscape(root: Path, seed: int = 20260827) -> dict[str, object]:
    """Write descriptive PCA tables and control-vs-drug direction diagnostics."""
    run_metadata_audit(root)
    annotations = pd.read_csv(root / "cell_line_annotations.csv").sort_values(
        "cell_line", ignore_index=True
    )
    cell_lines = annotations["cell_line"].astype(str).tolist()

    primary_meta, primary = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_24h.parquet", "log1p_cpm"
    )
    response_meta, response_values = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet", "delta_log1p_cpm"
    )
    control_time_meta, control_time = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_control_time.parquet", "log1p_cpm"
    )
    control_24h = _aligned_rows(
        primary_meta, primary, cell_lines, condition="control", time_hours=24
    )
    response = _aligned_rows(response_meta, response_values, cell_lines)
    control_6h = _aligned_rows(
        control_time_meta, control_time, cell_lines, condition="control", time_hours=6
    )
    control_24h_again = _aligned_rows(
        control_time_meta, control_time, cell_lines, condition="control", time_hours=24
    )
    if not np.allclose(control_24h, control_24h_again, rtol=0, atol=0):
        raise ValueError("DMSO 24h values differ between processed artifacts")

    control_scores, control_pca, control_genes = _descriptive_pca(control_24h, seed)
    response_scores, response_pca, response_genes = _descriptive_pca(response, seed)
    landscape = annotations.copy()
    landscape["control_pc1"] = control_scores[:, 0]
    landscape["control_pc2"] = control_scores[:, 1]
    landscape["response_pc1"] = response_scores[:, 0]
    landscape["response_pc2"] = response_scores[:, 1]
    landscape.to_csv(root / "results" / "tables" / "response_landscape.csv", index=False)

    control_change = control_24h - control_6h
    comparison_rows: list[dict[str, float | str]] = []
    for index, cell_line in enumerate(cell_lines):
        control_vector = control_change[index]
        drug_vector = response[index]
        denominator = float(np.linalg.norm(control_vector) * np.linalg.norm(drug_vector))
        cosine = (
            float(np.dot(control_vector, drug_vector) / denominator)
            if denominator
            else float("nan")
        )
        comparison_rows.append(
            {
                "cell_line": cell_line,
                "control_6h_to_24h_rmse": float(np.sqrt(np.mean(np.square(control_vector)))),
                "drug_24h_response_rmse": float(np.sqrt(np.mean(np.square(drug_vector)))),
                "control_vs_drug_pcc": pearson_or_nan(control_vector, drug_vector),
                "control_vs_drug_cosine": cosine,
            }
        )
    directions = pd.DataFrame(comparison_rows).merge(
        annotations[["cell_line", "lineage"]], on="cell_line", validate="one_to_one"
    )
    directions.to_csv(root / "results" / "tables" / "direction_comparison.csv", index=False)

    sensitivity = annotations["trametinib_sensitivity"].to_numpy(dtype=float)
    pc1_pearson = pearson_or_nan(sensitivity, response_scores[:, 0])
    pc1_spearman = spearman_or_nan(sensitivity, response_scores[:, 0])
    report = {
        "status": "exploratory_full_cohort_not_used_for_model_selection",
        "cell_lines": len(cell_lines),
        "control_pca_genes": control_genes,
        "response_pca_genes": response_genes,
        "control_pc_explained_variance_ratio": control_pca.explained_variance_ratio_.tolist(),
        "response_pc_explained_variance_ratio": response_pca.explained_variance_ratio_.tolist(),
        "response_pc1_sensitivity_pearson": pc1_pearson,
        "response_pc1_sensitivity_spearman": pc1_spearman,
        "median_control_vs_drug_pcc": float(directions["control_vs_drug_pcc"].median()),
        "median_control_vs_drug_cosine": float(directions["control_vs_drug_cosine"].median()),
        "fraction_positive_cosine": float(directions["control_vs_drug_cosine"].gt(0).mean()),
        "shared_anchor_warning": (
            "control change and drug response share DMSO 24h with opposite signs; "
            "direction comparison is diagnostic, not an independent causal contrast"
        ),
    }
    log_path = root / "results" / "logs" / "landscape_summary.json"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_landscape_figures(root, landscape, directions, report)
    return report
