"""Gated W11 single-cell distribution translation in a fixed PCA space."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.stats import wasserstein_distance
from sklearn.decomposition import PCA

from yakseopdong.metrics import bootstrap_mean_interval
from yakseopdong.pseudobulk import read_vector_parquet

SEED = 20260827
MODELS = ("NO_CHANGE", "B1", "B4_FIXED_D20_A100", "CCLR_FIXED_D20_R20_A100")


def extract_h5ad_csc_rows(
    path: Path,
    selected_indices: NDArray[np.int64],
    chunk_size: int = 256,
) -> sparse.csr_matrix:
    """Extract selected H5AD CSC rows without materializing the full sparse matrix."""
    selected = np.asarray(selected_indices, dtype=np.int64)
    with h5py.File(path, "r") as handle:
        group = handle["X"]
        if str(group.attrs.get("encoding-type")) != "csc_matrix":
            raise ValueError("single-cell H5AD X must be CSC")
        n_rows, n_genes = (int(value) for value in group.attrs["shape"])
        row_to_local = np.full(n_rows, -1, dtype=np.int32)
        row_to_local[selected] = np.arange(len(selected), dtype=np.int32)
        indptr = np.asarray(group["indptr"], dtype=np.int64)
        row_parts: list[NDArray[np.int32]] = []
        column_parts: list[NDArray[np.int32]] = []
        data_parts: list[NDArray[np.float32]] = []
        for start in range(0, n_genes, chunk_size):
            stop = min(start + chunk_size, n_genes)
            pointer = indptr[start : stop + 1]
            first, last = int(pointer[0]), int(pointer[-1])
            source_rows = np.asarray(group["indices"][first:last], dtype=np.int32)
            local_rows = row_to_local[source_rows]
            keep = local_rows >= 0
            if not keep.any():
                continue
            lengths = np.diff(pointer)
            columns = np.repeat(
                np.arange(start, stop, dtype=np.int32), lengths
            )[keep]
            values = np.asarray(group["data"][first:last], dtype=np.float32)[keep]
            if (values < 0).any() or not np.allclose(values, np.round(values)):
                raise ValueError("single-cell X is not non-negative integer-like counts")
            row_parts.append(local_rows[keep])
            column_parts.append(columns)
            data_parts.append(values)
    matrix = sparse.coo_matrix(
        (
            np.concatenate(data_parts),
            (np.concatenate(row_parts), np.concatenate(column_parts)),
        ),
        shape=(len(selected), n_genes),
    )
    return matrix.tocsr()


def sparse_log1p_cpm(matrix: sparse.csr_matrix) -> sparse.csr_matrix:
    """Normalize single-cell counts while retaining sparse zero entries."""
    output = matrix.astype(np.float32, copy=True)
    library_sizes = np.asarray(output.sum(axis=1)).ravel().astype(float)
    if (library_sizes <= 0).any():
        raise ValueError("selected single cells include zero libraries")
    output = sparse.diags((1_000_000 / library_sizes).astype(np.float32)) @ output
    np.log1p(output.data, out=output.data)
    return output.tocsr()


def energy_distance_multivariate(
    left: NDArray[np.floating], right: NDArray[np.floating]
) -> float:
    """Biased non-negative multivariate energy distance."""
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    value = 2 * cdist(x, y).mean() - cdist(x, x).mean() - cdist(y, y).mean()
    return float(max(value, 0.0))


def sliced_wasserstein_distance(
    left: NDArray[np.floating],
    right: NDArray[np.floating],
    projections: NDArray[np.floating],
) -> float:
    """Mean 1D Wasserstein distance over fixed unit projections."""
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    directions = np.asarray(projections, dtype=float)
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    return float(
        np.mean(
            [
                wasserstein_distance(x @ direction, y @ direction)
                for direction in directions
            ]
        )
    )


def _load_cells(root: Path, external_depmap: set[str]) -> tuple[pd.DataFrame, sparse.csr_matrix]:
    path = root / "data" / "mcfarland_2020.h5ad"
    adata = ad.read_h5ad(path, backed="r")
    try:
        obs = adata.obs.reset_index(drop=False).rename(columns={"index": "cell_barcode"})
        condition = obs["hash_tag"].astype(str).map(
            {"DMSO_24hr": "control", "Tram_24hr": "trametinib"}
        )
        selected = (
            obs["cell_quality"].eq("normal")
            & obs["DepMap_ID"].astype(str).isin(external_depmap)
            & condition.notna()
        )
        metadata = obs.loc[
            selected,
            ["cell_barcode", "singlet_ID", "DepMap_ID", "disease", "ncounts"],
        ].copy()
        metadata["condition"] = condition.loc[selected].to_numpy()
        metadata = metadata.rename(
            columns={"singlet_ID": "cell_line", "DepMap_ID": "depmap_id"}
        ).reset_index(drop=True)
        indices = np.flatnonzero(selected.to_numpy()).astype(np.int64)
    finally:
        adata.file.close()
    matrix = extract_h5ad_csc_rows(path, indices)
    return metadata, sparse_log1p_cpm(matrix)


def _write_figure(
    root: Path,
    metrics: pd.DataFrame,
    scores: pd.DataFrame,
    shifts: pd.DataFrame,
) -> None:
    figure_dir = root / "results" / "figures"
    b1 = metrics.loc[metrics["model"].eq("B1")].set_index("cell_line")
    b4 = metrics.loc[metrics["model"].eq("B4_FIXED_D20_A100")].set_index("cell_line")
    order = (b1["energy_distance"] - b4["energy_distance"]).abs().sort_values()
    representative = str(order.index[len(order) // 2])
    line_scores = scores.loc[scores["cell_line"].eq(representative)]
    control = line_scores.loc[line_scores["condition"].eq("control"), ["pc1", "pc2"]].to_numpy()
    treated = line_scores.loc[
        line_scores["condition"].eq("trametinib"), ["pc1", "pc2"]
    ].to_numpy()
    shift_lookup = shifts.set_index(["cell_line", "model"])
    b1_shift = shift_lookup.loc[
        (representative, "B1"), ["pc1_shift", "pc2_shift"]
    ].to_numpy(dtype=float)
    b4_shift = shift_lookup.loc[
        (representative, "B4_FIXED_D20_A100"), ["pc1_shift", "pc2_shift"]
    ].to_numpy(dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    axes[0].scatter(
        control[:, 0],
        control[:, 1],
        s=18,
        alpha=0.45,
        label="control",
        color="#6B7280",
    )
    axes[0].scatter(
        treated[:, 0],
        treated[:, 1],
        s=18,
        alpha=0.55,
        label="observed treated",
        color="#136F63",
    )
    axes[0].scatter(
        control[:, 0] + b1_shift[0],
        control[:, 1] + b1_shift[1],
        s=18,
        alpha=0.45,
        label="B1 shift",
        color="#EDAE49",
    )
    axes[0].scatter(
        control[:, 0] + b4_shift[0],
        control[:, 1] + b4_shift[1],
        s=18,
        alpha=0.45,
        label="B4 shift",
        color="#D1495B",
    )
    axes[0].set(
        title=f"Representative external line: {representative.split('_')[0]}",
        xlabel="Single-cell baseline PCA PC1",
        ylabel="PC2",
    )
    axes[0].legend(frameon=False, fontsize=8)

    comparison = metrics.loc[metrics["model"].isin(["B1", "B4_FIXED_D20_A100"])]
    wide = comparison.pivot(index="cell_line", columns="model", values="energy_distance")
    for _, row in wide.iterrows():
        axes[1].plot(
            [0, 1],
            [row["B1"], row["B4_FIXED_D20_A100"]],
            color="#CBD5E1",
            linewidth=0.8,
        )
    axes[1].scatter(np.zeros(len(wide)), wide["B1"], color="#EDAE49", s=34, label="B1")
    axes[1].scatter(
        np.ones(len(wide)),
        wide["B4_FIXED_D20_A100"],
        color="#D1495B",
        s=34,
        label="B4",
    )
    axes[1].set_xticks([0, 1], ["B1", "B4"])
    axes[1].set(title="Paired external-line energy distance", ylabel="Lower is better")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(color="#E5E7EB", linewidth=0.7)
    fig.suptitle("Gated single-cell distribution extension", x=0.02, ha="left", fontweight="bold")
    fig.savefig(figure_dir / "single_cell_distribution.png", dpi=180)
    plt.close(fig)


def run_distribution(root: Path) -> dict[str, object]:
    """Execute the fixed W11 gate without tuning on distribution metrics."""
    prediction_meta, predictions = read_vector_parquet(
        root / "results" / "predictions" / "temporal_transfer_predictions.parquet",
        "predicted_delta_log1p_cpm",
    )
    prediction_meta = prediction_meta.loc[prediction_meta["time_hours"].eq(24)].reset_index(
        drop=True
    )
    predictions = predictions[
        read_vector_parquet(
            root / "results" / "predictions" / "temporal_transfer_predictions.parquet",
            "predicted_delta_log1p_cpm",
        )[0]["time_hours"].eq(24).to_numpy()
    ]
    external_depmap = set(prediction_meta["depmap_id"].astype(str))
    metadata, expression = _load_cells(root, external_depmap)
    control_mask = metadata["condition"].eq("control").to_numpy()
    control_expression = expression[control_mask]
    means = np.asarray(control_expression.mean(axis=0)).ravel()
    second = np.asarray(control_expression.power(2).mean(axis=0)).ravel()
    variances = np.maximum(second - np.square(means), 0.0)
    eligible = np.flatnonzero(variances > 1e-10)
    selected = eligible[np.argsort(-variances[eligible], kind="stable")[:2_000]]
    pca = PCA(n_components=10, svd_solver="randomized", random_state=SEED)
    pca.fit(control_expression[:, selected].toarray())
    all_scores = pca.transform(expression[:, selected].toarray())
    score_table = metadata.copy()
    for index in range(10):
        score_table[f"pc{index + 1}"] = all_scores[:, index]
    score_table.to_parquet(
        root / "results" / "tables" / "single_cell_pca_scores.parquet", index=False
    )

    prediction_lookup = {
        (str(row.cell_line), str(row.model)): predictions[index]
        for index, row in prediction_meta.iterrows()
    }
    shift_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    projections = np.random.default_rng(SEED).normal(size=(100, 10))
    for cell_line in sorted(score_table["cell_line"].astype(str).unique()):
        block = score_table.loc[score_table["cell_line"].eq(cell_line)]
        control_scores = block.loc[
            block["condition"].eq("control"), [f"pc{i}" for i in range(1, 11)]
        ].to_numpy()
        treated_scores = block.loc[
            block["condition"].eq("trametinib"), [f"pc{i}" for i in range(1, 11)]
        ].to_numpy()
        for model in MODELS:
            if model == "NO_CHANGE":
                shift = np.zeros(10, dtype=float)
            else:
                delta = prediction_lookup[(cell_line, model)]
                shift = np.asarray(delta[selected], dtype=float) @ pca.components_.T
            predicted_scores = control_scores + shift
            shift_rows.append(
                {
                    "cell_line": cell_line,
                    "model": model,
                    **{f"pc{i + 1}_shift": float(value) for i, value in enumerate(shift)},
                }
            )
            metric_rows.append(
                {
                    "cell_line": cell_line,
                    "model": model,
                    "n_control_cells": len(control_scores),
                    "n_treated_cells": len(treated_scores),
                    "energy_distance": energy_distance_multivariate(
                        predicted_scores, treated_scores
                    ),
                    "sliced_wasserstein": sliced_wasserstein_distance(
                        predicted_scores, treated_scores, projections
                    ),
                    "distribution_metric_used_for_tuning": False,
                    "paired_cell_trajectory_claim": False,
                }
            )
    shifts = pd.DataFrame(shift_rows)
    metrics = pd.DataFrame(metric_rows)
    shifts.to_csv(root / "results" / "tables" / "single_cell_shift_scores.csv", index=False)
    metrics.to_csv(
        root / "results" / "tables" / "single_cell_distribution_metrics.csv", index=False
    )

    summary_rows: list[dict[str, object]] = []
    for model in MODELS:
        block = metrics.loc[metrics["model"].eq(model)].set_index("cell_line")
        b1 = metrics.loc[metrics["model"].eq("B1")].set_index("cell_line")
        for metric in ("energy_distance", "sliced_wasserstein"):
            mean, low, high = bootstrap_mean_interval(block[metric].to_numpy(), seed=SEED)
            gain = b1[metric] - block[metric]
            gain_mean, gain_low, gain_high = bootstrap_mean_interval(
                gain.to_numpy(), seed=SEED + 1
            )
            summary_rows.append(
                {
                    "model": model,
                    "metric": metric,
                    "macro_mean": mean,
                    "ci95_low": low,
                    "ci95_high": high,
                    "gain_vs_b1": gain_mean,
                    "gain_ci95_low": gain_low,
                    "gain_ci95_high": gain_high,
                    "n_cell_lines": len(block),
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        root / "results" / "tables" / "single_cell_distribution_summary.csv", index=False
    )
    b4_gate = summary.loc[summary["model"].eq("B4_FIXED_D20_A100")]
    promoted = bool(b4_gate["gain_ci95_low"].gt(0).all())
    decision = "promote_to_main_text" if promoted else "not_promoted_future_work"
    gate = {
        "stage": "W11_single_cell_distribution",
        "decision": decision,
        "rule": "B4 paired gain vs B1 CI lower bound > 0 for both metrics",
        "passed": promoted,
        "b4_results": b4_gate.to_dict(orient="records"),
        "external_lines": int(metrics["cell_line"].nunique()),
        "cells": int(len(score_table)),
        "pca_variable_genes": len(selected),
        "pca_components": 10,
        "distribution_metric_used_for_tuning": False,
        "paired_cell_trajectory_claim": False,
    }
    (root / "results" / "logs" / "distribution_gate.json").write_text(
        json.dumps(gate, indent=2), encoding="utf-8"
    )
    decision_text = (
        "# W11 single-cell distribution gate\n\n"
        f"Decision: **{decision}**.\n\n"
        "The extension translates control-cell PCA scores by a fixed predicted pseudobulk "
        "delta. It does not infer paired trajectories and was not tuned on distribution metrics.\n"
    )
    (root / "report" / "w11_gate_decision.md").write_text(decision_text, encoding="utf-8")
    _write_figure(root, metrics, score_table, shifts)
    return gate
