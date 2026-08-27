"""Nested held-out-cell-line evaluation for the CCLR main model."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import yaml

from yakseopdong.benchmark import (
    METRICS,
    _comparison_table,
    _line_metrics,
    _load_matrices,
    _sha256,
    _summarize_metrics,
    run_baselines,
)
from yakseopdong.metrics import pearson_or_nan, rmse
from yakseopdong.models import CCLRArtifact, fit_predict_cclr
from yakseopdong.plots import write_cclr_figures
from yakseopdong.pseudobulk import write_vector_parquet
from yakseopdong.splits import validate_outer_assignments

CCLR_METRICS = [*METRICS, "rmse_gain_vs_b4", "nrmse_gain_vs_b4"]


def _write_fold_artifact(path: Path, artifact: CCLRArtifact) -> None:
    """Serialize the numerical state needed to reconstruct one fold model."""
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        control_gene_indices=artifact.control_embedding.gene_indices,
        control_pca_mean=artifact.control_embedding.pca.mean_.astype(np.float32),
        control_pca_components=artifact.control_embedding.pca.components_.astype(np.float32),
        control_pca_explained_variance=(
            artifact.control_embedding.pca.explained_variance_.astype(np.float32)
        ),
        control_pca_explained_variance_ratio=(
            artifact.control_embedding.pca.explained_variance_ratio_.astype(np.float32)
        ),
        response_pca_mean=artifact.response_embedding.pca.mean_.astype(np.float32),
        response_pca_components=artifact.response_embedding.pca.components_.astype(np.float32),
        response_pca_explained_variance_ratio=(
            artifact.response_embedding.pca.explained_variance_ratio_.astype(np.float32)
        ),
        ridge_coef=np.asarray(artifact.ridge.coef_, dtype=np.float32),
        ridge_intercept=np.asarray(artifact.ridge.intercept_, dtype=np.float32),
        control_dimension=np.asarray(artifact.control_dimension, dtype=np.int64),
        response_rank=np.asarray(artifact.response_rank, dtype=np.int64),
        alpha=np.asarray(artifact.alpha, dtype=np.float64),
    )


def _component_loading_rows(
    artifact: CCLRArtifact,
    genes: pd.DataFrame,
    outer_fold: int,
    top_k: int = 50,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for component_index, loadings in enumerate(
        artifact.response_embedding.pca.components_, start=1
    ):
        orders = {
            "negative": np.argsort(loadings, kind="stable")[:top_k],
            "positive": np.argsort(-loadings, kind="stable")[:top_k],
        }
        for direction, indices in orders.items():
            for rank, gene_index in enumerate(indices, start=1):
                gene = genes.iloc[int(gene_index)]
                rows.append(
                    {
                        "outer_fold": outer_fold,
                        "component": component_index,
                        "direction": direction,
                        "loading_rank": rank,
                        "gene_index": int(gene_index),
                        "gene_id": gene["gene_id"],
                        "gene_symbol": gene["gene_symbol"],
                        "loading": float(loadings[gene_index]),
                    }
                )
    return rows


def _component_summary(scores: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (outer_fold, component), group in scores.groupby(["outer_fold", "component"]):
        rows.append(
            {
                "outer_fold": int(outer_fold),
                "component": int(component),
                "n_test_lines": int(len(group)),
                "observed_vs_predicted_pcc": pearson_or_nan(
                    group["observed_score"], group["predicted_score"]
                ),
                "score_rmse": rmse(group["observed_score"], group["predicted_score"]),
                "explained_variance_ratio": float(
                    group["explained_variance_ratio"].iloc[0]
                ),
            }
        )
    return pd.DataFrame(rows)


def run_cclr(root: Path, seed: int = 20260827) -> dict[str, object]:
    """Run leakage-safe nested CV for CCLR and compare it with B0-B4."""
    run_baselines(root, seed=seed)
    annotations = pd.read_csv(root / "cell_line_annotations.csv").sort_values(
        "cell_line", ignore_index=True
    )
    outer = pd.read_csv(root / "split_assignments.csv")
    inner = pd.read_csv(root / "inner_split_assignments.csv")
    validate_outer_assignments(outer, n_splits=5)
    annotations = annotations.merge(
        outer[["cell_line", "outer_fold"]], on="cell_line", validate="one_to_one"
    ).sort_values("cell_line", ignore_index=True)
    cell_lines = annotations["cell_line"].astype(str).tolist()
    control, response = _load_matrices(root, cell_lines)
    genes = pd.read_parquet(root / "data" / "processed" / "gene_metadata.parquet")
    if len(genes) != response.shape[1]:
        raise ValueError("gene metadata and response matrix do not align")

    with (root / "config" / "models.yaml").open(encoding="utf-8") as handle:
        model_config = yaml.safe_load(handle)
    cclr_config = model_config["main_model"]
    feature_config = model_config["baseline_protocol"]["feature_selection"]
    bootstrap_replicates = int(cclr_config["bootstrap_replicates"])

    baseline_metrics = pd.read_csv(root / "results" / "tables" / "baseline_metrics_by_line.csv")
    b4_metrics = baseline_metrics.loc[baseline_metrics["model"].eq("B4")].set_index(
        "cell_line"
    )
    baseline_comparison = pd.read_csv(
        root / "results" / "tables" / "baseline_comparison.csv"
    )

    prediction_rows: list[dict[str, object]] = []
    prediction_vectors: list[np.ndarray] = []
    metric_rows: list[dict[str, object]] = []
    runtime_rows: list[dict[str, object]] = []
    hyperparameter_rows: list[dict[str, object]] = []
    inner_trace_rows: list[dict[str, object]] = []
    loading_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    artifact_rows: list[dict[str, object]] = []

    started = perf_counter()
    for outer_fold in range(5):
        test_mask = annotations["outer_fold"].eq(outer_fold).to_numpy()
        train_mask = ~test_mask
        train_lines = annotations.loc[train_mask, "cell_line"].astype(str).tolist()
        inner_map = (
            inner.loc[inner["outer_fold"].eq(outer_fold)]
            .set_index("cell_line")["inner_fold"]
        )
        if set(inner_map.index) != set(train_lines):
            raise ValueError(f"inner assignment coverage differs in outer fold {outer_fold}")
        inner_folds = inner_map.loc[train_lines].to_numpy(dtype=int)

        fold_started = perf_counter()
        prediction, info, artifact = fit_predict_cclr(
            train_control=control[train_mask],
            train_response=response[train_mask],
            test_control=control[test_mask],
            inner_fold_ids=inner_folds,
            config=cclr_config,
            feature_config=feature_config,
            seed=seed + 1_000 * outer_fold,
        )
        fold_seconds = perf_counter() - fold_started
        artifact_path = root / "results" / "models" / f"cclr_outer_fold_{outer_fold}.npz"
        _write_fold_artifact(artifact_path, artifact)
        artifact_rows.append(
            {
                "outer_fold": outer_fold,
                "path": str(artifact_path.relative_to(root)),
                "bytes": artifact_path.stat().st_size,
                "sha256": _sha256(artifact_path),
            }
        )

        test_indices = np.flatnonzero(test_mask)
        mean_response = artifact.response_embedding.pca.mean_.astype(np.float32)
        for test_position, global_index in enumerate(test_indices):
            cell_line = str(annotations.loc[global_index, "cell_line"])
            observed = response[global_index]
            predicted = prediction[test_position]
            row = {
                "cell_line": cell_line,
                "depmap_id": annotations.loc[global_index, "depmap_id"],
                "lineage": annotations.loc[global_index, "lineage"],
                "outer_fold": outer_fold,
                "model": "CCLR",
            }
            prediction_rows.append(row)
            prediction_vectors.append(predicted)
            line_metrics = _line_metrics(observed, predicted, mean_response)
            b4_line = b4_metrics.loc[cell_line]
            metric_rows.append(
                {
                    **row,
                    **line_metrics,
                    "rmse_gain_vs_b4": float(
                        b4_line["rmse_delta"] - line_metrics["rmse_delta"]
                    ),
                    "nrmse_gain_vs_b4": float(
                        b4_line["nrmse_delta"] - line_metrics["nrmse_delta"]
                    ),
                }
            )

        test_control_scores = artifact.control_embedding.transform(control[test_mask])[
            :, : artifact.control_dimension
        ]
        predicted_scores = np.asarray(
            artifact.ridge.predict(test_control_scores), dtype=np.float64
        )
        if predicted_scores.ndim == 1:
            predicted_scores = predicted_scores[:, None]
        observed_scores = artifact.response_embedding.transform(response[test_mask])
        variance = artifact.response_embedding.pca.explained_variance_ratio_
        for test_position, global_index in enumerate(test_indices):
            for component_index in range(artifact.response_rank):
                score_rows.append(
                    {
                        "cell_line": annotations.loc[global_index, "cell_line"],
                        "outer_fold": outer_fold,
                        "component": component_index + 1,
                        "observed_score": float(
                            observed_scores[test_position, component_index]
                        ),
                        "predicted_score": float(
                            predicted_scores[test_position, component_index]
                        ),
                        "explained_variance_ratio": float(variance[component_index]),
                        "observed_test_response_used_for_fit": False,
                    }
                )
        loading_rows.extend(_component_loading_rows(artifact, genes, outer_fold))
        hyperparameter_rows.append(
            {
                "outer_fold": outer_fold,
                "control_dimension": int(info["control_dimension"]),
                "response_rank": int(info["response_rank"]),
                "alpha": float(info["alpha"]),
                "selected_control_genes": int(info["selected_control_genes"]),
                "response_variance_explained": float(
                    np.sum(info["response_explained_variance_ratio"])
                ),
                "parameter_count": int(info["parameter_count"]),
            }
        )
        runtime_rows.append(
            {
                "outer_fold": outer_fold,
                "model": "CCLR",
                "fold_pipeline_seconds": fold_seconds,
                "parameter_count": int(info["parameter_count"]),
                "selected_control_genes": int(info["selected_control_genes"]),
            }
        )
        for row in info["inner_cv_trace"]:
            inner_trace_rows.append({"outer_fold": outer_fold, **row})

    table_dir = root / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    prediction_metadata = pd.DataFrame(prediction_rows)
    prediction_matrix = np.vstack(prediction_vectors).astype(np.float32)
    prediction_path = root / "results" / "predictions" / "cclr_predictions.parquet"
    prior_prediction_sha256 = _sha256(prediction_path) if prediction_path.exists() else None
    write_vector_parquet(
        prediction_path,
        prediction_metadata,
        prediction_matrix,
        "predicted_delta_log1p_cpm",
    )
    prediction_sha256 = _sha256(prediction_path)

    metrics = pd.DataFrame(metric_rows)
    summary = _summarize_metrics(
        metrics,
        seed=seed + 50_000,
        models=["CCLR"],
        metric_names=CCLR_METRICS,
        bootstrap_replicates=bootstrap_replicates,
    )
    cclr_comparison = _comparison_table(
        summary, models=["CCLR"], metric_names=CCLR_METRICS
    )
    comparison = pd.concat([baseline_comparison, cclr_comparison], ignore_index=True)
    hyperparameters = pd.DataFrame(hyperparameter_rows)
    scores = pd.DataFrame(score_rows)
    component_summary = _component_summary(scores)
    artifact_manifest = pd.DataFrame(artifact_rows)

    metrics.to_csv(table_dir / "cclr_metrics_by_line.csv", index=False)
    summary.to_csv(table_dir / "cclr_metrics_summary.csv", index=False)
    comparison.to_csv(table_dir / "model_comparison_w6.csv", index=False)
    hyperparameters.to_csv(table_dir / "cclr_hyperparameters.csv", index=False)
    pd.DataFrame(inner_trace_rows).to_csv(table_dir / "cclr_inner_cv.csv", index=False)
    pd.DataFrame(runtime_rows).to_csv(table_dir / "cclr_runtime.csv", index=False)
    pd.DataFrame(loading_rows).to_csv(
        table_dir / "cclr_component_top_loadings.csv", index=False
    )
    scores.to_csv(table_dir / "cclr_component_scores.csv", index=False)
    component_summary.to_csv(table_dir / "cclr_component_summary.csv", index=False)
    artifact_manifest.to_csv(table_dir / "cclr_model_artifacts.csv", index=False)
    write_cclr_figures(root, comparison, metrics, hyperparameters, component_summary)

    cclr_row = cclr_comparison.set_index("model").loc["CCLR"]
    report = {
        "status": "complete",
        "device": "cpu",
        "cell_lines": len(annotations),
        "genes": int(response.shape[1]),
        "outer_folds": 5,
        "inner_folds": 4,
        "model": "CCLR",
        "total_seconds": perf_counter() - started,
        "predictions_rows": int(len(prediction_metadata)),
        "each_cell_line_predicted_once": bool(
            prediction_metadata.groupby("cell_line").size().eq(1).all()
        ),
        "selected_hyperparameters": hyperparameters.to_dict(orient="records"),
        "primary_metrics": {
            "rmse_delta_mean": float(cclr_row["rmse_delta_mean"]),
            "pcc_delta_mean": float(cclr_row["pcc_delta_mean"]),
            "pcc_context_mean": float(cclr_row["pcc_context_mean"]),
            "rmse_gain_vs_b1_mean": float(cclr_row["rmse_gain_vs_b1_mean"]),
            "rmse_gain_vs_b1_ci95_low": float(
                cclr_row["rmse_gain_vs_b1_ci95_low"]
            ),
            "rmse_gain_vs_b1_ci95_high": float(
                cclr_row["rmse_gain_vs_b1_ci95_high"]
            ),
            "rmse_gain_vs_b4_mean": float(cclr_row["rmse_gain_vs_b4_mean"]),
            "rmse_gain_vs_b4_ci95_low": float(
                cclr_row["rmse_gain_vs_b4_ci95_low"]
            ),
            "rmse_gain_vs_b4_ci95_high": float(
                cclr_row["rmse_gain_vs_b4_ci95_high"]
            ),
        },
        "prediction_sha256": prediction_sha256,
        "repeat_run_prediction_sha256_match": bool(
            prior_prediction_sha256 is not None
            and prior_prediction_sha256 == prediction_sha256
        ),
        "models_config_sha256": _sha256(root / "config" / "models.yaml"),
        "splits_config_sha256": _sha256(root / "config" / "splits.yaml"),
        "artifact_manifest_sha256": hashlib.sha256(
            artifact_manifest.to_csv(index=False).encode("utf-8")
        ).hexdigest(),
        "outer_test_response_used_for_fit": False,
        "outer_test_response_used_for_component_evaluation_only": True,
        "sensitivity_used_as_predictor": False,
        "mutations_used_as_predictor": False,
    }
    log_path = root / "results" / "logs" / "cclr_summary.json"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
