"""Nested held-out-cell-line evaluation for the B0-B4 baselines."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
import yaml

from yakseopdong.landscape import _aligned_rows
from yakseopdong.metrics import (
    bootstrap_mean_interval,
    context_residual,
    nrmse,
    pearson_or_nan,
    rmse,
    signed_topk_overlap,
    spearman_or_nan,
)
from yakseopdong.models import fit_predict_baselines
from yakseopdong.plots import write_baseline_figure
from yakseopdong.pseudobulk import read_vector_parquet, write_vector_parquet
from yakseopdong.splits import run_splits, validate_outer_assignments

MODELS = ["B0", "B1", "B2", "B3", "B4"]
METRICS = [
    "rmse_delta",
    "nrmse_delta",
    "pcc_delta",
    "spearman_delta",
    "signed_top50_overlap",
    "pcc_context",
    "rmse_gain_vs_b1",
    "nrmse_gain_vs_b1",
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_matrices(root: Path, cell_lines: list[str]) -> tuple[np.ndarray, np.ndarray]:
    primary_meta, primary = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_24h.parquet", "log1p_cpm"
    )
    response_meta, response_values = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet", "delta_log1p_cpm"
    )
    control = _aligned_rows(
        primary_meta, primary, cell_lines, condition="control", time_hours=24
    )
    response = _aligned_rows(response_meta, response_values, cell_lines)
    return control, response


def _line_metrics(
    observed: np.ndarray,
    predicted: np.ndarray,
    b1_prediction: np.ndarray,
) -> dict[str, float]:
    b1_rmse = rmse(observed, b1_prediction)
    b1_nrmse = nrmse(observed, b1_prediction)
    return {
        "rmse_delta": rmse(observed, predicted),
        "nrmse_delta": nrmse(observed, predicted),
        "pcc_delta": pearson_or_nan(observed, predicted),
        "spearman_delta": spearman_or_nan(observed, predicted),
        "signed_top50_overlap": signed_topk_overlap(observed, predicted, k=50),
        "pcc_context": pearson_or_nan(
            context_residual(observed, b1_prediction),
            context_residual(predicted, b1_prediction),
        ),
        "rmse_gain_vs_b1": b1_rmse - rmse(observed, predicted),
        "nrmse_gain_vs_b1": b1_nrmse - nrmse(observed, predicted),
    }


def _summarize_metrics(metrics: pd.DataFrame, seed: int) -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for model_index, model in enumerate(MODELS):
        model_rows = metrics.loc[metrics["model"].eq(model)]
        for metric_index, metric in enumerate(METRICS):
            values = model_rows[metric].to_numpy(dtype=float)
            mean, low, high = bootstrap_mean_interval(
                values,
                n_bootstrap=2_000,
                seed=seed + 100 * model_index + metric_index,
            )
            rows.append(
                {
                    "model": model,
                    "metric": metric,
                    "macro_mean": mean,
                    "ci95_low": low,
                    "ci95_high": high,
                    "n_cell_lines": int(np.isfinite(values).sum()),
                    "bootstrap_replicates": 2_000,
                }
            )
    return pd.DataFrame(rows)


def _comparison_table(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    for model in MODELS:
        selected = summary.loc[summary["model"].eq(model)].set_index("metric")
        row: dict[str, float | str] = {"model": model}
        for metric in METRICS:
            row[f"{metric}_mean"] = float(selected.loc[metric, "macro_mean"])
            row[f"{metric}_ci95_low"] = float(selected.loc[metric, "ci95_low"])
            row[f"{metric}_ci95_high"] = float(selected.loc[metric, "ci95_high"])
        rows.append(row)
    return pd.DataFrame(rows)


def run_baselines(root: Path, seed: int = 20260827) -> dict[str, object]:
    """Run B0-B4 in nested CV and write predictions, metrics, and audit logs."""
    run_splits(root, seed=seed)
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

    with (root / "config" / "models.yaml").open(encoding="utf-8") as handle:
        model_config = yaml.safe_load(handle)
    baseline_config = model_config["baseline_protocol"]

    prediction_rows: list[dict[str, object]] = []
    prediction_vectors: list[np.ndarray] = []
    metric_rows: list[dict[str, float | int | str | bool]] = []
    runtime_rows: list[dict[str, float | int | str]] = []
    hyperparameter_rows: list[dict[str, float | int | str]] = []
    inner_trace_rows: list[dict[str, float | int | str]] = []

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
        predictions, info = fit_predict_baselines(
            train_control=control[train_mask],
            train_response=response[train_mask],
            train_lineages=annotations.loc[train_mask, "lineage"].to_numpy(dtype=str),
            test_control=control[test_mask],
            test_lineages=annotations.loc[test_mask, "lineage"].to_numpy(dtype=str),
            inner_fold_ids=inner_folds,
            config=baseline_config,
            seed=seed + 1_000 * outer_fold,
        )
        fold_seconds = perf_counter() - fold_started

        train_indices = np.flatnonzero(train_mask)
        test_indices = np.flatnonzero(test_mask)
        mean_response = predictions["B1"][0]
        for test_position, global_index in enumerate(test_indices):
            observed = response[global_index]
            for model in MODELS:
                predicted = predictions[model][test_position]
                row = {
                    "cell_line": annotations.loc[global_index, "cell_line"],
                    "depmap_id": annotations.loc[global_index, "depmap_id"],
                    "lineage": annotations.loc[global_index, "lineage"],
                    "outer_fold": outer_fold,
                    "model": model,
                    "b2_fallback": (
                        bool(info["b2_fallback"][test_position]) if model == "B2" else False
                    ),
                    "b3_neighbor_cell_line": (
                        annotations.loc[
                            train_indices[info["b3_neighbor_train_position"][test_position]],
                            "cell_line",
                        ]
                        if model == "B3"
                        else ""
                    ),
                }
                prediction_rows.append(row)
                prediction_vectors.append(predicted)
                metric_rows.append(
                    {
                        **row,
                        **_line_metrics(observed, predicted, mean_response),
                    }
                )

        for model in MODELS:
            runtime_rows.append(
                {
                    "outer_fold": outer_fold,
                    "model": model,
                    "fold_pipeline_seconds": fold_seconds,
                    "parameter_count": int(info["parameter_count"][model]),
                    "selected_control_genes": int(info["selected_control_genes"]),
                }
            )
        hyperparameter_rows.append(
            {
                "outer_fold": outer_fold,
                "b3_dimension": int(info["b3_dimension"]),
                "b4_dimension": int(info["b4_dimension"]),
                "b4_alpha": float(info["b4_alpha"]),
                "b2_fallback_test_lines": int(sum(info["b2_fallback"])),
            }
        )
        for row in info["inner_cv_trace"]:
            inner_trace_rows.append({"outer_fold": outer_fold, **row})

    prediction_metadata = pd.DataFrame(prediction_rows)
    prediction_matrix = np.vstack(prediction_vectors).astype(np.float32)
    prediction_path = root / "results" / "predictions" / "baseline_predictions.parquet"
    write_vector_parquet(
        prediction_path,
        prediction_metadata,
        prediction_matrix,
        "predicted_delta_log1p_cpm",
    )
    metrics = pd.DataFrame(metric_rows)
    summary = _summarize_metrics(metrics, seed)
    comparison = _comparison_table(summary)
    table_dir = root / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(table_dir / "baseline_metrics_by_line.csv", index=False)
    summary.to_csv(table_dir / "baseline_metrics_summary.csv", index=False)
    comparison.to_csv(table_dir / "baseline_comparison.csv", index=False)
    pd.DataFrame(runtime_rows).to_csv(table_dir / "baseline_runtime.csv", index=False)
    pd.DataFrame(hyperparameter_rows).to_csv(
        table_dir / "baseline_hyperparameters.csv", index=False
    )
    pd.DataFrame(inner_trace_rows).to_csv(table_dir / "baseline_inner_cv.csv", index=False)
    write_baseline_figure(root, comparison)

    total_seconds = perf_counter() - started
    raw_primary_metrics = comparison.set_index("model")[[
        "rmse_delta_mean",
        "pcc_delta_mean",
        "pcc_context_mean",
        "rmse_gain_vs_b1_mean",
        "rmse_gain_vs_b1_ci95_low",
        "rmse_gain_vs_b1_ci95_high",
    ]].to_dict(orient="index")
    primary_metrics = {
        model: {
            metric: (float(value) if np.isfinite(value) else None)
            for metric, value in values.items()
        }
        for model, values in raw_primary_metrics.items()
    }
    report = {
        "status": "complete",
        "device": "cpu",
        "cell_lines": len(annotations),
        "genes": int(response.shape[1]),
        "outer_folds": 5,
        "inner_folds": 4,
        "models": MODELS,
        "total_seconds": total_seconds,
        "predictions_rows": int(len(prediction_metadata)),
        "each_cell_line_predicted_once_per_model": bool(
            prediction_metadata.groupby(["model", "cell_line"]).size().eq(1).all()
        ),
        "primary_metrics": primary_metrics,
        "prediction_sha256": _sha256(prediction_path),
        "models_config_sha256": _sha256(root / "config" / "models.yaml"),
        "splits_config_sha256": _sha256(root / "config" / "splits.yaml"),
        "sensitivity_used_as_predictor": False,
        "mutations_used_as_predictor": False,
    }
    log_path = root / "results" / "logs" / "baseline_summary.json"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
