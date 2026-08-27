"""Independent validation of W6 metrics, predictions, and fold artifacts."""

from __future__ import annotations

import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd

from yakseopdong.benchmark import _sha256
from yakseopdong.models import fit_predict_cclr
from yakseopdong.pseudobulk import read_vector_parquet


def _row_rmse(observed: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(observed.astype(float) - predicted.astype(float)))))


def validate_cclr(root: Path, seed: int = 20260827) -> dict[str, object]:
    """Recalculate W6 evidence without using the benchmark metric helpers."""
    response_meta, response = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet", "delta_log1p_cpm"
    )
    control_meta, control = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_24h.parquet", "log1p_cpm"
    )
    baseline_meta, baseline = read_vector_parquet(
        root / "results" / "predictions" / "baseline_predictions.parquet",
        "predicted_delta_log1p_cpm",
    )
    cclr_meta, cclr = read_vector_parquet(
        root / "results" / "predictions" / "cclr_predictions.parquet",
        "predicted_delta_log1p_cpm",
    )
    response_by_line = {
        str(row.cell_line): response[position]
        for position, row in response_meta.iterrows()
    }
    control_mask = control_meta["condition"].eq("control") & control_meta[
        "time_hours"
    ].eq(24)
    control_by_line = {
        str(control_meta.iloc[position].cell_line): control[position]
        for position in np.flatnonzero(control_mask.to_numpy())
    }
    baseline_by_line_model = {
        (str(row.cell_line), str(row.model)): baseline[position]
        for position, row in baseline_meta.iterrows()
    }
    stored_metrics = pd.read_csv(
        root / "results" / "tables" / "cclr_metrics_by_line.csv"
    ).set_index("cell_line")

    recalculated: list[dict[str, float | str]] = []
    maximum_metric_difference = 0.0
    for position, row in cclr_meta.iterrows():
        cell_line = str(row.cell_line)
        observed = response_by_line[cell_line]
        cclr_prediction = cclr[position]
        cclr_rmse = _row_rmse(observed, cclr_prediction)
        b1_rmse = _row_rmse(observed, baseline_by_line_model[(cell_line, "B1")])
        b4_rmse = _row_rmse(observed, baseline_by_line_model[(cell_line, "B4")])
        values = {
            "cell_line": cell_line,
            "rmse_delta": cclr_rmse,
            "rmse_gain_vs_b1": b1_rmse - cclr_rmse,
            "rmse_gain_vs_b4": b4_rmse - cclr_rmse,
        }
        recalculated.append(values)
        for metric in ["rmse_delta", "rmse_gain_vs_b1", "rmse_gain_vs_b4"]:
            maximum_metric_difference = max(
                maximum_metric_difference,
                abs(float(values[metric]) - float(stored_metrics.loc[cell_line, metric])),
            )
    recalculated_frame = pd.DataFrame(recalculated)
    gain_vs_b4 = recalculated_frame["rmse_gain_vs_b4"].to_numpy(dtype=float)
    rng = np.random.default_rng(seed + 50_000 + 8)
    bootstrap_draws = rng.choice(
        gain_vs_b4, size=(2_000, len(gain_vs_b4)), replace=True
    ).mean(axis=1)
    independent_ci = np.quantile(bootstrap_draws, [0.025, 0.975])
    summary = json.loads((root / "results" / "logs" / "cclr_summary.json").read_text())
    reported_ci = np.asarray(
        [
            summary["primary_metrics"]["rmse_gain_vs_b4_ci95_low"],
            summary["primary_metrics"]["rmse_gain_vs_b4_ci95_high"],
        ],
        dtype=float,
    )

    maximum_artifact_prediction_difference = 0.0
    artifact_manifest = pd.read_csv(
        root / "results" / "tables" / "cclr_model_artifacts.csv"
    )
    artifact_hashes_match = True
    for fold in range(5):
        manifest_row = artifact_manifest.loc[artifact_manifest["outer_fold"].eq(fold)].iloc[0]
        artifact_path = root / str(manifest_row["path"])
        artifact_hashes_match &= _sha256(artifact_path) == str(manifest_row["sha256"])
        artifact = np.load(artifact_path)
        prediction_positions = cclr_meta.index[cclr_meta["outer_fold"].eq(fold)].to_numpy()
        fold_control = np.vstack(
            [
                control_by_line[str(cclr_meta.loc[position, "cell_line"])]
                for position in prediction_positions
            ]
        ).astype(float)
        selected = artifact["control_gene_indices"]
        control_scores = (
            (fold_control[:, selected] - artifact["control_pca_mean"])
            @ artifact["control_pca_components"].T
        ) / np.sqrt(artifact["control_pca_explained_variance"])
        ridge_coef = artifact["ridge_coef"]
        if ridge_coef.ndim == 1:
            ridge_coef = ridge_coef[None, :]
        predicted_scores = (
            control_scores[:, : int(artifact["control_dimension"])] @ ridge_coef.T
            + np.atleast_1d(artifact["ridge_intercept"])
        )
        reconstructed = (
            predicted_scores @ artifact["response_pca_components"]
            + artifact["response_pca_mean"]
        )
        maximum_artifact_prediction_difference = max(
            maximum_artifact_prediction_difference,
            float(np.max(np.abs(reconstructed - cclr[prediction_positions]))),
        )

    api_parameters = set(inspect.signature(fit_predict_cclr).parameters)
    checks = {
        "prediction_rows_are_94": len(cclr_meta) == 94,
        "each_cell_line_predicted_once": cclr_meta["cell_line"].is_unique,
        "fit_api_has_no_test_response": "test_response" not in api_parameters,
        "stored_metrics_match_independent_calculation": maximum_metric_difference < 1e-12,
        "paired_bootstrap_ci_matches": bool(
            np.allclose(independent_ci, reported_ci, atol=1e-15)
        ),
        "five_artifacts_present": len(artifact_manifest) == 5,
        "artifact_hashes_match_manifest": bool(artifact_hashes_match),
        "artifact_predictions_reconstruct": maximum_artifact_prediction_difference < 2e-5,
        "repeat_prediction_hash_matches": bool(
            summary["repeat_run_prediction_sha256_match"]
        ),
        "prediction_hash_matches_summary": (
            _sha256(root / "results" / "predictions" / "cclr_predictions.parquet")
            == summary["prediction_sha256"]
        ),
    }
    report = {
        "status": "pass" if all(checks.values()) else "fail",
        "checks": checks,
        "independent_macro_rmse": float(recalculated_frame["rmse_delta"].mean()),
        "independent_macro_rmse_gain_vs_b1": float(
            recalculated_frame["rmse_gain_vs_b1"].mean()
        ),
        "independent_macro_rmse_gain_vs_b4": float(gain_vs_b4.mean()),
        "independent_rmse_gain_vs_b4_ci95": independent_ci.tolist(),
        "maximum_metric_absolute_difference": maximum_metric_difference,
        "maximum_artifact_prediction_absolute_difference": (
            maximum_artifact_prediction_difference
        ),
    }
    log_path = root / "results" / "logs" / "cclr_validation.json"
    log_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if report["status"] != "pass":
        raise RuntimeError(f"CCLR validation failed: {checks}")
    return report
