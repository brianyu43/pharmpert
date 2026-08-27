"""Independent validation of the gated W11 single-cell distribution extension."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from yakseopdong.distribution import (
    energy_distance_multivariate,
    sliced_wasserstein_distance,
)
from yakseopdong.metrics import bootstrap_mean_interval


def validate_distribution(root: Path) -> dict[str, object]:
    """Recompute W11 distances, summary gains, and the promotion gate."""
    metrics = pd.read_csv(
        root / "results" / "tables" / "single_cell_distribution_metrics.csv"
    )
    summary = pd.read_csv(
        root / "results" / "tables" / "single_cell_distribution_summary.csv"
    )
    scores = pd.read_parquet(root / "results" / "tables" / "single_cell_pca_scores.parquet")
    shifts = pd.read_csv(root / "results" / "tables" / "single_cell_shift_scores.csv")
    gate = json.loads((root / "results" / "logs" / "distribution_gate.json").read_text())
    if len(metrics) != 68 or metrics["cell_line"].nunique() != 17:
        raise ValueError("W11 metric grid must contain 17 lines x 4 models")
    if len(shifts) != 68 or shifts[["cell_line", "model"]].duplicated().any():
        raise ValueError("W11 latent-shift grid is incomplete")
    if metrics["distribution_metric_used_for_tuning"].any():
        raise ValueError("W11 distribution metrics were marked as tuning inputs")
    if metrics["paired_cell_trajectory_claim"].any():
        raise ValueError("W11 incorrectly claims paired cell trajectories")

    score_columns = [f"pc{index}" for index in range(1, 11)]
    shift_columns = [f"pc{index}_shift" for index in range(1, 11)]
    lookup = metrics.set_index(["cell_line", "model"])
    projections = np.random.default_rng(20260827).normal(size=(100, 10))
    errors: list[float] = []
    for _, row in shifts.iterrows():
        block = scores.loc[scores["cell_line"].eq(row.cell_line)]
        control = block.loc[block["condition"].eq("control"), score_columns].to_numpy()
        treated = block.loc[
            block["condition"].eq("trametinib"), score_columns
        ].to_numpy()
        predicted = control + row[shift_columns].to_numpy(dtype=float)
        reported = lookup.loc[(row.cell_line, row.model)]
        errors.extend(
            [
                abs(
                    energy_distance_multivariate(predicted, treated)
                    - float(reported.energy_distance)
                ),
                abs(
                    sliced_wasserstein_distance(predicted, treated, projections)
                    - float(reported.sliced_wasserstein)
                ),
            ]
        )
    max_metric_error = float(max(errors))
    if max_metric_error > 1e-9:
        raise ValueError("W11 distribution metrics do not recompute from latent scores")

    b1 = metrics.loc[metrics["model"].eq("B1")].set_index("cell_line")
    b4 = metrics.loc[metrics["model"].eq("B4_FIXED_D20_A100")].set_index("cell_line")
    summary_errors: list[float] = []
    for metric in ("energy_distance", "sliced_wasserstein"):
        gain = b1[metric] - b4[metric]
        mean, low, high = bootstrap_mean_interval(gain.to_numpy(), seed=20260828)
        row = summary.loc[
            summary["model"].eq("B4_FIXED_D20_A100")
            & summary["metric"].eq(metric)
        ].iloc[0]
        summary_errors.extend(
            [
                abs(mean - row.gain_vs_b1),
                abs(low - row.gain_ci95_low),
                abs(high - row.gain_ci95_high),
            ]
        )
    max_summary_error = float(max(summary_errors))
    if max_summary_error > 1e-12:
        raise ValueError("W11 B4 paired-gain summary does not recompute")
    expected_pass = bool(
        summary.loc[
            summary["model"].eq("B4_FIXED_D20_A100"), "gain_ci95_low"
        ].gt(0).all()
    )
    if bool(gate["passed"]) != expected_pass:
        raise ValueError("W11 promotion decision differs from the frozen gate")
    figure = root / "results" / "figures" / "single_cell_distribution.png"
    if not figure.exists() or figure.stat().st_size < 20_000:
        raise ValueError("the W11 distribution figure is absent or unexpectedly small")
    report = {
        "stage": "W11_single_cell_distribution",
        "status": "passed",
        "decision": gate["decision"],
        "gate_passed": expected_pass,
        "cell_lines": metrics["cell_line"].nunique(),
        "cells": len(scores),
        "metric_rows": len(metrics),
        "max_metric_recomputation_error": max_metric_error,
        "max_summary_recomputation_error": max_summary_error,
        "distribution_metric_used_for_tuning": False,
        "paired_cell_trajectory_claim": False,
    }
    output = root / "results" / "logs" / "distribution_validation.json"
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
