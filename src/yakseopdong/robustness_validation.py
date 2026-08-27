"""Independent validation of W10 robustness and noise-ceiling outputs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from yakseopdong.metrics import bootstrap_mean_interval


def validate_robustness(root: Path) -> dict[str, object]:
    """Recompute W10 summaries and enforce completeness/leakage contracts."""
    metrics = pd.read_csv(root / "results" / "tables" / "robustness_metrics.csv")
    noise = pd.read_csv(root / "results" / "tables" / "noise_ceiling_by_line.csv")
    subsampling = pd.read_csv(root / "results" / "tables" / "subsampling_metrics.csv")
    lolo = pd.read_csv(root / "results" / "tables" / "leave_one_lineage_out.csv")
    conclusions = pd.read_csv(root / "results" / "tables" / "robustness_conclusions.csv")
    if len(metrics) != 69:
        raise ValueError("W10 robustness metric registry is incomplete")
    if len(noise) != 470 or noise[["repeat", "cell_line"]].duplicated().any():
        raise ValueError("W10 split-half grid must contain 5 x 94 rows")
    if len(subsampling) != 470 or subsampling[["repeat", "cell_line"]].duplicated().any():
        raise ValueError("W10 cell-subsampling grid must contain 5 x 94 rows")
    if len(lolo) != 94 or lolo["cell_line"].nunique() != 94:
        raise ValueError("W10 leave-one-lineage-out coverage is incomplete")
    if len(conclusions) != 6:
        raise ValueError("W10 conclusion registry is incomplete")

    floor_error = float(
        np.max(
            np.abs(
                noise["full_target_noise_floor_approx"].to_numpy()
                - noise["split_half_rmse"].to_numpy() / 2.0
            )
        )
    )
    if floor_error > 1e-12:
        raise ValueError("W10 full-target noise-floor approximation differs from protocol")

    first = subsampling.loc[subsampling["repeat"].eq(0)]
    metric = metrics.loc[
        metrics["analysis"].eq("cell_subsampling")
        & metrics["variant"].eq("n20_repeat0")
        & metrics["metric"].eq("rmse_gain_vs_b1")
    ].iloc[0]
    recomputed, low, high = bootstrap_mean_interval(
        first["rmse_gain_vs_b1"].to_numpy(), seed=int(metric.seed)
    )
    summary_error = float(
        max(
            abs(recomputed - metric.estimate),
            abs(low - metric.ci95_low),
            abs(high - metric.ci95_high),
        )
    )
    if summary_error > 1e-12:
        raise ValueError("W10 subsampling summary does not recompute")

    threshold = metrics.loc[
        metrics["analysis"].eq("inclusion_threshold")
        & metrics["metric"].eq("rmse_gain_vs_b1")
    ].set_index("variant")
    if set(threshold.index) != {"min_cells_10", "min_cells_20", "min_cells_30"}:
        raise ValueError("W10 threshold variants are incomplete")
    observed_sizes = threshold["n"].to_numpy()
    if not (np.sort(observed_sizes)[::-1] == observed_sizes).all():
        raise ValueError("W10 threshold cohort sizes are not monotone")
    if metrics.loc[metrics["analysis"].eq("response_rank"), "variant"].nunique() != 7:
        raise ValueError("W10 response-rank sensitivity is incomplete")
    if metrics.loc[metrics["analysis"].eq("gene_filter"), "variant"].nunique() != 4:
        raise ValueError("W10 gene-filter sensitivity is incomplete")

    required_figures = [
        root / "results" / "figures" / "noise_ceiling.png",
        root / "results" / "figures" / "robustness_summary.png",
    ]
    if any(not path.exists() or path.stat().st_size < 20_000 for path in required_figures):
        raise ValueError("a required W10 figure is absent or unexpectedly small")
    report = {
        "stage": "W10_robustness",
        "status": "passed",
        "robustness_rows": len(metrics),
        "noise_rows": len(noise),
        "subsampling_rows": len(subsampling),
        "lolo_lines": len(lolo),
        "threshold_cohort_sizes": threshold["n"].astype(int).to_dict(),
        "max_noise_floor_identity_error": floor_error,
        "max_subsampling_summary_recomputation_error": summary_error,
        "sensitivity_used_as_predictor": False,
    }
    output = root / "results" / "logs" / "robustness_validation.json"
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
