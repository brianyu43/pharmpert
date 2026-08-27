"""Independent consistency and leakage checks for the W7 ablation outputs."""

from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from yakseopdong.ablation import ABLATION_METRICS, VARIANT_ORDER, predict_fixed_low_rank
from yakseopdong.benchmark import _sha256
from yakseopdong.splits import validate_outer_assignments


def _bootstrap_interval(
    values: np.ndarray, replicates: int, seed: int
) -> tuple[float, float, float]:
    """Recalculate a percentile interval without the production metric helper."""
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    sampled = rng.choice(finite, size=(replicates, len(finite)), replace=True)
    bootstrap_means = sampled.mean(axis=1)
    low, high = np.quantile(bootstrap_means, [0.025, 0.975])
    return float(finite.mean()), float(low), float(high)


def _assert_close(actual: float, expected: float, label: str) -> None:
    if np.isnan(actual) and np.isnan(expected):
        return
    if not np.isclose(actual, expected, rtol=0, atol=1e-12):
        raise ValueError(f"{label} differs: {actual} != {expected}")


def validate_ablation(root: Path) -> dict[str, object]:
    """Validate W7 tables, paired metrics, fixed protocol, and output hashes."""
    with (root / "config" / "models.yaml").open(encoding="utf-8") as handle:
        models_config = yaml.safe_load(handle)
    config = models_config["w7_ablation"]
    if config["comparisons"] != VARIANT_ORDER:
        raise ValueError("configured W7 comparison list is not frozen as implemented")
    if config["selection_rule"] != "report_all_fixed_variants_without_outer_test_selection":
        raise ValueError("outer-test selection prohibition is not explicit")
    if bool(config["sensitivity_as_predictor"]):
        raise ValueError("sensitivity must not be a W7 predictor")

    outer = pd.read_csv(root / "split_assignments.csv")
    validate_outer_assignments(outer, n_splits=5)
    annotations = pd.read_csv(root / "cell_line_annotations.csv")
    if annotations[["lineage", "braf_mut", "kras_mut"]].isna().any().any():
        raise ValueError("W7 metadata fields are incomplete")
    if annotations["cell_line"].duplicated().any():
        raise ValueError("cell-line annotations contain aliases or duplicates")

    table_dir = root / "results" / "tables"
    metrics = pd.read_csv(table_dir / "ablation_metrics_by_line.csv")
    summary = pd.read_csv(table_dir / "ablation_metrics_summary.csv")
    comparison = pd.read_csv(table_dir / "ablation_metrics.csv")
    runtime = pd.read_csv(table_dir / "ablation_runtime.csv")
    coverage = pd.read_csv(table_dir / "pathway_panel_coverage.csv")
    enrichment = pd.read_csv(table_dir / "component_pathway_enrichment.csv")
    stability = pd.read_csv(table_dir / "cclr_subspace_stability.csv")
    report = json.loads((root / "results" / "logs" / "ablation_summary.json").read_text())

    if metrics.groupby(["model", "cell_line"]).size().ne(1).any():
        raise ValueError("each variant must predict each cell line exactly once")
    if set(metrics["model"]) != set(VARIANT_ORDER):
        raise ValueError("ablation metrics do not contain the exact frozen variants")
    counts = metrics.groupby("model")["cell_line"].nunique()
    if not counts.eq(94).all() or len(metrics) != 94 * len(VARIANT_ORDER):
        raise ValueError("W7 metrics do not have complete 94-line coverage")
    if set(metrics["cell_line"]) != set(outer["cell_line"]):
        raise ValueError("W7 metrics and frozen outer assignments differ")
    if metrics["outer_test_response_used_for_fit"].astype(bool).any():
        raise ValueError("a W7 row claims outer-test response use during fitting")
    signature = inspect.signature(predict_fixed_low_rank)
    if "test_response" in signature.parameters:
        raise ValueError("fixed prediction API must not accept outer-test responses")

    b4 = metrics.loc[metrics["model"].eq("B4_DIRECT_RIDGE_W5")].set_index("cell_line")
    cclr = metrics.loc[metrics["model"].eq("CCLR_NESTED_W6")].set_index("cell_line")
    for row in metrics.itertuples():
        _assert_close(
            row.rmse_gain_vs_b4,
            float(b4.loc[row.cell_line, "rmse_delta"] - row.rmse_delta),
            f"{row.model}/{row.cell_line} gain vs B4",
        )
        _assert_close(
            row.rmse_gain_vs_cclr,
            float(cclr.loc[row.cell_line, "rmse_delta"] - row.rmse_delta),
            f"{row.model}/{row.cell_line} gain vs CCLR",
        )

    indexed_summary = summary.set_index(["model", "metric"])
    replicates = int(config["bootstrap_replicates"])
    base_seed = int(config["bootstrap_seed"]) + 70_000
    for model_index, model in enumerate(VARIANT_ORDER):
        rows = metrics.loc[metrics["model"].eq(model)]
        for metric_index, metric in enumerate(ABLATION_METRICS):
            expected = _bootstrap_interval(
                rows[metric].to_numpy(dtype=float),
                replicates,
                base_seed + 100 * model_index + metric_index,
            )
            actual = indexed_summary.loc[(model, metric)]
            for value, column in zip(
                expected, ["macro_mean", "ci95_low", "ci95_high"], strict=True
            ):
                _assert_close(float(actual[column]), value, f"{model}/{metric}/{column}")
    if set(comparison["model"]) != set(VARIANT_ORDER):
        raise ValueError("wide comparison table has incomplete variants")

    if runtime.groupby(["model", "outer_fold"]).size().ne(1).any():
        raise ValueError("runtime/complexity table must have one row per variant-fold")
    if set(runtime["model"]) != set(VARIANT_ORDER) or len(runtime) != 5 * len(VARIANT_ORDER):
        raise ValueError("runtime/complexity table is incomplete")
    generated = runtime["source"].eq("w7_fixed_outer_fit")
    expected_parameters = (
        32_738 * (runtime.loc[generated, "response_rank"] + 1)
        + runtime.loc[generated, "response_rank"]
        * (runtime.loc[generated, "predictor_dimension"] + 1)
    )
    if not np.array_equal(
        runtime.loc[generated, "parameter_count"].to_numpy(dtype=int),
        expected_parameters.to_numpy(dtype=int),
    ):
        raise ValueError("fixed-variant parameter counts are inconsistent")

    union = coverage.loc[coverage["row_type"].eq("panel_union")]
    if len(union) != 1 or int(union.iloc[0]["defined_symbols"]) != 1_039:
        raise ValueError("frozen pathway union is not the expected 1,039 symbols")
    if int(union.iloc[0]["mapped_symbols"]) < 950:
        raise ValueError("pathway mapping coverage is unexpectedly low")
    fold_coverage = coverage.loc[coverage["row_type"].eq("outer_training_selection")]
    if fold_coverage["outer_fold"].nunique() != 5:
        raise ValueError("pathway training-only selection is missing a fold")
    if len(enrichment) != 5 * 20 * 2 * 6:
        raise ValueError("component enrichment does not cover all fold/rank/directions")
    if not enrichment[["p_value", "fdr_bh"]].apply(
        lambda column: column.between(0, 1).all()
    ).all():
        raise ValueError("component enrichment probabilities are outside [0, 1]")
    if len(stability) != 10:
        raise ValueError("five response subspaces require ten pairwise comparisons")
    if not stability[["mean_squared_cosine", "minimum_cosine", "maximum_cosine"]].apply(
        lambda column: column.between(0, 1 + 1e-6).all()
    ).all():
        raise ValueError("subspace stability values are outside [0, 1]")

    required_artifacts = [
        root / "results" / "figures" / "ablation_performance.png",
        root / "results" / "figures" / "complexity_vs_performance.png",
        root / "results" / "figures" / "component_diagnostics.png",
        root / "notebooks" / "05_ablation.ipynb",
    ]
    missing = [str(path.relative_to(root)) for path in required_artifacts if not path.exists()]
    if missing:
        raise ValueError(f"W7 presentation artifacts are missing: {missing}")
    comparison_sha = _sha256(table_dir / "ablation_metrics.csv")
    if comparison_sha != report["comparison_sha256"]:
        raise ValueError("ablation comparison hash differs from its run report")
    if not report["repeat_run_comparison_sha256_match"]:
        raise ValueError("W7 deterministic comparison hash has not matched on a repeat run")
    if report["outer_test_response_used_for_fit"]:
        raise ValueError("run report claims outer-test response use")
    if report["outer_test_used_for_variant_selection"]:
        raise ValueError("run report claims outer-test variant selection")

    validation = {
        "status": "pass",
        "cell_lines": 94,
        "variants": len(VARIANT_ORDER),
        "metrics_rows": len(metrics),
        "paired_metrics_recalculated": True,
        "bootstrap_intervals_recalculated": True,
        "fixed_parameter_counts_recalculated": True,
        "outer_test_response_api_absent": True,
        "pathway_panel_and_component_diagnostics_checked": True,
        "repeat_run_comparison_sha256_match": True,
        "comparison_sha256": comparison_sha,
        "validation_inputs_sha256": hashlib.sha256(
            (
                _sha256(table_dir / "ablation_metrics_by_line.csv")
                + _sha256(table_dir / "ablation_metrics_summary.csv")
                + comparison_sha
            ).encode()
        ).hexdigest(),
    }
    output = root / "results" / "logs" / "ablation_validation.json"
    output.write_text(json.dumps(validation, indent=2), encoding="utf-8")
    return validation
