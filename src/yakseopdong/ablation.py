"""W7 fixed diagnostic ablations on the frozen held-out-cell-line splits."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any

import numpy as np
import pandas as pd
import yaml
from numpy.typing import NDArray
from sklearn.linear_model import Ridge

from yakseopdong.benchmark import (
    METRICS,
    _comparison_table,
    _line_metrics,
    _load_matrices,
    _sha256,
    _summarize_metrics,
)
from yakseopdong.models import (
    ControlEmbedding,
    ResponseEmbedding,
    fit_binary_covariate_encoder,
    fit_control_embedding,
    fit_lineage_encoder,
    fit_response_embedding,
)
from yakseopdong.pathways import (
    cclr_subspace_stability,
    component_pathway_enrichment,
    load_geneset_config,
    pathway_coverage_table,
    pathway_panel_indices,
)
from yakseopdong.splits import validate_outer_assignments

ABLATION_METRICS = [*METRICS, "rmse_gain_vs_b4", "rmse_gain_vs_cclr"]
VARIANT_ORDER = [
    "B1_MEAN_W5",
    "B4_DIRECT_RIDGE_W5",
    "CCLR_NESTED_W6",
    "FIXED_FULL_D20_R20",
    "CONTROL_D05_R20",
    "CONTROL_D10_R20",
    "CONTROL_D30_R20",
    "RANK_D20_R02",
    "RANK_D20_R05",
    "RANK_D20_R10",
    "RANK_D20_R30",
    "RANK_D20_R40",
    "RANK_D20_R50",
    "PATHWAY_D20_R20",
    "LINEAGE_D20_R20",
    "BRAF_KRAS_D20_R20",
]


@dataclass(frozen=True)
class FixedVariant:
    """One predeclared fixed W7 model variant."""

    variant: str
    family: str
    control_source: str
    control_dimension: int
    response_rank: int
    covariates: str = "none"


def fixed_variants() -> list[FixedVariant]:
    """Return all de-duplicated, generated W7 variants in report order."""
    return [
        FixedVariant("FIXED_FULL_D20_R20", "fixed_reference", "full", 20, 20),
        FixedVariant("CONTROL_D05_R20", "control_dimension", "full", 5, 20),
        FixedVariant("CONTROL_D10_R20", "control_dimension", "full", 10, 20),
        FixedVariant("CONTROL_D30_R20", "control_dimension", "full", 30, 20),
        FixedVariant("RANK_D20_R02", "response_rank", "full", 20, 2),
        FixedVariant("RANK_D20_R05", "response_rank", "full", 20, 5),
        FixedVariant("RANK_D20_R10", "response_rank", "full", 20, 10),
        FixedVariant("RANK_D20_R30", "response_rank", "full", 20, 30),
        FixedVariant("RANK_D20_R40", "response_rank", "full", 20, 40),
        FixedVariant("RANK_D20_R50", "response_rank", "full", 20, 50),
        FixedVariant("PATHWAY_D20_R20", "feature_input", "pathway", 20, 20),
        FixedVariant(
            "LINEAGE_D20_R20", "metadata_input", "full", 20, 20, "lineage"
        ),
        FixedVariant(
            "BRAF_KRAS_D20_R20",
            "metadata_input",
            "full",
            20,
            20,
            "braf_kras",
        ),
    ]


def predict_fixed_low_rank(
    train_response: NDArray[np.floating],
    train_predictors: NDArray[np.floating],
    test_predictors: NDArray[np.floating],
    response_embedding: ResponseEmbedding,
    response_rank: int,
    alpha: float,
) -> NDArray[np.float32]:
    """Fit one fixed low-rank ridge without accepting an outer-test response."""
    response = np.asarray(train_response, dtype=np.float64)
    train_x = np.asarray(train_predictors, dtype=np.float64)
    test_x = np.asarray(test_predictors, dtype=np.float64)
    if len(response) != len(train_x) or train_x.shape[1] != test_x.shape[1]:
        raise ValueError("training rows or predictor columns do not align")
    if not 0 < response_rank <= response_embedding.pca.n_components_:
        raise ValueError("response_rank is not feasible for the fitted response basis")
    response_scores = response_embedding.transform(response)[:, :response_rank]
    ridge = Ridge(alpha=float(alpha), fit_intercept=True, solver="svd")
    ridge.fit(train_x, response_scores)
    predicted_scores = np.asarray(ridge.predict(test_x), dtype=np.float64)
    if predicted_scores.ndim == 1:
        predicted_scores = predicted_scores[:, None]
    padded = np.zeros(
        (len(test_x), response_embedding.pca.n_components_), dtype=np.float64
    )
    padded[:, :response_rank] = predicted_scores
    return response_embedding.inverse_transform(padded).astype(np.float32)


def _predict_variant(
    variant: FixedVariant,
    train_response: NDArray[np.floating],
    full_embedding: ControlEmbedding,
    pathway_embedding: ControlEmbedding,
    response_embedding: ResponseEmbedding,
    train_control: NDArray[np.floating],
    test_control: NDArray[np.floating],
    train_lineages: NDArray[np.str_],
    test_lineages: NDArray[np.str_],
    train_mutations: NDArray[np.floating],
    test_mutations: NDArray[np.floating],
    alpha: float,
) -> tuple[NDArray[np.float32], int, int]:
    embedding = pathway_embedding if variant.control_source == "pathway" else full_embedding
    train_x = embedding.transform(train_control)[:, : variant.control_dimension]
    test_x = embedding.transform(test_control)[:, : variant.control_dimension]
    if variant.covariates == "lineage":
        encoder = fit_lineage_encoder(train_lineages)
        train_x = np.column_stack([train_x, encoder.transform(train_lineages)])
        test_x = np.column_stack([test_x, encoder.transform(test_lineages)])
    elif variant.covariates == "braf_kras":
        encoder = fit_binary_covariate_encoder(train_mutations)
        train_x = np.column_stack([train_x, encoder.transform(train_mutations)])
        test_x = np.column_stack([test_x, encoder.transform(test_mutations)])
    elif variant.covariates != "none":
        raise ValueError(f"unknown covariate mode: {variant.covariates}")
    prediction = predict_fixed_low_rank(
        train_response,
        train_x,
        test_x,
        response_embedding,
        response_rank=variant.response_rank,
        alpha=alpha,
    )
    return prediction, train_x.shape[1], len(embedding.gene_indices)


def _snapshot_metrics(root: Path) -> pd.DataFrame:
    baseline = pd.read_csv(root / "results" / "tables" / "baseline_metrics_by_line.csv")
    cclr = pd.read_csv(root / "results" / "tables" / "cclr_metrics_by_line.csv")
    snapshots = [
        baseline.loc[baseline["model"].eq("B1")].assign(model="B1_MEAN_W5"),
        baseline.loc[baseline["model"].eq("B4")].assign(model="B4_DIRECT_RIDGE_W5"),
        cclr.assign(model="CCLR_NESTED_W6"),
    ]
    combined = pd.concat(snapshots, ignore_index=True)
    for variant in ("B1_MEAN_W5", "B4_DIRECT_RIDGE_W5", "CCLR_NESTED_W6"):
        if combined.loc[combined["model"].eq(variant), "cell_line"].nunique() != 94:
            raise ValueError(f"existing benchmark snapshot is incomplete for {variant}")
    b4 = combined.loc[combined["model"].eq("B4_DIRECT_RIDGE_W5")].set_index(
        "cell_line"
    )
    cclr_rows = combined.loc[combined["model"].eq("CCLR_NESTED_W6")].set_index(
        "cell_line"
    )
    combined["rmse_gain_vs_b4"] = [
        float(b4.loc[line, "rmse_delta"] - value)
        for line, value in zip(combined["cell_line"], combined["rmse_delta"], strict=True)
    ]
    combined["rmse_gain_vs_cclr"] = [
        float(cclr_rows.loc[line, "rmse_delta"] - value)
        for line, value in zip(combined["cell_line"], combined["rmse_delta"], strict=True)
    ]
    combined["outer_test_response_used_for_fit"] = False
    return combined


def _snapshot_runtime(root: Path, genes: int) -> list[dict[str, object]]:
    baseline_runtime = pd.read_csv(root / "results" / "tables" / "baseline_runtime.csv")
    cclr_runtime = pd.read_csv(root / "results" / "tables" / "cclr_runtime.csv")
    rows: list[dict[str, object]] = []
    for fold in range(5):
        rows.append(
            {
                "outer_fold": fold,
                "model": "B1_MEAN_W5",
                "family": "benchmark_snapshot",
                "control_dimension": 0,
                "response_rank": 0,
                "predictor_dimension": 0,
                "selected_control_genes": 0,
                "parameter_count": genes,
                "fit_runtime_seconds": np.nan,
                "source": "existing_w5",
            }
        )
        b4 = baseline_runtime.loc[
            baseline_runtime["outer_fold"].eq(fold)
            & baseline_runtime["model"].eq("B4")
        ].iloc[0]
        b4_hyper = pd.read_csv(
            root / "results" / "tables" / "baseline_hyperparameters.csv"
        ).set_index("outer_fold").loc[fold]
        rows.append(
            {
                "outer_fold": fold,
                "model": "B4_DIRECT_RIDGE_W5",
                "family": "benchmark_snapshot",
                "control_dimension": int(b4_hyper["b4_dimension"]),
                "response_rank": genes,
                "predictor_dimension": int(b4_hyper["b4_dimension"]),
                "selected_control_genes": int(b4["selected_control_genes"]),
                "parameter_count": int(b4["parameter_count"]),
                "fit_runtime_seconds": np.nan,
                "source": "existing_w5",
            }
        )
        cclr = cclr_runtime.loc[cclr_runtime["outer_fold"].eq(fold)].iloc[0]
        cclr_hyper = pd.read_csv(
            root / "results" / "tables" / "cclr_hyperparameters.csv"
        ).set_index("outer_fold").loc[fold]
        rows.append(
            {
                "outer_fold": fold,
                "model": "CCLR_NESTED_W6",
                "family": "benchmark_snapshot",
                "control_dimension": int(cclr_hyper["control_dimension"]),
                "response_rank": int(cclr_hyper["response_rank"]),
                "predictor_dimension": int(cclr_hyper["control_dimension"]),
                "selected_control_genes": int(cclr["selected_control_genes"]),
                "parameter_count": int(cclr["parameter_count"]),
                "fit_runtime_seconds": np.nan,
                "source": "existing_w6",
            }
        )
    return rows


def _variant_table(runtime: pd.DataFrame, comparison: pd.DataFrame) -> pd.DataFrame:
    aggregates = (
        runtime.groupby(["model", "family"], sort=False)
        .agg(
            folds=("outer_fold", "nunique"),
            control_dimension_min=("control_dimension", "min"),
            control_dimension_max=("control_dimension", "max"),
            response_rank_min=("response_rank", "min"),
            response_rank_max=("response_rank", "max"),
            predictor_dimension_mean=("predictor_dimension", "mean"),
            selected_control_genes_mean=("selected_control_genes", "mean"),
            parameter_count_mean=("parameter_count", "mean"),
        )
        .reset_index()
    )
    metrics = comparison[
        ["model", "rmse_delta_mean", "rmse_gain_vs_b1_mean", "rmse_gain_vs_b4_mean"]
    ]
    output = aggregates.merge(metrics, on="model", validate="one_to_one")
    order = {name: index for index, name in enumerate(VARIANT_ORDER)}
    output["report_order"] = output["model"].map(order)
    return output.sort_values("report_order", ignore_index=True)


def _config_hash(config: dict[str, Any]) -> str:
    return hashlib.sha256(
        yaml.safe_dump(config, sort_keys=True).encode("utf-8")
    ).hexdigest()


def run_ablation(root: Path, seed: int = 20260827) -> dict[str, object]:
    """Run all frozen W7 ablations without outer-test model selection."""
    with (root / "config" / "models.yaml").open(encoding="utf-8") as handle:
        model_config = yaml.safe_load(handle)
    config = model_config["w7_ablation"]
    if config["comparisons"] != VARIANT_ORDER:
        raise ValueError("W7 config comparison order differs from the implementation")
    feature_config = model_config["baseline_protocol"]["feature_selection"]
    geneset_config = load_geneset_config(root / "config" / "genesets.yaml")
    alpha = float(config["fixed_reference"]["ridge_alpha"])

    annotations = pd.read_csv(root / "cell_line_annotations.csv").sort_values(
        "cell_line", ignore_index=True
    )
    outer = pd.read_csv(root / "split_assignments.csv")
    validate_outer_assignments(outer, n_splits=5)
    annotations = annotations.merge(
        outer[["cell_line", "outer_fold"]], on="cell_line", validate="one_to_one"
    ).sort_values("cell_line", ignore_index=True)
    if annotations[["lineage", "braf_mut", "kras_mut"]].isna().any().any():
        raise ValueError("W7 metadata covariates must be complete")
    cell_lines = annotations["cell_line"].astype(str).tolist()
    control, response = _load_matrices(root, cell_lines)
    genes = pd.read_parquet(root / "data" / "processed" / "gene_metadata.parquet")
    if len(genes) != response.shape[1]:
        raise ValueError("gene metadata and response matrix do not align")
    pathway_indices = pathway_panel_indices(genes, geneset_config)

    metric_frames = [_snapshot_metrics(root)]
    runtime_rows = _snapshot_runtime(root, response.shape[1])
    selected_pathway_counts: dict[int, int] = {}
    started = perf_counter()
    for outer_fold in range(5):
        test_mask = annotations["outer_fold"].eq(outer_fold).to_numpy()
        train_mask = ~test_mask
        fold_seed = seed + 1_000 * outer_fold
        full_embedding = fit_control_embedding(
            control[train_mask],
            max_components=max(config["control_dimensions"]),
            max_variable_genes=int(feature_config["max_variable_genes"]),
            min_mean_log1p_cpm=float(feature_config["min_mean_log1p_cpm"]),
            seed=fold_seed,
        )
        pathway_embedding = fit_control_embedding(
            control[train_mask],
            max_components=int(config["fixed_reference"]["control_dimension"]),
            max_variable_genes=int(feature_config["max_variable_genes"]),
            min_mean_log1p_cpm=float(feature_config["min_mean_log1p_cpm"]),
            seed=fold_seed + 1,
            candidate_gene_indices=pathway_indices,
        )
        response_embedding = fit_response_embedding(
            response[train_mask],
            max_components=max(config["response_ranks"]),
            seed=fold_seed + 2,
        )
        selected_pathway_counts[outer_fold] = len(pathway_embedding.gene_indices)
        mean_response = response[train_mask].mean(axis=0, dtype=np.float64)
        test_indices = np.flatnonzero(test_mask)
        for variant in fixed_variants():
            variant_started = perf_counter()
            prediction, predictor_dimension, selected_control_genes = _predict_variant(
                variant=variant,
                train_response=response[train_mask],
                full_embedding=full_embedding,
                pathway_embedding=pathway_embedding,
                response_embedding=response_embedding,
                train_control=control[train_mask],
                test_control=control[test_mask],
                train_lineages=annotations.loc[train_mask, "lineage"].to_numpy(dtype=str),
                test_lineages=annotations.loc[test_mask, "lineage"].to_numpy(dtype=str),
                train_mutations=annotations.loc[
                    train_mask, ["braf_mut", "kras_mut"]
                ].to_numpy(dtype=float),
                test_mutations=annotations.loc[
                    test_mask, ["braf_mut", "kras_mut"]
                ].to_numpy(dtype=float),
                alpha=alpha,
            )
            variant_seconds = perf_counter() - variant_started
            rows: list[dict[str, object]] = []
            for test_position, global_index in enumerate(test_indices):
                rows.append(
                    {
                        "cell_line": annotations.loc[global_index, "cell_line"],
                        "depmap_id": annotations.loc[global_index, "depmap_id"],
                        "lineage": annotations.loc[global_index, "lineage"],
                        "outer_fold": outer_fold,
                        "model": variant.variant,
                        **_line_metrics(
                            response[global_index],
                            prediction[test_position],
                            mean_response,
                        ),
                        "outer_test_response_used_for_fit": False,
                    }
                )
            metric_frames.append(pd.DataFrame(rows))
            parameter_count = response.shape[1] * (variant.response_rank + 1) + (
                variant.response_rank * (predictor_dimension + 1)
            )
            runtime_rows.append(
                {
                    "outer_fold": outer_fold,
                    "model": variant.variant,
                    "family": variant.family,
                    "control_dimension": variant.control_dimension,
                    "response_rank": variant.response_rank,
                    "predictor_dimension": predictor_dimension,
                    "selected_control_genes": selected_control_genes,
                    "parameter_count": parameter_count,
                    "fit_runtime_seconds": variant_seconds,
                    "source": "w7_fixed_outer_fit",
                }
            )

    metrics = pd.concat(metric_frames, ignore_index=True)
    b4 = metrics.loc[metrics["model"].eq("B4_DIRECT_RIDGE_W5")].set_index("cell_line")
    cclr = metrics.loc[metrics["model"].eq("CCLR_NESTED_W6")].set_index("cell_line")
    generated = ~metrics["model"].isin(
        ["B1_MEAN_W5", "B4_DIRECT_RIDGE_W5", "CCLR_NESTED_W6"]
    )
    metrics.loc[generated, "rmse_gain_vs_b4"] = [
        float(b4.loc[line, "rmse_delta"] - value)
        for line, value in zip(
            metrics.loc[generated, "cell_line"],
            metrics.loc[generated, "rmse_delta"],
            strict=True,
        )
    ]
    metrics.loc[generated, "rmse_gain_vs_cclr"] = [
        float(cclr.loc[line, "rmse_delta"] - value)
        for line, value in zip(
            metrics.loc[generated, "cell_line"],
            metrics.loc[generated, "rmse_delta"],
            strict=True,
        )
    ]
    order = {name: index for index, name in enumerate(VARIANT_ORDER)}
    metrics["report_order"] = metrics["model"].map(order)
    metrics = metrics.sort_values(["report_order", "cell_line"], ignore_index=True)
    summary = _summarize_metrics(
        metrics,
        seed=int(config["bootstrap_seed"]) + 70_000,
        models=VARIANT_ORDER,
        metric_names=ABLATION_METRICS,
        bootstrap_replicates=int(config["bootstrap_replicates"]),
    )
    comparison = _comparison_table(
        summary, models=VARIANT_ORDER, metric_names=ABLATION_METRICS
    )
    runtime = pd.DataFrame(runtime_rows)
    variants = _variant_table(runtime, comparison)
    coverage = pathway_coverage_table(genes, geneset_config, selected_pathway_counts)
    loadings = pd.read_csv(
        root / "results" / "tables" / "cclr_component_top_loadings.csv"
    )
    enrichment = component_pathway_enrichment(loadings, genes, geneset_config)
    stability = cclr_subspace_stability(
        sorted((root / "results" / "models").glob("cclr_outer_fold_*.npz"))
    )

    table_dir = root / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = table_dir / "ablation_metrics_by_line.csv"
    comparison_path = table_dir / "ablation_metrics.csv"
    prior_sha = _sha256(comparison_path) if comparison_path.exists() else None
    metrics.to_csv(metrics_path, index=False)
    summary.to_csv(table_dir / "ablation_metrics_summary.csv", index=False)
    comparison.to_csv(comparison_path, index=False)
    runtime.to_csv(table_dir / "ablation_runtime.csv", index=False)
    variants.to_csv(table_dir / "ablation_variants.csv", index=False)
    coverage.to_csv(table_dir / "pathway_panel_coverage.csv", index=False)
    enrichment.to_csv(table_dir / "component_pathway_enrichment.csv", index=False)
    stability.to_csv(table_dir / "cclr_subspace_stability.csv", index=False)

    from yakseopdong.plots import write_ablation_figures

    write_ablation_figures(root, comparison, variants, enrichment, stability)
    comparison_sha = _sha256(comparison_path)
    generated_names = [variant.variant for variant in fixed_variants()]
    best_fixed = (
        comparison.loc[comparison["model"].isin(generated_names)]
        .sort_values("rmse_delta_mean")
        .iloc[0]
    )
    report = {
        "status": "complete",
        "protocol_version": str(config["protocol_version"]),
        "mode": str(config["mode"]),
        "cell_lines": len(annotations),
        "genes": int(response.shape[1]),
        "variants": len(VARIANT_ORDER),
        "outer_folds": 5,
        "metrics_rows": len(metrics),
        "pathway_panel_defined_symbols": int(
            coverage.loc[coverage["row_type"].eq("panel_union"), "defined_symbols"].iloc[0]
        ),
        "pathway_panel_mapped_columns": len(pathway_indices),
        "best_fixed_diagnostic_variant_not_promoted": str(best_fixed["model"]),
        "best_fixed_diagnostic_rmse": float(best_fixed["rmse_delta_mean"]),
        "total_seconds": perf_counter() - started,
        "comparison_sha256": comparison_sha,
        "repeat_run_comparison_sha256_match": bool(
            prior_sha is not None and prior_sha == comparison_sha
        ),
        "w7_config_sha256": _config_hash(config),
        "genesets_config_sha256": _sha256(root / "config" / "genesets.yaml"),
        "outer_test_response_used_for_fit": False,
        "outer_test_used_for_variant_selection": False,
        "sensitivity_used_as_predictor": False,
        "mutation_columns_used": ["braf_mut", "kras_mut"],
    }
    log_path = root / "results" / "logs" / "ablation_summary.json"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
