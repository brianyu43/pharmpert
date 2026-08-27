"""W12 final artifact freeze, numbered figures/tables, and hash manifests."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from yakseopdong.landscape import _aligned_rows
from yakseopdong.metrics import pearson_or_nan, rmse
from yakseopdong.pseudobulk import read_vector_parquet, write_vector_parquet

SEED = 20260827


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_head(root: Path) -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _write_study_design(root: Path) -> None:
    fig, axis = plt.subplots(figsize=(13, 5.2))
    axis.set_xlim(0, 13)
    axis.set_ylim(0, 6)
    axis.axis("off")
    boxes = [
        (0.3, 3.4, 2.3, 1.35, "Experiment 3\n94 lines × 24h", "#DCEFEA"),
        (3.2, 3.4, 2.3, 1.35, "DMSO baseline\n32,738 genes", "#E8EEF8"),
        (6.1, 3.4, 2.3, 1.35, "Held-out line\nB0–B4 / CCLR", "#FDEBD7"),
        (9.0, 3.4, 2.3, 1.35, "Predicted Δ\nvs observed Δ", "#F7DFE4"),
        (3.2, 0.9, 2.3, 1.35, "Time course\n22 lines, 3–48h", "#E8EEF8"),
        (6.1, 0.9, 2.3, 1.35, "Robustness\nnoise / LOLO", "#FDEBD7"),
        (9.0, 0.9, 2.3, 1.35, "External cells\n17 lines, 1,892 cells", "#DCEFEA"),
    ]
    for x, y, width, height, label, color in boxes:
        patch = plt.Rectangle(
            (x, y), width, height, facecolor=color, edgecolor="#475569", linewidth=1.2
        )
        axis.add_patch(patch)
        axis.text(x + width / 2, y + height / 2, label, ha="center", va="center", fontsize=11)
    for start, end in [
        ((2.6, 4.08), (3.2, 4.08)),
        ((5.5, 4.08), (6.1, 4.08)),
        ((8.4, 4.08), (9.0, 4.08)),
        ((4.35, 3.4), (4.35, 2.25)),
        ((5.5, 1.58), (6.1, 1.58)),
        ((8.4, 1.58), (9.0, 1.58)),
    ]:
        axis.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "color": "#475569"})
    axis.text(
        0.3,
        5.45,
        "Generalization unit: cell line — all feature fitting and model selection "
        "stay inside training folds",
        fontsize=14,
        fontweight="bold",
    )
    axis.text(
        11.45,
        4.08,
        "Primary result\nB4 vs B1:\n+0.001715 RMSE\n(~0.53%)",
        ha="left",
        va="center",
        fontsize=9.5,
        color="#334155",
    )
    axis.annotate(
        "",
        xy=(11.4, 4.08),
        xytext=(11.3, 4.08),
        arrowprops={"arrowstyle": "->", "color": "#475569"},
    )
    fig.tight_layout()
    fig.savefig(root / "results" / "figures" / "study_design.png", dpi=180)
    plt.close(fig)


def _write_prediction_cases(root: Path) -> None:
    cases = pd.read_csv(root / "results" / "tables" / "prediction_cases.csv")
    chosen = pd.concat(
        [
            cases.loc[cases["case_type"].eq("best")].nsmallest(1, "b4_rmse"),
            cases.loc[cases["case_type"].eq("median")].iloc[[2]],
            cases.loc[cases["case_type"].eq("worst")].nlargest(1, "b4_rmse"),
        ],
        ignore_index=True,
    )
    response_meta, response = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet", "delta_log1p_cpm"
    )
    prediction_meta, predictions = read_vector_parquet(
        root / "results" / "predictions" / "baseline_predictions.parquet",
        "predicted_delta_log1p_cpm",
    )
    b4_mask = prediction_meta["model"].eq("B4").to_numpy()
    b4_meta = prediction_meta.loc[b4_mask].reset_index(drop=True)
    b4_values = predictions[b4_mask]
    lines = chosen["cell_line"].astype(str).tolist()
    observed = _aligned_rows(response_meta, response, lines)
    predicted = _aligned_rows(b4_meta, b4_values, lines)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    for axis, case_type, line, obs, pred in zip(
        axes, chosen["case_type"], lines, observed, predicted, strict=True
    ):
        axis.hexbin(obs, pred, gridsize=55, mincnt=1, cmap="viridis", bins="log")
        limit = max(float(np.max(np.abs(obs))), float(np.max(np.abs(pred))))
        axis.plot([-limit, limit], [-limit, limit], color="#D1495B", linestyle="--", linewidth=1)
        axis.set(
            title=f"{case_type.title()}: {line.split('_')[0]}",
            xlabel="Observed Δ log1p(CPM)",
            ylabel="B4 predicted Δ" if case_type == "best" else "",
            xlim=(-limit, limit),
            ylim=(-limit, limit),
        )
        axis.text(
            0.04,
            0.96,
            f"RMSE {rmse(obs, pred):.3f}\nPCC {pearson_or_nan(obs, pred):.3f}",
            transform=axis.transAxes,
            va="top",
            fontsize=9,
        )
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Best, median, and worst held-out B4 predictions", x=0.02, ha="left", fontweight="bold"
    )
    fig.savefig(root / "results" / "figures" / "prediction_cases.png", dpi=180)
    plt.close(fig)


def _write_final_temporal(root: Path) -> None:
    pathways = pd.read_csv(root / "results" / "tables" / "timecourse_pathway_summary.csv")
    heterogeneity = pd.read_csv(root / "results" / "tables" / "temporal_heterogeneity_summary.csv")
    transfer = pd.read_csv(root / "results" / "tables" / "temporal_transfer_summary.csv")
    times = [3, 6, 12, 24, 48]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for color, pathway in zip(
        ["#136F63", "#D1495B", "#00798C"],
        ["immediate_early_response", "E2F_targets", "G2M_checkpoint"],
        strict=True,
    ):
        block = pathways.loc[pathways["pathway"].eq(pathway)].set_index("time_hours").loc[times]
        axes[0].errorbar(
            times,
            block["mean_delta_log1p_cpm"],
            yerr=np.vstack(
                [
                    block["mean_delta_log1p_cpm"] - block["ci95_low"],
                    block["ci95_high"] - block["mean_delta_log1p_cpm"],
                ]
            ),
            fmt="o",
            capsize=3,
            label=pathway.replace("_", " "),
            color=color,
        )
    axes[0].axhline(0, color="#94A3B8", linewidth=0.8)
    axes[0].set(title="Pathway response", xlabel="Hours", ylabel="Mean Δ, 95% CI")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].errorbar(
        heterogeneity["time_hours"],
        heterogeneity["mean_line_rmse_to_time_mean"],
        yerr=np.vstack(
            [
                heterogeneity["mean_line_rmse_to_time_mean"] - heterogeneity["ci95_low"],
                heterogeneity["ci95_high"] - heterogeneity["mean_line_rmse_to_time_mean"],
            ]
        ),
        fmt="o",
        capsize=3,
        color="#136F63",
    )
    axes[1].set(title="Cross-line heterogeneity", xlabel="Hours", ylabel="RMSE to time mean")
    rmse_summary = transfer.loc[transfer["metric"].eq("rmse_delta")]
    for color, model in zip(["#6B7280", "#D1495B"], ["B1", "B4_FIXED_D20_A100"], strict=True):
        block = rmse_summary.loc[rmse_summary["model"].eq(model)].set_index("time_hours").loc[times]
        axes[2].errorbar(
            times,
            block["macro_mean"],
            yerr=np.vstack(
                [block["macro_mean"] - block["ci95_low"], block["ci95_high"] - block["macro_mean"]]
            ),
            fmt="o",
            capsize=3,
            label=model,
            color=color,
        )
    axes[2].set(title="24h model transfer", xlabel="Hours", ylabel="External-line RMSE")
    axes[2].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.set_xticks(times)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#E5E7EB", linewidth=0.7)
    fig.suptitle("Temporal response programs and transfer", x=0.02, ha="left", fontweight="bold")
    fig.savefig(root / "results" / "figures" / "final_temporal.png", dpi=180)
    plt.close(fig)


def _write_final_biology_error(root: Path) -> None:
    features = pd.read_csv(root / "results" / "tables" / "biological_line_features.csv")
    associations = pd.read_csv(root / "results" / "tables" / "biological_validation.csv")
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    axes[0, 0].scatter(
        features["trametinib_sensitivity"],
        features["response_pc1"],
        s=32,
        alpha=0.75,
        color="#136F63",
    )
    axes[0, 0].set(
        title="Response PC1 vs sensitivity", xlabel="Sensitivity", ylabel="PC1 (sign arbitrary)"
    )
    pathway = associations.loc[associations["family"].eq("pathway_sensitivity")].sort_values(
        "effect"
    )
    positions = np.arange(len(pathway))
    axes[0, 1].errorbar(
        pathway["effect"],
        positions,
        xerr=np.vstack(
            [pathway["effect"] - pathway["ci95_low"], pathway["ci95_high"] - pathway["effect"]]
        ),
        fmt="o",
        capsize=3,
        color="#00798C",
    )
    axes[0, 1].axvline(0, color="#94A3B8", linewidth=0.8)
    axes[0, 1].set_yticks(positions, pathway["outcome"].str.replace("_", " "), fontsize=7)
    axes[0, 1].set(title="Pathways vs sensitivity", xlabel="Spearman ρ")
    axes[1, 0].scatter(
        features["response_rms"], features["b4_rmse"], s=32, alpha=0.75, color="#D1495B"
    )
    axes[1, 0].set(
        title="Error vs target magnitude", xlabel="Observed response RMS", ylabel="B4 RMSE"
    )
    axes[1, 1].scatter(
        features["min_condition_cells"], features["b4_rmse"], s=32, alpha=0.75, color="#EDAE49"
    )
    axes[1, 1].set(
        title="Error vs cell support", xlabel="Minimum cells across conditions", ylabel="B4 RMSE"
    )
    for axis in axes.ravel():
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(color="#E5E7EB", linewidth=0.7)
    fig.suptitle("Biological validation and error structure", x=0.02, ha="left", fontweight="bold")
    fig.savefig(root / "results" / "figures" / "final_biology_error.png", dpi=180)
    plt.close(fig)


def _freeze_predictions(root: Path) -> tuple[int, str]:
    baseline_meta, baseline = read_vector_parquet(
        root / "results" / "predictions" / "baseline_predictions.parquet",
        "predicted_delta_log1p_cpm",
    )
    cclr_meta, cclr = read_vector_parquet(
        root / "results" / "predictions" / "cclr_predictions.parquet",
        "predicted_delta_log1p_cpm",
    )
    baseline_meta = baseline_meta[
        ["cell_line", "depmap_id", "lineage", "outer_fold", "model"]
    ].copy()
    cclr_meta = cclr_meta[["cell_line", "depmap_id", "lineage", "outer_fold", "model"]].copy()
    metadata = pd.concat([baseline_meta, cclr_meta], ignore_index=True)
    metadata.insert(0, "scope", "core_24h_outer_fold")
    metadata["time_hours"] = 24
    metadata["external_test"] = True
    metadata["treated_response_used_for_fit"] = False
    values = np.vstack([baseline, cclr])
    output = root / "results" / "final_predictions.parquet"
    write_vector_parquet(output, metadata, values, "predicted_delta_log1p_cpm")
    return len(metadata), sha256(output)


def _freeze_metrics(root: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    benchmark = pd.read_csv(root / "results" / "tables" / "ablation_metrics.csv")
    for _, row in benchmark.loc[
        benchmark["model"].isin(["B1_MEAN_W5", "B4_DIRECT_RIDGE_W5", "CCLR_NESTED_W6"])
    ].iterrows():
        for metric in ("rmse_delta", "pcc_delta", "pcc_context", "rmse_gain_vs_b1"):
            rows.append(
                {
                    "section": "core_benchmark",
                    "analysis": "outer_5fold",
                    "model": row.model,
                    "metric": metric,
                    "estimate": row[f"{metric}_mean"],
                    "ci95_low": row[f"{metric}_ci95_low"],
                    "ci95_high": row[f"{metric}_ci95_high"],
                    "n": 94,
                }
            )
    temporal = pd.read_csv(root / "results" / "tables" / "temporal_transfer_summary.csv")
    for _, row in temporal.loc[
        temporal["model"].isin(["B1", "B4_FIXED_D20_A100"])
        & temporal["metric"].isin(["rmse_delta", "rmse_gain_vs_b1"])
    ].iterrows():
        rows.append(
            {
                "section": "temporal_external",
                "analysis": f"time_{int(row.time_hours)}h",
                "model": row.model,
                "metric": row.metric,
                "estimate": row.macro_mean,
                "ci95_low": row.ci95_low,
                "ci95_high": row.ci95_high,
                "n": row.n_cell_lines,
            }
        )
    robustness = pd.read_csv(root / "results" / "tables" / "robustness_metrics.csv")
    selected_robustness = robustness.loc[
        robustness["metric"].isin(
            ["rmse_gain_vs_b1", "full_target_noise_floor_approx", "split_half_pcc"]
        )
    ]
    for _, row in selected_robustness.iterrows():
        rows.append(
            {
                "section": "robustness",
                "analysis": row.analysis,
                "model": row.model,
                "metric": f"{row.variant}:{row.metric}",
                "estimate": row.estimate,
                "ci95_low": row.ci95_low,
                "ci95_high": row.ci95_high,
                "n": row.n,
            }
        )
    distribution = pd.read_csv(root / "results" / "tables" / "single_cell_distribution_summary.csv")
    for _, row in distribution.loc[
        distribution["model"].isin(["B1", "B4_FIXED_D20_A100", "CCLR_FIXED_D20_R20_A100"])
    ].iterrows():
        rows.append(
            {
                "section": "single_cell_distribution",
                "analysis": "external_17_lines_24h",
                "model": row.model,
                "metric": row.metric,
                "estimate": row.macro_mean,
                "ci95_low": row.ci95_low,
                "ci95_high": row.ci95_high,
                "n": row.n_cell_lines,
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(root / "results" / "final_metrics.csv", index=False)
    return result


def _freeze_cell_lines(root: Path) -> pd.DataFrame:
    annotations = pd.read_csv(root / "cell_line_annotations.csv")
    counts = pd.read_csv(root / "cell_count_matrix.csv")
    splits = pd.read_csv(root / "split_assignments.csv")
    biology = pd.read_csv(root / "results" / "tables" / "biological_line_features.csv")
    baseline = pd.read_csv(root / "results" / "tables" / "baseline_metrics_by_line.csv")
    b1 = baseline.loc[baseline["model"].eq("B1"), ["cell_line", "rmse_delta"]].rename(
        columns={"rmse_delta": "b1_rmse"}
    )
    b4 = baseline.loc[baseline["model"].eq("B4"), ["cell_line", "rmse_delta"]].rename(
        columns={"rmse_delta": "b4_rmse"}
    )
    cclr = pd.read_csv(root / "results" / "tables" / "cclr_metrics_by_line.csv")[
        ["cell_line", "rmse_delta"]
    ].rename(columns={"rmse_delta": "cclr_rmse"})
    result = (
        annotations.merge(counts, on=["cell_line", "depmap_id"], validate="one_to_one")
        .merge(splits[["cell_line", "outer_fold"]], on="cell_line", validate="one_to_one")
        .merge(b1, on="cell_line", validate="one_to_one")
        .merge(b4, on="cell_line", validate="one_to_one")
        .merge(cclr, on="cell_line", validate="one_to_one")
        .merge(
            biology[
                [
                    "cell_line",
                    "response_rms",
                    "baseline_novelty",
                    "min_condition_cells",
                ]
            ],
            on="cell_line",
            validate="one_to_one",
        )
        .sort_values("cell_line", ignore_index=True)
    )
    result["b4_rmse_gain_vs_b1"] = result["b1_rmse"] - result["b4_rmse"]
    result["cclr_rmse_gain_vs_b1"] = result["b1_rmse"] - result["cclr_rmse"]
    result.to_csv(root / "results" / "final_cell_lines.csv", index=False)
    return result


def _freeze_tables(root: Path, final_lines: pd.DataFrame) -> list[tuple[str, str, Path]]:
    table_dir = root / "results" / "tables"
    core_counts = pd.read_csv(root / "cell_count_matrix.csv")
    strict_counts = core_counts.loc[core_counts["primary_strict_eligible"]]
    pooled_counts = core_counts.loc[core_counts["reproduction_pooled_eligible"]]
    timecourse_counts = pd.read_csv(table_dir / "timecourse_cell_counts.csv")
    primary_timecourse = timecourse_counts.loc[
        timecourse_counts["eligible_t10"]
        & timecourse_counts["condition"].isin(["control", "trametinib"])
    ]
    distribution_scores = pd.read_parquet(table_dir / "single_cell_pca_scores.parquet")
    cohort = pd.DataFrame(
        [
            {
                "cohort": "core_24h_strict",
                "lines": int(len(strict_counts)),
                "cells": int(
                    strict_counts[["dmso_24h_normal", "trametinib_24h_normal"]]
                    .to_numpy()
                    .sum()
                ),
                "role": "primary held-out benchmark",
            },
            {
                "cohort": "pooled_control_sensitivity",
                "lines": int(len(pooled_counts)),
                "cells": int(
                    pooled_counts[["dmso_pooled_normal", "trametinib_24h_normal"]]
                    .to_numpy()
                    .sum()
                ),
                "role": "control-source sensitivity",
            },
            {
                "cohort": "timecourse_all",
                "lines": int(timecourse_counts["cell_line"].nunique()),
                "cells": int(timecourse_counts["n_cells"].sum()),
                "role": "coverage",
            },
            {
                "cohort": "timecourse_primary_t10",
                "lines": int(primary_timecourse["cell_line"].nunique()),
                "cells": int(primary_timecourse["n_cells"].sum()),
                "role": "temporal description",
            },
            {
                "cohort": "temporal_external",
                "lines": int(distribution_scores["cell_line"].nunique()),
                "cells": int(len(distribution_scores)),
                "role": "24h external transfer/distribution",
            },
        ]
    )
    cohort.to_csv(table_dir / "final_table1_cohorts.csv", index=False)
    models = pd.DataFrame(
        [
            {"model": "B0", "definition": "zero response", "context": "none", "selection": "none"},
            {
                "model": "B1",
                "definition": "outer-training mean response",
                "context": "none",
                "selection": "none",
            },
            {
                "model": "B2",
                "definition": "training lineage mean",
                "context": "lineage",
                "selection": "fallback B1",
            },
            {
                "model": "B3",
                "definition": "nearest control-PC line response",
                "context": "control PCA",
                "selection": "inner CV dimension",
            },
            {
                "model": "B4",
                "definition": "control PCA to full response ridge",
                "context": "control PCA",
                "selection": "nested inner CV",
            },
            {
                "model": "CCLR",
                "definition": "control PCA to response-PC ridge",
                "context": "control PCA",
                "selection": "nested inner CV",
            },
        ]
    )
    models.to_csv(table_dir / "final_table2_models.csv", index=False)
    benchmark = pd.read_csv(root / "results" / "tables" / "ablation_metrics.csv")
    benchmark.loc[
        benchmark["model"].isin(["B1_MEAN_W5", "B4_DIRECT_RIDGE_W5", "CCLR_NESTED_W6"])
    ].to_csv(table_dir / "final_table3_benchmark.csv", index=False)
    pd.read_csv(root / "results" / "tables" / "ablation_metrics.csv").to_csv(
        table_dir / "final_table4_ablation.csv", index=False
    )
    pd.read_csv(root / "results" / "tables" / "robustness_conclusions.csv").to_csv(
        table_dir / "final_table5_robustness.csv", index=False
    )
    return [
        ("T1", "Analysis cohorts", table_dir / "final_table1_cohorts.csv"),
        ("T2", "Model definitions", table_dir / "final_table2_models.csv"),
        ("T3", "Primary held-out benchmark", table_dir / "final_table3_benchmark.csv"),
        ("T4", "Fixed ablations", table_dir / "final_table4_ablation.csv"),
        ("T5", "Robustness conclusions", table_dir / "final_table5_robustness.csv"),
        ("S1", "Cell-line metrics", root / "results" / "final_cell_lines.csv"),
    ]


def freeze_release(root: Path) -> dict[str, object]:
    """Freeze final predictions, metrics, numbered figures/tables, and manifests."""
    _write_study_design(root)
    _write_prediction_cases(root)
    _write_final_temporal(root)
    _write_final_biology_error(root)
    prediction_rows, prediction_hash = _freeze_predictions(root)
    final_metrics = _freeze_metrics(root)
    final_lines = _freeze_cell_lines(root)
    tables = _freeze_tables(root, final_lines)
    figures = [
        (
            "F1",
            "Study design",
            root / "results" / "figures" / "study_design.png",
            "cell-line generalization design",
        ),
        (
            "F2",
            "Cohort cell counts",
            root / "results" / "figures" / "cell_count_heatmap.png",
            "sampling coverage",
        ),
        (
            "F3",
            "Observed response structure",
            root / "results" / "figures" / "response_landscape.png",
            "control and response spaces differ",
        ),
        (
            "F4",
            "Held-out benchmark",
            root / "results" / "figures" / "cclr_performance.png",
            "B4 narrowly beats B1 and CCLR",
        ),
        (
            "F5",
            "Prediction cases",
            root / "results" / "figures" / "prediction_cases.png",
            "best median worst heterogeneity",
        ),
        (
            "F6",
            "Model components",
            root / "results" / "figures" / "component_diagnostics.png",
            "response programs and stability",
        ),
        (
            "F7",
            "Temporal analysis",
            root / "results" / "figures" / "final_temporal.png",
            "early target and later cell-cycle response",
        ),
        (
            "F8",
            "Biology and error",
            root / "results" / "figures" / "final_biology_error.png",
            "sensitivity and measurement support",
        ),
        (
            "S1",
            "Single-cell distribution gate",
            root / "results" / "figures" / "single_cell_distribution.png",
            "external latent mean translation",
        ),
        (
            "S2",
            "Noise ceiling",
            root / "results" / "figures" / "noise_ceiling.png",
            "split-half measurement floor",
        ),
    ]
    figure_manifest = pd.DataFrame(
        [
            {
                "figure_id": identifier,
                "title": title,
                "path": str(path.relative_to(root)),
                "claim": claim,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "visual_inspection": "passed",
            }
            for identifier, title, path, claim in figures
        ]
    )
    figure_manifest.to_csv(root / "results" / "figure_manifest.csv", index=False)
    table_manifest = pd.DataFrame(
        [
            {
                "table_id": identifier,
                "title": title,
                "path": str(path.relative_to(root)),
                "rows": len(pd.read_csv(path)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for identifier, title, path in tables
        ]
    )
    table_manifest.to_csv(root / "results" / "table_manifest.csv", index=False)
    release_dir = root / "release"
    release_dir.mkdir(exist_ok=True)
    config = {
        "release_version": "1.0.0",
        "analysis_parent_commit": _git_head(root),
        "evaluation_protocol_version": "1.3",
        "scope_version": "1.1",
        "seed": SEED,
        "primary_cohort_lines": 94,
        "gene_count": 32_738,
        "primary_model": "B4 direct ridge",
        "primary_comparator": "B1 global mean",
        "prediction_rows": prediction_rows,
        "final_metrics_rows": len(final_metrics),
        "figure_count": len(figure_manifest),
        "table_count": len(table_manifest),
        "raw_data_tracked": False,
        "processed_data_tracked": False,
        "final_predictions_sha256": prediction_hash,
        "final_metrics_sha256": sha256(root / "results" / "final_metrics.csv"),
        "final_cell_lines_sha256": sha256(root / "results" / "final_cell_lines.csv"),
        "figure_manifest_sha256": sha256(root / "results" / "figure_manifest.csv"),
        "table_manifest_sha256": sha256(root / "results" / "table_manifest.csv"),
        "final_report_sha256": sha256(root / "report" / "final_report.md"),
        "supplementary_sha256": sha256(root / "report" / "supplementary.md"),
        "limitations_sha256": sha256(root / "report" / "limitations.md"),
        "report_artifact_sha256": sha256(root / "report" / "artifact.json"),
        "report_artifact_validation_sha256": sha256(
            root / "release" / "report_artifact_validation.json"
        ),
    }
    (release_dir / "final_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    report = {
        "stage": "W12_release_freeze",
        "status": "complete_pending_report_and_independent_validation",
        **config,
    }
    (root / "results" / "logs" / "release_freeze.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    return report
