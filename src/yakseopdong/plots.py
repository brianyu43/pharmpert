"""Static, reproducible research figures for pseudobulk QC."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BLUE = "#2F6B9A"
GOLD = "#C7922B"
INK = "#222831"
GRID = "#D9DEE5"


def _finish(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_qc_figures(root: Path, control_qc: pd.DataFrame, markers: pd.DataFrame) -> None:
    """Write the control-time and marker-response figures."""
    figure_dir = root / "results" / "figures"

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    axes[0].hist(
        control_qc["dmso_6h_vs_24h_pcc"], bins=14, color=BLUE, edgecolor="white"
    )
    axes[0].set_title("DMSO 6h–24h pseudobulk correlation", loc="left", color=INK)
    axes[0].set_xlabel("Pearson correlation across QC genes")
    axes[0].set_ylabel("Cell lines")

    strict = control_qc.loc[control_qc["primary_strict_eligible"]]
    axes[1].scatter(
        strict["dmso_6h_vs_24h_rmse"],
        strict["trametinib_24h_vs_dmso_24h_rmse"],
        s=28,
        color=BLUE,
        alpha=0.78,
        edgecolor="white",
        linewidth=0.4,
    )
    low = min(
        strict["dmso_6h_vs_24h_rmse"].min(),
        strict["trametinib_24h_vs_dmso_24h_rmse"].min(),
    )
    high = max(
        strict["dmso_6h_vs_24h_rmse"].max(),
        strict["trametinib_24h_vs_dmso_24h_rmse"].max(),
    )
    axes[1].plot([low, high], [low, high], linestyle="--", color=INK, linewidth=1)
    axes[1].set_title("Control-time vs treatment response RMSE", loc="left", color=INK)
    axes[1].set_xlabel("DMSO 6h vs 24h RMSE")
    axes[1].set_ylabel("Trametinib 24h vs DMSO 24h RMSE")

    for axis in axes:
        axis.grid(axis="y", color=GRID, linewidth=0.7, alpha=0.8)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Control-time quality check · 97 lines; treatment comparison · strict 94 lines",
        x=0.06,
        y=1.02,
        ha="left",
        fontsize=10,
        color="#59636E",
    )
    fig.tight_layout()
    _finish(fig, figure_dir / "control_time_qc.png")

    ordered = markers.sort_values("mean_delta_log1p_cpm")
    positions = np.arange(len(ordered))
    means = ordered["mean_delta_log1p_cpm"].to_numpy()
    errors = np.vstack(
        [
            means - ordered["ci95_low"].to_numpy(),
            ordered["ci95_high"].to_numpy() - means,
        ]
    )
    colors = [BLUE if value < 0 else GOLD for value in means]
    fig, axis = plt.subplots(figsize=(7.4, 4.8))
    axis.errorbar(
        means,
        positions,
        xerr=errors,
        fmt="none",
        ecolor=INK,
        elinewidth=1,
        capsize=3,
        zorder=1,
    )
    axis.scatter(means, positions, s=55, c=colors, edgecolor=INK, linewidth=0.5, zorder=2)
    axis.axvline(0, color=INK, linewidth=1)
    axis.set_yticks(positions, ordered["gene_symbol"])
    axis.set_xlabel("Mean Δ log1p(CPM), trametinib 24h − DMSO 24h")
    fig.suptitle(
        "Immediate-early MAPK marker response",
        x=0.125,
        y=0.98,
        ha="left",
        color=INK,
        fontsize=14,
    )
    fig.text(
        0.125,
        0.925,
        "Strict 94-line cohort; bars are 95% cell-line bootstrap intervals",
        color="#59636E",
        fontsize=9,
    )
    axis.grid(axis="x", color=GRID, linewidth=0.7, alpha=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=[0, 0, 1, 0.89])
    _finish(fig, figure_dir / "marker_response.png")


def write_cell_count_heatmap(root: Path, counts: pd.DataFrame) -> None:
    """Write the first-stage cell-line by condition count heatmap."""
    columns = ["dmso_6h_normal", "dmso_24h_normal", "trametinib_24h_normal"]
    matrix = np.log1p(counts.set_index("cell_line")[columns].to_numpy(dtype=float))
    fig, axis = plt.subplots(figsize=(6.3, 11))
    image = axis.imshow(matrix, aspect="auto", cmap="Blues", interpolation="nearest")
    axis.set_xticks(range(3), ["DMSO 6h", "DMSO 24h", "Trametinib 24h"])
    tick_positions = np.arange(0, len(counts), 8)
    axis.set_yticks(tick_positions, counts.iloc[tick_positions]["cell_line"], fontsize=7)
    fig.suptitle(
        "Normal-cell counts by experiment 3 condition",
        x=0.17,
        y=0.995,
        ha="left",
        color=INK,
        fontsize=14,
    )
    fig.text(
        0.17,
        0.974,
        "97 cell lines; color scale is log1p(cell count)",
        color="#59636E",
        fontsize=9,
    )
    colorbar = fig.colorbar(image, ax=axis, shrink=0.45, pad=0.03)
    colorbar.set_label("log1p(normal cells)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _finish(fig, root / "results" / "figures" / "cell_count_heatmap.png")


def _display_lineages(lineages: pd.Series, minimum_count: int = 4) -> pd.Series:
    counts = lineages.value_counts()
    common = set(counts[counts >= minimum_count].index)
    return lineages.where(lineages.isin(common), "other / rare")


def write_landscape_figures(
    root: Path,
    landscape: pd.DataFrame,
    directions: pd.DataFrame,
    report: dict[str, object],
) -> None:
    """Write descriptive PCA and perturbation-direction figures."""
    figure_dir = root / "results" / "figures"
    display_lineage = _display_lineages(landscape["lineage"])
    categories = sorted(display_lineage.unique())
    palette = plt.get_cmap("tab10")
    colors = {category: palette(index % 10) for index, category in enumerate(categories)}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    panels = [
        (
            "control_pc1",
            "control_pc2",
            "Baseline DMSO 24h PCA",
            "control_pc_explained_variance_ratio",
        ),
        (
            "response_pc1",
            "response_pc2",
            "Trametinib response PCA",
            "response_pc_explained_variance_ratio",
        ),
    ]
    for axis, (x_column, y_column, title, variance_key) in zip(axes, panels, strict=True):
        variance = report[variance_key]
        for category in categories:
            mask = display_lineage.eq(category)
            axis.scatter(
                landscape.loc[mask, x_column],
                landscape.loc[mask, y_column],
                s=34,
                alpha=0.82,
                color=colors[category],
                edgecolor="white",
                linewidth=0.4,
                label=category,
            )
        axis.set_title(title, loc="left", color=INK)
        axis.set_xlabel(f"PC1 ({100 * variance[0]:.1f}%)")
        axis.set_ylabel(f"PC2 ({100 * variance[1]:.1f}%)")
        axis.grid(color=GRID, linewidth=0.7, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.suptitle(
        "Control context and transcriptional response occupy different low-dimensional spaces",
        x=0.07,
        y=1.02,
        ha="left",
        fontsize=13,
        color=INK,
    )
    fig.tight_layout()
    _finish(fig, figure_dir / "response_landscape.png")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].hist(
        directions["control_vs_drug_cosine"].dropna(),
        bins=16,
        color=BLUE,
        edgecolor="white",
    )
    axes[0].axvline(0, color=INK, linewidth=1)
    axes[0].set_title("Direction similarity", loc="left", color=INK)
    axes[0].set_xlabel("Cosine: DMSO 6→24h change vs drug response")
    axes[0].set_ylabel("Cell lines")

    axes[1].scatter(
        directions["control_6h_to_24h_rmse"],
        directions["drug_24h_response_rmse"],
        s=34,
        color=BLUE,
        alpha=0.8,
        edgecolor="white",
        linewidth=0.4,
    )
    low = float(
        min(
            directions["control_6h_to_24h_rmse"].min(),
            directions["drug_24h_response_rmse"].min(),
        )
    )
    high = float(
        max(
            directions["control_6h_to_24h_rmse"].max(),
            directions["drug_24h_response_rmse"].max(),
        )
    )
    axes[1].plot([low, high], [low, high], linestyle="--", color=INK, linewidth=1)
    axes[1].set_title("Magnitude comparison", loc="left", color=INK)
    axes[1].set_xlabel("DMSO 6→24h RMSE")
    axes[1].set_ylabel("Trametinib response RMSE")
    for axis in axes:
        axis.grid(color=GRID, linewidth=0.7, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Shared-anchor diagnostic: control-time/source and drug-response directions oppose",
        x=0.08,
        y=1.02,
        ha="left",
        fontsize=13,
        color=INK,
    )
    fig.text(
        0.08,
        -0.02,
        "Caution: both contrasts contain DMSO 24h with opposite signs; "
        "this is not an independent causal comparison.",
        color="#59636E",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    _finish(fig, figure_dir / "direction_comparison.png")

    fig, axis = plt.subplots(figsize=(6.4, 4.8))
    axis.scatter(
        landscape["trametinib_sensitivity"],
        landscape["response_pc1"],
        s=38,
        color=BLUE,
        alpha=0.82,
        edgecolor="white",
        linewidth=0.4,
    )
    axis.set_xlabel("Trametinib sensitivity (1 − author-combined AUC)")
    axis.set_ylabel("Response PC1 score")
    axis.set_title("Response PC1 and external drug sensitivity", loc="left", color=INK)
    axis.text(
        0.02,
        0.98,
        f"Pearson r = {report['response_pc1_sensitivity_pearson']:.3f}\n"
        f"Spearman ρ = {report['response_pc1_sensitivity_spearman']:.3f}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        color="#59636E",
    )
    axis.grid(color=GRID, linewidth=0.7, alpha=0.7)
    axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    _finish(fig, figure_dir / "response_pc1_sensitivity.png")


def write_baseline_figure(root: Path, comparison: pd.DataFrame) -> None:
    """Write macro performance and B1-gain panels for B0-B4."""
    models = comparison["model"].tolist()
    positions = np.arange(len(models))
    colors = ["#9AA3AD", GOLD, "#5A8F62", "#7A6FAC", BLUE]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5))

    for axis, metric, title, ylabel in [
        (axes[0], "rmse_delta", "Held-out error", "RMSE-Δ (lower is better)"),
        (axes[1], "rmse_gain_vs_b1", "Improvement over B1", "RMSE gain (higher is better)"),
    ]:
        means = comparison[f"{metric}_mean"].to_numpy(dtype=float)
        lower = means - comparison[f"{metric}_ci95_low"].to_numpy(dtype=float)
        upper = comparison[f"{metric}_ci95_high"].to_numpy(dtype=float) - means
        axis.bar(positions, means, color=colors, width=0.72)
        axis.errorbar(
            positions,
            means,
            yerr=np.vstack([lower, upper]),
            fmt="none",
            ecolor=INK,
            capsize=3,
            linewidth=1,
        )
        if metric.endswith("gain_vs_b1"):
            axis.axhline(0, color=INK, linewidth=1)
        axis.set_title(title, loc="left", color=INK)
        axis.set_ylabel(ylabel)
        axis.set_xticks(positions, models)

    context = comparison["pcc_context_mean"].to_numpy(dtype=float)
    finite = np.isfinite(context)
    axes[2].bar(positions[finite], context[finite], color=np.asarray(colors)[finite], width=0.72)
    axes[2].axhline(0, color=INK, linewidth=1)
    axes[2].set_title("Cell-line-specific signal", loc="left", color=INK)
    axes[2].set_ylabel("PCC-context (higher is better)")
    axes[2].set_xticks(positions, models)
    for axis in axes:
        axis.grid(axis="y", color=GRID, linewidth=0.7, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Five-fold held-out cell-line baseline benchmark",
        x=0.07,
        y=1.02,
        ha="left",
        fontsize=13,
        color=INK,
    )
    fig.tight_layout()
    _finish(fig, root / "results" / "figures" / "baseline_performance.png")


def write_cclr_figures(
    root: Path,
    comparison: pd.DataFrame,
    cclr_metrics: pd.DataFrame,
    hyperparameters: pd.DataFrame,
    component_summary: pd.DataFrame,
) -> None:
    """Write honest CCLR benchmark and response-component diagnostic figures."""
    selected_models = ["B1", "B4", "CCLR"]
    selected = comparison.set_index("model").loc[selected_models]
    positions = np.arange(len(selected_models))
    colors = [GOLD, BLUE, "#7A6FAC"]
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.6))

    rmse_means = selected["rmse_delta_mean"].to_numpy(dtype=float)
    rmse_errors = np.vstack(
        [
            rmse_means - selected["rmse_delta_ci95_low"].to_numpy(dtype=float),
            selected["rmse_delta_ci95_high"].to_numpy(dtype=float) - rmse_means,
        ]
    )
    axes[0].errorbar(
        positions,
        rmse_means,
        yerr=rmse_errors,
        fmt="none",
        ecolor=INK,
        capsize=4,
        linewidth=1,
    )
    axes[0].scatter(positions, rmse_means, s=70, c=colors, edgecolor=INK, linewidth=0.5)
    for position, value in zip(positions, rmse_means, strict=True):
        axes[0].text(position, value + 0.0008, f"{value:.4f}", ha="center", fontsize=8)
    axes[0].set_title("Held-out response error", loc="left", color=INK)
    axes[0].set_ylabel("Macro RMSE-Δ · focused scale")
    axes[0].set_xticks(positions, selected_models)
    padding = max(float(np.ptp(rmse_means)) * 0.8, 0.002)
    axes[0].set_ylim(float(rmse_means.min() - padding), float(rmse_means.max() + padding))

    gain_models = ["B4", "CCLR"]
    gain_selected = selected.loc[gain_models]
    gain_positions = np.arange(len(gain_models))
    gains = gain_selected["rmse_gain_vs_b1_mean"].to_numpy(dtype=float)
    gain_errors = np.vstack(
        [
            gains
            - gain_selected["rmse_gain_vs_b1_ci95_low"].to_numpy(dtype=float),
            gain_selected["rmse_gain_vs_b1_ci95_high"].to_numpy(dtype=float) - gains,
        ]
    )
    axes[1].axhline(0, color=INK, linewidth=1)
    axes[1].errorbar(
        gain_positions,
        gains,
        yerr=gain_errors,
        fmt="none",
        ecolor=INK,
        capsize=4,
        linewidth=1,
    )
    axes[1].scatter(
        gain_positions,
        gains,
        s=70,
        c=[BLUE, "#7A6FAC"],
        edgecolor=INK,
        linewidth=0.5,
    )
    for position, value in zip(gain_positions, gains, strict=True):
        axes[1].text(position, value + 0.00035, f"{value:+.4f}", ha="center", fontsize=8)
    axes[1].set_title("Improvement over B1", loc="left", color=INK)
    axes[1].set_ylabel("Paired RMSE gain (higher is better)")
    axes[1].set_xticks(gain_positions, gain_models)

    paired_gain = cclr_metrics["rmse_gain_vs_b4"].to_numpy(dtype=float)
    axes[2].hist(paired_gain, bins=16, color="#7A6FAC", edgecolor="white")
    axes[2].axvline(0, color=INK, linewidth=1)
    axes[2].axvline(
        paired_gain.mean(), color=GOLD, linewidth=1.5, linestyle="--", label="macro mean"
    )
    axes[2].set_title("CCLR paired difference vs B4", loc="left", color=INK)
    axes[2].set_xlabel("Per-cell-line RMSE gain")
    axes[2].set_ylabel("Cell lines")
    axes[2].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.grid(axis="y", color=GRID, linewidth=0.7, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Context-conditioned low-rank response benchmark",
        x=0.06,
        y=1.02,
        ha="left",
        fontsize=13,
        color=INK,
    )
    fig.text(
        0.06,
        -0.015,
        "Strict 94-line cohort · lineage-aware outer 5-fold · 95% paired cell-line bootstrap CI",
        color="#59636E",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    _finish(fig, root / "results" / "figures" / "cclr_performance.png")

    folds = hyperparameters["outer_fold"].to_numpy(dtype=int)
    component_pivot = component_summary.pivot(
        index="outer_fold", columns="component", values="observed_vs_predicted_pcc"
    ).sort_index()
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.5))
    width = 0.36
    axes[0].bar(
        folds - width / 2,
        hyperparameters["control_dimension"],
        width=width,
        color=BLUE,
        label="control PCs",
    )
    axes[0].bar(
        folds + width / 2,
        hyperparameters["response_rank"],
        width=width,
        color="#7A6FAC",
        label="response rank",
    )
    axes[0].set_title("Selected model dimensions", loc="left", color=INK)
    axes[0].set_xlabel("Outer fold")
    axes[0].set_ylabel("Dimensions")
    axes[0].set_xticks(folds)
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].bar(
        folds,
        100 * hyperparameters["response_variance_explained"],
        color="#7A6FAC",
        edgecolor=INK,
        linewidth=0.4,
    )
    axes[1].set_title("Selected response subspace", loc="left", color=INK)
    axes[1].set_xlabel("Outer fold")
    axes[1].set_ylabel("Training-response variance explained (%)")
    axes[1].set_xticks(folds)

    matrix = component_pivot.to_numpy(dtype=float)
    image = axes[2].imshow(
        matrix,
        aspect="auto",
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        interpolation="nearest",
    )
    axes[2].set_title("Held-out component-score PCC", loc="left", color=INK)
    axes[2].set_xlabel("Response component")
    axes[2].set_ylabel("Outer fold")
    axes[2].set_xticks(np.arange(len(component_pivot.columns)), component_pivot.columns)
    axes[2].set_yticks(np.arange(len(component_pivot.index)), component_pivot.index)
    colorbar = fig.colorbar(image, ax=axes[2], shrink=0.75, pad=0.03)
    colorbar.set_label("Pearson r")
    for axis in axes[:2]:
        axis.grid(axis="y", color=GRID, linewidth=0.7, alpha=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Fold-specific CCLR response programs",
        x=0.06,
        y=1.02,
        ha="left",
        fontsize=13,
        color=INK,
    )
    fig.text(
        0.06,
        -0.015,
        "Response bases are fit independently on each outer-training fold; "
        "component numbers are not aligned across folds.",
        color="#59636E",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    _finish(fig, root / "results" / "figures" / "cclr_components.png")
