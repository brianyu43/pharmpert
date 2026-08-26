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
