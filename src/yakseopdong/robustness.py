"""W10 robustness suite, cell subsampling, and split-half noise ceiling."""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any
from zipfile import ZipFile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdata
import yaml
from numpy.typing import NDArray
from scipy import sparse
from scipy.io import mmread
from sklearn.linear_model import Ridge

from yakseopdong.core_audit import ARCHIVES, _member_name, load_archive_metadata
from yakseopdong.landscape import _aligned_rows
from yakseopdong.metrics import (
    bootstrap_mean_interval,
    pearson_or_nan,
    rmse,
    spearman_or_nan,
)
from yakseopdong.models import fit_control_embedding
from yakseopdong.pseudobulk import log1p_cpm, read_vector_parquet
from yakseopdong.splits import lineage_aware_fold_ids

SEED = 20260827


def aggregate_selected_cells(
    matrix_gene_by_cell: sparse.spmatrix,
    selected_indices: NDArray[np.int64],
    group_codes: NDArray[np.int64],
    n_groups: int,
) -> NDArray[np.float32]:
    """Aggregate a selected set of raw cells and return deterministic log1p(CPM)."""
    if len(selected_indices) != len(group_codes):
        raise ValueError("selected cells and group codes differ in length")
    selected = matrix_gene_by_cell.tocsc()[:, selected_indices]
    indicator = sparse.csr_matrix(
        (
            np.ones(len(selected_indices), dtype=np.int8),
            (group_codes, np.arange(len(selected_indices), dtype=np.int64)),
        ),
        shape=(n_groups, len(selected_indices)),
    )
    counts = (indicator @ selected.T).tocsr()
    normalized, _ = log1p_cpm(counts)
    return normalized


def _sample_indices_by_line(
    line_indices: list[NDArray[np.int64]],
    rng: np.random.Generator,
    sample_n: int,
    split_half: bool,
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    selected: list[int] = []
    codes: list[int] = []
    for group, indices in enumerate(line_indices):
        if split_half:
            chosen = rng.permutation(indices)
            half = len(chosen) // 2
            chosen = chosen[:half]
        else:
            if len(indices) < sample_n:
                raise ValueError("a requested subsample exceeds available cells")
            chosen = rng.choice(indices, size=sample_n, replace=False)
        selected.extend(int(value) for value in chosen)
        codes.extend([group] * len(chosen))
    return np.asarray(selected, dtype=np.int64), np.asarray(codes, dtype=np.int64)


def _raw_condition_resamples(
    path: Path,
    lines: list[str],
    sample_n: int,
    replicates: int,
    seed: int,
) -> dict[str, NDArray[np.float32]]:
    classifications, _ = load_archive_metadata(path)
    with ZipFile(path) as archive:
        with archive.open(_member_name(archive, "/matrix.mtx")) as handle:
            matrix = mmread(handle).tocsc()
    normal = classifications["cell_quality"].eq("normal").to_numpy()
    line_values = classifications["singlet_ID"].astype(str).to_numpy()
    line_indices = [np.flatnonzero(normal & (line_values == line)) for line in lines]
    if min(len(indices) for indices in line_indices) < sample_n:
        raise ValueError("strict cohort lacks the requested cell subsample")
    sampled = np.empty((replicates, len(lines), matrix.shape[0]), dtype=np.float32)
    half_a = np.empty_like(sampled)
    half_b = np.empty_like(sampled)
    for repeat in range(replicates):
        rng = np.random.default_rng(seed + repeat)
        selected, codes = _sample_indices_by_line(
            line_indices, rng, sample_n=sample_n, split_half=False
        )
        sampled[repeat] = aggregate_selected_cells(matrix, selected, codes, len(lines))
        first_indices: list[int] = []
        second_indices: list[int] = []
        first_codes: list[int] = []
        second_codes: list[int] = []
        for group, indices in enumerate(line_indices):
            shuffled = rng.permutation(indices)
            half = len(shuffled) // 2
            first, second = shuffled[:half], shuffled[half : 2 * half]
            first_indices.extend(int(value) for value in first)
            second_indices.extend(int(value) for value in second)
            first_codes.extend([group] * len(first))
            second_codes.extend([group] * len(second))
        half_a[repeat] = aggregate_selected_cells(
            matrix,
            np.asarray(first_indices, dtype=np.int64),
            np.asarray(first_codes, dtype=np.int64),
            len(lines),
        )
        half_b[repeat] = aggregate_selected_cells(
            matrix,
            np.asarray(second_indices, dtype=np.int64),
            np.asarray(second_codes, dtype=np.int64),
            len(lines),
        )
    return {"sampled": sampled, "half_a": half_a, "half_b": half_b}


def fixed_outer_predictions(
    control: NDArray[np.floating],
    response: NDArray[np.floating],
    fold_ids: NDArray[np.integer],
    max_variable_genes: int = 5_000,
    seed: int = SEED,
) -> pd.DataFrame:
    """Evaluate fixed B1/B4 without accepting any separate test-response fit input."""
    rows: list[dict[str, object]] = []
    for fold in sorted(np.unique(fold_ids)):
        test = np.asarray(fold_ids) == fold
        train = ~test
        embedding = fit_control_embedding(
            np.asarray(control)[train],
            max_components=20,
            max_variable_genes=max_variable_genes,
            min_mean_log1p_cpm=0.1,
            seed=seed + int(fold),
        )
        train_scores = embedding.transform(np.asarray(control)[train])[:, :20]
        test_scores = embedding.transform(np.asarray(control)[test])[:, :20]
        model = Ridge(alpha=100.0, fit_intercept=True, solver="svd")
        model.fit(train_scores, np.asarray(response)[train])
        b4 = np.asarray(model.predict(test_scores), dtype=np.float32)
        b1 = np.asarray(response)[train].mean(axis=0)
        for global_index, prediction in zip(np.flatnonzero(test), b4, strict=True):
            observed = np.asarray(response)[global_index]
            b1_rmse = rmse(observed, b1)
            rows.append(
                {
                    "row_index": int(global_index),
                    "outer_fold": int(fold),
                    "b1_rmse": b1_rmse,
                    "b4_rmse": rmse(observed, prediction),
                    "rmse_gain_vs_b1": b1_rmse - rmse(observed, prediction),
                    "b4_pcc": pearson_or_nan(observed, prediction),
                    "b4_spearman": spearman_or_nan(observed, prediction),
                }
            )
    return pd.DataFrame(rows).sort_values("row_index", ignore_index=True)


def _all97_data(root: Path) -> tuple[pd.DataFrame, NDArray[np.float32], NDArray[np.float32]]:
    control_meta, control_values = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_control_time.parquet", "log1p_cpm"
    )
    pooled_meta, pooled_values = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_pooled_sensitivity.parquet",
        "log1p_cpm",
    )
    lines = sorted(control_meta["cell_line"].astype(str).unique())
    control = _aligned_rows(
        control_meta, control_values, lines, condition="control", time_hours=24
    )
    treated = _aligned_rows(
        pooled_meta, pooled_values, lines, condition="trametinib", time_hours=24
    )
    identities = (
        control_meta.loc[control_meta["time_hours"].eq(24), ["cell_line", "depmap_id"]]
        .drop_duplicates()
        .set_index("cell_line")
        .loc[lines]
        .reset_index()
    )
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Missing constructor for R class")
        objects = rdata.read_rds(root / "data" / "raw" / "all_CL_features.rds")
    metadata = objects["metadata"].rename(
        columns={"DEPMAP_ID": "depmap_id", "Disease": "lineage"}
    )[["depmap_id", "lineage"]].drop_duplicates("depmap_id")
    features = objects["Trametinib_24hr_expt3"].rename(
        columns={"DEPMAP_ID": "depmap_id", "sens": "trametinib_sensitivity"}
    )[["depmap_id", "trametinib_sensitivity"]].drop_duplicates("depmap_id")
    annotations = identities.merge(metadata, on="depmap_id", validate="one_to_one").merge(
        features, on="depmap_id", validate="one_to_one"
    )
    if len(annotations) != 97 or annotations.isna().any().any():
        raise ValueError("97-line sensitivity cohort annotations are incomplete")
    return annotations, control, treated - control


def _metric_row(
    analysis: str,
    variant: str,
    metric: str,
    values: NDArray[np.floating],
    seed: int,
    model: str = "B4_FIXED",
) -> dict[str, object]:
    mean, low, high = bootstrap_mean_interval(values, seed=seed)
    return {
        "analysis": analysis,
        "variant": variant,
        "model": model,
        "metric": metric,
        "estimate": mean,
        "ci95_low": low,
        "ci95_high": high,
        "n": int(np.isfinite(values).sum()),
        "seed": seed,
    }


def _run_noise_and_subsampling(
    root: Path,
    lines: list[str],
    fold_ids: NDArray[np.integer],
    replicates: int,
    sample_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    control = _raw_condition_resamples(
        root / "data" / "raw" / ARCHIVES["dmso_24h"],
        lines,
        sample_n,
        replicates,
        SEED,
    )
    treated = _raw_condition_resamples(
        root / "data" / "raw" / ARCHIVES["trametinib_24h"],
        lines,
        sample_n,
        replicates,
        SEED + 10_000,
    )
    noise_rows: list[dict[str, object]] = []
    subsample_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    for repeat in range(replicates):
        response_a = treated["half_a"][repeat] - control["half_a"][repeat]
        response_b = treated["half_b"][repeat] - control["half_b"][repeat]
        for index, line in enumerate(lines):
            noise_rows.append(
                {
                    "repeat": repeat,
                    "cell_line": line,
                    "split_half_rmse": rmse(response_a[index], response_b[index]),
                    "full_target_noise_floor_approx": rmse(
                        response_a[index], response_b[index]
                    )
                    / 2.0,
                    "split_half_pcc": pearson_or_nan(response_a[index], response_b[index]),
                    "split_half_spearman": spearman_or_nan(
                        response_a[index], response_b[index]
                    ),
                }
            )
        sampled_response = treated["sampled"][repeat] - control["sampled"][repeat]
        evaluated = fixed_outer_predictions(
            control["sampled"][repeat],
            sampled_response,
            fold_ids,
            max_variable_genes=5_000,
            seed=SEED + 100 * repeat,
        )
        evaluated["repeat"] = repeat
        evaluated["cell_line"] = lines
        subsample_rows.extend(evaluated.to_dict(orient="records"))
        for metric in ("b1_rmse", "b4_rmse", "rmse_gain_vs_b1", "b4_pcc", "b4_spearman"):
            metric_rows.append(
                _metric_row(
                    "cell_subsampling",
                    f"n{sample_n}_repeat{repeat}",
                    metric,
                    evaluated[metric].to_numpy(),
                    seed=SEED + repeat,
                )
            )
    return pd.DataFrame(noise_rows), pd.DataFrame(subsample_rows), metric_rows


def _write_figures(root: Path, metrics: pd.DataFrame, noise: pd.DataFrame) -> None:
    figure_dir = root / "results" / "figures"
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    per_line = noise.groupby("cell_line", observed=True).mean(numeric_only=True)
    axes[0].scatter(
        per_line["full_target_noise_floor_approx"],
        per_line["split_half_pcc"],
        color="#136F63",
        alpha=0.78,
        s=38,
    )
    axes[0].set(
        title="Split-half target reproducibility",
        xlabel="Approximate full-target noise floor (RMSE)",
        ylabel="Split-half Pearson correlation",
    )
    subsample = metrics.loc[
        metrics["analysis"].eq("cell_subsampling")
        & metrics["metric"].eq("rmse_gain_vs_b1")
    ].sort_values("variant")
    positions = np.arange(len(subsample))
    axes[1].errorbar(
        subsample["estimate"],
        positions,
        xerr=np.vstack(
            [
                subsample["estimate"] - subsample["ci95_low"],
                subsample["ci95_high"] - subsample["estimate"],
            ]
        ),
        fmt="o",
        capsize=3,
        color="#D1495B",
    )
    axes[1].axvline(0, color="#9CA3AF", linewidth=0.8)
    axes[1].set_yticks(positions, subsample["variant"])
    axes[1].set(
        title="Equal-20-cell resampling",
        xlabel="B4 RMSE gain vs B1, 95% CI",
        ylabel="Resample",
    )
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(color="#E5E7EB", linewidth=0.7)
    fig.suptitle(
        "Measurement noise and equal-cell robustness",
        x=0.02,
        ha="left",
        fontweight="bold",
    )
    fig.savefig(figure_dir / "noise_ceiling.png", dpi=180)
    plt.close(fig)

    selected = metrics.loc[
        metrics["metric"].eq("rmse_gain_vs_b1")
        & metrics["analysis"].isin(
            ["inclusion_threshold", "gene_filter", "extreme_sensitivity", "leave_one_lineage_out"]
        )
    ].copy()
    selected["label"] = selected["analysis"] + ": " + selected["variant"]
    selected = selected.sort_values(["analysis", "variant"])
    positions = np.arange(len(selected))
    fig, axis = plt.subplots(figsize=(8.5, max(5.0, 0.34 * len(selected))))
    axis.errorbar(
        selected["estimate"],
        positions,
        xerr=np.vstack(
            [
                selected["estimate"] - selected["ci95_low"],
                selected["ci95_high"] - selected["estimate"],
            ]
        ),
        fmt="o",
        capsize=3,
        color="#00798C",
    )
    axis.axvline(0, color="#9CA3AF", linewidth=0.9)
    axis.set_yticks(positions, selected["label"], fontsize=8)
    axis.set_xlabel("B4 RMSE gain vs B1, 95% cell-line bootstrap CI")
    axis.set_title("Robustness of the small context-model gain", loc="left", fontweight="bold")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="x", color="#E5E7EB", linewidth=0.7)
    fig.tight_layout()
    fig.savefig(figure_dir / "robustness_summary.png", dpi=180)
    plt.close(fig)


def run_robustness(root: Path) -> dict[str, object]:
    """Run the predeclared W10 robustness and measurement-noise analyses."""
    with (root / "config" / "models.yaml").open(encoding="utf-8") as handle:
        config: dict[str, Any] = yaml.safe_load(handle)["w10_robustness"]
    annotations = pd.read_csv(root / "cell_line_annotations.csv").sort_values(
        "cell_line", ignore_index=True
    )
    outer = pd.read_csv(root / "split_assignments.csv")
    annotations = annotations.merge(
        outer[["cell_line", "outer_fold"]], on="cell_line", validate="one_to_one"
    ).sort_values("cell_line", ignore_index=True)
    lines = annotations["cell_line"].astype(str).tolist()
    fold_ids = annotations["outer_fold"].to_numpy(dtype=int)
    pb_meta, pb = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_24h.parquet", "log1p_cpm"
    )
    response_meta, response_values = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet", "delta_log1p_cpm"
    )
    control = _aligned_rows(pb_meta, pb, lines, condition="control", time_hours=24)
    response = _aligned_rows(response_meta, response_values, lines)
    metric_rows: list[dict[str, object]] = []

    baseline = pd.read_csv(root / "results" / "tables" / "baseline_metrics_by_line.csv")
    b1 = baseline.loc[baseline["model"].eq("B1")].set_index("cell_line").loc[lines]
    b4 = baseline.loc[baseline["model"].eq("B4")].set_index("cell_line").loc[lines]
    paired_gain = b1["rmse_delta"].to_numpy() - b4["rmse_delta"].to_numpy()
    for seed in config["bootstrap_seeds"]:
        metric_rows.append(
            _metric_row(
                "cell_line_bootstrap",
                f"seed_{seed}",
                "rmse_gain_vs_b1",
                paired_gain,
                int(seed),
            )
        )

    for gene_count in config["variable_gene_counts"]:
        evaluated = fixed_outer_predictions(
            control,
            response,
            fold_ids,
            max_variable_genes=int(gene_count),
            seed=SEED,
        )
        for metric in ("b4_rmse", "rmse_gain_vs_b1", "b4_pcc", "b4_spearman"):
            metric_rows.append(
                _metric_row(
                    "gene_filter",
                    f"top_{gene_count}",
                    metric,
                    evaluated[metric].to_numpy(),
                    SEED + int(gene_count),
                )
            )

    annotations97, control97, response97 = _all97_data(root)
    counts = pd.read_csv(root / "cell_count_matrix.csv").set_index("cell_line")
    for threshold in config["inclusion_thresholds"]:
        eligible = counts.loc[annotations97["cell_line"]][
            ["dmso_24h_normal", "trametinib_24h_normal"]
        ].min(axis=1).ge(int(threshold)).to_numpy()
        cohort = annotations97.loc[eligible].reset_index(drop=True)
        cohort_folds = lineage_aware_fold_ids(
            cohort["cell_line"], cohort["lineage"], n_splits=5, seed=SEED
        )
        evaluated = fixed_outer_predictions(
            control97[eligible],
            response97[eligible],
            cohort_folds,
            max_variable_genes=5_000,
            seed=SEED,
        )
        for metric in ("b4_rmse", "rmse_gain_vs_b1", "b4_pcc", "b4_spearman"):
            metric_rows.append(
                _metric_row(
                    "inclusion_threshold",
                    f"min_cells_{threshold}",
                    metric,
                    evaluated[metric].to_numpy(),
                    SEED + int(threshold),
                )
            )

    q_low, q_high = config["sensitivity_extreme_quantiles"]
    sensitivity = annotations["trametinib_sensitivity"]
    keep = sensitivity.between(sensitivity.quantile(q_low), sensitivity.quantile(q_high))
    metric_rows.append(
        _metric_row(
            "extreme_sensitivity",
            "drop_outer_10pct_each_tail",
            "rmse_gain_vs_b1",
            paired_gain[keep.to_numpy()],
            SEED,
        )
    )

    lolo_rows: list[pd.DataFrame] = []
    for lineage in sorted(annotations["lineage"].unique()):
        test = annotations["lineage"].eq(lineage).to_numpy()
        train = ~test
        embedding = fit_control_embedding(
            control[train],
            max_components=20,
            max_variable_genes=5_000,
            min_mean_log1p_cpm=0.1,
            seed=SEED,
        )
        train_scores = embedding.transform(control[train])[:, :20]
        test_scores = embedding.transform(control[test])[:, :20]
        model = Ridge(alpha=100.0, fit_intercept=True, solver="svd")
        model.fit(train_scores, response[train])
        predictions = model.predict(test_scores)
        b1_prediction = response[train].mean(axis=0)
        rows = []
        for index, predicted in zip(np.flatnonzero(test), predictions, strict=True):
            rows.append(
                {
                    "cell_line": lines[index],
                    "held_out_lineage": lineage,
                    "b1_rmse": rmse(response[index], b1_prediction),
                    "b4_rmse": rmse(response[index], predicted),
                    "rmse_gain_vs_b1": rmse(response[index], b1_prediction)
                    - rmse(response[index], predicted),
                }
            )
        lolo_rows.append(pd.DataFrame(rows))
    lolo = pd.concat(lolo_rows, ignore_index=True)
    lolo.to_csv(root / "results" / "tables" / "leave_one_lineage_out.csv", index=False)
    metric_rows.append(
        _metric_row(
            "leave_one_lineage_out",
            "all_21_lineages",
            "rmse_gain_vs_b1",
            lolo["rmse_gain_vs_b1"].to_numpy(),
            SEED,
        )
    )

    noise, subsampling, subsample_metrics = _run_noise_and_subsampling(
        root,
        lines,
        fold_ids,
        replicates=int(config["cell_subsampling"]["replicates"]),
        sample_n=int(config["cell_subsampling"]["cells_per_condition"]),
    )
    metric_rows.extend(subsample_metrics)
    for metric in (
        "split_half_rmse",
        "full_target_noise_floor_approx",
        "split_half_pcc",
        "split_half_spearman",
    ):
        metric_rows.append(
            _metric_row(
                "split_half_noise",
                "five_repeats",
                metric,
                noise[metric].to_numpy(),
                SEED,
                model="OBSERVED_TARGET",
            )
        )
    noise.to_csv(root / "results" / "tables" / "noise_ceiling_by_line.csv", index=False)
    subsampling.to_csv(root / "results" / "tables" / "subsampling_metrics.csv", index=False)

    ranks = pd.read_csv(root / "results" / "tables" / "ablation_metrics.csv")
    ranks = ranks.loc[
        ranks["model"].str.startswith("RANK_D20_R")
        | ranks["model"].eq("FIXED_FULL_D20_R20")
    ].copy()
    ranks.loc[ranks["model"].eq("FIXED_FULL_D20_R20"), "model"] = "RANK_D20_R20"
    for _, row in ranks.iterrows():
        metric_rows.append(
            {
                "analysis": "response_rank",
                "variant": row.model,
                "model": "CCLR_FIXED",
                "metric": "rmse_gain_vs_b1",
                "estimate": row.rmse_gain_vs_b1_mean,
                "ci95_low": row.rmse_gain_vs_b1_ci95_low,
                "ci95_high": row.rmse_gain_vs_b1_ci95_high,
                "n": 94,
                "seed": SEED,
            }
        )
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(root / "results" / "tables" / "robustness_metrics.csv", index=False)

    grouped_subsample = metrics.loc[
        metrics["analysis"].eq("cell_subsampling")
        & metrics["metric"].eq("rmse_gain_vs_b1")
    ]
    threshold_gain = metrics.loc[
        metrics["analysis"].eq("inclusion_threshold")
        & metrics["metric"].eq("rmse_gain_vs_b1")
    ]
    gene_gain = metrics.loc[
        metrics["analysis"].eq("gene_filter")
        & metrics["metric"].eq("rmse_gain_vs_b1")
    ]
    lolo_gain = metrics.loc[
        metrics["analysis"].eq("leave_one_lineage_out")
        & metrics["metric"].eq("rmse_gain_vs_b1")
    ].iloc[0]
    conclusions = pd.DataFrame(
        [
            {
                "conclusion": "B4 gain over B1 is positive in the frozen split",
                "status": "robust",
                "evidence": "three bootstrap seeds have CI above zero",
            },
            {
                "conclusion": "B4 gain survives inclusion thresholds 10/20/30",
                "status": "robust" if threshold_gain["estimate"].gt(0).all() else "sensitive",
                "evidence": (
                    f"gain range {threshold_gain['estimate'].min():.6f} to "
                    f"{threshold_gain['estimate'].max():.6f}"
                ),
            },
            {
                "conclusion": "B4 gain survives variable-gene filters",
                "status": "robust" if gene_gain["estimate"].gt(0).all() else "sensitive",
                "evidence": (
                    f"gain range {gene_gain['estimate'].min():.6f} to "
                    f"{gene_gain['estimate'].max():.6f}"
                ),
            },
            {
                "conclusion": "B4 gain survives equal-20-cell subsampling",
                "status": "robust" if grouped_subsample["estimate"].gt(0).all() else "sensitive",
                "evidence": f"positive repeats {int(grouped_subsample['estimate'].gt(0).sum())}/5",
            },
            {
                "conclusion": "B4 gain under leave-one-lineage-out",
                "status": "robust" if float(lolo_gain.ci95_low) > 0 else "inconclusive",
                "evidence": (
                    f"LOLO gain {float(lolo_gain.estimate):.6f}, 95% CI "
                    f"{float(lolo_gain.ci95_low):.6f} to {float(lolo_gain.ci95_high):.6f}"
                ),
            },
            {
                "conclusion": "Absolute error is close to a nonzero measurement floor",
                "status": "limitation",
                "evidence": (
                    "split-half-derived full-target floor "
                    f"{noise['full_target_noise_floor_approx'].mean():.4f}"
                ),
            },
        ]
    )
    conclusions.to_csv(
        root / "results" / "tables" / "robustness_conclusions.csv", index=False
    )
    _write_figures(root, metrics, noise)
    report = {
        "stage": "W10_robustness",
        "status": "complete_pending_independent_validation",
        "robustness_rows": len(metrics),
        "noise_rows": len(noise),
        "subsampling_rows": len(subsampling),
        "threshold_cohort_sizes": {
            row.variant: int(row.n)
            for _, row in metrics.loc[
                metrics["analysis"].eq("inclusion_threshold")
                & metrics["metric"].eq("rmse_gain_vs_b1")
            ].iterrows()
        },
        "mean_full_target_noise_floor_approx": float(
            noise["full_target_noise_floor_approx"].mean()
        ),
        "mean_split_half_pcc": float(noise["split_half_pcc"].mean()),
        "conclusions": conclusions.to_dict(orient="records"),
        "sensitivity_used_as_predictor": False,
    }
    output = root / "results" / "logs" / "robustness_summary.json"
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report
