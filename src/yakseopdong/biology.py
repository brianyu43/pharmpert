"""W9 biological interpretation and held-out prediction error diagnostics."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from numpy.typing import NDArray
from scipy.stats import hypergeom, rankdata
from sklearn.decomposition import PCA

from yakseopdong.landscape import _aligned_rows
from yakseopdong.metrics import bootstrap_mean_interval, spearman_or_nan
from yakseopdong.models import fit_control_embedding
from yakseopdong.pathways import benjamini_hochberg, load_geneset_config
from yakseopdong.pseudobulk import read_vector_parquet

SEED = 20260827
TARGET_GENES = (
    "EGR1",
    "ETV4",
    "ETV5",
    "DUSP4",
    "DUSP5",
    "DUSP6",
    "SPRY2",
    "SPRY4",
    "MCM2",
    "MCM3",
    "MCM4",
    "MCM5",
    "MCM6",
    "MCM7",
    "MCM10",
)


def spearman_bootstrap_permutation(
    x: NDArray[np.floating],
    y: NDArray[np.floating],
    seed: int,
    n_bootstrap: int = 2_000,
    n_permutation: int = 10_000,
) -> tuple[float, float, float, float]:
    """Spearman effect, paired bootstrap interval, and two-sided permutation p."""
    left = np.asarray(x, dtype=float)
    right = np.asarray(y, dtype=float)
    finite = np.isfinite(left) & np.isfinite(right)
    left, right = left[finite], right[finite]
    if len(left) < 5:
        return float("nan"), float("nan"), float("nan"), float("nan")
    observed = spearman_or_nan(left, right)
    rng = np.random.default_rng(seed)
    bootstrap = np.empty(n_bootstrap, dtype=float)
    for index in range(n_bootstrap):
        sample = rng.integers(0, len(left), size=len(left))
        bootstrap[index] = spearman_or_nan(left[sample], right[sample])
    low, high = np.nanquantile(bootstrap, [0.025, 0.975])
    x_rank = rankdata(left).astype(float)
    y_rank = rankdata(right).astype(float)
    x_rank -= x_rank.mean()
    y_rank -= y_rank.mean()
    denominator = np.sqrt(np.square(x_rank).sum() * np.square(y_rank).sum())
    permuted = np.vstack([rng.permutation(y_rank) for _ in range(n_permutation)])
    null = permuted @ x_rank / denominator
    p_value = (1 + int((np.abs(null) >= abs(observed)).sum())) / (n_permutation + 1)
    return float(observed), float(low), float(high), float(p_value)


def partial_spearman_bootstrap_permutation(
    x: NDArray[np.floating],
    y: NDArray[np.floating],
    covariate: NDArray[np.floating],
    seed: int,
    n_bootstrap: int = 2_000,
    n_permutation: int = 10_000,
) -> tuple[float, float, float, float]:
    """Partial rank correlation with bootstrap CI and residual permutation p."""
    matrix = np.column_stack([x, y, covariate]).astype(float)
    matrix = matrix[np.isfinite(matrix).all(axis=1)]

    def effect(current: NDArray[np.float64]) -> float:
        ranked = np.column_stack([rankdata(current[:, index]) for index in range(3)])
        design = np.column_stack([np.ones(len(ranked)), ranked[:, 2]])
        left = ranked[:, 0] - design @ np.linalg.lstsq(
            design, ranked[:, 0], rcond=None
        )[0]
        right = ranked[:, 1] - design @ np.linalg.lstsq(
            design, ranked[:, 1], rcond=None
        )[0]
        return float(np.corrcoef(left, right)[0, 1])

    observed = effect(matrix)
    rng = np.random.default_rng(seed)
    bootstrap = np.asarray(
        [effect(matrix[rng.integers(0, len(matrix), size=len(matrix))]) for _ in range(n_bootstrap)]
    )
    low, high = np.nanquantile(bootstrap, [0.025, 0.975])
    ranked = np.column_stack([rankdata(matrix[:, index]) for index in range(3)])
    design = np.column_stack([np.ones(len(ranked)), ranked[:, 2]])
    left = ranked[:, 0] - design @ np.linalg.lstsq(design, ranked[:, 0], rcond=None)[0]
    right = ranked[:, 1] - design @ np.linalg.lstsq(design, ranked[:, 1], rcond=None)[0]
    denominator = np.sqrt(np.square(left).sum() * np.square(right).sum())
    null = np.asarray(
        [float(left @ rng.permutation(right) / denominator) for _ in range(n_permutation)]
    )
    p_value = (1 + int((np.abs(null) >= abs(observed)).sum())) / (n_permutation + 1)
    return observed, float(low), float(high), float(p_value)


def cliffs_delta_bootstrap_permutation(
    values: NDArray[np.floating],
    binary: NDArray[np.bool_],
    seed: int,
    n_bootstrap: int = 2_000,
    n_permutation: int = 10_000,
) -> tuple[float, float, float, float]:
    """Cliff's delta with stratified bootstrap CI and label-permutation p-value."""
    array = np.asarray(values, dtype=float)
    labels = np.asarray(binary, dtype=bool)
    positive, negative = array[labels], array[~labels]
    if min(len(positive), len(negative)) < 2:
        return float("nan"), float("nan"), float("nan"), float("nan")

    def delta(left: NDArray[np.floating], right: NDArray[np.floating]) -> float:
        return float(np.sign(left[:, None] - right[None, :]).mean())

    observed = delta(positive, negative)
    rng = np.random.default_rng(seed)
    bootstrap = np.asarray(
        [
            delta(
                rng.choice(positive, size=len(positive), replace=True),
                rng.choice(negative, size=len(negative), replace=True),
            )
            for _ in range(n_bootstrap)
        ]
    )
    low, high = np.quantile(bootstrap, [0.025, 0.975])
    observed_mean_difference = abs(float(positive.mean() - negative.mean()))
    null = np.empty(n_permutation, dtype=float)
    for index in range(n_permutation):
        permuted = rng.permutation(labels)
        null[index] = abs(float(array[permuted].mean() - array[~permuted].mean()))
    p_value = (1 + int((null >= observed_mean_difference).sum())) / (n_permutation + 1)
    return observed, float(low), float(high), float(p_value)


def eta_squared_permutation(
    values: NDArray[np.floating],
    groups: NDArray[np.str_],
    seed: int,
    n_permutation: int = 10_000,
) -> tuple[float, float]:
    """Global lineage eta-squared and a group-label permutation p-value."""
    array = np.asarray(values, dtype=float)
    labels = np.asarray(groups, dtype=str)

    def eta(current: NDArray[np.str_]) -> float:
        grand = array.mean()
        between = sum(
            int(mask.sum()) * float(np.square(array[mask].mean() - grand))
            for label in np.unique(current)
            if (mask := current == label).any()
        )
        total = float(np.square(array - grand).sum())
        return between / total if total else float("nan")

    observed = eta(labels)
    rng = np.random.default_rng(seed)
    null = np.asarray([eta(rng.permutation(labels)) for _ in range(n_permutation)])
    p_value = (1 + int((null >= observed).sum())) / (n_permutation + 1)
    return float(observed), float(p_value)


def _response_pca(
    response: NDArray[np.floating], genes: pd.DataFrame
) -> tuple[NDArray[np.float64], pd.DataFrame, NDArray[np.int64], PCA]:
    variances = np.asarray(response, dtype=float).var(axis=0)
    eligible = np.flatnonzero(variances > 1e-10)
    selected = eligible[np.argsort(-variances[eligible], kind="stable")[:5_000]]
    pca = PCA(n_components=5, svd_solver="randomized", random_state=SEED)
    scores = pca.fit_transform(np.asarray(response)[:, selected])
    rows: list[dict[str, object]] = []
    for component, loading in enumerate(pca.components_, start=1):
        order = np.argsort(loading, kind="stable")
        for direction, positions in (("negative", order[:50]), ("positive", order[-50:])):
            for position in positions:
                gene_index = int(selected[position])
                rows.append(
                    {
                        "outer_fold": 0,
                        "component": component,
                        "direction": direction,
                        "gene_id": genes.loc[gene_index, "gene_id"],
                        "gene_symbol": genes.loc[gene_index, "gene_symbol"],
                        "loading": float(loading[position]),
                    }
                )
    return scores, pd.DataFrame(rows), selected, pca


def _pathway_enrichment(
    loadings: pd.DataFrame, genes: pd.DataFrame, genesets: dict[str, Any]
) -> pd.DataFrame:
    universe = set(genes["gene_symbol"].astype(str))
    rows: list[dict[str, object]] = []
    for (component, direction), block in loadings.groupby(
        ["component", "direction"], observed=True
    ):
        top = set(block["gene_symbol"].astype(str)) & universe
        for name, collection in genesets["collections"].items():
            pathway = set(collection["genes"]) & universe
            overlap = top & pathway
            rows.append(
                {
                    "component": int(component),
                    "direction": direction,
                    "pathway": name,
                    "top_gene_count": len(top),
                    "pathway_gene_count": len(pathway),
                    "overlap_count": len(overlap),
                    "overlap_genes": ";".join(sorted(overlap)),
                    "p_value": float(
                        hypergeom.sf(
                            len(overlap) - 1,
                            len(universe),
                            len(pathway),
                            len(top),
                        )
                    ),
                }
            )
    result = pd.DataFrame(rows)
    result["fdr_bh"] = benjamini_hochberg(result["p_value"].to_numpy())
    return result


def _baseline_novelty(
    control: NDArray[np.floating], annotations: pd.DataFrame
) -> NDArray[np.float64]:
    novelty = np.zeros(len(annotations), dtype=float)
    for outer_fold in sorted(annotations["outer_fold"].unique()):
        test = annotations["outer_fold"].eq(outer_fold).to_numpy()
        train = ~test
        embedding = fit_control_embedding(
            control[train],
            max_components=20,
            max_variable_genes=5_000,
            min_mean_log1p_cpm=0.1,
            seed=SEED + 1_000 * int(outer_fold),
        )
        train_scores = embedding.transform(control[train])[:, :20]
        test_scores = embedding.transform(control[test])[:, :20]
        distances = np.sqrt(
            np.square(test_scores[:, None, :] - train_scores[None, :, :]).sum(axis=2)
        )
        novelty[test] = distances.min(axis=1)
    return novelty


def _association_row(
    family: str,
    outcome: str,
    predictor: str,
    effect: float,
    low: float,
    high: float,
    p_value: float,
    n: int,
    method: str,
) -> dict[str, object]:
    return {
        "family": family,
        "outcome": outcome,
        "predictor": predictor,
        "effect": effect,
        "ci95_low": low,
        "ci95_high": high,
        "p_value": p_value,
        "n_cell_lines": n,
        "method": method,
        "interpretation_only": True,
    }


def _write_figures(root: Path, analysis: pd.DataFrame, scores: pd.DataFrame) -> None:
    figure_dir = root / "results" / "figures"
    colors = {False: "#136F63", True: "#D1495B"}
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    for status in (False, True):
        block = scores.loc[scores["braf_or_kras_mut"].eq(status)]
        axes[0].scatter(
            block["trametinib_sensitivity"],
            block["response_pc1"],
            s=38,
            alpha=0.78,
            color=colors[status],
            label="BRAF/KRAS mutant" if status else "BRAF/KRAS wild type",
        )
    axes[0].set(
        title="Observed response PC1 and sensitivity",
        xlabel="Trametinib sensitivity",
        ylabel="Response PC1 (sign arbitrary)",
    )
    axes[0].legend(frameon=False, fontsize=8)
    pathway = analysis.loc[analysis["family"].eq("pathway_sensitivity")].sort_values(
        "effect"
    )
    positions = np.arange(len(pathway))
    axes[1].errorbar(
        pathway["effect"],
        positions,
        xerr=np.vstack(
            [pathway["effect"] - pathway["ci95_low"], pathway["ci95_high"] - pathway["effect"]]
        ),
        fmt="o",
        capsize=3,
        color="#00798C",
    )
    axes[1].axvline(0, color="#9CA3AF", linewidth=0.8)
    axes[1].set_yticks(positions, pathway["outcome"].str.replace("_", " "), fontsize=8)
    axes[1].set(title="24h pathway score vs sensitivity", xlabel="Spearman ρ, 95% CI")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(color="#E5E7EB", linewidth=0.7)
    fig.suptitle("Biological correlates of the 24h response", x=0.02, ha="left", fontweight="bold")
    fig.savefig(figure_dir / "biological_correlates.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    axes[0].scatter(
        scores["response_rms"], scores["b4_rmse"], color="#136F63", alpha=0.78, s=38
    )
    axes[0].set(
        title="Error increases with target magnitude",
        xlabel="Observed response RMS",
        ylabel="B4 held-out RMSE",
    )
    axes[1].scatter(
        scores["baseline_novelty"], scores["b4_rmse"], color="#D1495B", alpha=0.78, s=38
    )
    axes[1].set(
        title="Error and baseline-state novelty",
        xlabel="Nearest training distance in fold-fitted PC20",
        ylabel="B4 held-out RMSE",
    )
    worst = scores.nlargest(3, "b4_rmse")
    for axis, x_column in zip(axes, ("response_rms", "baseline_novelty"), strict=True):
        for _, row in worst.iterrows():
            axis.annotate(
                str(row.cell_line).split("_")[0],
                (row[x_column], row.b4_rmse),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
            )
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(color="#E5E7EB", linewidth=0.7)
    fig.suptitle("Held-out error diagnostics", x=0.02, ha="left", fontweight="bold")
    fig.savefig(figure_dir / "prediction_error_drivers.png", dpi=180)
    plt.close(fig)


def run_biology(root: Path) -> dict[str, object]:
    """Execute W9 associations, enrichment, cases, error taxonomy, and figures."""
    annotations = pd.read_csv(root / "cell_line_annotations.csv").sort_values(
        "cell_line", ignore_index=True
    )
    outer = pd.read_csv(root / "split_assignments.csv")
    annotations = annotations.merge(
        outer[["cell_line", "outer_fold"]], on="cell_line", validate="one_to_one"
    ).sort_values("cell_line", ignore_index=True)
    lines = annotations["cell_line"].astype(str).tolist()
    pb_meta, pb = read_vector_parquet(
        root / "data" / "processed" / "pseudobulk_24h.parquet", "log1p_cpm"
    )
    response_meta, response_values = read_vector_parquet(
        root / "data" / "processed" / "response_24h.parquet", "delta_log1p_cpm"
    )
    control = _aligned_rows(pb_meta, pb, lines, condition="control", time_hours=24)
    response = _aligned_rows(response_meta, response_values, lines)
    genes = pq.read_table(root / "data" / "processed" / "gene_metadata.parquet").to_pandas()
    scores, loadings, _, pca = _response_pca(response, genes)

    b4 = (
        pd.read_csv(root / "results" / "tables" / "baseline_metrics_by_line.csv")
        .loc[lambda x: x["model"].eq("B4")]
        .set_index("cell_line")
        .loc[lines]
    )
    cclr = (
        pd.read_csv(root / "results" / "tables" / "cclr_metrics_by_line.csv")
        .set_index("cell_line")
        .loc[lines]
    )
    counts = pd.read_csv(root / "cell_count_matrix.csv").set_index("cell_line").loc[lines]
    features = annotations.copy()
    for component in range(5):
        features[f"response_pc{component + 1}"] = scores[:, component]
    features["braf_or_kras_mut"] = features["braf_mut"] | features["kras_mut"]
    features["response_rms"] = np.sqrt(np.mean(np.square(response), axis=1))
    features["min_condition_cells"] = counts[
        ["dmso_24h_normal", "trametinib_24h_normal"]
    ].min(axis=1).to_numpy()
    features["baseline_novelty"] = _baseline_novelty(control, annotations)
    features["b4_rmse"] = b4["rmse_delta"].to_numpy()
    features["cclr_rmse"] = cclr["rmse_delta"].to_numpy()
    features["b4_pcc_context"] = b4["pcc_context"].to_numpy()

    association_rows: list[dict[str, object]] = []
    test_index = 0
    for component in range(1, 6):
        outcome = f"response_pc{component}"
        effect, low, high, p_value = spearman_bootstrap_permutation(
            features["trametinib_sensitivity"].to_numpy(),
            features[outcome].to_numpy(),
            seed=SEED + test_index,
        )
        association_rows.append(
            _association_row(
                "component_sensitivity",
                outcome,
                "trametinib_sensitivity",
                effect,
                low,
                high,
                p_value,
                len(features),
                "spearman_bootstrap_permutation",
            )
        )
        test_index += 1
        for mutation in ("braf_mut", "kras_mut"):
            effect, low, high, p_value = cliffs_delta_bootstrap_permutation(
                features[outcome].to_numpy(),
                features[mutation].to_numpy(dtype=bool),
                seed=SEED + test_index,
            )
            association_rows.append(
                _association_row(
                    "component_mutation",
                    outcome,
                    mutation,
                    effect,
                    low,
                    high,
                    p_value,
                    len(features),
                    "cliffs_delta_mutant_minus_wildtype",
                )
            )
            test_index += 1
        eta, p_value = eta_squared_permutation(
            features[outcome].to_numpy(),
            features["lineage"].to_numpy(dtype=str),
            seed=SEED + test_index,
        )
        association_rows.append(
            _association_row(
                "component_lineage",
                outcome,
                "lineage_global",
                eta,
                float("nan"),
                float("nan"),
                p_value,
                len(features),
                "eta_squared_label_permutation",
            )
        )
        test_index += 1

    genesets = load_geneset_config(root / "config" / "genesets.yaml")
    symbols = genes["gene_symbol"].astype(str)
    for pathway, collection in genesets["collections"].items():
        indices = np.flatnonzero(symbols.isin(set(collection["genes"])).to_numpy())
        pathway_score = response[:, indices].mean(axis=1)
        features[f"pathway__{pathway}"] = pathway_score
        effect, low, high, p_value = spearman_bootstrap_permutation(
            features["trametinib_sensitivity"].to_numpy(),
            pathway_score,
            seed=SEED + test_index,
        )
        association_rows.append(
            _association_row(
                "pathway_sensitivity",
                pathway,
                "trametinib_sensitivity",
                effect,
                low,
                high,
                p_value,
                len(features),
                "spearman_bootstrap_permutation",
            )
        )
        test_index += 1

    for outcome in ("b4_rmse", "cclr_rmse"):
        for predictor in (
            "response_rms",
            "baseline_novelty",
            "min_condition_cells",
            "trametinib_sensitivity",
        ):
            effect, low, high, p_value = spearman_bootstrap_permutation(
                features[predictor].to_numpy(),
                features[outcome].to_numpy(),
                seed=SEED + test_index,
            )
            association_rows.append(
                _association_row(
                    "error_factor",
                    outcome,
                    predictor,
                    effect,
                    low,
                    high,
                    p_value,
                    len(features),
                    "spearman_bootstrap_permutation",
                )
            )
            test_index += 1
        for predictor, covariate in (
            ("response_rms", "min_condition_cells"),
            ("min_condition_cells", "response_rms"),
        ):
            effect, low, high, p_value = partial_spearman_bootstrap_permutation(
                features[predictor].to_numpy(),
                features[outcome].to_numpy(),
                features[covariate].to_numpy(),
                seed=SEED + test_index,
            )
            association_rows.append(
                _association_row(
                    "error_factor_partial",
                    outcome,
                    f"{predictor}|{covariate}",
                    effect,
                    low,
                    high,
                    p_value,
                    len(features),
                    "partial_spearman_bootstrap_residual_permutation",
                )
            )
            test_index += 1

    associations = pd.DataFrame(association_rows)
    finite_p = associations["p_value"].notna()
    associations["fdr_bh"] = np.nan
    associations.loc[finite_p, "fdr_bh"] = benjamini_hochberg(
        associations.loc[finite_p, "p_value"].to_numpy()
    )
    associations.to_csv(root / "results" / "tables" / "biological_validation.csv", index=False)

    enrichment = _pathway_enrichment(loadings, genes, genesets)
    enrichment.to_csv(
        root / "results" / "tables" / "biological_pathway_enrichment.csv", index=False
    )
    loadings.to_csv(root / "results" / "tables" / "biological_component_loadings.csv", index=False)

    target_rows: list[dict[str, object]] = []
    for gene_index, symbol in enumerate(TARGET_GENES):
        positions = np.flatnonzero(symbols.eq(symbol).to_numpy())
        if not len(positions):
            continue
        values = response[:, positions].mean(axis=1)
        mean, low, high = bootstrap_mean_interval(values, seed=SEED + gene_index)
        rho, rho_low, rho_high, p_value = spearman_bootstrap_permutation(
            features["trametinib_sensitivity"].to_numpy(),
            values,
            seed=SEED + 10_000 + gene_index,
        )
        target_rows.append(
            {
                "gene_symbol": symbol,
                "mapped_columns": len(positions),
                "mean_delta_log1p_cpm": mean,
                "mean_ci95_low": low,
                "mean_ci95_high": high,
                "sensitivity_spearman": rho,
                "spearman_ci95_low": rho_low,
                "spearman_ci95_high": rho_high,
                "p_value": p_value,
            }
        )
    targets = pd.DataFrame(target_rows)
    targets["fdr_bh"] = benjamini_hochberg(targets["p_value"].to_numpy())
    targets.to_csv(root / "results" / "tables" / "target_gene_validation.csv", index=False)

    ordered = features.sort_values("b4_rmse").reset_index(drop=True)
    case_indices = {
        "best": list(range(5)),
        "median": list(range(len(ordered) // 2 - 2, len(ordered) // 2 + 3)),
        "worst": list(range(len(ordered) - 5, len(ordered))),
    }
    case_rows = []
    for case_type, indices in case_indices.items():
        block = ordered.iloc[indices].copy()
        block["case_type"] = case_type
        case_rows.append(block)
    cases = pd.concat(case_rows, ignore_index=True)
    case_columns = [
        "case_type",
        "cell_line",
        "depmap_id",
        "lineage",
        "trametinib_sensitivity",
        "braf_mut",
        "kras_mut",
        "b4_rmse",
        "cclr_rmse",
        "b4_pcc_context",
        "response_rms",
        "baseline_novelty",
        "min_condition_cells",
    ]
    cases[case_columns].to_csv(
        root / "results" / "tables" / "prediction_cases.csv", index=False
    )

    high_error = features["b4_rmse"].ge(features["b4_rmse"].quantile(0.75))
    high_magnitude = features["response_rms"].ge(features["response_rms"].quantile(0.75))
    high_novelty = features["baseline_novelty"].ge(features["baseline_novelty"].quantile(0.75))
    low_cells = features["min_condition_cells"].le(features["min_condition_cells"].quantile(0.25))
    classification = features[case_columns[1:]].copy()
    classification["high_error"] = high_error
    classification["high_response_magnitude"] = high_magnitude
    classification["high_baseline_novelty"] = high_novelty
    classification["low_cell_support"] = low_cells
    classification["error_type"] = np.select(
        [
            high_error & high_magnitude,
            high_error & high_novelty,
            high_error & low_cells,
            high_error,
        ],
        [
            "high_magnitude_target",
            "baseline_novelty",
            "low_cell_support",
            "mixed_unexplained",
        ],
        default="not_high_error",
    )
    classification.to_csv(root / "results" / "tables" / "error_classification.csv", index=False)
    features.to_csv(root / "results" / "tables" / "biological_line_features.csv", index=False)

    _write_figures(root, associations, features)
    significant = associations.loc[associations["fdr_bh"].lt(0.05)]
    error_factors = associations.loc[associations["family"].eq("error_factor")]
    report = {
        "stage": "W9_biological_interpretation",
        "status": "complete_pending_independent_validation",
        "cell_lines": len(features),
        "response_pca_explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "association_tests": len(associations),
        "significant_fdr05": int(len(significant)),
        "significant_associations": significant[
            ["family", "outcome", "predictor", "effect", "fdr_bh"]
        ].to_dict(orient="records"),
        "error_factor_effects": error_factors[
            ["outcome", "predictor", "effect", "ci95_low", "ci95_high", "fdr_bh"]
        ].to_dict(orient="records"),
        "pathway_enrichment_fdr05": int(enrichment["fdr_bh"].lt(0.05).sum()),
        "target_genes": len(targets),
        "high_error_lines": int(high_error.sum()),
        "interpretation_only": True,
        "sensitivity_used_as_predictor": False,
    }
    output = root / "results" / "logs" / "biology_summary.json"
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    return report
