"""Quality-control calculations for pseudobulk matrices."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

QC_MIN_TOTAL_CONTROL_COUNTS = 100
QC_MIN_LINES_PER_TIME = 10
BOOTSTRAP_SEED = 20260827
BOOTSTRAP_DRAWS = 2_000
MANUAL_MARKERS = ("EGR1", "ETV4", "ETV5", "DUSP4", "DUSP5", "DUSP6", "SPRY2", "SPRY4")


def row_pearson(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Pearson correlation for matching rows, returning NaN for constants."""
    left_centered = left - left.mean(axis=1, keepdims=True)
    right_centered = right - right.mean(axis=1, keepdims=True)
    numerator = np.sum(left_centered * right_centered, axis=1)
    denominator = np.sqrt(
        np.sum(left_centered**2, axis=1) * np.sum(right_centered**2, axis=1)
    )
    return np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan, dtype=np.float64),
        where=denominator > 0,
    )


def control_qc_gene_mask(
    dmso6_counts: sparse.spmatrix,
    dmso24_counts: sparse.spmatrix,
) -> np.ndarray:
    """Predeclared control-only descriptive gene universe for the time QC."""
    total = np.asarray((dmso6_counts + dmso24_counts).sum(axis=0)).ravel()
    expressed_6h = np.asarray((dmso6_counts > 0).sum(axis=0)).ravel()
    expressed_24h = np.asarray((dmso24_counts > 0).sum(axis=0)).ravel()
    return (
        (total >= QC_MIN_TOTAL_CONTROL_COUNTS)
        & (expressed_6h >= QC_MIN_LINES_PER_TIME)
        & (expressed_24h >= QC_MIN_LINES_PER_TIME)
    )


def summarize_markers(
    genes: pd.DataFrame,
    response: np.ndarray,
    markers: tuple[str, ...] = MANUAL_MARKERS,
) -> pd.DataFrame:
    """Summarize cell-line response distributions for frozen MAPK markers."""
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    rows: list[dict[str, float | int | str]] = []
    symbols = genes["gene_symbol"].astype(str)
    for marker in markers:
        matches = np.flatnonzero(symbols.eq(marker).to_numpy())
        if len(matches) != 1:
            raise ValueError(f"expected one row for marker {marker}, found {len(matches)}")
        values = response[:, matches[0]].astype(np.float64)
        boot_indices = rng.integers(0, len(values), size=(BOOTSTRAP_DRAWS, len(values)))
        boot_means = values[boot_indices].mean(axis=1)
        rows.append(
            {
                "gene_symbol": marker,
                "gene_id": str(genes.iloc[matches[0]]["gene_id"]),
                "n_cell_lines": len(values),
                "mean_delta_log1p_cpm": float(values.mean()),
                "median_delta_log1p_cpm": float(np.median(values)),
                "ci95_low": float(np.quantile(boot_means, 0.025)),
                "ci95_high": float(np.quantile(boot_means, 0.975)),
                "fraction_negative": float((values < 0).mean()),
            }
        )
    return pd.DataFrame(rows)


def run_qc(
    root: Path,
    genes: pd.DataFrame,
    primary_lines: list[str],
    sensitivity_lines: list[str],
    dmso6_counts: sparse.spmatrix,
    dmso24_counts: sparse.spmatrix,
    dmso6_log: np.ndarray,
    dmso24_log: np.ndarray,
    tram24_log: np.ndarray,
    cell_qc: pd.DataFrame,
) -> dict[str, object]:
    """Calculate control-time, marker, and source-cell QC tables and figures."""
    gene_mask = control_qc_gene_mask(dmso6_counts, dmso24_counts)
    if int(gene_mask.sum()) == 0:
        raise ValueError("control-time QC gene universe is empty")

    control_6h = dmso6_log[:, gene_mask].astype(np.float64)
    control_24h = dmso24_log[:, gene_mask].astype(np.float64)
    pcc = row_pearson(control_6h, control_24h)
    spearman = row_pearson(rankdata(control_6h, axis=1), rankdata(control_24h, axis=1))
    control_rmse = np.sqrt(np.mean((control_24h - control_6h) ** 2, axis=1))
    tram_rmse = np.sqrt(
        np.mean((tram24_log[:, gene_mask] - control_24h) ** 2, axis=1)
    )
    qc_table = pd.DataFrame(
        {
            "cell_line": sensitivity_lines,
            "dmso_6h_vs_24h_pcc": pcc,
            "dmso_6h_vs_24h_spearman": spearman,
            "dmso_6h_vs_24h_rmse": control_rmse,
            "trametinib_24h_vs_dmso_24h_rmse": tram_rmse,
            "control_to_treatment_rmse_ratio": control_rmse / tram_rmse,
            "primary_strict_eligible": [line in set(primary_lines) for line in sensitivity_lines],
        }
    )

    primary_positions = [sensitivity_lines.index(line) for line in primary_lines]
    primary_response = tram24_log[primary_positions, :] - dmso24_log[primary_positions, :]
    marker_table = summarize_markers(genes, primary_response)

    table_dir = root / "results" / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    qc_table.to_csv(table_dir / "control_time_qc.csv", index=False)
    marker_table.to_csv(table_dir / "marker_response_summary.csv", index=False)
    cell_qc.to_csv(table_dir / "source_cell_qc.csv", index=False)

    from yakseopdong.plots import write_qc_figures

    write_qc_figures(root, qc_table, marker_table)

    strict_qc = qc_table.loc[qc_table["primary_strict_eligible"]]
    summary: dict[str, object] = {
        "qc_gene_rule": {
            "min_total_control_counts": QC_MIN_TOTAL_CONTROL_COUNTS,
            "min_expressing_lines_per_time": QC_MIN_LINES_PER_TIME,
            "gene_count": int(gene_mask.sum()),
            "model_feature_selection": False,
        },
        "control_time_all_97": {
            "median_pcc": float(np.nanmedian(qc_table["dmso_6h_vs_24h_pcc"])),
            "p05_pcc": float(np.nanquantile(qc_table["dmso_6h_vs_24h_pcc"], 0.05)),
            "median_spearman": float(np.nanmedian(qc_table["dmso_6h_vs_24h_spearman"])),
            "median_rmse": float(np.median(qc_table["dmso_6h_vs_24h_rmse"])),
        },
        "strict_94_comparison": {
            "median_control_time_rmse": float(
                np.median(strict_qc["dmso_6h_vs_24h_rmse"])
            ),
            "median_trametinib_response_rmse": float(
                np.median(strict_qc["trametinib_24h_vs_dmso_24h_rmse"])
            ),
            "median_control_to_treatment_rmse_ratio": float(
                np.median(strict_qc["control_to_treatment_rmse_ratio"])
            ),
        },
        "marker_response": marker_table.to_dict(orient="records"),
        "pooling_decision": "sensitivity_only",
        "pooling_used_for_primary": False,
    }
    pd.DataFrame(
        [
            {"metric": "qc_gene_count", "value": int(gene_mask.sum()), "cohort": "97"},
            {
                "metric": "dmso_6h_vs_24h_median_pcc",
                "value": summary["control_time_all_97"]["median_pcc"],
                "cohort": "97",
            },
            {
                "metric": "dmso_6h_vs_24h_p05_pcc",
                "value": summary["control_time_all_97"]["p05_pcc"],
                "cohort": "97",
            },
            {
                "metric": "dmso_6h_vs_24h_median_spearman",
                "value": summary["control_time_all_97"]["median_spearman"],
                "cohort": "97",
            },
            {
                "metric": "control_time_source_median_rmse",
                "value": summary["strict_94_comparison"]["median_control_time_rmse"],
                "cohort": "94",
            },
            {
                "metric": "trametinib_response_median_rmse",
                "value": summary["strict_94_comparison"]["median_trametinib_response_rmse"],
                "cohort": "94",
            },
            {
                "metric": "control_to_treatment_median_rmse_ratio",
                "value": summary["strict_94_comparison"][
                    "median_control_to_treatment_rmse_ratio"
                ],
                "cohort": "94",
            },
        ]
    ).to_csv(table_dir / "qc_summary.csv", index=False)
    log_path = root / "results" / "logs" / "qc_summary.json"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
