"""Frozen pathway-panel mapping and W6 component diagnostics for W7."""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from numpy.typing import NDArray
from scipy.stats import hypergeom


def load_geneset_config(path: Path) -> dict[str, Any]:
    """Load a frozen gene-set snapshot and reject incomplete configurations."""
    with path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if config.get("status") != "frozen":
        raise ValueError("gene-set configuration must be frozen before W7 evaluation")
    collections = config.get("collections")
    if not isinstance(collections, dict) or not collections:
        raise ValueError("gene-set configuration has no collections")
    for name, collection in collections.items():
        genes = collection.get("genes", [])
        if not genes or len(genes) != len(set(genes)):
            raise ValueError(f"collection {name} must contain unique genes")
    return config


def gene_symbol_indices(
    genes: pd.DataFrame, symbols: set[str]
) -> NDArray[np.int64]:
    """Map a fixed symbol set to every aligned gene-matrix column carrying it."""
    if "gene_symbol" not in genes:
        raise ValueError("gene metadata lacks gene_symbol")
    mask = genes["gene_symbol"].astype(str).isin(symbols).to_numpy()
    indices = np.flatnonzero(mask).astype(np.int64, copy=False)
    if not len(indices):
        raise ValueError("the fixed gene panel has no overlap with the matrix")
    return indices


def pathway_panel_indices(
    genes: pd.DataFrame, config: dict[str, Any]
) -> NDArray[np.int64]:
    """Return the union of all frozen input-panel collections in matrix order."""
    collection_names = config["input_panel_collections"]
    symbols = set().union(
        *(set(config["collections"][name]["genes"]) for name in collection_names)
    )
    return gene_symbol_indices(genes, symbols)


def pathway_coverage_table(
    genes: pd.DataFrame,
    config: dict[str, Any],
    fold_selected_counts: dict[int, int] | None = None,
) -> pd.DataFrame:
    """Summarize fixed-set symbol coverage and per-fold input selection."""
    dataset_symbols = set(genes["gene_symbol"].astype(str))
    rows: list[dict[str, object]] = []
    panel_symbols: set[str] = set()
    for name, collection in config["collections"].items():
        symbols = set(collection["genes"])
        panel_symbols.update(symbols)
        rows.append(
            {
                "row_type": "collection",
                "collection": name,
                "outer_fold": pd.NA,
                "defined_symbols": len(symbols),
                "mapped_symbols": len(symbols & dataset_symbols),
                "mapped_matrix_columns": len(gene_symbol_indices(genes, symbols)),
                "selected_training_columns": pd.NA,
            }
        )
    rows.append(
        {
            "row_type": "panel_union",
            "collection": "all_input_collections",
            "outer_fold": pd.NA,
            "defined_symbols": len(panel_symbols),
            "mapped_symbols": len(panel_symbols & dataset_symbols),
            "mapped_matrix_columns": len(gene_symbol_indices(genes, panel_symbols)),
            "selected_training_columns": pd.NA,
        }
    )
    for fold, selected in sorted((fold_selected_counts or {}).items()):
        rows.append(
            {
                "row_type": "outer_training_selection",
                "collection": "all_input_collections",
                "outer_fold": fold,
                "defined_symbols": len(panel_symbols),
                "mapped_symbols": len(panel_symbols & dataset_symbols),
                "mapped_matrix_columns": len(gene_symbol_indices(genes, panel_symbols)),
                "selected_training_columns": selected,
            }
        )
    return pd.DataFrame(rows)


def benjamini_hochberg(p_values: NDArray[np.floating]) -> NDArray[np.float64]:
    """Return monotone Benjamini-Hochberg adjusted p-values."""
    values = np.asarray(p_values, dtype=np.float64)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("p-values must be a finite vector")
    if ((values < 0) | (values > 1)).any():
        raise ValueError("p-values must lie in [0, 1]")
    order = np.argsort(values, kind="stable")
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return output


def component_pathway_enrichment(
    loadings: pd.DataFrame,
    genes: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Test W6 top-loading genes against each prespecified collection."""
    required = {"outer_fold", "component", "direction", "gene_symbol"}
    if not required.issubset(loadings.columns):
        raise ValueError("component loading table is missing required columns")
    universe = set(genes["gene_symbol"].astype(str))
    rows: list[dict[str, object]] = []
    for (outer_fold, component, direction), group in loadings.groupby(
        ["outer_fold", "component", "direction"], sort=True
    ):
        top_symbols = set(group["gene_symbol"].astype(str)) & universe
        for name, collection in config["collections"].items():
            pathway = set(collection["genes"]) & universe
            overlap = top_symbols & pathway
            p_value = float(
                hypergeom.sf(
                    len(overlap) - 1,
                    len(universe),
                    len(pathway),
                    len(top_symbols),
                )
            )
            rows.append(
                {
                    "outer_fold": int(outer_fold),
                    "component": int(component),
                    "direction": str(direction),
                    "collection": name,
                    "top_gene_count": len(top_symbols),
                    "pathway_gene_count": len(pathway),
                    "overlap_count": len(overlap),
                    "overlap_genes": ";".join(sorted(overlap)),
                    "p_value": p_value,
                }
            )
    result = pd.DataFrame(rows)
    result["fdr_bh"] = benjamini_hochberg(result["p_value"].to_numpy())
    return result


def cclr_subspace_stability(model_paths: list[Path]) -> pd.DataFrame:
    """Compare fold-specific W6 response subspaces via principal-angle overlap."""
    components: dict[int, NDArray[np.float64]] = {}
    for path in model_paths:
        fold = int(path.stem.rsplit("_", maxsplit=1)[-1])
        with np.load(path) as artifact:
            matrix = np.asarray(artifact["response_pca_components"], dtype=np.float64)
        components[fold] = matrix
    rows: list[dict[str, object]] = []
    for fold_a, fold_b in combinations(sorted(components), 2):
        left = components[fold_a]
        right = components[fold_b]
        rank = min(len(left), len(right))
        singular_values = np.linalg.svd(left[:rank] @ right[:rank].T, compute_uv=False)
        rows.append(
            {
                "outer_fold_a": fold_a,
                "outer_fold_b": fold_b,
                "common_rank": rank,
                "mean_squared_cosine": float(np.mean(np.square(singular_values))),
                "minimum_cosine": float(np.min(singular_values)),
                "maximum_cosine": float(np.max(singular_values)),
            }
        )
    return pd.DataFrame(rows)
