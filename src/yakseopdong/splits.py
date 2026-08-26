"""Deterministic, lineage-aware cell-line split utilities."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from yakseopdong.metadata import run_metadata_audit


def assert_disjoint_groups(train_groups: Iterable[str], test_groups: Iterable[str]) -> None:
    """Raise when a generalization group appears in both train and test."""
    train = set(train_groups)
    test = set(test_groups)
    overlap = train & test
    if overlap:
        raise ValueError(f"train/test group leakage: {sorted(overlap)}")


def lineage_aware_fold_ids(
    groups: Iterable[str],
    lineages: Iterable[str],
    n_splits: int,
    seed: int,
) -> np.ndarray:
    """Assign unique groups to balanced folds while spreading each lineage."""
    group_array = np.asarray(list(groups), dtype=str)
    lineage_array = np.asarray(list(lineages), dtype=str)
    if len(group_array) != len(lineage_array):
        raise ValueError("groups and lineages must have the same length")
    if len(np.unique(group_array)) != len(group_array):
        raise ValueError("generalization groups must be unique")
    if not 2 <= n_splits <= len(group_array):
        raise ValueError("n_splits must be between 2 and the number of groups")

    rng = np.random.default_rng(seed)
    fold_sizes = np.zeros(n_splits, dtype=int)
    assignments = np.full(len(group_array), -1, dtype=int)
    unique_lineages = sorted(
        np.unique(lineage_array),
        key=lambda name: (-int(np.sum(lineage_array == name)), name),
    )
    for lineage in unique_lineages:
        indices = np.flatnonzero(lineage_array == lineage)
        indices = indices[rng.permutation(len(indices))]
        lineage_fold_counts = np.zeros(n_splits, dtype=int)
        tie_order = rng.permutation(n_splits)
        tie_rank = np.empty(n_splits, dtype=int)
        tie_rank[tie_order] = np.arange(n_splits)
        for index in indices:
            fold = min(
                range(n_splits),
                key=lambda candidate: (
                    lineage_fold_counts[candidate],
                    fold_sizes[candidate],
                    tie_rank[candidate],
                ),
            )
            assignments[index] = fold
            lineage_fold_counts[fold] += 1
            fold_sizes[fold] += 1

    if (assignments < 0).any():
        raise RuntimeError("split assignment is incomplete")
    if fold_sizes.max() - fold_sizes.min() > 1:
        raise RuntimeError(f"fold sizes are unexpectedly imbalanced: {fold_sizes.tolist()}")
    return assignments


def validate_outer_assignments(assignments: pd.DataFrame, n_splits: int) -> None:
    """Validate the one-test-fold-per-cell-line outer split contract."""
    required = {"cell_line", "depmap_id", "lineage", "outer_fold"}
    missing = required - set(assignments.columns)
    if missing:
        raise ValueError(f"missing split columns: {sorted(missing)}")
    if not assignments["cell_line"].is_unique or not assignments["depmap_id"].is_unique:
        raise ValueError("cell-line and DepMap identifiers must be unique")
    if set(assignments["outer_fold"]) != set(range(n_splits)):
        raise ValueError("outer fold identifiers are incomplete")
    sizes = assignments.groupby("outer_fold").size()
    if sizes.max() - sizes.min() > 1:
        raise ValueError(f"outer fold sizes are imbalanced: {sizes.to_dict()}")
    for fold in range(n_splits):
        test = assignments.loc[assignments["outer_fold"].eq(fold), "cell_line"]
        train = assignments.loc[~assignments["outer_fold"].eq(fold), "cell_line"]
        assert_disjoint_groups(train, test)


def run_splits(root: Path, seed: int = 20260827) -> dict[str, object]:
    """Create and validate the frozen outer and nested-inner assignments."""
    run_metadata_audit(root)
    annotations = pd.read_csv(root / "cell_line_annotations.csv").sort_values(
        "cell_line", ignore_index=True
    )
    annotations["outer_fold"] = lineage_aware_fold_ids(
        annotations["cell_line"], annotations["lineage"], n_splits=5, seed=seed
    )
    outer = annotations[["cell_line", "depmap_id", "lineage", "outer_fold"]].copy()
    validate_outer_assignments(outer, n_splits=5)
    outer.to_csv(root / "split_assignments.csv", index=False)

    inner_parts: list[pd.DataFrame] = []
    for outer_fold in range(5):
        train = outer.loc[~outer["outer_fold"].eq(outer_fold)].copy()
        train["inner_fold"] = lineage_aware_fold_ids(
            train["cell_line"],
            train["lineage"],
            n_splits=4,
            seed=seed + 101 * (outer_fold + 1),
        )
        train.insert(0, "outer_fold", train.pop("outer_fold"))
        train["outer_fold"] = outer_fold
        inner_parts.append(train)
    inner = pd.concat(inner_parts, ignore_index=True)
    if inner.duplicated(["outer_fold", "cell_line"]).any():
        raise ValueError("duplicate inner split assignment")
    inner.to_csv(root / "inner_split_assignments.csv", index=False)

    outer_sizes = outer.groupby("outer_fold").size().astype(int).to_dict()
    inner_sizes = {
        str(outer_fold): (
            inner.loc[inner["outer_fold"].eq(outer_fold)]
            .groupby("inner_fold")
            .size()
            .astype(int)
            .to_dict()
        )
        for outer_fold in range(5)
    }
    lineage_counts = annotations["lineage"].value_counts().sort_index().astype(int).to_dict()
    report = {
        "status": "frozen_v1",
        "seed": seed,
        "generalization_unit": "cell_line",
        "strategy": "greedy lineage-spread with balanced fold sizes",
        "outer_fold_sizes": outer_sizes,
        "inner_fold_sizes_by_outer": inner_sizes,
        "lineage_counts": lineage_counts,
        "all_outer_train_test_disjoint": True,
        "each_cell_line_tested_once": True,
    }
    path = root / "results" / "logs" / "split_audit.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
