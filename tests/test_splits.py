import numpy as np
import pytest

from yakseopdong.splits import assert_disjoint_groups, lineage_aware_fold_ids


def test_disjoint_cell_lines_pass() -> None:
    assert_disjoint_groups(["A", "B"], ["C"])


def test_overlapping_cell_line_fails() -> None:
    with pytest.raises(ValueError, match="group leakage"):
        assert_disjoint_groups(["A", "B"], ["B", "C"])


def test_lineage_aware_folds_are_deterministic_balanced_and_spread() -> None:
    groups = [f"CL{index:02d}" for index in range(23)]
    lineages = ["lung"] * 9 + ["skin"] * 6 + ["brain"] * 5 + ["rare"] * 3
    first = lineage_aware_fold_ids(groups, lineages, n_splits=5, seed=11)
    second = lineage_aware_fold_ids(groups, lineages, n_splits=5, seed=11)
    np.testing.assert_array_equal(first, second)
    sizes = np.bincount(first, minlength=5)
    assert sizes.max() - sizes.min() <= 1
    assert len(set(first[:9])) == 5


def test_lineage_aware_folds_reject_duplicate_groups() -> None:
    with pytest.raises(ValueError, match="unique"):
        lineage_aware_fold_ids(["A", "A", "B"], ["x", "x", "y"], 2, 1)
