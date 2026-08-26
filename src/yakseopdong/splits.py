"""Split validation utilities."""

from __future__ import annotations

from collections.abc import Iterable


def assert_disjoint_groups(train_groups: Iterable[str], test_groups: Iterable[str]) -> None:
    """Raise when a generalization group appears in both train and test."""
    train = set(train_groups)
    test = set(test_groups)
    overlap = train & test
    if overlap:
        raise ValueError(f"train/test group leakage: {sorted(overlap)}")
