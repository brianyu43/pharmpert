import inspect

import numpy as np
from scipy import sparse

from yakseopdong.robustness import (
    _sample_indices_by_line,
    aggregate_selected_cells,
    fixed_outer_predictions,
)


def test_aggregate_selected_cells_matches_manual_sums() -> None:
    matrix = sparse.csc_matrix(
        np.asarray(
            [
                [1, 2, 3, 4],
                [0, 5, 0, 6],
                [7, 0, 8, 0],
            ],
            dtype=np.int32,
        )
    )
    result = aggregate_selected_cells(
        matrix,
        np.asarray([0, 1, 2, 3], dtype=np.int64),
        np.asarray([0, 0, 1, 1], dtype=np.int64),
        n_groups=2,
    )
    expected_counts = np.asarray([[3, 5, 7], [7, 6, 8]], dtype=float)
    expected = np.log1p(expected_counts * 1_000_000 / expected_counts.sum(axis=1)[:, None])
    np.testing.assert_allclose(result, expected, rtol=1e-6)


def test_sampling_is_deterministic_and_without_replacement() -> None:
    indices = [np.arange(0, 10), np.arange(10, 20)]
    first, codes = _sample_indices_by_line(
        indices, np.random.default_rng(4), sample_n=5, split_half=False
    )
    second, _ = _sample_indices_by_line(
        indices, np.random.default_rng(4), sample_n=5, split_half=False
    )
    np.testing.assert_array_equal(first, second)
    assert len(np.unique(first)) == 10
    np.testing.assert_array_equal(np.bincount(codes), [5, 5])


def test_fixed_outer_predictions_is_deterministic_and_has_no_test_response_api() -> None:
    assert "test_response" not in inspect.signature(fixed_outer_predictions).parameters
    rng = np.random.default_rng(5)
    control = rng.normal(loc=2.0, size=(30, 40))
    response = rng.normal(size=(30, 40))
    folds = np.arange(30) % 5
    first = fixed_outer_predictions(
        control, response, folds, max_variable_genes=30, seed=8
    )
    second = fixed_outer_predictions(
        control, response, folds, max_variable_genes=30, seed=8
    )
    np.testing.assert_allclose(first["b4_rmse"], second["b4_rmse"])
    assert len(first) == 30
