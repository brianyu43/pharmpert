from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

from yakseopdong.distribution import (
    energy_distance_multivariate,
    extract_h5ad_csc_rows,
    sliced_wasserstein_distance,
    sparse_log1p_cpm,
)


def _write_sparse_h5ad(path: Path, matrix: sparse.csc_matrix) -> None:
    with h5py.File(path, "w") as handle:
        group = handle.create_group("X")
        group.attrs["encoding-type"] = "csc_matrix"
        group.attrs["shape"] = matrix.shape
        group.create_dataset("data", data=matrix.data)
        group.create_dataset("indices", data=matrix.indices)
        group.create_dataset("indptr", data=matrix.indptr)


def test_extract_h5ad_csc_rows_preserves_requested_order(tmp_path: Path) -> None:
    matrix = sparse.csc_matrix(np.arange(20, dtype=np.float32).reshape(4, 5))
    path = tmp_path / "small.h5ad"
    _write_sparse_h5ad(path, matrix)
    extracted = extract_h5ad_csc_rows(
        path, np.asarray([3, 1], dtype=np.int64), chunk_size=2
    )
    np.testing.assert_array_equal(extracted.toarray(), matrix.toarray()[[3, 1]])


def test_sparse_log1p_cpm_keeps_sparse_shape() -> None:
    matrix = sparse.csr_matrix([[1, 0, 3], [0, 2, 2]], dtype=np.float32)
    normalized = sparse_log1p_cpm(matrix)
    assert normalized.shape == matrix.shape
    assert normalized.nnz == matrix.nnz
    np.testing.assert_allclose(
        normalized.toarray().sum(axis=1),
        np.log1p(np.asarray([[250_000, 0, 750_000], [0, 500_000, 500_000]])).sum(axis=1),
    )


def test_distribution_distances_are_zero_for_identical_samples() -> None:
    rng = np.random.default_rng(4)
    values = rng.normal(size=(12, 3))
    projections = rng.normal(size=(20, 3))
    assert energy_distance_multivariate(values, values) == 0.0
    assert sliced_wasserstein_distance(values, values, projections) == 0.0
    assert energy_distance_multivariate(values, values + 1.0) > 0
    assert sliced_wasserstein_distance(values, values + 1.0, projections) > 0
