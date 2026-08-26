from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from yakseopdong.pseudobulk import (
    aggregate_normal_cells,
    log1p_cpm,
    read_vector_parquet,
    write_vector_parquet,
)


def test_aggregate_normal_cells_filters_and_groups() -> None:
    matrix = sparse.coo_matrix(
        np.array(
            [
                [1, 2, 100, 4],
                [0, 3, 100, 5],
            ],
            dtype=np.int32,
        )
    )
    classifications = pd.DataFrame(
        {
            "singlet_ID": ["B", "A", "A", "A"],
            "DepMap_ID": ["ACH-B", "ACH-A", "ACH-A", "ACH-A"],
            "cell_quality": ["normal", "normal", "doublet", "normal"],
        }
    )

    lines, depmap, counts, n_cells = aggregate_normal_cells(matrix, classifications)
    assert lines == ("A", "B")
    assert depmap == ("ACH-A", "ACH-B")
    np.testing.assert_array_equal(counts.toarray(), [[6, 8], [1, 0]])
    np.testing.assert_array_equal(n_cells, [2, 1])


def test_log1p_cpm_has_fixed_row_sum_before_log() -> None:
    counts = sparse.csr_matrix([[1, 1], [3, 1]], dtype=np.int64)
    values, library_sizes = log1p_cpm(counts)
    np.testing.assert_array_equal(library_sizes, [2, 4])
    np.testing.assert_allclose(values[0], np.log1p([500_000, 500_000]), rtol=1e-6)
    np.testing.assert_allclose(values[1], np.log1p([750_000, 250_000]), rtol=1e-6)


def test_vector_parquet_round_trip(tmp_path: Path) -> None:
    metadata = pd.DataFrame({"cell_line": ["A", "B"], "n_cells": [20, 30]})
    values = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32)
    path = tmp_path / "vectors.parquet"
    write_vector_parquet(path, metadata, values, "expression")
    observed_metadata, observed_values = read_vector_parquet(path, "expression")
    pd.testing.assert_frame_equal(observed_metadata, metadata)
    np.testing.assert_array_equal(observed_values, values)


def test_tracked_processed_manifest_matches_frozen_contract() -> None:
    root = Path(__file__).resolve().parents[1]
    manifest = pd.read_csv(root / "processed_manifest.csv").set_index("path")
    assert manifest.loc["data/processed/pseudobulk_24h.parquet", "row_count"] == 188
    assert manifest.loc["data/processed/response_24h.parquet", "row_count"] == 94
    assert manifest.loc["data/processed/response_pooled_sensitivity.parquet", "row_count"] == 97
    assert set(manifest["vector_length"]) == {0, 32_738}
    assert manifest["sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
