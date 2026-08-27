import inspect

import numpy as np
import pandas as pd
from scipy import sparse

from yakseopdong.temporal import (
    _fixed_transfer_predictions,
    aggregate_grouped_rows,
    eligibility_from_counts,
    parse_timecourse_tag,
)


def test_parse_timecourse_tag_uses_hash_assignment_contract() -> None:
    assert parse_timecourse_tag("DMSO_3hr") == ("control", 3)
    assert parse_timecourse_tag("Tram_48hr") == ("trametinib", 48)
    assert parse_timecourse_tag("Untreated_48hr") == ("untreated", 48)
    assert parse_timecourse_tag("unknown") is None
    assert parse_timecourse_tag("multiplet") is None


def test_grouped_row_aggregation_preserves_raw_counts() -> None:
    matrix = sparse.csc_matrix(
        np.asarray(
            [
                [1, 0, 2, 0],
                [0, 3, 0, 4],
                [5, 0, 0, 6],
                [0, 7, 8, 0],
            ],
            dtype=np.int32,
        )
    )
    result = aggregate_grouped_rows(
        matrix,
        np.asarray([0, 1, 2, 3], dtype=np.int64),
        np.asarray([0, 0, 1, 1], dtype=np.int64),
        n_groups=2,
        chunk_size=2,
    )
    np.testing.assert_array_equal(result, [[1, 3, 2, 4], [5, 7, 8, 6]])


def test_eligibility_requires_all_ten_matched_groups() -> None:
    rows = []
    for line, minimum in (("A", 10), ("B", 9)):
        for condition in ("control", "trametinib"):
            for hour in (3, 6, 12, 24, 48):
                rows.append(
                    {
                        "cell_line": line,
                        "depmap_id": f"D-{line}",
                        "disease": "x",
                        "condition": condition,
                        "time_hours": hour,
                        "n_cells": minimum,
                    }
                )
    result = eligibility_from_counts(pd.DataFrame(rows)).set_index("cell_line")
    assert bool(result.loc["A", "eligible_t10"])
    assert not bool(result.loc["A", "eligible_t20"])
    assert not bool(result.loc["B", "eligible_t10"])


def test_external_transfer_has_no_test_response_argument() -> None:
    assert "test_response" not in inspect.signature(_fixed_transfer_predictions).parameters
