from pathlib import Path

import pandas as pd

from yakseopdong.core_audit import build_cell_count_matrix


def _cells(cell_line: str, depmap_id: str, n_normal: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "singlet_ID": [cell_line] * n_normal,
            "DepMap_ID": [depmap_id] * n_normal,
            "cell_quality": ["normal"] * n_normal,
        }
    )


def test_strict_and_pooled_cohorts_are_separate() -> None:
    dmso_6h = _cells("LINE_A", "ACH-1", 15)
    dmso_24h = _cells("LINE_A", "ACH-1", 10)
    trametinib_24h = _cells("LINE_A", "ACH-1", 20)

    counts = build_cell_count_matrix(dmso_6h, dmso_24h, trametinib_24h).iloc[0]
    assert counts["dmso_pooled_normal"] == 25
    assert bool(counts["primary_strict_eligible"]) is False
    assert bool(counts["reproduction_pooled_eligible"]) is True


def test_tracked_core_count_matrix_matches_frozen_cohorts() -> None:
    root = Path(__file__).resolve().parents[1]
    counts = pd.read_csv(root / "cell_count_matrix.csv")
    assert len(counts) == 97
    assert counts["cell_line"].is_unique
    assert counts["depmap_id"].is_unique
    assert int(counts["primary_strict_eligible"].sum()) == 94
    assert int(counts["reproduction_pooled_eligible"].sum()) == 97
