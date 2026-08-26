import pytest

from yakseopdong.splits import assert_disjoint_groups


def test_disjoint_cell_lines_pass() -> None:
    assert_disjoint_groups(["A", "B"], ["C"])


def test_overlapping_cell_line_fails() -> None:
    with pytest.raises(ValueError, match="group leakage"):
        assert_disjoint_groups(["A", "B"], ["B", "C"])
