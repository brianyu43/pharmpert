import numpy as np
import pandas as pd
from scipy import sparse

from yakseopdong.qc import control_qc_gene_mask, row_pearson, summarize_markers


def test_row_pearson_handles_matching_and_constant_rows() -> None:
    left = np.array([[1.0, 2.0, 3.0], [1.0, 1.0, 1.0]])
    right = np.array([[2.0, 4.0, 6.0], [2.0, 3.0, 4.0]])
    observed = row_pearson(left, right)
    assert observed[0] == 1.0
    assert np.isnan(observed[1])


def test_control_qc_mask_uses_only_control_counts() -> None:
    dmso6 = sparse.csr_matrix(np.tile([10, 1], (10, 1)))
    dmso24 = sparse.csr_matrix(np.tile([10, 1], (10, 1)))
    np.testing.assert_array_equal(control_qc_gene_mask(dmso6, dmso24), [True, False])


def test_marker_summary_requires_and_reports_one_gene() -> None:
    genes = pd.DataFrame(
        {"gene_id": ["ENSG1", "ENSG2"], "gene_symbol": ["EGR1", "DUSP6"]}
    )
    response = np.array([[-1.0, -0.5], [-2.0, 0.5]], dtype=np.float32)
    summary = summarize_markers(genes, response, markers=("EGR1", "DUSP6"))
    assert summary["gene_symbol"].tolist() == ["EGR1", "DUSP6"]
    assert summary.loc[0, "fraction_negative"] == 1.0
    assert summary.loc[1, "mean_delta_log1p_cpm"] == 0.0
