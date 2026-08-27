from pathlib import Path

import numpy as np

from yakseopdong.pathways import benjamini_hochberg, cclr_subspace_stability


def test_benjamini_hochberg_is_bounded_and_monotone_in_rank_order() -> None:
    values = np.asarray([0.04, 0.001, 0.2, 0.01])
    adjusted = benjamini_hochberg(values)
    assert np.all((0 <= adjusted) & (adjusted <= 1))
    order = np.argsort(values)
    assert np.all(np.diff(adjusted[order]) >= 0)


def test_identical_response_subspaces_have_unit_overlap(tmp_path: Path) -> None:
    components = np.eye(3, 6)
    paths = []
    for fold in range(3):
        path = tmp_path / f"cclr_outer_fold_{fold}.npz"
        np.savez_compressed(path, response_pca_components=components)
        paths.append(path)
    stability = cclr_subspace_stability(paths)
    assert len(stability) == 3
    np.testing.assert_allclose(stability["mean_squared_cosine"], 1.0)
