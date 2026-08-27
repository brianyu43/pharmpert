import numpy as np

from yakseopdong.biology import (
    cliffs_delta_bootstrap_permutation,
    eta_squared_permutation,
    partial_spearman_bootstrap_permutation,
    spearman_bootstrap_permutation,
)


def test_spearman_bootstrap_permutation_detects_monotone_signal() -> None:
    x = np.arange(20, dtype=float)
    effect, low, high, p_value = spearman_bootstrap_permutation(
        x, x, seed=7, n_bootstrap=100, n_permutation=200
    )
    assert effect == 1.0
    assert low > 0.99 and high <= 1.0
    assert p_value < 0.02


def test_cliffs_delta_direction_is_mutant_minus_wildtype() -> None:
    values = np.asarray([5, 6, 7, 0, 1, 2], dtype=float)
    binary = np.asarray([True, True, True, False, False, False])
    effect, low, _, p_value = cliffs_delta_bootstrap_permutation(
        values, binary, seed=8, n_bootstrap=100, n_permutation=200
    )
    assert effect == 1.0
    assert low == 1.0
    assert p_value < 0.2


def test_eta_squared_permutation_is_bounded() -> None:
    values = np.asarray([0, 0.1, 4, 4.1, 8, 8.1])
    groups = np.asarray(["a", "a", "b", "b", "c", "c"])
    effect, p_value = eta_squared_permutation(
        values, groups, seed=9, n_permutation=200
    )
    assert 0.9 < effect <= 1.0
    assert 0.0 < p_value < 0.1


def test_partial_spearman_removes_shared_monotone_covariate() -> None:
    rng = np.random.default_rng(12)
    covariate = np.arange(30, dtype=float)
    x = covariate + rng.normal(scale=4.0, size=30)
    y = covariate + rng.normal(scale=4.0, size=30)
    effect, _, _, _ = partial_spearman_bootstrap_permutation(
        x,
        y,
        covariate,
        seed=10,
        n_bootstrap=100,
        n_permutation=200,
    )
    assert abs(effect) < 0.5
