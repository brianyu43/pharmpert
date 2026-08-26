import math

import numpy as np
import pytest

from yakseopdong.metrics import (
    bootstrap_mean_interval,
    context_residual,
    nrmse,
    pearson_or_nan,
    rmse,
    rmse_gain_vs_b1,
    signed_topk_overlap,
    spearman_or_nan,
)


def test_rmse_and_nrmse() -> None:
    observed = np.array([1.0, -1.0])
    predicted = np.array([0.0, 0.0])
    assert rmse(observed, predicted) == pytest.approx(1.0)
    assert nrmse(observed, predicted) == pytest.approx(1.0)


def test_b0_constant_correlation_is_nan() -> None:
    assert math.isnan(pearson_or_nan([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]))
    assert math.isnan(spearman_or_nan([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]))


def test_spearman_uses_rank_order() -> None:
    assert spearman_or_nan([1.0, 2.0, 3.0], [10.0, 20.0, 40.0]) == pytest.approx(1.0)


def test_context_residual_uses_training_mean() -> None:
    residual = context_residual([2.0, -1.0], [1.0, -2.0])
    np.testing.assert_allclose(residual, [1.0, 1.0])


def test_positive_rmse_gain_means_better_than_b1() -> None:
    observed = [2.0, -2.0]
    prediction = [1.5, -1.5]
    training_mean = [0.0, 0.0]
    assert rmse_gain_vs_b1(observed, prediction, training_mean) > 0


def test_signed_topk_overlap_is_direction_aware() -> None:
    observed = [-3.0, -2.0, -1.0, 1.0, 2.0, 3.0]
    assert signed_topk_overlap(observed, observed, k=2) == 1.0
    assert signed_topk_overlap(observed, observed[::-1], k=2) == 0.0


def test_metric_shape_mismatch_fails() -> None:
    with pytest.raises(ValueError, match="shape mismatch"):
        rmse([1.0], [1.0, 2.0])


def test_bootstrap_interval_is_deterministic() -> None:
    first = bootstrap_mean_interval([1.0, 2.0, 3.0], n_bootstrap=100, seed=7)
    second = bootstrap_mean_interval([1.0, 2.0, 3.0], n_bootstrap=100, seed=7)
    assert first == second
    assert first[0] == pytest.approx(2.0)
