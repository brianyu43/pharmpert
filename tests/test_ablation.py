import inspect

import numpy as np

from yakseopdong.ablation import predict_fixed_low_rank
from yakseopdong.models import fit_response_embedding


def test_fixed_low_rank_prediction_has_no_test_response_api_and_is_deterministic() -> None:
    assert "test_response" not in inspect.signature(predict_fixed_low_rank).parameters
    rng = np.random.default_rng(19)
    train_response = rng.normal(size=(18, 25))
    train_predictors = rng.normal(size=(18, 4))
    test_predictors = rng.normal(size=(3, 4))
    embedding = fit_response_embedding(train_response, max_components=6, seed=11)
    first = predict_fixed_low_rank(
        train_response,
        train_predictors,
        test_predictors,
        embedding,
        response_rank=3,
        alpha=100.0,
    )
    second = predict_fixed_low_rank(
        train_response,
        train_predictors,
        test_predictors,
        embedding,
        response_rank=3,
        alpha=100.0,
    )
    assert first.shape == (3, 25)
    np.testing.assert_allclose(first, second)
