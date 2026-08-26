"""Leakage-safe evaluation primitives for response prediction."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.stats import rankdata

EPSILON = 1e-8


def _paired_vectors(observed: ArrayLike, predicted: ArrayLike) -> tuple[NDArray, NDArray]:
    obs = np.asarray(observed, dtype=float).ravel()
    pred = np.asarray(predicted, dtype=float).ravel()
    if obs.shape != pred.shape:
        raise ValueError(f"shape mismatch: observed={obs.shape}, predicted={pred.shape}")
    if obs.size == 0:
        raise ValueError("metric inputs must not be empty")
    if not (np.isfinite(obs).all() and np.isfinite(pred).all()):
        raise ValueError("metric inputs must contain only finite values")
    return obs, pred


def rmse(observed: ArrayLike, predicted: ArrayLike) -> float:
    """Root mean squared error over a single response vector."""
    obs, pred = _paired_vectors(observed, predicted)
    return float(np.sqrt(np.mean(np.square(pred - obs))))


def nrmse(observed: ArrayLike, predicted: ArrayLike, epsilon: float = EPSILON) -> float:
    """RMSE normalized by the root mean square magnitude of the observation."""
    obs, pred = _paired_vectors(observed, predicted)
    denominator = max(float(np.sqrt(np.mean(np.square(obs)))), epsilon)
    return rmse(obs, pred) / denominator


def pearson_or_nan(observed: ArrayLike, predicted: ArrayLike) -> float:
    """Pearson correlation, returning NaN for either constant vector."""
    obs, pred = _paired_vectors(observed, predicted)
    if np.ptp(obs) == 0 or np.ptp(pred) == 0:
        return float("nan")
    return float(np.corrcoef(obs, pred)[0, 1])


def spearman_or_nan(observed: ArrayLike, predicted: ArrayLike) -> float:
    """Spearman correlation, returning NaN for either constant vector."""
    obs, pred = _paired_vectors(observed, predicted)
    if np.ptp(obs) == 0 or np.ptp(pred) == 0:
        return float("nan")
    return pearson_or_nan(rankdata(obs), rankdata(pred))


def context_residual(response: ArrayLike, training_mean_response: ArrayLike) -> NDArray:
    """Subtract the outer-training mean response from a response vector."""
    response_vector, training_mean = _paired_vectors(response, training_mean_response)
    return response_vector - training_mean


def rmse_gain_vs_b1(
    observed: ArrayLike, predicted: ArrayLike, training_mean_response: ArrayLike
) -> float:
    """Positive values mean the prediction improves on B1 global mean."""
    return rmse(observed, training_mean_response) - rmse(observed, predicted)


def signed_topk_overlap(observed: ArrayLike, predicted: ArrayLike, k: int = 50) -> float:
    """Mean overlap fraction for top upregulated and downregulated genes."""
    obs, pred = _paired_vectors(observed, predicted)
    if k <= 0:
        raise ValueError("k must be positive")
    if 2 * k > obs.size:
        raise ValueError("2 * k must not exceed the gene universe size")

    obs_order = np.argsort(obs, kind="stable")
    pred_order = np.argsort(pred, kind="stable")
    down_overlap = len(set(obs_order[:k]) & set(pred_order[:k])) / k
    up_overlap = len(set(obs_order[-k:]) & set(pred_order[-k:])) / k
    return float((down_overlap + up_overlap) / 2)


def bootstrap_mean_interval(
    values: ArrayLike,
    n_bootstrap: int = 2_000,
    seed: int = 20260827,
) -> tuple[float, float, float]:
    """Cell-line bootstrap mean and percentile 95% interval."""
    array = np.asarray(values, dtype=float).ravel()
    array = array[np.isfinite(array)]
    if array.size == 0:
        return float("nan"), float("nan"), float("nan")
    if n_bootstrap <= 0:
        raise ValueError("n_bootstrap must be positive")
    rng = np.random.default_rng(seed)
    draws = rng.choice(array, size=(n_bootstrap, len(array)), replace=True).mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(array.mean()), float(low), float(high)
