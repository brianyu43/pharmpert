"""Leakage-safe baseline response models for held-out cell lines."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge

from yakseopdong.metrics import rmse


@dataclass(frozen=True)
class ControlEmbedding:
    """Training-fitted gene selector and whitened control PCA."""

    gene_indices: NDArray[np.int64]
    pca: PCA

    def transform(self, control: NDArray[np.floating]) -> NDArray[np.float64]:
        return self.pca.transform(np.asarray(control)[:, self.gene_indices])


@dataclass(frozen=True)
class ResponseEmbedding:
    """Training-fitted low-rank response basis."""

    pca: PCA

    def transform(self, response: NDArray[np.floating]) -> NDArray[np.float64]:
        return self.pca.transform(np.asarray(response, dtype=np.float64))

    def inverse_transform(self, scores: NDArray[np.floating]) -> NDArray[np.float64]:
        return self.pca.inverse_transform(np.asarray(scores, dtype=np.float64))


@dataclass(frozen=True)
class CCLRArtifact:
    """All training-fitted objects needed for one outer-fold CCLR prediction."""

    control_embedding: ControlEmbedding
    response_embedding: ResponseEmbedding
    ridge: Ridge
    control_dimension: int
    response_rank: int
    alpha: float

    def predict(self, control: NDArray[np.floating]) -> NDArray[np.float32]:
        control_scores = self.control_embedding.transform(control)[:, : self.control_dimension]
        response_scores = np.asarray(self.ridge.predict(control_scores), dtype=np.float64)
        if response_scores.ndim == 1:
            response_scores = response_scores[:, None]
        return self.response_embedding.inverse_transform(response_scores).astype(np.float32)


def fit_control_embedding(
    train_control: NDArray[np.floating],
    max_components: int,
    max_variable_genes: int,
    min_mean_log1p_cpm: float,
    seed: int,
) -> ControlEmbedding:
    """Fit control-only feature filtering and PCA on training lines."""
    train = np.asarray(train_control, dtype=np.float64)
    if train.ndim != 2 or len(train) < 3:
        raise ValueError("train_control must be a 2D matrix with at least three rows")
    means = train.mean(axis=0)
    variances = train.var(axis=0)
    eligible = np.flatnonzero((means >= min_mean_log1p_cpm) & (variances > 1e-10))
    if len(eligible) < max_components:
        raise ValueError("too few eligible genes for the requested PCA dimension")
    order = np.argsort(-variances[eligible], kind="stable")
    selected = eligible[order[:max_variable_genes]].astype(np.int64, copy=False)
    n_components = min(max_components, len(train) - 1, len(selected))
    pca = PCA(
        n_components=n_components,
        whiten=True,
        svd_solver="randomized",
        random_state=seed,
    )
    pca.fit(train[:, selected])
    return ControlEmbedding(gene_indices=selected, pca=pca)


def fit_response_embedding(
    train_response: NDArray[np.floating],
    max_components: int,
    seed: int,
) -> ResponseEmbedding:
    """Fit a response PCA basis using training responses only."""
    train = np.asarray(train_response, dtype=np.float64)
    if train.ndim != 2 or len(train) < 3:
        raise ValueError("train_response must be a 2D matrix with at least three rows")
    n_components = min(int(max_components), len(train) - 1, train.shape[1])
    if n_components < 1:
        raise ValueError("response PCA requires at least one feasible component")
    pca = PCA(
        n_components=n_components,
        whiten=False,
        svd_solver="randomized",
        random_state=seed,
    )
    pca.fit(train)
    return ResponseEmbedding(pca=pca)


def nearest_neighbor_predictions(
    train_scores: NDArray[np.floating],
    train_response: NDArray[np.floating],
    query_scores: NDArray[np.floating],
    dimension: int,
) -> tuple[NDArray[np.float32], NDArray[np.int64]]:
    """Copy the response of the nearest training line in whitened PC space."""
    train = np.asarray(train_scores, dtype=np.float64)[:, :dimension]
    query = np.asarray(query_scores, dtype=np.float64)[:, :dimension]
    squared_distances = np.square(query[:, None, :] - train[None, :, :]).sum(axis=2)
    neighbors = np.argmin(squared_distances, axis=1).astype(np.int64)
    predictions = np.asarray(train_response, dtype=np.float32)[neighbors]
    return predictions, neighbors


def _candidate_dimensions(dimensions: list[int], available: int) -> list[int]:
    candidates = sorted({int(value) for value in dimensions if 0 < int(value) <= available})
    if not candidates:
        raise ValueError("no PCA dimension is feasible")
    return candidates


def select_baseline_hyperparameters(
    train_control: NDArray[np.floating],
    train_response: NDArray[np.floating],
    inner_fold_ids: NDArray[np.integer],
    dimensions: list[int],
    ridge_alphas: list[float],
    max_variable_genes: int,
    min_mean_log1p_cpm: float,
    seed: int,
) -> tuple[int, int, float, list[dict[str, float | int | str]]]:
    """Choose B3/B4 settings using only nested training folds."""
    control = np.asarray(train_control)
    response = np.asarray(train_response)
    fold_ids = np.asarray(inner_fold_ids, dtype=int)
    if len(control) != len(response) or len(control) != len(fold_ids):
        raise ValueError("inner-CV arrays have inconsistent row counts")
    unique_folds = sorted(np.unique(fold_ids).tolist())
    if len(unique_folds) < 2:
        raise ValueError("inner CV requires at least two folds")

    b3_losses: dict[int, list[float]] = {dimension: [] for dimension in dimensions}
    b4_losses: dict[tuple[int, float], list[float]] = {
        (dimension, float(alpha)): [] for dimension in dimensions for alpha in ridge_alphas
    }
    trace: list[dict[str, float | int | str]] = []
    for inner_fold in unique_folds:
        validation_mask = fold_ids == inner_fold
        inner_train_mask = ~validation_mask
        max_available = min(max(dimensions), int(inner_train_mask.sum()) - 1)
        embedding = fit_control_embedding(
            control[inner_train_mask],
            max_components=max_available,
            max_variable_genes=max_variable_genes,
            min_mean_log1p_cpm=min_mean_log1p_cpm,
            seed=seed + inner_fold,
        )
        feasible = _candidate_dimensions(dimensions, embedding.pca.n_components_)
        inner_train_scores = embedding.transform(control[inner_train_mask])
        validation_scores = embedding.transform(control[validation_mask])
        y_train = response[inner_train_mask]
        y_validation = response[validation_mask]

        for dimension in feasible:
            b3_prediction, _ = nearest_neighbor_predictions(
                inner_train_scores, y_train, validation_scores, dimension
            )
            line_losses = [
                rmse(observed, predicted)
                for observed, predicted in zip(y_validation, b3_prediction, strict=True)
            ]
            b3_losses[dimension].extend(line_losses)
            trace.append(
                {
                    "model": "B3",
                    "inner_fold": int(inner_fold),
                    "dimension": dimension,
                    "alpha": float("nan"),
                    "mean_rmse": float(np.mean(line_losses)),
                }
            )

            for alpha in ridge_alphas:
                model = Ridge(alpha=float(alpha), fit_intercept=True, solver="svd")
                model.fit(inner_train_scores[:, :dimension], y_train)
                prediction = model.predict(validation_scores[:, :dimension])
                line_losses = [
                    rmse(observed, predicted)
                    for observed, predicted in zip(y_validation, prediction, strict=True)
                ]
                b4_losses[(dimension, float(alpha))].extend(line_losses)
                trace.append(
                    {
                        "model": "B4",
                        "inner_fold": int(inner_fold),
                        "dimension": dimension,
                        "alpha": float(alpha),
                        "mean_rmse": float(np.mean(line_losses)),
                    }
                )

    b3_scores = [
        (float(np.mean(losses)), dimension)
        for dimension, losses in b3_losses.items()
        if losses
    ]
    b4_scores = [
        (float(np.mean(losses)), dimension, alpha)
        for (dimension, alpha), losses in b4_losses.items()
        if losses
    ]
    if not b3_scores or not b4_scores:
        raise RuntimeError("inner CV produced no candidate scores")
    _, best_b3_dimension = min(b3_scores, key=lambda item: (item[0], item[1]))
    _, best_b4_dimension, best_alpha = min(
        b4_scores, key=lambda item: (item[0], item[1], item[2])
    )
    return best_b3_dimension, best_b4_dimension, best_alpha, trace


def select_cclr_hyperparameters(
    train_control: NDArray[np.floating],
    train_response: NDArray[np.floating],
    inner_fold_ids: NDArray[np.integer],
    control_dimensions: list[int],
    response_ranks: list[int],
    ridge_alphas: list[float],
    max_variable_genes: int,
    min_mean_log1p_cpm: float,
    seed: int,
) -> tuple[int, int, float, list[dict[str, float | int | str]]]:
    """Select CCLR hyperparameters entirely inside the outer training fold."""
    control = np.asarray(train_control, dtype=np.float64)
    response = np.asarray(train_response, dtype=np.float64)
    fold_ids = np.asarray(inner_fold_ids, dtype=int)
    if len(control) != len(response) or len(control) != len(fold_ids):
        raise ValueError("inner-CV arrays have inconsistent row counts")
    unique_folds = sorted(np.unique(fold_ids).tolist())
    if len(unique_folds) < 2:
        raise ValueError("inner CV requires at least two folds")

    candidates = {
        (int(dimension), int(rank), float(alpha)): []
        for dimension in control_dimensions
        for rank in response_ranks
        for alpha in ridge_alphas
    }
    trace: list[dict[str, float | int | str]] = []
    for inner_fold in unique_folds:
        validation_mask = fold_ids == inner_fold
        inner_train_mask = ~validation_mask
        fold_seed = seed + 10_000 * int(inner_fold)
        max_control = min(max(control_dimensions), int(inner_train_mask.sum()) - 1)
        max_response = min(max(response_ranks), int(inner_train_mask.sum()) - 1)
        control_embedding = fit_control_embedding(
            control[inner_train_mask],
            max_components=max_control,
            max_variable_genes=max_variable_genes,
            min_mean_log1p_cpm=min_mean_log1p_cpm,
            seed=fold_seed,
        )
        response_embedding = fit_response_embedding(
            response[inner_train_mask],
            max_components=max_response,
            seed=fold_seed + 1,
        )
        feasible_dimensions = _candidate_dimensions(
            control_dimensions, control_embedding.pca.n_components_
        )
        feasible_ranks = _candidate_dimensions(
            response_ranks, response_embedding.pca.n_components_
        )
        train_control_scores = control_embedding.transform(control[inner_train_mask])
        validation_control_scores = control_embedding.transform(control[validation_mask])
        train_response_scores = response_embedding.transform(response[inner_train_mask])
        validation_response = response[validation_mask]

        for dimension in feasible_dimensions:
            for rank in feasible_ranks:
                for alpha in ridge_alphas:
                    ridge = Ridge(alpha=float(alpha), fit_intercept=True, solver="svd")
                    ridge.fit(
                        train_control_scores[:, :dimension],
                        train_response_scores[:, :rank],
                    )
                    predicted_scores = np.asarray(
                        ridge.predict(validation_control_scores[:, :dimension]),
                        dtype=np.float64,
                    )
                    if predicted_scores.ndim == 1:
                        predicted_scores = predicted_scores[:, None]
                    padded_scores = np.zeros(
                        (len(predicted_scores), response_embedding.pca.n_components_),
                        dtype=np.float64,
                    )
                    padded_scores[:, :rank] = predicted_scores
                    prediction = response_embedding.inverse_transform(padded_scores)
                    line_losses = [
                        rmse(observed, predicted)
                        for observed, predicted in zip(
                            validation_response, prediction, strict=True
                        )
                    ]
                    key = (dimension, rank, float(alpha))
                    candidates[key].extend(line_losses)
                    trace.append(
                        {
                            "model": "CCLR",
                            "inner_fold": int(inner_fold),
                            "control_dimension": dimension,
                            "response_rank": rank,
                            "alpha": float(alpha),
                            "mean_rmse": float(np.mean(line_losses)),
                            "n_validation_lines": int(validation_mask.sum()),
                        }
                    )

    scored = [
        (float(np.mean(losses)), rank, dimension, -alpha, dimension, rank, alpha)
        for (dimension, rank, alpha), losses in candidates.items()
        if losses
    ]
    if not scored:
        raise RuntimeError("inner CV produced no feasible CCLR candidate scores")
    _, _, _, _, best_dimension, best_rank, best_alpha = min(scored)
    return best_dimension, best_rank, best_alpha, trace


def fit_predict_cclr(
    train_control: NDArray[np.floating],
    train_response: NDArray[np.floating],
    test_control: NDArray[np.floating],
    inner_fold_ids: NDArray[np.integer],
    config: dict[str, Any],
    feature_config: dict[str, Any],
    seed: int,
) -> tuple[NDArray[np.float32], dict[str, Any], CCLRArtifact]:
    """Fit CCLR without accepting or inspecting outer-test responses."""
    control_train = np.asarray(train_control, dtype=np.float64)
    response_train = np.asarray(train_response, dtype=np.float64)
    control_test = np.asarray(test_control, dtype=np.float64)
    if len(control_train) != len(response_train):
        raise ValueError("training control and response rows must align")
    if control_train.shape[1] != control_test.shape[1]:
        raise ValueError("training and test control genes must align")

    control_dimensions = [int(value) for value in config["control_pca_dimensions"]]
    response_ranks = [int(value) for value in config["response_ranks"]]
    ridge_alphas = [float(value) for value in config["ridge_alphas"]]
    best_dimension, best_rank, best_alpha, trace = select_cclr_hyperparameters(
        control_train,
        response_train,
        inner_fold_ids,
        control_dimensions=control_dimensions,
        response_ranks=response_ranks,
        ridge_alphas=ridge_alphas,
        max_variable_genes=int(feature_config["max_variable_genes"]),
        min_mean_log1p_cpm=float(feature_config["min_mean_log1p_cpm"]),
        seed=seed,
    )
    control_embedding = fit_control_embedding(
        control_train,
        max_components=best_dimension,
        max_variable_genes=int(feature_config["max_variable_genes"]),
        min_mean_log1p_cpm=float(feature_config["min_mean_log1p_cpm"]),
        seed=seed,
    )
    response_embedding = fit_response_embedding(
        response_train,
        max_components=best_rank,
        seed=seed + 1,
    )
    train_control_scores = control_embedding.transform(control_train)[:, :best_dimension]
    train_response_scores = response_embedding.transform(response_train)[:, :best_rank]
    ridge = Ridge(alpha=best_alpha, fit_intercept=True, solver="svd")
    ridge.fit(train_control_scores, train_response_scores)
    artifact = CCLRArtifact(
        control_embedding=control_embedding,
        response_embedding=response_embedding,
        ridge=ridge,
        control_dimension=best_dimension,
        response_rank=best_rank,
        alpha=best_alpha,
    )
    prediction = artifact.predict(control_test)
    genes = int(response_train.shape[1])
    parameter_count = genes * (best_rank + 1) + best_rank * (best_dimension + 1)
    info = {
        "control_dimension": best_dimension,
        "response_rank": best_rank,
        "alpha": best_alpha,
        "selected_control_genes": int(len(control_embedding.gene_indices)),
        "parameter_count": int(parameter_count),
        "response_explained_variance_ratio": (
            response_embedding.pca.explained_variance_ratio_.tolist()
        ),
        "inner_cv_trace": trace,
    }
    return prediction, info, artifact


def fit_predict_baselines(
    train_control: NDArray[np.floating],
    train_response: NDArray[np.floating],
    train_lineages: NDArray[np.str_],
    test_control: NDArray[np.floating],
    test_lineages: NDArray[np.str_],
    inner_fold_ids: NDArray[np.integer],
    config: dict[str, Any],
    seed: int,
) -> tuple[dict[str, NDArray[np.float32]], dict[str, Any]]:
    """Fit B0-B4 without accepting or inspecting outer-test responses."""
    control_train = np.asarray(train_control)
    response_train = np.asarray(train_response, dtype=np.float32)
    control_test = np.asarray(test_control)
    lineage_train = np.asarray(train_lineages, dtype=str)
    lineage_test = np.asarray(test_lineages, dtype=str)
    if len(control_train) != len(response_train) or len(control_test) != len(lineage_test):
        raise ValueError("training or test metadata do not align to matrices")

    feature_config = config["feature_selection"]
    dimensions = [int(value) for value in config["control_pca_dimensions"]]
    alphas = [float(value) for value in config["ridge_alphas"]]
    best_b3_dim, best_b4_dim, best_alpha, trace = select_baseline_hyperparameters(
        control_train,
        response_train,
        inner_fold_ids,
        dimensions=dimensions,
        ridge_alphas=alphas,
        max_variable_genes=int(feature_config["max_variable_genes"]),
        min_mean_log1p_cpm=float(feature_config["min_mean_log1p_cpm"]),
        seed=seed,
    )
    max_dimension = max(best_b3_dim, best_b4_dim)
    embedding = fit_control_embedding(
        control_train,
        max_components=max_dimension,
        max_variable_genes=int(feature_config["max_variable_genes"]),
        min_mean_log1p_cpm=float(feature_config["min_mean_log1p_cpm"]),
        seed=seed,
    )
    train_scores = embedding.transform(control_train)
    test_scores = embedding.transform(control_test)

    mean_response = response_train.mean(axis=0, dtype=np.float64).astype(np.float32)
    predictions: dict[str, NDArray[np.float32]] = {
        "B0": np.zeros((len(control_test), response_train.shape[1]), dtype=np.float32),
        "B1": np.repeat(mean_response[None, :], len(control_test), axis=0),
    }
    minimum_lineage_lines = int(config["lineage_mean_min_train_lines"])
    b2_rows: list[NDArray[np.float32]] = []
    b2_fallback: list[bool] = []
    for lineage in lineage_test:
        same_lineage = lineage_train == lineage
        if int(same_lineage.sum()) >= minimum_lineage_lines:
            b2_rows.append(
                response_train[same_lineage].mean(axis=0, dtype=np.float64).astype(np.float32)
            )
            b2_fallback.append(False)
        else:
            b2_rows.append(mean_response)
            b2_fallback.append(True)
    predictions["B2"] = np.vstack(b2_rows)

    b3_prediction, neighbor_positions = nearest_neighbor_predictions(
        train_scores, response_train, test_scores, best_b3_dim
    )
    predictions["B3"] = b3_prediction
    ridge = Ridge(alpha=best_alpha, fit_intercept=True, solver="svd")
    ridge.fit(train_scores[:, :best_b4_dim], response_train)
    predictions["B4"] = ridge.predict(test_scores[:, :best_b4_dim]).astype(np.float32)

    info = {
        "b3_dimension": best_b3_dim,
        "b4_dimension": best_b4_dim,
        "b4_alpha": best_alpha,
        "selected_control_genes": int(len(embedding.gene_indices)),
        "b2_fallback": b2_fallback,
        "b3_neighbor_train_position": neighbor_positions.tolist(),
        "parameter_count": {
            "B0": 0,
            "B1": int(response_train.shape[1]),
            "B2": int(response_train.shape[1]),
            "B3": 0,
            "B4": int((best_b4_dim + 1) * response_train.shape[1]),
        },
        "inner_cv_trace": trace,
    }
    return predictions, info
