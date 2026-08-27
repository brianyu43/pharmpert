import numpy as np

from yakseopdong.models import (
    fit_binary_covariate_encoder,
    fit_control_embedding,
    fit_lineage_encoder,
    fit_predict_baselines,
    fit_predict_cclr,
    fit_response_embedding,
)


def test_control_embedding_is_fit_only_from_training_matrix() -> None:
    rng = np.random.default_rng(3)
    train = rng.normal(size=(12, 40))
    first = fit_control_embedding(train, 3, 20, -10.0, seed=5)
    second = fit_control_embedding(train.copy(), 3, 20, -10.0, seed=5)
    np.testing.assert_array_equal(first.gene_indices, second.gene_indices)
    np.testing.assert_allclose(first.pca.components_, second.pca.components_)


def test_control_embedding_respects_predefined_candidate_panel() -> None:
    rng = np.random.default_rng(33)
    train = rng.normal(loc=2.0, size=(15, 40))
    candidates = np.asarray([2, 5, 8, 13, 21, 34])
    embedding = fit_control_embedding(
        train,
        max_components=3,
        max_variable_genes=5,
        min_mean_log1p_cpm=0.1,
        seed=8,
        candidate_gene_indices=candidates,
    )
    assert set(embedding.gene_indices).issubset(set(candidates))
    assert len(embedding.gene_indices) == 5


def test_training_fitted_metadata_encoders_handle_unknown_and_constant_columns() -> None:
    lineage = fit_lineage_encoder(np.asarray(["lung", "skin", "lung"]))
    encoded = lineage.transform(np.asarray(["skin", "unknown"]))
    assert encoded.shape == (2, 2)
    np.testing.assert_array_equal(encoded[1], np.zeros(2))

    binary = fit_binary_covariate_encoder(
        np.asarray([[0.0, 1.0], [1.0, 1.0], [0.0, 1.0]])
    )
    assert binary.kept_indices.tolist() == [0]
    assert binary.transform(np.asarray([[1.0, 1.0]])).shape == (1, 1)


def test_all_baselines_predict_without_outer_test_response() -> None:
    rng = np.random.default_rng(4)
    train_control = rng.normal(loc=2.0, scale=0.5, size=(16, 50))
    train_response = rng.normal(scale=0.2, size=(16, 50))
    test_control = rng.normal(loc=2.0, scale=0.5, size=(4, 50))
    train_lineages = np.asarray(["a"] * 8 + ["b"] * 6 + ["rare"] * 2)
    test_lineages = np.asarray(["a", "b", "unseen", "rare"])
    inner_folds = np.tile(np.arange(4), 4)
    config = {
        "lineage_mean_min_train_lines": 2,
        "feature_selection": {
            "max_variable_genes": 30,
            "min_mean_log1p_cpm": 0.1,
        },
        "control_pca_dimensions": [2, 3],
        "ridge_alphas": [0.1, 1.0],
    }
    predictions, info = fit_predict_baselines(
        train_control,
        train_response,
        train_lineages,
        test_control,
        test_lineages,
        inner_folds,
        config,
        seed=9,
    )
    assert set(predictions) == {"B0", "B1", "B2", "B3", "B4"}
    assert all(prediction.shape == (4, 50) for prediction in predictions.values())
    assert np.count_nonzero(predictions["B0"]) == 0
    assert info["b2_fallback"] == [False, False, True, False]


def test_response_embedding_reconstructs_training_mean_at_zero_score() -> None:
    rng = np.random.default_rng(12)
    train_response = rng.normal(size=(14, 30))
    embedding = fit_response_embedding(train_response, max_components=4, seed=7)
    reconstructed_mean = embedding.inverse_transform(np.zeros((1, 4)))[0]
    np.testing.assert_allclose(reconstructed_mean, train_response.mean(axis=0))


def test_cclr_predicts_without_outer_test_response_and_is_deterministic() -> None:
    rng = np.random.default_rng(22)
    train_control = rng.normal(loc=2.0, scale=0.5, size=(20, 60))
    latent = train_control[:, :3] @ rng.normal(size=(3, 2))
    response_basis = rng.normal(scale=0.2, size=(2, 60))
    train_response = latent @ response_basis + rng.normal(scale=0.01, size=(20, 60))
    test_control = rng.normal(loc=2.0, scale=0.5, size=(5, 60))
    inner_folds = np.tile(np.arange(4), 5)
    config = {
        "control_pca_dimensions": [2, 3],
        "response_ranks": [1, 2],
        "ridge_alphas": [0.1, 1.0],
    }
    feature_config = {
        "max_variable_genes": 40,
        "min_mean_log1p_cpm": 0.1,
    }
    first, first_info, first_artifact = fit_predict_cclr(
        train_control,
        train_response,
        test_control,
        inner_folds,
        config,
        feature_config,
        seed=31,
    )
    second, second_info, _ = fit_predict_cclr(
        train_control,
        train_response,
        test_control,
        inner_folds,
        config,
        feature_config,
        seed=31,
    )
    assert first.shape == (5, 60)
    assert first_info["response_rank"] in {1, 2}
    assert first_info["control_dimension"] in {2, 3}
    assert first_info == second_info
    np.testing.assert_allclose(first, second)
    np.testing.assert_allclose(first, first_artifact.predict(test_control))
