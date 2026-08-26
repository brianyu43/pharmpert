import numpy as np

from yakseopdong.models import fit_control_embedding, fit_predict_baselines


def test_control_embedding_is_fit_only_from_training_matrix() -> None:
    rng = np.random.default_rng(3)
    train = rng.normal(size=(12, 40))
    first = fit_control_embedding(train, 3, 20, -10.0, seed=5)
    second = fit_control_embedding(train.copy(), 3, 20, -10.0, seed=5)
    np.testing.assert_array_equal(first.gene_indices, second.gene_indices)
    np.testing.assert_allclose(first.pca.components_, second.pca.components_)


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
