import numpy as np

from scripts.evaluate_dnabert2_strategies import aggregate_embeddings, make_model


def test_aggregate_embeddings_caps_evenly_spaced_windows_per_isolate() -> None:
    isolate_ids = np.array(["a", "a", "a", "b", "b"])
    embeddings = np.array([[1.0], [3.0], [100.0], [4.0], [8.0]])
    frame = aggregate_embeddings(isolate_ids, embeddings, cap=2, aggregation="mean")
    assert frame.to_dict(orient="records") == [{"isolate_id": "a", 0: 50.5}, {"isolate_id": "b", 0: 6.0}]


def test_strategy_grid_supports_three_classifiers() -> None:
    assert make_model("logistic_regression", 42)["classifier"].random_state == 42
    assert make_model("linear_svm", 42)["classifier"].random_state == 42
    assert make_model("random_forest", 42).random_state == 42
