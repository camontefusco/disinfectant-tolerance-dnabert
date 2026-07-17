import numpy as np
import pandas as pd

from scripts.evaluate_embedding_repeated_grouped_splits import best_threshold


def test_best_threshold_uses_balanced_accuracy() -> None:
    labels = pd.Series([0, 0, 1, 1])
    probabilities = np.array([0.1, 0.4, 0.6, 0.9])
    assert best_threshold(labels, probabilities) == 0.6

