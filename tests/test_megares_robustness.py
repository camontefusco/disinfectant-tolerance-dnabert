import numpy as np

from scripts.evaluate_megares_similarity_robustness import brier_score, expected_calibration_error


def test_perfect_probabilities_have_zero_brier_and_calibration_error() -> None:
    labels = np.array([0, 1, 2])
    probabilities = np.eye(3)
    assert brier_score(labels, probabilities, 3) == 0
    assert expected_calibration_error(labels, probabilities) == 0
