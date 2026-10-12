"""Baseline B3, input size: a model that sees how large a perturbation's input is and nothing else.

Why it exists. A drug enters the model with magnitude 1 at each of its mechanism targets, so a drug with three targets
carries a total input of 3 and the size of its input follows the number of targets curated for it (the user's decision
of 11 October 2026, docs/confirmatory_architecture_review.md question 6, and docs/label_source_review.md section 13.1).
The objection that answers is: a model could score well by reading input size rather than what the targets do. This
baseline is the control for that objection. It is given the two numbers that describe input size, the number of seeds
and the total absolute magnitude, and no graph, no node identity and no target identity. Whatever it scores is the
part of a model's score that input size alone explains.

Per symptom it fits a logistic regression on (1, number of seeds, total absolute magnitude) by Newton's method, with a
small fixed ridge on the two slopes. The ridge is there because the two features are the same number under the input
rule as built (total magnitude = number of targets for every drug but one), which leaves the design singular; keeping
both features and a ridge means the same baseline still reads correctly under --normalise-drug-input (total 1, count
varying) and under an occupancy input (both varying, not proportional). A gene perturbation has one seed of magnitude
1, so for genes this baseline is the popularity baseline; its reading is about the drugs.

Interface as the other baselines: fit(...) then predict(...) -> [num_perturbations, num_symptoms].
"""
from __future__ import annotations

import numpy as np

RIDGE = 1e-3  # on the two slopes only, never on the intercept, so an unused feature goes to zero instead of diverging
NEWTON_STEPS = 50
CONVERGED = 1e-10


def input_size_features(perturbation_seeds, perturbation_magnitudes) -> np.ndarray:
    """[num_perturbations, 2]: the number of seeds and the total absolute input magnitude of each perturbation."""
    counts = np.array([len(seeds) for seeds in perturbation_seeds], dtype=float)
    totals = np.array([float(np.abs(np.asarray(magnitudes, dtype=float)).sum()) for magnitudes in perturbation_magnitudes])
    return np.stack([counts, totals], axis=1)


def _fit_one_symptom(design: np.ndarray, outcomes: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted logistic regression by Newton's method; returns the coefficients of the design's columns."""
    coefficients = np.zeros(design.shape[1])
    penalty = np.diag([0.0] + [RIDGE] * (design.shape[1] - 1))
    for _ in range(NEWTON_STEPS):
        probabilities = 1.0 / (1.0 + np.exp(-np.clip(design @ coefficients, -30.0, 30.0)))
        gradient = design.T @ (weights * (outcomes - probabilities)) - penalty @ coefficients
        curvature = design.T @ (design * (weights * probabilities * (1.0 - probabilities))[:, None]) + penalty
        try:
            step = np.linalg.solve(curvature, gradient)
        except np.linalg.LinAlgError:
            break
        coefficients = coefficients + step
        if np.abs(step).max() < CONVERGED:
            break
    return coefficients


class InputSizeBaseline:
    def __init__(self) -> None:
        self.coefficients: np.ndarray | None = None
        self.feature_mean: np.ndarray | None = None
        self.feature_scale: np.ndarray | None = None

    def fit(self, training_features: np.ndarray, training_outcomes: np.ndarray,
            training_label_mask: np.ndarray | None = None) -> InputSizeBaseline:
        """training_features: [num_training_perturbations, 2] from input_size_features. training_outcomes in {0, 1}.

        With training_label_mask (a label selection) a pair the selection set aside carries weight 0, so it is neither
        a positive nor a negative, as in the trainer's loss and in the popularity baseline."""
        features = np.asarray(training_features, dtype=float)
        self.feature_mean = features.mean(axis=0)
        self.feature_scale = np.maximum(features.std(axis=0), 1e-8)  # standardised, so one fixed ridge suits any input rule
        design = np.concatenate([np.ones((len(features), 1)), (features - self.feature_mean) / self.feature_scale], axis=1)
        mask = np.ones_like(training_outcomes, dtype=float) if training_label_mask is None else np.asarray(training_label_mask, dtype=float)
        self.coefficients = np.stack([_fit_one_symptom(design, training_outcomes[:, symptom].astype(float), mask[:, symptom])
                                      for symptom in range(training_outcomes.shape[1])], axis=1)
        return self

    def predict(self, features: np.ndarray) -> np.ndarray:
        if self.coefficients is None:
            raise RuntimeError("fit before predict")
        features = np.asarray(features, dtype=float)
        design = np.concatenate([np.ones((len(features), 1)), (features - self.feature_mean) / self.feature_scale], axis=1)
        return 1.0 / (1.0 + np.exp(-np.clip(design @ self.coefficients, -30.0, 30.0)))
