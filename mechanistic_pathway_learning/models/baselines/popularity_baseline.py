"""Baseline B0, popularity (design section 5.6).

Predicts each symptom at its training base rate, optionally scaled by the
perturbation's graph degree so hub bias has a baseline of its own. Any model
that cannot beat this is reading symptom frequency, not biology.

Every baseline exposes fit(training_observations, graph) and
predict(perturbations) -> array [num_perturbations, num_symptoms] so the
evaluation harness treats baselines and the proposed model identically.
"""
from __future__ import annotations

import numpy as np


class PopularityBaseline:
    def __init__(self, scale_by_degree: bool = False) -> None:
        self.scale_by_degree = scale_by_degree
        self.symptom_base_rate: np.ndarray | None = None

    def fit(self, training_outcomes: np.ndarray, training_perturbation_degrees: np.ndarray | None = None) -> "PopularityBaseline":
        """training_outcomes: [num_training_perturbations, num_symptoms] in {0, 1}."""
        self.symptom_base_rate = training_outcomes.mean(axis=0)
        return self

    def predict(self, num_perturbations: int, perturbation_degrees: np.ndarray | None = None) -> np.ndarray:
        if self.symptom_base_rate is None:
            raise RuntimeError("fit before predict")
        scores = np.tile(self.symptom_base_rate, (num_perturbations, 1))
        if self.scale_by_degree and perturbation_degrees is not None:
            scores = scores * np.log1p(perturbation_degrees)[:, None]
        return scores
