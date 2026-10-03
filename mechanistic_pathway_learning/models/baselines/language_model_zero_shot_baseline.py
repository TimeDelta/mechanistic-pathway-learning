"""Baseline: language_model_zero_shot_baseline (design section 5.6). Not implemented in version 0.1.

Every baseline exposes fit(training_observations, graph) and
predict(perturbations) -> array [num_perturbations, num_symptoms] so the
evaluation harness treats baselines and the proposed model identically.
"""


class LanguageModelZeroShotBaseline:
    def fit(self, training_observations, graph):
        raise NotImplementedError

    def predict(self, perturbations):
        raise NotImplementedError
