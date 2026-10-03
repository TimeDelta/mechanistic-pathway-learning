"""Baseline B1, network proximity by random walk with restart (design section 5.6).

Score(perturbation, symptom) = stationary random-walk mass that starts at the
perturbation's nodes and lands on the symptom's anchor set, where the anchor set
is the union of the training perturbations' nodes for that symptom. This is the
network-propagation baseline of Cowen et al. (2017): no learning beyond the
restart probability, so it measures how much of the task is graph proximity.

The walk runs on the undirected graph with currency metabolites removed. Sources
are propagated in batches to keep memory bounded on a 34,000-node graph.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sparse


def build_normalized_adjacency(num_nodes: int, source_index: np.ndarray, target_index: np.ndarray, excluded_nodes: np.ndarray | None = None) -> sparse.csr_matrix:
    """Column-normalized undirected adjacency W so that W @ p redistributes probability mass."""
    keep = np.ones(len(source_index), dtype=bool)
    if excluded_nodes is not None and len(excluded_nodes):
        excluded = np.zeros(num_nodes, dtype=bool)
        excluded[excluded_nodes] = True
        keep = ~(excluded[source_index] | excluded[target_index])
    rows = np.concatenate([source_index[keep], target_index[keep]])
    cols = np.concatenate([target_index[keep], source_index[keep]])
    adjacency = sparse.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(num_nodes, num_nodes)).tocsr()
    adjacency.data[:] = 1.0
    column_sums = np.asarray(adjacency.sum(axis=0)).ravel()
    column_sums[column_sums == 0] = 1.0
    return adjacency @ sparse.diags(1.0 / column_sums)


class RandomWalkWithRestartBaseline:
    def __init__(self, restart_probability: float = 0.3, num_iterations: int = 50, batch_size: int = 256) -> None:
        self.restart_probability = restart_probability
        self.num_iterations = num_iterations
        self.batch_size = batch_size
        self.normalized_adjacency: sparse.csr_matrix | None = None
        self.anchor_mask_by_symptom: np.ndarray | None = None  # [num_symptoms, num_nodes]

    def stationary_distributions(self, seed_matrix: np.ndarray) -> np.ndarray:
        """seed_matrix: [num_nodes, num_sources] column-stochastic restart vectors -> stationary [num_nodes, num_sources]."""
        state = seed_matrix.copy()
        for _ in range(self.num_iterations):
            state = (1.0 - self.restart_probability) * (self.normalized_adjacency @ state) + self.restart_probability * seed_matrix
        return state

    def fit(self, normalized_adjacency: sparse.csr_matrix, training_perturbation_seeds: list[np.ndarray], training_outcomes: np.ndarray) -> "RandomWalkWithRestartBaseline":
        """training_perturbation_seeds: per training perturbation, the node indices it perturbs."""
        self.normalized_adjacency = normalized_adjacency
        num_nodes = normalized_adjacency.shape[0]
        num_symptoms = training_outcomes.shape[1]
        anchor_mask = np.zeros((num_symptoms, num_nodes))
        for seeds, outcomes in zip(training_perturbation_seeds, training_outcomes):
            for symptom_index in np.where(outcomes > 0)[0]:
                anchor_mask[symptom_index, seeds] = 1.0
        anchor_sizes = anchor_mask.sum(axis=1, keepdims=True)
        anchor_sizes[anchor_sizes == 0] = 1.0
        self.anchor_mask_by_symptom = anchor_mask / anchor_sizes
        return self

    def predict(self, perturbation_seeds: list[np.ndarray]) -> np.ndarray:
        if self.normalized_adjacency is None or self.anchor_mask_by_symptom is None:
            raise RuntimeError("fit before predict")
        num_nodes = self.normalized_adjacency.shape[0]
        scores = np.zeros((len(perturbation_seeds), self.anchor_mask_by_symptom.shape[0]))
        for start in range(0, len(perturbation_seeds), self.batch_size):
            batch = perturbation_seeds[start : start + self.batch_size]
            seed_matrix = np.zeros((num_nodes, len(batch)))
            for column, seeds in enumerate(batch):
                if len(seeds):
                    seed_matrix[seeds, column] = 1.0 / len(seeds)
            stationary = self.stationary_distributions(seed_matrix)  # [num_nodes, batch]
            scores[start : start + len(batch)] = (self.anchor_mask_by_symptom @ stationary).T
        return scores
