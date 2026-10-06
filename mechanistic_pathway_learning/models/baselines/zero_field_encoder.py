"""An encoder whose field is zero everywhere: with a head's degree offset it makes the degree-only control.

Trained through the same harness, splits, loss and early stopping as the real encoders, a head on this field can use
nothing but its biases and the per-perturbation covariate (the standardised log degree of the perturbed nodes): the
sigmoid head becomes a per-symptom logistic regression on log degree, and the noisy-OR head reduces to its
degree-dependent leak. A real encoder with the same offset beats it only by what the graph adds beyond degree.
"""
from __future__ import annotations

import torch
from torch import Tensor, nn


class ZeroFieldEncoder(nn.Module):
    def __init__(self, num_graph_nodes: int, node_state_dim: int) -> None:
        super().__init__()
        self.num_graph_nodes = num_graph_nodes
        self.node_state_dim = node_state_dim
        self.register_buffer("device_anchor", torch.zeros(()), persistent=False)

    def forward(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:  # noqa: ARG002
        return torch.zeros(perturbation_node_index.shape[0], self.num_graph_nodes, self.node_state_dim, device=self.device_anchor.device)

    def perturbation_difference_field(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:
        return self.forward(perturbation_node_index, perturbation_sign_and_magnitude, relation_adjacencies)
