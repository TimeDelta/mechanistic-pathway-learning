"""The descriptors-only control: the perturbed node's own fixed features and nothing from the graph.

Rich node descriptors invite a shortcut that skips the mechanism: similar proteins, similar symptoms (Gillis and
Pavlidis, PLoS ONE 6, e17258, 2011, found that "high-quality gene function predictions can be made using data that
possesses no information on which gene interacts with which"). This encoder writes, at each perturbed node only,
the difference a one-layer network makes there, relu(W x_v + injection) - relu(W x_v), with x_v the node's features
(structural features plus descriptors), and leaves every other node at zero. Trained through the same harness and
head as the real encoders, it measures what the descriptors give on their own, so an encoder with propagation is
credited only with what the graph adds beyond them.
"""
from __future__ import annotations

import torch
from torch import Tensor, nn

PERTURBATION_FEATURE_DIM = 2  # sign and magnitude, as in RelationalMessagePassingEncoder


class LocalDescriptorEncoder(nn.Module):
    def __init__(self, num_graph_nodes: int, node_features: Tensor, node_state_dim: int) -> None:
        super().__init__()
        if node_features.shape[0] != num_graph_nodes:
            raise ValueError("node_features must have one row per graph node")
        self.num_graph_nodes = num_graph_nodes
        self.node_state_dim = node_state_dim
        self.register_buffer("node_features", node_features.to(torch.float32))
        self.feature_projection = nn.Linear(node_features.shape[1], node_state_dim)
        self.perturbation_injection = nn.Linear(PERTURBATION_FEATURE_DIM, node_state_dim)

    def forward(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:  # noqa: ARG002
        batch_size = perturbation_node_index.shape[0]
        valid = perturbation_node_index >= 0
        safe_index = perturbation_node_index.clamp_min(0)
        base_state = self.feature_projection(self.node_features[safe_index])  # [B, P, D]
        injected = self.perturbation_injection(perturbation_sign_and_magnitude)
        local_difference = (torch.relu(base_state + injected) - torch.relu(base_state)) * valid[:, :, None]
        field = torch.zeros(batch_size, self.num_graph_nodes, self.node_state_dim, device=local_difference.device, dtype=local_difference.dtype)
        return field.scatter_add(1, safe_index[:, :, None].expand(-1, -1, self.node_state_dim), local_difference)

    def perturbation_difference_field(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:
        return self.forward(perturbation_node_index, perturbation_sign_and_magnitude, relation_adjacencies)
