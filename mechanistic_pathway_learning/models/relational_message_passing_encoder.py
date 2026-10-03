"""Relational message passing encoder (experiment design, section 5.2, route 1).

A perturbation is written onto its nodes as a (sign, magnitude) injection added
to a learned base state for every graph node. L rounds of typed message passing
over the hard edges of the physiology graph produce a node-state field of shape
[batch_size, num_graph_nodes, node_state_dim], which the noisy-OR pathway module
head reads locally. The encoder uses only torch so the scaffold runs without
PyTorch Geometric; swapping in a library implementation later only has to keep
the forward signature.

Edges are shared across the batch and given once as edge_index [2, num_edges]
(source row, destination row) with edge_relation_type [num_edges]. Messages are
mean-aggregated per relation type at the destination node. Compartment changes
are only possible through transport edges because the graph contains no other
cross-compartment edges (design section 4.1, hard constraint 4).
"""
from __future__ import annotations

import math

import torch
from torch import Tensor, nn

PERTURBATION_FEATURE_DIM = 2  # (sign, magnitude)


class RelationalMessagePassingEncoder(nn.Module):
    def __init__(
        self,
        num_graph_nodes: int,
        num_relation_types: int,
        node_state_dim: int,
        num_message_passing_layers: int = 3,
    ) -> None:
        super().__init__()
        self.num_graph_nodes = num_graph_nodes
        self.num_relation_types = num_relation_types
        self.node_state_dim = node_state_dim
        self.num_message_passing_layers = num_message_passing_layers
        self.base_node_state = nn.Embedding(num_graph_nodes, node_state_dim)
        self.perturbation_injection = nn.Linear(PERTURBATION_FEATURE_DIM, node_state_dim)
        scale = 1.0 / math.sqrt(node_state_dim)
        self.relation_weight = nn.Parameter(
            torch.randn(num_message_passing_layers, num_relation_types, node_state_dim, node_state_dim) * scale
        )
        self.self_weight = nn.Parameter(torch.randn(num_message_passing_layers, node_state_dim, node_state_dim) * scale)
        self.layer_bias = nn.Parameter(torch.zeros(num_message_passing_layers, node_state_dim))

    def initial_node_state_field(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor) -> Tensor:
        """Base state for every node plus the injected perturbation at the perturbed nodes.

        perturbation_node_index: [batch_size, max_perturbed_nodes], padded with -1.
        perturbation_sign_and_magnitude: [batch_size, max_perturbed_nodes, 2].
        """
        batch_size = perturbation_node_index.shape[0]
        node_state_field = self.base_node_state.weight[None, :, :].expand(batch_size, -1, -1).clone()
        valid_mask = perturbation_node_index >= 0
        injected_state = self.perturbation_injection(perturbation_sign_and_magnitude) * valid_mask[:, :, None]
        safe_index = perturbation_node_index.clamp_min(0)
        node_state_field = node_state_field.scatter_add(
            dim=1,
            index=safe_index[:, :, None].expand(-1, -1, self.node_state_dim),
            src=injected_state,
        )
        return node_state_field

    def forward(
        self,
        perturbation_node_index: Tensor,
        perturbation_sign_and_magnitude: Tensor,
        edge_index: Tensor,
        edge_relation_type: Tensor,
    ) -> Tensor:
        node_state_field = self.initial_node_state_field(perturbation_node_index, perturbation_sign_and_magnitude)
        source_index, destination_index = edge_index[0], edge_index[1]
        for layer_index in range(self.num_message_passing_layers):
            aggregated_messages = torch.zeros_like(node_state_field)
            for relation_index in range(self.num_relation_types):
                edges_of_relation = edge_relation_type == relation_index
                if not torch.any(edges_of_relation):
                    continue
                relation_source = source_index[edges_of_relation]
                relation_destination = destination_index[edges_of_relation]
                source_states = node_state_field[:, relation_source, :]
                messages = source_states @ self.relation_weight[layer_index, relation_index]
                relation_sum = torch.zeros_like(node_state_field).index_add(1, relation_destination, messages)
                in_degree = torch.zeros(self.num_graph_nodes, device=node_state_field.device).index_add(
                    0, relation_destination, torch.ones_like(relation_destination, dtype=node_state_field.dtype)
                )
                aggregated_messages = aggregated_messages + relation_sum / in_degree.clamp_min(1.0)[None, :, None]
            self_messages = node_state_field @ self.self_weight[layer_index]
            node_state_field = torch.relu(self_messages + aggregated_messages + self.layer_bias[layer_index])
        return node_state_field
