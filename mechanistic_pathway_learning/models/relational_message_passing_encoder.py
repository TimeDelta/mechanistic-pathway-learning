"""Relational message passing encoder (experiment design, section 5.2, route 1).

A perturbation is written onto its nodes as a (sign, magnitude) injection added
to a learned base state for every graph node. L rounds of typed message passing
over the hard edges of the physiology graph produce a node-state field of shape
[batch_size, num_graph_nodes, node_state_dim], which the noisy-OR pathway module
head reads locally. The encoder uses only torch so the scaffold runs without
PyTorch Geometric; swapping in a library implementation later only has to keep
the forward signature.

Two readings of the field are available. The absolute field carries the learned base state of
every node plus the propagated perturbation, so a readout can memorize node identity. The
difference field, forward(perturbed) minus forward(unperturbed), carries only what the
perturbation changed; a module that reads it responds to the propagated perturbation and not to
which nodes exist, which is the mechanism reading the design asks for (section 5.2) and the
default of the training harness. The unperturbed field is one extra batch-of-one pass.

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

    @staticmethod
    def build_relation_adjacencies(edge_index: Tensor, edge_relation_type: Tensor, num_graph_nodes: int, num_relation_types: int) -> list[Tensor | None]:
        """One sparse, in-degree-normalized adjacency per relation: entry [destination, source] = 1 / in_degree(destination)."""
        adjacencies: list[Tensor | None] = []
        for relation_index in range(num_relation_types):
            edges_of_relation = edge_relation_type == relation_index
            if not torch.any(edges_of_relation):
                adjacencies.append(None)
                continue
            source = edge_index[0][edges_of_relation]
            destination = edge_index[1][edges_of_relation]
            in_degree = torch.zeros(num_graph_nodes, dtype=torch.float32).index_add(0, destination, torch.ones(len(destination)))
            values = 1.0 / in_degree[destination]
            adjacency = torch.sparse_coo_tensor(torch.stack([destination, source]), values, (num_graph_nodes, num_graph_nodes)).coalesce()
            adjacencies.append(adjacency)
        return adjacencies

    def forward(
        self,
        perturbation_node_index: Tensor,
        perturbation_sign_and_magnitude: Tensor,
        edge_index: Tensor | None = None,
        edge_relation_type: Tensor | None = None,
        relation_adjacencies: list[Tensor | None] | None = None,
    ) -> Tensor:
        """Either pass edge arrays (adjacencies are built and cached) or precomputed relation_adjacencies.

        Aggregation uses one sparse matrix product per relation and layer over the whole batch
        (states reshaped to [num_graph_nodes, batch_size * node_state_dim]), which avoids the
        per-edge intermediate of size batch x edges x dim.
        """
        if relation_adjacencies is None:
            if edge_index is None or edge_relation_type is None:
                raise ValueError("pass edge_index and edge_relation_type or relation_adjacencies")
            cache_key = (int(edge_index.shape[1]), int(edge_relation_type.sum()))
            if getattr(self, "_adjacency_cache_key", None) != cache_key:
                self._cached_adjacencies = self.build_relation_adjacencies(edge_index, edge_relation_type, self.num_graph_nodes, self.num_relation_types)
                self._adjacency_cache_key = cache_key
            relation_adjacencies = self._cached_adjacencies
        node_state_field = self.initial_node_state_field(perturbation_node_index, perturbation_sign_and_magnitude)
        return self.propagate(node_state_field, relation_adjacencies)

    def unperturbed_node_state_field(self, relation_adjacencies: list[Tensor | None]) -> Tensor:
        """Field of shape [1, num_graph_nodes, node_state_dim] with no perturbation written on it."""
        device = self.base_node_state.weight.device
        empty_index = torch.full((1, 1), -1, dtype=torch.long, device=device)
        empty_injection = torch.zeros((1, 1, PERTURBATION_FEATURE_DIM), device=device)
        return self.forward(empty_index, empty_injection, relation_adjacencies=relation_adjacencies)

    def perturbation_difference_field(
        self,
        perturbation_node_index: Tensor,
        perturbation_sign_and_magnitude: Tensor,
        relation_adjacencies: list[Tensor | None],
    ) -> Tensor:
        """forward(perturbed) - forward(unperturbed): the propagated change caused by the perturbation."""
        perturbed = self.forward(perturbation_node_index, perturbation_sign_and_magnitude, relation_adjacencies=relation_adjacencies)
        return perturbed - self.unperturbed_node_state_field(relation_adjacencies)

    def propagate(self, node_state_field: Tensor, relation_adjacencies: list[Tensor | None]) -> Tensor:
        batch_size = node_state_field.shape[0]
        for layer_index in range(self.num_message_passing_layers):
            aggregated_messages = torch.zeros_like(node_state_field)
            for relation_index, adjacency in enumerate(relation_adjacencies):
                if adjacency is None:
                    continue
                messages = node_state_field @ self.relation_weight[layer_index, relation_index]  # [B, N, D]
                flattened = messages.permute(1, 0, 2).reshape(self.num_graph_nodes, batch_size * self.node_state_dim)
                aggregated = torch.sparse.mm(adjacency.to(flattened.device), flattened)
                aggregated_messages = aggregated_messages + aggregated.reshape(self.num_graph_nodes, batch_size, self.node_state_dim).permute(1, 0, 2)
            self_messages = node_state_field @ self.self_weight[layer_index]
            node_state_field = torch.relu(self_messages + aggregated_messages + self.layer_bias[layer_index])
        return node_state_field
