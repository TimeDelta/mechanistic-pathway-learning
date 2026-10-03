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

Node states start from a learned embedding per node (identity; the default) or, when a fixed
feature matrix is given, from a linear map of structural node features (type, compartment, degree,
flags). The second form has no node identity at all: whatever it predicts for a held-out gene comes
from the structure around it, which is the inductive reading the mechanism claim needs and the
ablation that separates structure from memorized identity.

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
        node_features: Tensor | None = None,
    ) -> None:
        super().__init__()
        self.num_graph_nodes = num_graph_nodes
        self.num_relation_types = num_relation_types
        self.node_state_dim = node_state_dim
        self.num_message_passing_layers = num_message_passing_layers
        if node_features is None:
            self.base_node_state = nn.Embedding(num_graph_nodes, node_state_dim)
            self.register_buffer("node_features", None)
            self.feature_projection = None
        else:
            if node_features.shape[0] != num_graph_nodes:
                raise ValueError("node_features must have one row per graph node")
            self.base_node_state = None
            self.register_buffer("node_features", node_features.to(torch.float32))
            self.feature_projection = nn.Linear(node_features.shape[1], node_state_dim)
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
        base_state = self.base_node_state.weight if self.base_node_state is not None else self.feature_projection(self.node_features)
        node_state_field = base_state[None, :, :].expand(batch_size, -1, -1).clone()
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
        device = next(self.parameters()).device
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

    def stacked_relation_adjacency(self, relation_adjacencies: list[Tensor | None]) -> tuple[Tensor, Tensor]:
        """All present relation adjacencies stacked vertically into one sparse matrix of shape
        [num_present_relations * num_graph_nodes, num_graph_nodes], plus the present relation indices.

        One sparse product with the stacked matrix aggregates every relation at once; the result is
        viewed as [num_present_relations, num_graph_nodes, ...] with no copy. The stack is cached per
        adjacency list (keyed by the identity of its tensors) because the harness reuses one list.
        """
        cache_key = tuple(id(adjacency) for adjacency in relation_adjacencies)
        cached = getattr(self, "_stacked_adjacency_cache", None)
        if cached is not None and cached[0] == cache_key:
            return cached[1], cached[2]
        present_relations = [index for index, adjacency in enumerate(relation_adjacencies) if adjacency is not None]
        indices, values = [], []
        for stack_position, relation_index in enumerate(present_relations):
            adjacency = relation_adjacencies[relation_index].coalesce()
            shifted = adjacency.indices().clone()
            shifted[0] += stack_position * self.num_graph_nodes
            indices.append(shifted)
            values.append(adjacency.values())
        device = relation_adjacencies[present_relations[0]].device if present_relations else next(self.parameters()).device
        if present_relations:
            stacked = torch.sparse_coo_tensor(torch.cat(indices, dim=1), torch.cat(values), (len(present_relations) * self.num_graph_nodes, self.num_graph_nodes)).coalesce()
        else:
            stacked = torch.sparse_coo_tensor(torch.zeros((2, 0), dtype=torch.long, device=device), torch.zeros(0, device=device), (0, self.num_graph_nodes))
        present_index = torch.tensor(present_relations, dtype=torch.long, device=device)
        self._stacked_adjacency_cache = (cache_key, stacked, present_index)
        return stacked, present_index

    def propagate(self, node_state_field: Tensor, relation_adjacencies: list[Tensor | None]) -> Tensor:
        """L rounds of typed mean aggregation; returns the field as [batch_size, num_graph_nodes, node_state_dim].

        Per layer the new state is relu(X W_self + sum_r A_r X W_r + b). The states are kept node-major,
        [num_graph_nodes, batch_size, node_state_dim], so the sparse product reads them as a
        [num_graph_nodes, batch_size * node_state_dim] view without the permute-and-copy that a batch-major
        layout needs; all relations are aggregated by one sparse product with the stacked adjacency and
        transformed by one batched matrix product, using A_r (X W_r) = (A_r X) W_r. On the metabolic slice
        this halves the memory traffic of the batch-major, per-relation loop it replaces (profiling in
        docs/experiment_design.md section 8); the output is identical up to floating point rounding.
        """
        batch_size = node_state_field.shape[0]
        stacked_adjacency, present_relations = self.stacked_relation_adjacency(relation_adjacencies)
        num_present = int(present_relations.numel())
        node_major_state = node_state_field.permute(1, 0, 2).contiguous()  # [N, B, D]
        for layer_index in range(self.num_message_passing_layers):
            self_messages = node_major_state @ self.self_weight[layer_index]
            if num_present:
                flattened = node_major_state.reshape(self.num_graph_nodes, batch_size * self.node_state_dim)
                aggregated = torch.sparse.mm(stacked_adjacency.to(flattened.device), flattened)  # [R_present * N, B * D]
                aggregated = aggregated.view(num_present, self.num_graph_nodes * batch_size, self.node_state_dim)
                relation_messages = torch.bmm(aggregated, self.relation_weight[layer_index, present_relations]).sum(dim=0)
                relation_messages = relation_messages.view(self.num_graph_nodes, batch_size, self.node_state_dim)
                pre_activation = self_messages + relation_messages + self.layer_bias[layer_index]
            else:
                pre_activation = self_messages + self.layer_bias[layer_index]
            node_major_state = torch.relu(pre_activation)
        return node_major_state.permute(1, 0, 2).contiguous()
