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

With node descriptors among the features, descriptor_treatment (mechanistic_pathway_learning/models/
descriptor_treatments.py) decides how they enter: plain (one map of all the columns), seed_masked (each perturbation's
perturbed nodes start from their structural columns only, in the perturbed pass and in its own unperturbed reference,
so the reference is one pass per perturbation instead of one per batch and a forward costs about twice as much), or
zero_init_slow (a separate descriptor map that starts at zero and trains at its own learning rate).

Edges are shared across the batch and given once as edge_index [2, num_edges]
(source row, destination row) with edge_relation_type [num_edges]. Messages are
mean-aggregated per relation type at the destination node. Compartment changes
are only possible through transport edges because the graph contains no other
cross-compartment edges (design section 4.1, hard constraint 4).

One relation is not a mean. A protein complex is present only to the extent that every subunit is
present, so the members of a complex combine by a conjunction rather than an average: a mean lets a
plentiful subunit stand in for a missing one, which is the opposite of what an obligate complex does.
conjunction_relations (by default none, for the message-passing encoder "member_of") are therefore
aggregated by an element-wise minimum over the member nodes, one minimum per node-state channel, so
each preserved quantity is set by its own scarcest subunit and no member has to be summarised into a
single number first (the user, 9 October 2026: "i thought it would be multiple mins (one for each
preserved quantity)"). conjunction_aggregation picks the form: "minimum" is the hard minimum,
"soft_minimum" is -T log sum exp(-x/T) over the members, which approaches the hard minimum as the
temperature falls and is the control that separates "a conjunction helps" from "a hard minimum
helps", and "mean" is the in-degree mean every other relation uses.

Entry nodes are not always there. A drug node (mechanistic_pathway_learning/graph/drug_entry_nodes.py)
exists only in the perturbation that takes the drug, although the graph table holds a node for
every drug of the study. Given entry_node_mask, an entry node that is not among a perturbation's
seeds is held at zero in that perturbation's field, before the first layer and after every layer,
and in the unperturbed reference every entry node is absent. Its edges to the rest of the graph
(entry_edge_index, the mechanism edges of a drug) are added as a sum rather than a mean:

    message into target = weight * (x_source W_unsigned) + sign * weight * (x_source W_signed)

with one pair of maps per layer, which is the edge form of the injection Linear([sign, magnitude])
that a seed on the target itself receives. A mean would divide a drug's signal at its target by the
number of study drugs sharing that target, and with every drug node present a gene knockout's field
would depend on which drugs the study contains; both are measured in
docs/drug_entry_nodes_measured.md. The relations of entry_relation_indices are left out of the
stacked product because these edges are handled here.

An entry node sits one edge behind the nodes a seed on the targets starts from, so with L layers its
signal reaches L - 1 edges beyond the targets instead of L. On the confirmatory graph at three
layers that leaves a drug a median 0.11 of the nodes its targets reach
(docs/drug_entry_nodes_measured.md). entry_before_first_layer delivers the entry edges once before
the first layer as well, with maps of their own, so the targets start the first layer already
carrying the drug's signal and the depth beyond them is the L of a seed on the targets. The entry
edges still send at every layer, which makes a drug a sustained input where a seed is only an
initial state.
"""
from __future__ import annotations

import math

import torch
from torch import Tensor, nn

from mechanistic_pathway_learning.models.descriptor_treatments import check_descriptor_treatment, features_without_descriptors, seed_unit_mask

PERTURBATION_FEATURE_DIM = 2  # (sign, magnitude)
CONJUNCTION_AGGREGATIONS = ("mean", "minimum", "soft_minimum")
DEFAULT_CONJUNCTION_RELATIONS = ("member_of",)  # the obligate-complex relation of the physiology graph
DEFAULT_SOFT_MINIMUM_TEMPERATURE = 0.1


class RelationalMessagePassingEncoder(nn.Module):
    def __init__(
        self,
        num_graph_nodes: int,
        num_relation_types: int,
        node_state_dim: int,
        num_message_passing_layers: int = 3,
        node_features: Tensor | None = None,
        num_descriptor_columns: int = 0,
        descriptor_treatment: str = "plain",
        seed_mask_partner_index: Tensor | None = None,
        conjunction_relation_indices: tuple[int, ...] = (),
        conjunction_aggregation: str = "minimum",
        soft_minimum_temperature: float = DEFAULT_SOFT_MINIMUM_TEMPERATURE,
        entry_node_mask: Tensor | None = None,
        entry_edge_index: Tensor | None = None,
        entry_edge_sign: Tensor | None = None,
        entry_edge_weight: Tensor | None = None,
        entry_relation_indices: tuple[int, ...] = (),
        entry_before_first_layer: bool = False,
    ) -> None:
        super().__init__()
        if conjunction_aggregation not in CONJUNCTION_AGGREGATIONS:
            raise ValueError(f"conjunction_aggregation must be one of {CONJUNCTION_AGGREGATIONS}")
        if soft_minimum_temperature <= 0:
            raise ValueError("soft_minimum_temperature must be positive")
        if entry_node_mask is None and (entry_edge_index is not None or entry_relation_indices or entry_before_first_layer):
            raise ValueError("entry edges need entry_node_mask, the nodes that exist only where they are seeded")
        if entry_node_mask is not None:
            if entry_node_mask.shape != (num_graph_nodes,):
                raise ValueError("entry_node_mask must have one entry per graph node")
            # under seed_masked the trainer passes each entry node's targets as its mask partners (seed_mask_partner_pairs),
            # so a drug seeded on its node hides the descriptors it hides when it is seeded on its targets
            if entry_edge_index is None:
                entry_edge_index = torch.zeros((2, 0), dtype=torch.long)
            num_entry_edges = entry_edge_index.shape[1]
            entry_edge_sign = torch.ones(num_entry_edges) if entry_edge_sign is None else entry_edge_sign.to(torch.float32)
            entry_edge_weight = torch.ones(num_entry_edges) if entry_edge_weight is None else entry_edge_weight.to(torch.float32)
            if entry_edge_sign.shape != (num_entry_edges,) or entry_edge_weight.shape != (num_entry_edges,):
                raise ValueError("entry_edge_sign and entry_edge_weight must have one entry per entry edge")
            if num_entry_edges and not bool(entry_node_mask.bool()[entry_edge_index[0]].all()):
                raise ValueError("every entry edge must leave an entry node")
        self.register_buffer("entry_node_mask", None if entry_node_mask is None else entry_node_mask.bool(), persistent=False)
        self.register_buffer("entry_edge_source", None if entry_node_mask is None else entry_edge_index[0].long(), persistent=False)
        self.register_buffer("entry_edge_target", None if entry_node_mask is None else entry_edge_index[1].long(), persistent=False)
        self.register_buffer("entry_edge_weight", None if entry_node_mask is None else entry_edge_weight, persistent=False)
        self.register_buffer("entry_edge_signed_weight", None if entry_node_mask is None else entry_edge_sign * entry_edge_weight, persistent=False)
        self.entry_relation_indices = tuple(sorted(set(entry_relation_indices)))
        # with "mean" the conjunction relations are left in the stacked product, so the encoder is the one it was
        # before this option existed and the member edges are averaged like any other relation
        self.conjunction_relation_indices = () if conjunction_aggregation == "mean" else tuple(sorted(set(conjunction_relation_indices)))
        self.conjunction_aggregation = conjunction_aggregation
        self.soft_minimum_temperature = soft_minimum_temperature
        self.register_buffer("seed_mask_partner_index", None if seed_mask_partner_index is None else seed_mask_partner_index.long(), persistent=False)
        self.num_graph_nodes = num_graph_nodes
        self.num_relation_types = num_relation_types
        self.node_state_dim = node_state_dim
        self.num_message_passing_layers = num_message_passing_layers
        self.descriptor_treatment = descriptor_treatment
        self.num_descriptor_columns = num_descriptor_columns
        self.descriptor_projection = None
        if node_features is None:
            if descriptor_treatment != "plain":
                raise ValueError("a descriptor treatment needs node features")
            self.base_node_state = nn.Embedding(num_graph_nodes, node_state_dim)
            self.register_buffer("node_features", None)
            self.feature_projection = None
        else:
            if node_features.shape[0] != num_graph_nodes:
                raise ValueError("node_features must have one row per graph node")
            check_descriptor_treatment(descriptor_treatment, num_descriptor_columns, node_features.shape[1])
            self.base_node_state = None
            self.register_buffer("node_features", node_features.to(torch.float32))
            if descriptor_treatment == "zero_init_slow":
                self.num_structural_columns = node_features.shape[1] - num_descriptor_columns
                self.feature_projection = nn.Linear(self.num_structural_columns, node_state_dim)
                self.descriptor_projection = nn.Linear(num_descriptor_columns, node_state_dim, bias=False)
                nn.init.zeros_(self.descriptor_projection.weight)
            else:
                self.feature_projection = nn.Linear(node_features.shape[1], node_state_dim)
            if descriptor_treatment == "seed_masked":
                self.register_buffer("node_features_without_descriptors", features_without_descriptors(self.node_features, num_descriptor_columns), persistent=False)
        self.perturbation_injection = nn.Linear(PERTURBATION_FEATURE_DIM, node_state_dim)
        scale = 1.0 / math.sqrt(node_state_dim)
        self.relation_weight = nn.Parameter(
            torch.randn(num_message_passing_layers, num_relation_types, node_state_dim, node_state_dim) * scale
        )
        self.self_weight = nn.Parameter(torch.randn(num_message_passing_layers, node_state_dim, node_state_dim) * scale)
        self.layer_bias = nn.Parameter(torch.zeros(num_message_passing_layers, node_state_dim))
        # created last and only with entry nodes, so an encoder without them draws the parameters it always drew
        self.entry_unsigned_weight, self.entry_signed_weight = None, None
        self.entry_first_unsigned_weight, self.entry_first_signed_weight = None, None
        if entry_node_mask is not None:
            self.entry_unsigned_weight = nn.Parameter(torch.randn(num_message_passing_layers, node_state_dim, node_state_dim) * scale)
            self.entry_signed_weight = nn.Parameter(torch.randn(num_message_passing_layers, node_state_dim, node_state_dim) * scale)
            if entry_before_first_layer:
                self.entry_first_unsigned_weight = nn.Parameter(torch.randn(node_state_dim, node_state_dim) * scale)
                self.entry_first_signed_weight = nn.Parameter(torch.randn(node_state_dim, node_state_dim) * scale)

    def descriptor_parameters(self) -> list[nn.Parameter]:
        """The parameters that read the descriptor columns on their own (zero_init_slow), for their own learning rate."""
        return [] if self.descriptor_projection is None else list(self.descriptor_projection.parameters())

    def base_state(self, without_descriptors: bool = False) -> Tensor:
        """[num_graph_nodes, node_state_dim]: the learned embedding, or the map of the node features (with the descriptor
        columns set to zero when without_descriptors, under seed_masked)."""
        if self.base_node_state is not None:
            return self.base_node_state.weight
        if self.descriptor_projection is not None:
            return (self.feature_projection(self.node_features[:, : self.num_structural_columns])
                    + self.descriptor_projection(self.node_features[:, self.num_structural_columns :]))
        return self.feature_projection(self.node_features_without_descriptors if without_descriptors else self.node_features)

    def initial_node_state_field(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, inject: bool = True) -> Tensor:
        """Base state for every node plus the injected perturbation at the perturbed nodes.

        perturbation_node_index: [batch_size, max_perturbed_nodes], padded with -1.
        perturbation_sign_and_magnitude: [batch_size, max_perturbed_nodes, 2].
        Under seed_masked the perturbed nodes take their base state without descriptors; inject=False gives that field
        with nothing injected, the unperturbed reference of each perturbation.
        """
        batch_size = perturbation_node_index.shape[0]
        base_state = self.base_state()
        node_state_field = base_state[None, :, :].expand(batch_size, -1, -1).clone()
        if self.descriptor_treatment == "seed_masked":
            seed_rows = seed_unit_mask(perturbation_node_index, self.num_graph_nodes, self.seed_mask_partner_index)
            node_state_field = torch.where(seed_rows[:, :, None], self.base_state(without_descriptors=True)[None, :, :], node_state_field)
        if not inject:
            return node_state_field
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
        return self.propagate(node_state_field, relation_adjacencies, self.present_nodes(perturbation_node_index))

    def present_nodes(self, perturbation_node_index: Tensor) -> Tensor | None:
        """[batch_size, num_graph_nodes], False for an entry node the perturbation is not seeded on; None without entry
        nodes. An entry node exists only where it is taken, so a padded index (-1) and every unseeded row leave it out,
        and the batch-of-one unperturbed reference, whose index is all padding, has no entry node at all."""
        if self.entry_node_mask is None:
            return None
        batch_size = perturbation_node_index.shape[0]
        seeded = torch.zeros((batch_size, self.num_graph_nodes + 1), dtype=torch.bool, device=perturbation_node_index.device)
        # padding goes to an extra last column that is then dropped, so index -1 marks no node
        seeded.scatter_(1, torch.where(perturbation_node_index >= 0, perturbation_node_index, torch.full_like(perturbation_node_index, self.num_graph_nodes)), True)
        return ~self.entry_node_mask[None, :] | seeded[:, : self.num_graph_nodes]

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
        """forward(perturbed) - forward(unperturbed): the propagated change caused by the perturbation. Under seed_masked the
        unperturbed reference differs per perturbation (its seeds lack descriptors in both passes), so the perturbed and
        reference fields propagate together as one batch of twice the size."""
        if self.descriptor_treatment == "seed_masked":
            batch_size = perturbation_node_index.shape[0]
            both = torch.cat([self.initial_node_state_field(perturbation_node_index, perturbation_sign_and_magnitude),
                              self.initial_node_state_field(perturbation_node_index, perturbation_sign_and_magnitude, inject=False)], dim=0)
            present = self.present_nodes(perturbation_node_index)
            if present is not None:  # the reference of a perturbation is the body without the drug: no entry node in it
                present = torch.cat([present, (~self.entry_node_mask)[None, :].expand(batch_size, -1)], dim=0)
            propagated = self.propagate(both, relation_adjacencies, present)
            return propagated[:batch_size] - propagated[batch_size:]
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
        present_relations = [index for index, adjacency in enumerate(relation_adjacencies)
                             if adjacency is not None and index not in self.conjunction_relation_indices and index not in self.entry_relation_indices]
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
        self._compressed_adjacency_cache = None
        return stacked, present_index

    def compressed_relation_adjacency(self, relation_adjacencies: list[Tensor | None]) -> tuple[Tensor, Tensor, list[tuple[int, int, int]]]:
        """The stacked rows that hold an edge, as their own sparse matrix [K, N], with the target node of each row and,
        per present relation, (relation index, first row, end row): rows stay in stack order, so a relation's rows are
        contiguous. On graph_full_neuronal K is 54,982 of the 700,435 stacked rows, so the product and the per-relation
        transform touch about a thirteenth of the memory."""
        stacked, present_index = self.stacked_relation_adjacency(relation_adjacencies)
        cached = getattr(self, "_compressed_adjacency_cache", None)
        if cached is not None:
            return cached
        indices = stacked.indices()
        nonempty_rows, compressed_row = torch.unique(indices[0], return_inverse=True)
        compressed = torch.sparse_coo_tensor(torch.stack([compressed_row, indices[1]]), stacked.values(), (len(nonempty_rows), self.num_graph_nodes)).coalesce()
        stack_position = torch.div(nonempty_rows, self.num_graph_nodes, rounding_mode="floor")
        segments = []
        for position, relation_index in enumerate(present_index.tolist()):
            rows_of_relation = torch.nonzero(stack_position == position).flatten()
            if len(rows_of_relation):
                segments.append((relation_index, int(rows_of_relation[0]), int(rows_of_relation[-1]) + 1))
        self._compressed_adjacency_cache = (compressed, nonempty_rows % self.num_graph_nodes, segments)
        return self._compressed_adjacency_cache

    def conjunction_member_edges(self, relation_adjacencies: list[Tensor | None]) -> list[tuple[int, Tensor, Tensor]]:
        """Per conjunction relation, (relation index, member source rows, the node each member feeds).

        Read off the relation's own adjacency, whose indices are [destination, source]; the in-degree values are
        dropped, because a minimum over the members is not an average over them and must not be scaled by how many
        there are. Cached per adjacency list, as the stacked product is.
        """
        cache_key = tuple(id(adjacency) for adjacency in relation_adjacencies)
        cached = getattr(self, "_conjunction_edge_cache", None)
        if cached is not None and cached[0] == cache_key:
            return cached[1]
        member_edges = []
        for relation_index in self.conjunction_relation_indices:
            adjacency = relation_adjacencies[relation_index] if relation_index < len(relation_adjacencies) else None
            if adjacency is None:
                continue
            indices = adjacency.coalesce().indices()
            member_edges.append((relation_index, indices[1].clone(), indices[0].clone()))
        self._conjunction_edge_cache = (cache_key, member_edges)
        return member_edges

    def conjunction_messages(self, node_major_state: Tensor, layer_index: int,
                             member_edges: list[tuple[int, Tensor, Tensor]]) -> Tensor:
        """The messages of the conjunction relations: one minimum per node-state channel over each node's members.

        node_major_state is [num_graph_nodes, batch_size, node_state_dim], so the minimum runs channel by channel and
        a complex is as present as its scarcest subunit in every preserved quantity separately. A node with no member
        gets zero, which the relation's weight maps to zero, so nothing is added where there is no complex.
        """
        messages = torch.zeros_like(node_major_state)
        for relation_index, member_source, member_target in member_edges:
            source = member_source.to(node_major_state.device)
            target = member_target.to(node_major_state.device)
            members = node_major_state[source]  # [num_member_edges, batch_size, node_state_dim]
            target_rows = target[:, None, None].expand_as(members)
            smallest = torch.zeros_like(node_major_state).scatter_reduce(0, target_rows, members, reduce="amin", include_self=False)
            if self.conjunction_aggregation == "soft_minimum":
                # -T log sum exp(-x/T), shifted by the hard minimum of each node so the exponential cannot overflow;
                # the shift is exact, not an approximation, because log sum exp(-(x - m)/T) = log sum exp(-x/T) + m/T
                shifted = torch.exp(-(members - smallest[target]) / self.soft_minimum_temperature)
                summed = torch.zeros_like(node_major_state).index_add(0, target, shifted)
                has_member = summed > 0
                smallest = torch.where(has_member, smallest - self.soft_minimum_temperature * torch.log(summed.clamp_min(1e-12)), smallest)
            messages = messages + smallest @ self.relation_weight[layer_index, relation_index]
        return messages

    def entry_messages(self, node_major_state: Tensor, layer_index: int | None) -> Tensor:
        """The messages of the entry edges, summed into their targets: weight * (x W_unsigned) + sign * weight * (x W_signed).

        node_major_state is [num_graph_nodes, batch_size, node_state_dim] with absent entry nodes already at zero, so an
        entry node sends only in the perturbations seeded on it and no mask is needed here. layer_index None takes the
        maps of the delivery before the first layer."""
        unsigned_weight = self.entry_first_unsigned_weight if layer_index is None else self.entry_unsigned_weight[layer_index]
        signed_weight = self.entry_first_signed_weight if layer_index is None else self.entry_signed_weight[layer_index]
        source_state = node_major_state[self.entry_edge_source]  # [num_entry_edges, batch_size, node_state_dim]
        per_edge = ((source_state * self.entry_edge_weight[:, None, None]) @ unsigned_weight
                    + (source_state * self.entry_edge_signed_weight[:, None, None]) @ signed_weight)
        return torch.zeros_like(node_major_state).index_add(0, self.entry_edge_target, per_edge)

    def propagate(self, node_state_field: Tensor, relation_adjacencies: list[Tensor | None], present_nodes: Tensor | None = None) -> Tensor:
        """L rounds of typed mean aggregation; returns the field as [batch_size, num_graph_nodes, node_state_dim].

        present_nodes ([batch_size, num_graph_nodes], or None when every node is always present) holds each row's
        absent entry nodes at zero before the first layer and after every layer, and adds the entry edges' messages.

        Per layer the new state is relu(X W_self + sum_r A_r X W_r + b). The states are kept node-major,
        [num_graph_nodes, batch_size, node_state_dim], so the sparse product reads them as a
        [num_graph_nodes, batch_size * node_state_dim] view without the permute-and-copy that a batch-major
        layout needs; all relations are aggregated by one sparse product with the stacked adjacency, using
        A_r (X W_r) = (A_r X) W_r. On the metabolic slice this halved the memory traffic of the batch-major,
        per-relation loop it replaced (profiling in docs/experiment_design.md section 8). Since 8 October 2026 the
        product keeps only the stacked rows that hold an edge (compressed_relation_adjacency), each relation's rows are
        transformed by its W_r, and the rows are added into their target nodes; on graph_full_neuronal that is 54,982
        of 700,435 rows. The output is identical to the full stacked product up to floating point rounding.
        """
        batch_size = node_state_field.shape[0]
        compressed_adjacency, row_target, segments = self.compressed_relation_adjacency(relation_adjacencies)
        member_edges = self.conjunction_member_edges(relation_adjacencies) if self.conjunction_relation_indices else []
        node_major_state = node_state_field.permute(1, 0, 2).contiguous()  # [N, B, D]
        node_major_presence = None
        if present_nodes is not None:
            node_major_presence = present_nodes.to(node_major_state.dtype).permute(1, 0)[:, :, None]  # [N, B, 1]
            node_major_state = node_major_state * node_major_presence
            if self.entry_first_unsigned_weight is not None and self.entry_edge_source.numel():
                # the targets enter the first layer already holding the drug's signal, as a seed on them would
                node_major_state = node_major_state + self.entry_messages(node_major_state, None)
        for layer_index in range(self.num_message_passing_layers):
            self_messages = node_major_state @ self.self_weight[layer_index]
            if member_edges:
                self_messages = self_messages + self.conjunction_messages(node_major_state, layer_index, member_edges)
            if node_major_presence is not None and self.entry_edge_source.numel():
                self_messages = self_messages + self.entry_messages(node_major_state, layer_index)
            if segments:
                flattened = node_major_state.reshape(self.num_graph_nodes, batch_size * self.node_state_dim)
                aggregated = torch.sparse.mm(compressed_adjacency.to(flattened.device), flattened).view(-1, batch_size, self.node_state_dim)  # [K, B, D]
                transformed = torch.cat([aggregated[start:end] @ self.relation_weight[layer_index, relation_index] for relation_index, start, end in segments], dim=0)
                relation_messages = torch.zeros_like(self_messages).index_add(0, row_target.to(flattened.device), transformed)
                pre_activation = self_messages + relation_messages + self.layer_bias[layer_index]
            else:
                pre_activation = self_messages + self.layer_bias[layer_index]
            node_major_state = torch.relu(pre_activation)
            if node_major_presence is not None:  # the layer bias alone would give an absent node a state
                node_major_state = node_major_state * node_major_presence
        return node_major_state.permute(1, 0, 2).contiguous()
