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
from mechanistic_pathway_learning.models.edge_gates import EdgeGates, GateLayout, flip_effect_estimates, relation_adjacency_parts, sparse_layout, sparse_values_product

PERTURBATION_FEATURE_DIM = 2  # (sign, magnitude)
CONJUNCTION_AGGREGATIONS = ("mean", "minimum", "soft_minimum")
DEFAULT_CONJUNCTION_RELATIONS = ("member_of",)  # the obligate-complex relation of the physiology graph
DEFAULT_SOFT_MINIMUM_TEMPERATURE = 0.1
# How a relation's transform is parameterised (preregistration amendment item 8 of 12 October 2026; the owner: "set
# the reduced rank bases with the example share flow through each edge type and 1 repeated layer regardless of depth").
#   "per_layer_dense", every run before 13 October 2026: a full node_state_dim by node_state_dim matrix per relation
#     per layer, 61,440 parameters at three layers, 20 relations and 32 channels, or 20.65 per kept training pair.
#   "shared_low_rank": one set applied at every round, so the depth is a number of rounds and not a parameter count.
#     Each relation owns a diagonal of node_state_dim parameters, and the relations that enough training examples
#     reach additionally own a rank-r factorisation, the transform being diag(d) + L R^T and applied as
#     x * d + (x L) R^T so the full matrix is never built. 1,792 parameters in all, 0.60 per kept pair.
# The rank follows the number of development training pairs that reach a relation, which is why the caller passes the
# indices rather than this module choosing them; the registered nine are encodes, catalyzed_by, substrate_of,
# activates, product_of, binds, regulates_transcription_of, inhibits and member_of, at 2,100 to 2,468 pairs against
# 451 or fewer for the other eleven. **Rank 4 is the only step up the registration permits**, named before any run so
# it cannot be chosen after a result is seen.
RELATION_PARAMETERISATIONS: tuple[str, ...] = ("per_layer_dense", "shared_low_rank")
DEFAULT_RELATION_FACTOR_RANK = 2
PERMITTED_RELATION_FACTOR_RANKS: tuple[int, ...] = (2, 4)
# The norm bound item 8 registers, for the same reason item 7 bounds the linear encoder: one set of weights applied at
# every round is a repeated map, and without a bound on its gain the iteration has no fixed point to run to. The
# diagonal is passed through tanh and scaled by the bound, and the factor is rescaled whenever the product of its two
# sides' spectral norms would exceed it, so the whole transform's induced norm stays at or below the bound.
DEFAULT_RELATION_NORM_BOUND = 0.9
RELATION_NORM_EPSILON = 1e-12
# Amendment item 7, for this encoder: the rounds run until the relative change of the state falls below a tolerance
# rather than to a fixed count (the owner, 12 October 2026: "Use the convergence instead of 8 steps. 8 steps is
# arbitrary, not physiological"; "Propagation window is settled. It goes to a fixed point to simulate steady state").
# This needs shared weights, which is why it waited on item 8: with one set per layer the depth IS the parameter
# count, and a run to convergence would be a different model at every round. With one shared set, bounded by
# DEFAULT_RELATION_NORM_BOUND, the update is the same contraction each round and the limit is what the encoder means.
# The perturbation is held on at every round for the same reason as in the linear encoder: given once, a stable
# shared layer forgets it.
MESSAGE_PASSING_FIXED_POINT_NORM_EPSILON = 1e-12
# What the bound has to cover for the update to contract, which the first attempt at the fixed point got wrong and the
# measurement caught: it ran to the cap of 60 rounds without settling. Each relation's adjacency is in-degree
# normalised, so one relation's message is bounded by its transform, but the relations are **summed** into a node, so a
# node fed by R relations receives up to R times that. And the self weight carried no bound at all. The update's
# Lipschitz constant is therefore ||self|| + R * relation bound, which is above 1 for any R worth having.
# Two things follow, and both are what the linear-response encoder already does. The summed relation messages are
# divided by the number of relations feeding the node, so the aggregate is bounded by one transform; and the self
# weight is spectrally bounded. The two then share relation_norm_bound, half each, so the whole update is bounded by
# it and the iteration contracts. Without shared weights none of this arises, because there is no repeated map, so the
# division and the self bound apply under "shared_low_rank" only and every earlier run is untouched.
SELF_AND_RELATION_SHARE_OF_THE_BOUND = 0.5


class RelationalMessagePassingEncoder(nn.Module):
    def __init__(
        self,
        num_graph_nodes: int,
        num_relation_types: int,
        node_state_dim: int,
        num_message_passing_layers: int = 3,
        node_features: Tensor | None = None,
        num_descriptor_columns: int = 0,
        descriptor_columns_kept: Tensor | None = None,
        relation_parameterisation: str = "per_layer_dense",
        fixed_point_tolerance: float | None = None,
        max_message_passing_rounds: int = 60,
        factorised_relation_indices: "tuple[int, ...] | None" = None,
        relation_factor_rank: int = DEFAULT_RELATION_FACTOR_RANK,
        relation_norm_bound: float = DEFAULT_RELATION_NORM_BOUND,
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
        change_only_edge_index: Tensor | None = None,
        change_only_relation_index: int | None = None,
        change_only_edge_relation: Tensor | None = None,
    ) -> None:
        super().__init__()
        self.edge_gates: GateLayout | None = None  # attach_edge_gates (decision 9.57); None is the registered model
        # the binary pass of decision 9.62 collects each forward's gated values and reads their gradients afterwards
        self.collect_flip_gradients = False
        self.flip_collection: list[Tensor] = []
        if change_only_edge_relation is not None and change_only_relation_index is not None:
            raise ValueError("give the relation of the change-only edges once: one index for all of them, or one index per edge")
        if (change_only_edge_index is None) != (change_only_relation_index is None and change_only_edge_relation is None):
            raise ValueError("change-only edges need both their edge index and the index of their relation")
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
        # edges that carry how far their source has moved from its unperturbed state, not the state itself (a carrier's
        # sequestration of a drug: silent while the carrier is at rest, so it cannot mark which drugs have a known carrier)
        # They may be of several relations (a carrier's sequestration and an enzyme's clearance): change_only_edge_relation
        # names each edge's relation, and each relation is averaged into its targets apart from the others.
        self.change_only_relation_index = change_only_relation_index
        if change_only_edge_index is None:
            self.register_buffer("change_only_edge_source", None, persistent=False)
            self.register_buffer("change_only_edge_target", None, persistent=False)
            self.register_buffer("change_only_edge_scale", None, persistent=False)
            self.register_buffer("change_only_edge_relation", None, persistent=False)
            self.change_only_relation_indices: tuple[int, ...] = ()
        else:
            change_only_target = change_only_edge_index[1].long()
            edge_relation = (torch.full((len(change_only_target),), int(change_only_relation_index), dtype=torch.long) if change_only_edge_relation is None
                             else change_only_edge_relation.long())
            if edge_relation.shape != change_only_target.shape:
                raise ValueError("change_only_edge_relation must have one entry per change-only edge")
            relation_indices = tuple(sorted({int(index) for index in edge_relation.tolist()})) or ((int(change_only_relation_index),) if change_only_relation_index is not None else ())
            if any(not 0 <= index < num_relation_types for index in relation_indices):
                raise ValueError("change_only_relation_index must name one of the relation types")
            # the in-degree is counted within a relation, so each relation's messages are a mean, as for every other relation
            relation_and_target = edge_relation * num_graph_nodes + change_only_target
            in_degree = torch.zeros(num_relation_types * num_graph_nodes).index_add(0, relation_and_target, torch.ones(len(change_only_target)))
            self.register_buffer("change_only_edge_source", change_only_edge_index[0].long(), persistent=False)
            self.register_buffer("change_only_edge_target", change_only_target, persistent=False)
            self.register_buffer("change_only_edge_scale", 1.0 / in_degree[relation_and_target], persistent=False)
            self.register_buffer("change_only_edge_relation", edge_relation, persistent=False)
            self.change_only_relation_indices = relation_indices
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
        if fixed_point_tolerance is not None:
            if fixed_point_tolerance <= 0:
                raise ValueError("fixed_point_tolerance must be positive, or None to run a fixed number of rounds")
            if relation_parameterisation != "shared_low_rank":
                raise ValueError("fixed_point_tolerance needs relation_parameterisation 'shared_low_rank': with one "
                                 "set of weights per layer the depth is the parameter count, so a run to convergence "
                                 "would be a different model at every round (amendment items 7 and 8)")
        if max_message_passing_rounds < 1:
            raise ValueError("max_message_passing_rounds must be at least 1")
        self.fixed_point_tolerance = fixed_point_tolerance
        self.max_message_passing_rounds = max_message_passing_rounds
        self.rounds_taken = 0  # what the last call took, for the run record
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
                self.register_buffer("node_features_without_descriptors", features_without_descriptors(self.node_features, num_descriptor_columns, descriptor_columns_kept), persistent=False)
        self.perturbation_injection = nn.Linear(PERTURBATION_FEATURE_DIM, node_state_dim)
        scale = 1.0 / math.sqrt(node_state_dim)
        if relation_parameterisation not in RELATION_PARAMETERISATIONS:
            raise ValueError(f"relation_parameterisation must be one of {RELATION_PARAMETERISATIONS}, not {relation_parameterisation!r}")
        self.relation_parameterisation = relation_parameterisation
        self.relation_factor_rank = relation_factor_rank
        self.relation_norm_bound = relation_norm_bound
        # the layers a parameter is held for: one shared set under "shared_low_rank" (amendment item 8), one per layer
        # under "per_layer_dense", which is what every run before 13 October 2026 built
        parameter_layers = 1 if relation_parameterisation == "shared_low_rank" else num_message_passing_layers
        self.parameter_layers = parameter_layers
        if relation_parameterisation == "shared_low_rank":
            if relation_factor_rank not in PERMITTED_RELATION_FACTOR_RANKS:
                raise ValueError(f"relation_factor_rank must be one of {PERMITTED_RELATION_FACTOR_RANKS} "
                                 f"(amendment item 8 names rank 4 as the only step up from 2), not {relation_factor_rank}")
            factorised = sorted(set(int(index) for index in (factorised_relation_indices or ())))
            if factorised and (factorised[0] < 0 or factorised[-1] >= num_relation_types):
                raise ValueError("factorised_relation_indices must index the relation types given")
            self.register_buffer("factorised_relation_indices", torch.tensor(factorised, dtype=torch.long), persistent=False)
            self.relation_diagonal = nn.Parameter(torch.zeros(num_relation_types, node_state_dim))
            self.relation_factor_left = nn.Parameter(torch.randn(len(factorised), node_state_dim, relation_factor_rank) * scale) \
                if factorised else None
            self.relation_factor_right = nn.Parameter(torch.randn(len(factorised), node_state_dim, relation_factor_rank) * scale) \
                if factorised else None
            self.relation_weight = None
        else:
            self.factorised_relation_indices = None
            self.relation_diagonal = self.relation_factor_left = self.relation_factor_right = None
            self.relation_weight = nn.Parameter(
                torch.randn(num_message_passing_layers, num_relation_types, node_state_dim, node_state_dim) * scale
            )
        self.self_weight = nn.Parameter(torch.randn(parameter_layers, node_state_dim, node_state_dim) * scale)
        self.bounds_the_update = relation_parameterisation == "shared_low_rank"
        self.layer_bias = nn.Parameter(torch.zeros(parameter_layers, node_state_dim))
        # created last and only with entry nodes, so an encoder without them draws the parameters it always drew
        self.entry_unsigned_weight, self.entry_signed_weight = None, None
        self.entry_first_unsigned_weight, self.entry_first_signed_weight = None, None
        if entry_node_mask is not None:
            self.entry_unsigned_weight = nn.Parameter(torch.randn(parameter_layers, node_state_dim, node_state_dim) * scale)
            self.entry_signed_weight = nn.Parameter(torch.randn(parameter_layers, node_state_dim, node_state_dim) * scale)
            if entry_before_first_layer:
                self.entry_first_unsigned_weight = nn.Parameter(torch.randn(node_state_dim, node_state_dim) * scale)
                self.entry_first_signed_weight = nn.Parameter(torch.randn(node_state_dim, node_state_dim) * scale)

    def parameter_layer(self, layer_index: int) -> int:
        """Which set of weights a round reads: its own under "per_layer_dense", the single shared set otherwise."""
        return 0 if self.parameter_layers == 1 else layer_index

    def relation_transform(self, relation_index: int, layer_index: int):
        """The transform of one relation at one round, as a function of [*, node_state_dim] messages.

        Under "per_layer_dense" this is that layer's matrix for the relation. Under "shared_low_rank" it is
        diag(d) + L R^T, bounded: the diagonal is tanh scaled by relation_norm_bound, and the factor is rescaled when
        the product of its sides' spectral norms would carry the transform past the bound. The matrix is never built,
        so the cost is two thin products rather than one square one."""
        if self.relation_parameterisation == "per_layer_dense":
            weight = self.relation_weight[layer_index, relation_index]
            return lambda messages: messages @ weight
        bound = self.relation_norm_bound * SELF_AND_RELATION_SHARE_OF_THE_BOUND
        diagonal = bound * torch.tanh(self.relation_diagonal[relation_index])
        position = (self.factorised_relation_indices == relation_index).nonzero()
        if position.numel() == 0:
            return lambda messages: messages * diagonal
        row = int(position[0])
        left, right = self.relation_factor_left[row], self.relation_factor_right[row]
        left_norm = torch.linalg.matrix_norm(left, ord=2).clamp_min(RELATION_NORM_EPSILON)
        right_norm = torch.linalg.matrix_norm(right, ord=2).clamp_min(RELATION_NORM_EPSILON)
        # the diagonal already spends up to relation_norm_bound, so the factor is held to what is left of it
        room = (bound - diagonal.abs().max()).clamp_min(0.0)
        rescale = torch.clamp(room / (left_norm * right_norm), max=1.0)
        return lambda messages: messages * diagonal + ((messages @ left) @ right.transpose(0, 1)) * rescale

    def bounded_self_weight(self, layer_index: int) -> Tensor:
        """The self transform, spectrally bounded to its share of relation_norm_bound under shared weights.

        Unbounded it is what stopped the first fixed-point attempt settling: a repeated map whose self term can grow
        without limit has no fixed point to run to."""
        weight = self.self_weight[layer_index]
        if not self.bounds_the_update:
            return weight
        bound = self.relation_norm_bound * SELF_AND_RELATION_SHARE_OF_THE_BOUND
        norm = torch.linalg.matrix_norm(weight, ord=2).clamp_min(RELATION_NORM_EPSILON)
        return weight * torch.clamp(bound / norm, max=1.0)

    def relation_transform_norm_bound(self) -> float:
        """The largest induced norm any relation's transform can have, which is what makes the repeated map a
        contraction and so gives the iteration of amendment item 7 a fixed point to reach."""
        return float(self.relation_norm_bound)

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
    def build_relation_adjacencies(edge_index: Tensor, edge_relation_type: Tensor, num_graph_nodes: int, num_relation_types: int,
                                   edge_sign: Tensor | None = None, edge_weight: Tensor | None = None) -> list[Tensor | None]:
        """One sparse, in-degree-normalised adjacency per relation: entry [destination, source] = sign x weight / in_degree(destination).

        The signs. Until 12 October 2026 this built an unsigned adjacency, so an inhibition and an activation of the same
        relation arrived alike and the message passing encoder read none of the graph's signs, while the linear-response
        encoder read all of them. The owner, 12 October 2026: "Why would MP not read the edge sign?! That seems like a
        bug"; and then "The edge signs are also settled: the MP has to adhere to the edge signs the same way the linear
        encoder does". Registered as item 9 of the amendment of 12 October. Twelve of the twenty relations have a
        constant sign that the relation index already carried, so nothing changes for them; the change is the eight whose
        sign varies within the relation, 42,903 edges, of which 7,321 are negative and were read as positive, 7,087 of
        them transcriptional repressions among `regulates_transcription_of`'s 41,638 edges.

        One consequence registered with it: with signs in the values, the in-degree mean can cancel, so a node with one
        activator and one repressor of equal weight receives exactly zero. That is the arithmetic of a mean over signed
        messages and not a defect of this function.

        edge_weight carries an edge's own strength where the graph gives it one, which today is the driving-force edges
        of the computed membrane potential (graph/membrane_potential.py). The normalisation stays the count of incoming
        edges, not the sum of their weights, so a relation's weights change what arrives without changing the divisor.
        """
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
            if edge_sign is not None:
                values = values * edge_sign[edges_of_relation].to(values.dtype)
            if edge_weight is not None:
                values = values * edge_weight[edges_of_relation].to(values.dtype)
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
            reference_row = torch.arange(batch_size, 2 * batch_size, device=both.device).repeat(2)  # each row and each reference reads its own reference
            propagated = self.propagate(both, relation_adjacencies, present, reference_row)
            return propagated[:batch_size] - propagated[batch_size:]
        if self.change_only_edge_source is not None:
            # one more row, the unperturbed reference, propagated beside the batch so that each layer can read how far a
            # source has moved from it; the row itself reads no change, so it is the reference it would be alone
            batch_size, width = perturbation_node_index.shape
            joint_index = torch.cat([perturbation_node_index, torch.full((1, width), -1, dtype=perturbation_node_index.dtype, device=perturbation_node_index.device)], dim=0)
            joint_injection = torch.cat([perturbation_sign_and_magnitude, torch.zeros_like(perturbation_sign_and_magnitude[:1])], dim=0)
            reference_row = torch.full((batch_size + 1,), batch_size, dtype=torch.long, device=perturbation_node_index.device)
            propagated = self.propagate(self.initial_node_state_field(joint_index, joint_injection), relation_adjacencies, self.present_nodes(joint_index), reference_row)
            return propagated[:batch_size] - propagated[batch_size:]
        perturbed = self.forward(perturbation_node_index, perturbation_sign_and_magnitude, relation_adjacencies=relation_adjacencies)
        return perturbed - self.unperturbed_node_state_field(relation_adjacencies)

    def change_only_messages(self, node_major_state: Tensor, layer_index: int, reference_row: Tensor) -> Tensor:
        """The messages of the change-only edges, averaged into their targets: ((x - x_reference) W_r) / in-degree.

        reference_row ([batch_size]) names, for each row of the batch, the row that holds its unperturbed reference; a
        reference names itself and so sends nothing. With every source at rest the message is exactly zero, whatever the
        source's own state is."""
        source_state = node_major_state[self.change_only_edge_source]  # [num_edges, batch_size, node_state_dim]
        change = (source_state - source_state[:, reference_row]) * self.change_only_edge_scale[:, None, None]
        if len(self.change_only_relation_indices) == 1:
            return torch.zeros_like(node_major_state).index_add(0, self.change_only_edge_target, self.relation_transform(self.change_only_relation_indices[0], self.parameter_layer(layer_index))(change))
        messages = torch.zeros_like(node_major_state)
        for relation_index in self.change_only_relation_indices:  # each relation through its own weight
            of_relation = self.change_only_edge_relation == relation_index
            messages = messages.index_add(0, self.change_only_edge_target[of_relation], self.relation_transform(relation_index, self.parameter_layer(layer_index))(change[of_relation]))
        return messages

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
                             if adjacency is not None and index not in self.conjunction_relation_indices and index not in self.entry_relation_indices
                             and index not in self.change_only_relation_indices]
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

    def attach_edge_gates(self, edge_gates: GateLayout, edge_index: Tensor, edge_relation_type: Tensor, relation_adjacencies: list[Tensor | None],
                          edge_sign: Tensor | None = None, edge_weight: Tensor | None = None) -> None:
        """Gate every edge of the stacked relations by its evidence prior (decision 9.57, models/edge_gates.py). The
        compressed adjacency's layout is read once from the static adjacencies; each propagate then rebuilds its values
        from the edges' values times their gates, so nothing is cached by tensor identity. The conjunction, entry and
        change-only relations are aggregated apart from the stack and are not gated; edge_sign and edge_weight must be
        the ones the static adjacencies were built with (run_main_model.message_passing_adjacencies), which this checks:
        with every gate at 1 the rebuilt values must equal the static ones."""
        if self.edge_gates is not None:
            raise ValueError("edge gates are attached once, to a freshly built encoder")
        if edge_gates.num_edges != edge_index.shape[1]:
            raise ValueError(f"edge_gates has {edge_gates.num_edges} edges and the graph {edge_index.shape[1]}")
        device = self.layer_bias.device
        static_compressed, row_target, segments = self.compressed_relation_adjacency(relation_adjacencies)
        _, present_index = self.stacked_relation_adjacency(relation_adjacencies)
        parts = relation_adjacency_parts(edge_index.long().cpu(), edge_relation_type.long().cpu(), self.num_graph_nodes, self.num_relation_types,
                                         None if edge_sign is None else edge_sign.cpu(), None if edge_weight is None else edge_weight.cpu())
        self.gate_present_relations = [int(index) for index in present_index.tolist()]
        for relation_index in self.gate_present_relations:
            part = parts[relation_index]
            for name in ("slot", "values", "numerator", "destination", "edges"):
                self.register_buffer(f"gate_part_{relation_index}_{name}", part[name].to(device), persistent=False)
            setattr(self, f"gate_part_{relation_index}_num_slots", part["num_slots"])
        self.register_buffer("gate_compressed_indices", static_compressed.indices().to(device), persistent=False)
        transposed_indices, transposed_order = sparse_layout(self.gate_compressed_indices, static_compressed.shape[0])
        self.register_buffer("gate_compressed_transposed_indices", transposed_indices, persistent=False)
        self.register_buffer("gate_compressed_transposed_order", transposed_order, persistent=False)
        self.register_buffer("gate_row_target", row_target.to(device), persistent=False)
        self.gate_segments = segments
        self.gate_compressed_shape = tuple(static_compressed.shape)
        self.gate_relation_adjacencies = list(relation_adjacencies)  # held, so their identities cannot be reused while the encoder lives
        with torch.no_grad():
            ungated = self.gated_compressed_values(torch.ones(edge_gates.num_edges, device=device))
        static_values = static_compressed.values().to(device)
        if ungated.shape != static_values.shape or not torch.allclose(ungated, static_values, rtol=1e-6, atol=1e-9):
            raise ValueError("with every gate at 1 the rebuilt adjacency differs from the static one: the edge signs or weights passed differ from the adjacencies'")
        read = [getattr(self, f"gate_part_{relation_index}_edges") for relation_index in self.gate_present_relations]
        self.edge_gates = edge_gates.to(device)
        self.edge_gates.set_read_edges(torch.cat(read) if read else torch.zeros(0, dtype=torch.long, device=device))

    def check_gated_adjacencies(self, relation_adjacencies: list[Tensor | None]) -> None:
        """A gated encoder rebuilds the stack from the adjacencies its gates were attached to; refuse any others, which
        would otherwise be silently replaced (a rewired graph needs a fresh encoder, as evaluate_on_rewired_graph.py builds)."""
        attached = self.gate_relation_adjacencies
        if relation_adjacencies is attached:
            return
        if len(relation_adjacencies) != len(attached) or any(given is not held for given, held in zip(relation_adjacencies, attached)):
            raise ValueError("these adjacencies are not the ones the edge gates were attached to; build a fresh encoder and attach its gates to them")

    def gated_compressed_values(self, edge_values: Tensor) -> Tensor:
        """The compressed adjacency's values in its stored order: per present relation, each coalesced entry's sum of its
        edges' sign x weight x gate over the gated in-degree, the sum of the gates of the relation's edges into the
        destination (decision 9.62). With every gate at 1 that is the count, so the static values come back; with gates
        of 0 and 1 it is the adjacency of the graph without the closed edges, as in the linear-response encoder, whose
        divisor is the gated sum of its weights."""
        values = []
        tiny = torch.finfo(edge_values.dtype).tiny
        for relation_index in self.gate_present_relations:
            gate = edge_values[getattr(self, f"gate_part_{relation_index}_edges")]
            destination = getattr(self, f"gate_part_{relation_index}_destination")
            gated_in_degree = torch.zeros(self.num_graph_nodes, dtype=gate.dtype, device=gate.device).index_add(0, destination, gate)
            edge_value = getattr(self, f"gate_part_{relation_index}_numerator") * gate / gated_in_degree[destination].clamp_min(tiny)
            slots = torch.zeros(getattr(self, f"gate_part_{relation_index}_num_slots"), dtype=edge_value.dtype, device=edge_value.device)
            values.append(slots.index_add(0, getattr(self, f"gate_part_{relation_index}_slot"), edge_value))
        return torch.cat(values) if values else torch.zeros(0, device=edge_values.device)

    def gated_row_is_open(self, edge_values: Tensor) -> Tensor:
        """[K] per compressed row (relation, target), whether any edge of it has a gate above 0: a row whose edges are all
        closed is a relation that no longer feeds its target, so the bounded update does not count it."""
        open_rows = []
        for relation_index, start, end in self.gate_segments:
            gate = edge_values[getattr(self, f"gate_part_{relation_index}_edges")]
            destination = getattr(self, f"gate_part_{relation_index}_destination")
            gated_in_degree = torch.zeros(self.num_graph_nodes, dtype=gate.dtype, device=gate.device).index_add(0, destination, gate)
            open_rows.append(gated_in_degree[self.gate_row_target[start:end]] > 0)
        return torch.cat(open_rows) if open_rows else torch.zeros(0, dtype=torch.bool, device=edge_values.device)

    def reset_flip_collection(self) -> None:
        self.flip_collection = []

    def flip_effect_estimates(self) -> Tensor:
        """[num_edges] the first-order change in the loss if each edge flipped from the gates' read-out (models/edge_gates.
        flip_effect_estimates), from the gradients of the gated values every forward since reset_flip_collection kept.
        A group is a relation's edges into one destination, and each edge counts 1 toward its divisor."""
        if self.edge_gates is None or not hasattr(self.edge_gates, "open_edges"):
            raise ValueError("flip estimates need hard edge gates (decision 9.62)")
        gradients = [values.grad for values in self.flip_collection if values.grad is not None]
        if not gradients:
            raise ValueError("no gated values with a gradient were collected; set collect_flip_gradients and run a backward pass first")
        slot_gradient = torch.stack(gradients).sum(dim=0)
        groups, origins, numerators, entry_gradients, offset = [], [], [], [], 0
        for relation_index in self.gate_present_relations:
            slot = getattr(self, f"gate_part_{relation_index}_slot")
            groups.append(relation_index * self.num_graph_nodes + getattr(self, f"gate_part_{relation_index}_destination"))
            origins.append(getattr(self, f"gate_part_{relation_index}_edges"))
            numerators.append(getattr(self, f"gate_part_{relation_index}_numerator"))
            entry_gradients.append(slot_gradient[offset + slot])
            offset += getattr(self, f"gate_part_{relation_index}_num_slots")
        group = torch.cat(groups)
        _, dense_group = torch.unique(group, return_inverse=True)
        numerator = torch.cat(numerators)
        return flip_effect_estimates(dense_group, torch.cat(origins), numerator, torch.ones_like(numerator), torch.cat(entry_gradients),
                                     self.edge_gates.open_edges(), self.edge_gates.num_edges)

    def gated_compressed_adjacency(self) -> tuple[Tensor, Tensor, list[tuple[int, int, int]]]:
        """The compressed adjacency with the gated values, as a sparse tensor for reading; propagate multiplies through
        edge_gates.sparse_values_product instead, whose value gradient forms no dense matrix."""
        values = self.gated_compressed_values(self.edge_gates.edge_values())
        compressed = torch.sparse_coo_tensor(self.gate_compressed_indices, values, self.gate_compressed_shape, is_coalesced=True)
        return compressed, self.gate_row_target, self.gate_segments

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
            messages = messages + self.relation_transform(relation_index, self.parameter_layer(layer_index))(smallest)
        return messages

    def entry_messages(self, node_major_state: Tensor, layer_index: int | None) -> Tensor:
        """The messages of the entry edges, summed into their targets: weight * (x W_unsigned) + sign * weight * (x W_signed).

        node_major_state is [num_graph_nodes, batch_size, node_state_dim] with absent entry nodes already at zero, so an
        entry node sends only in the perturbations seeded on it and no mask is needed here. layer_index None takes the
        maps of the delivery before the first layer."""
        unsigned_weight = self.entry_first_unsigned_weight if layer_index is None else self.entry_unsigned_weight[self.parameter_layer(layer_index)]
        signed_weight = self.entry_first_signed_weight if layer_index is None else self.entry_signed_weight[self.parameter_layer(layer_index)]
        source_state = node_major_state[self.entry_edge_source]  # [num_entry_edges, batch_size, node_state_dim]
        per_edge = ((source_state * self.entry_edge_weight[:, None, None]) @ unsigned_weight
                    + (source_state * self.entry_edge_signed_weight[:, None, None]) @ signed_weight)
        return torch.zeros_like(node_major_state).index_add(0, self.entry_edge_target, per_edge)

    def propagate(self, node_state_field: Tensor, relation_adjacencies: list[Tensor | None], present_nodes: Tensor | None = None,
                  reference_row: Tensor | None = None) -> Tensor:
        """L rounds of typed mean aggregation; returns the field as [batch_size, num_graph_nodes, node_state_dim].

        present_nodes ([batch_size, num_graph_nodes], or None when every node is always present) holds each row's
        absent entry nodes at zero before the first layer and after every layer, and adds the entry edges' messages.
        reference_row ([batch_size], or None) names each row's unperturbed reference in the same batch, which the
        change-only edges read; without it those edges send nothing, so a field computed alone (forward) treats every
        source as at rest.

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
        gated_values = None
        gated_row_open = None
        if self.edge_gates is not None:  # rebuilt each call from the gates (decision 9.57); never cached by tensor identity
            self.check_gated_adjacencies(relation_adjacencies)
            edge_values = self.edge_gates.edge_values()
            gated_values = self.gated_compressed_values(edge_values)
            if self.collect_flip_gradients and gated_values.requires_grad is False:
                gated_values = gated_values.detach().requires_grad_(True)
            if self.collect_flip_gradients:
                gated_values.retain_grad()
                self.flip_collection.append(gated_values)
            gated_row_open = self.gated_row_is_open(edge_values) if self.bounds_the_update else None
            compressed_adjacency, row_target, segments = None, self.gate_row_target, self.gate_segments
        else:
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
        rounds = self.max_message_passing_rounds if self.fixed_point_tolerance is not None else self.num_message_passing_layers
        previous_state, self.rounds_taken = node_major_state, 0
        for layer_index in range(rounds):
            parameter_layer = self.parameter_layer(layer_index)
            self_messages = node_major_state @ self.bounded_self_weight(parameter_layer)
            if member_edges:
                self_messages = self_messages + self.conjunction_messages(node_major_state, layer_index, member_edges)
            if node_major_presence is not None and self.entry_edge_source.numel():
                self_messages = self_messages + self.entry_messages(node_major_state, layer_index)
            if reference_row is not None and self.change_only_edge_source is not None and self.change_only_edge_source.numel():
                self_messages = self_messages + self.change_only_messages(node_major_state, layer_index, reference_row)
            if segments:
                flattened = node_major_state.reshape(self.num_graph_nodes, batch_size * self.node_state_dim)
                if gated_values is None:
                    aggregated = torch.sparse.mm(compressed_adjacency.to(flattened.device), flattened).view(-1, batch_size, self.node_state_dim)  # [K, B, D]
                else:
                    aggregated = sparse_values_product(self.gate_compressed_indices, gated_values, flattened, self.gate_compressed_shape[0], self.gate_compressed_transposed_indices,
                                                       self.gate_compressed_transposed_order).view(-1, batch_size, self.node_state_dim)
                transformed = torch.cat([self.relation_transform(relation_index, parameter_layer)(aggregated[start:end]) for relation_index, start, end in segments], dim=0)
                relation_messages = torch.zeros_like(self_messages).index_add(0, row_target.to(flattened.device), transformed)
                if self.bounds_the_update:  # the summed relations are a mean, so the aggregate is bounded by one transform
                    row_counts = torch.ones(row_target.numel(), device=relation_messages.device) if gated_row_open is None \
                        else gated_row_open.to(device=relation_messages.device, dtype=torch.float32)  # a relation whose gates are all closed no longer feeds
                    feeding = torch.zeros(self.num_graph_nodes, device=relation_messages.device).index_add(
                        0, row_target.to(relation_messages.device), row_counts)
                    relation_messages = relation_messages / feeding.clamp_min(1.0)[:, None, None]
                pre_activation = self_messages + relation_messages + self.layer_bias[parameter_layer]
            else:
                pre_activation = self_messages + self.layer_bias[parameter_layer]
            node_major_state = torch.relu(pre_activation)
            if node_major_presence is not None:  # the layer bias alone would give an absent node a state
                node_major_state = node_major_state * node_major_presence
            self.rounds_taken += 1
            if self.fixed_point_tolerance is not None:
                change = torch.linalg.vector_norm(node_major_state - previous_state) / \
                    torch.linalg.vector_norm(previous_state).clamp_min(MESSAGE_PASSING_FIXED_POINT_NORM_EPSILON)
                if float(change) < self.fixed_point_tolerance:
                    break
                previous_state = node_major_state
        return node_major_state.permute(1, 0, 2).contiguous()
