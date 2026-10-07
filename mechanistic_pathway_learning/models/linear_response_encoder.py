"""Linear-response encoder: a time-invariant, signed state space over the graph (design section 5.2, route 1).

The message passing encoder (relational_message_passing_encoder.py) emulates time badly for a chronic perturbation:
each of its L layers has its own weights, so the "dynamics" change from step to step; the perturbation is only the
initial state, so its signal is a transient rather than a sustained input; the difference field is exactly zero more
than L edges from the perturbed nodes (two on the default configuration, which on the metabolic layer reaches an
enzyme's reactions and their direct products and stops); and edges only run substrate -> reaction -> product, so a
blocked reaction cannot raise its substrate.

This encoder treats the field as the steady-state response of a linear system to a sustained input, the reading that
metabolic control analysis gives a small enzyme-level change (Kacser and Burns 1973; Heinrich and Rapoport 1974):

    h(t + 1) = (1 - damping) * h(t) + damping * (sum_r gain_r * (S_r h(t)) + u)

where u carries the perturbation at its nodes at every step, S_r is the signed, in-degree-normalised adjacency of
relation r (entry [destination, source] = sign of the edge / in-degree of the destination under r) and gain_r is a
learned per-channel gain. The same transition applies at every step, so after num_propagation_steps the field
approaches the fixed point h = sum_r gain_r * (S_r h) + u (damping sets the transient only). Gains of relations whose
edges carry a biological sign are positive, so the sign comes from the graph; relations without a sign (binds) get a
gain of either sign. Each node averages within each relation feeding it and then across those relations (entry =
sign / (in-degree under r * number of relations feeding the node)), so every row of sum_r |S_r| sums to at most one;
gains are bounded by contraction < 1, so the transition is a contraction in the max norm and the iteration converges.
(Bounding gains by contraction / R instead shrank the signal by about 0.2 per step, so changes six edges away were
of order 1e-9.)

Currency metabolites (water, protons, ATP, NAD(P)H and the like, flagged in the graph) receive the response but do
not pass it on: through them every reaction would reach every other in two steps, which is the hub problem of design
assumption A9 in another form. Carrier edges (cofactor_edges, mechanistic_pathway_learning/graph/cofactor_edges.py:
biopterins, quinones, folate carriers, transamination pairs and the like) can be given relations of their own,
cosubstrate_of (carrier -> reaction, supply) and coproduct_of (reaction -> carrier), whose gains are learned apart from
those of main substrates and products, and carriers get no depletion edge. A carrier pool is recycled (tetrahydrobiopterin
used by PAH is regenerated through PCBD1 and QDPR), so one consumer slowing down barely moves the pool, and the depletion
edge made a PAH loss raise tetrahydrobiopterin and through it tyrosine hydroxylase flux, L-dopa and dopamine, the
opposite of phenylketonuria. A synthesis defect does drain the pool, so the supply edges stay: dropping carrier edges
altogether left a GCH1 loss lowering tetrahydrobiopterin with nothing downstream able to see it, and 14 slice genes with
30 positive pairs make or recycle carriers (GCH1, PTS, SPR, QDPR, PCBD1, MTHFR, MTRR, SLC46A1, COQ2, COQ5, COQ7,
ALDH7A1, SLC19A3, HLCS). Untrained on the slice, with only the depletion edge dropped, a PAH loss raises phenylalanine
and lowers tyrosine, L-dopa and dopamine; a GCH1 loss lowers tetrahydrobiopterin, L-dopa and serotonin (dopamine still
rises, through a route not yet traced); an MTHFR loss raises homocysteine and lowers methionine. Carrier gains start low
(logit offset -3) and are learned.

The stoichiometric reading needs edges the graph does not store: a reaction's flux change depletes its substrates,
so every substrate_of edge (metabolite -> reaction) gets a reverse edge reaction -> metabolite with sign -1 under its
own relation, depletes_substrate. With catalyzed_by (gene -> reaction, +1), product_of (reaction -> metabolite, +1)
and substrate_of (metabolite -> reaction, +1, mass action), a loss of function at an enzyme lowers its reactions'
flux and its products and raises its substrates, and both changes travel on to the reactions that use them. The
in-degree normalisation reads isozymes as partial redundancy (losing one of k catalysts removes 1/k of the input) and
buffers metabolites with many producers. Transcriptional edges keep their own sign (CollecTRI activation +1,
repression -1), which the message passing encoder's single unsigned relation dropped.

The response propagates on a few channels (propagation_channels, default 4), each with its own per-relation gains,
so the channels settle at different time scales; a linear map then expands them to the node_state_dim channels the
heads read. Propagating all node_state_dim channels was memory-bound: each step wrote and reread a
[relations, nodes, batch, channels] tensor (about 200 MB on the metabolic slice at batch 16 and 32 channels), 19
seconds per batch forward and backward on one core.

The response decays steeply with distance under in-degree averaging (on the slice, the change at L-dopa is about 1e-5
of the change at the perturbed gene and at dopamine about 1e-8), so a linear readout sees little beyond the first few
edges. response_scale "signed_log" maps each channel through sign(h) * log(1 + |h| / s_c) with a learned scale s_c
(starting at 1e-4), which keeps the sign and the order of magnitude and makes distant changes readable; it is odd, so
a gain of function still mirrors a loss of function. response_scale "linear" leaves the response as it is.

The system is linear in the input, so the unperturbed response is zero and the field is a difference field by
construction: no base state, and no node identity, enters it. Node kinds enter only through an output gate,
field = (h W_expand) * sigmoid(node_features W + b), which lets a pooled readout tell where the response landed
without adding anything where it did not.
"""
from __future__ import annotations

import torch
from torch import Tensor, nn

SUBSTRATE_RELATION = "substrate_of"
PRODUCT_RELATION = "product_of"
DEPLETES_SUBSTRATE_RELATION = "depletes_substrate"
COSUBSTRATE_RELATION = "cosubstrate_of"
COPRODUCT_RELATION = "coproduct_of"
UNSIGNED_RELATIONS = ("binds",)
# how the signed adjacency is scaled so the iteration contracts: "in_degree" divides each entry by the
# destination's in-degree under its relation times the number of relations feeding it, which averages and so
# decays the response steeply with path length; "spectral" divides every entry by one global constant, the
# spectral radius of the unsigned aggregate, which contracts just as surely without the per-node divisor
# "in_degree_power" raises the in-degree to normalisation_exponent before dividing, then rescales the whole matrix
# once so the iteration still contracts. It exists because the decay is driven by the few hub destinations on a path
# (the in-degree divisor has median 2 but maximum 2048), so softening the divisor where it is large is what a long path
# needs; the exponent 1.0 reproduces "in_degree" exactly. "spectral" replaces the per-node divisor with one global
# constant and is kept as a measured negative result: that constant is 36.3 on the neuronal graph against an in-degree
# divisor below it at 98.8 percent of destinations, so it divides harder than what it replaces nearly everywhere and
# the field arriving at the membrane potential fell to 0.7 of its in-degree value (docs/membrane_potential_reach.md).
NORMALISATIONS = ("in_degree", "in_degree_power", "spectral")
CARRIER_INITIAL_GAIN_LOGIT_OFFSET = -3.0
INITIAL_LOG_RESPONSE_SCALE = -9.2  # log(1e-4)


class LinearResponseEncoder(nn.Module):
    def __init__(
        self,
        num_graph_nodes: int,
        relation_types: list[str],
        edge_source: Tensor,
        edge_target: Tensor,
        edge_relation: Tensor,
        edge_sign: Tensor,
        node_features: Tensor,
        node_state_dim: int,
        non_propagating_nodes: Tensor | None = None,
        cofactor_edges: Tensor | None = None,
        num_propagation_steps: int = 8,
        propagation_channels: int = 4,
        response_scale: str = "linear",
        damping: float = 0.5,
        contraction: float = 0.9,
        normalisation: str = "in_degree",
        normalisation_exponent: float = 0.5,
    ) -> None:
        super().__init__()
        if node_features.shape[0] != num_graph_nodes:
            raise ValueError("node_features must have one row per graph node")
        if not 0.0 < damping <= 1.0 or not 0.0 < contraction < 1.0:
            raise ValueError("damping must be in (0, 1] and contraction in (0, 1)")
        self.num_graph_nodes = num_graph_nodes
        self.node_state_dim = node_state_dim
        self.num_propagation_steps = num_propagation_steps
        self.damping = damping
        stacked_adjacency, relation_names = self.signed_stacked_adjacency(num_graph_nodes, relation_types, edge_source, edge_target, edge_relation, edge_sign,
                                                                          non_propagating_nodes, cofactor_edges, normalisation, normalisation_exponent)
        self.relation_names = relation_names
        self.register_buffer("stacked_adjacency", stacked_adjacency, persistent=False)
        self.register_buffer("node_features", node_features.to(torch.float32), persistent=False)
        self.register_buffer("relation_is_unsigned", torch.tensor([name in UNSIGNED_RELATIONS for name in relation_names]), persistent=False)
        self.maximum_gain = contraction
        self.propagation_channels = propagation_channels
        if response_scale not in ("linear", "signed_log"):
            raise ValueError("response_scale must be 'linear' or 'signed_log'")
        self.response_scale = response_scale
        self.log_response_scale = nn.Parameter(torch.full((propagation_channels,), INITIAL_LOG_RESPONSE_SCALE)) if response_scale == "signed_log" else None
        initial_gain_logit = torch.randn(len(relation_names), propagation_channels)  # spread so channels start at different time scales
        for carrier_relation in (COSUBSTRATE_RELATION, COPRODUCT_RELATION):
            if carrier_relation in relation_names:  # carrier pools are recycled and buffered: coupling through them starts weak (gain about 0.04)
                initial_gain_logit[relation_names.index(carrier_relation)] += CARRIER_INITIAL_GAIN_LOGIT_OFFSET
        self.gain_logit = nn.Parameter(initial_gain_logit)
        self.input_weight = nn.Parameter(torch.randn(propagation_channels) / propagation_channels**0.5)
        self.channel_expansion = nn.Parameter(torch.randn(propagation_channels, node_state_dim) / propagation_channels**0.5)
        self.output_gate = nn.Linear(node_features.shape[1], node_state_dim)

    @staticmethod
    def spectral_radius_of_unsigned_aggregate(num_graph_nodes: int, source: Tensor, target: Tensor, iterations: int = 200) -> float:
        """Largest eigenvalue magnitude of the unsigned aggregate adjacency M, where M[target, source] counts the edges
        from source to target over every relation. M is non-negative, so its Perron root is the largest eigenvalue and
        a power iteration from a positive start converges to it. This is what the spectral normalisation divides by:
        with entries scaled to give rho(sum_r |S_r|) = 1, the gain-weighted operator obeys
        rho(sum_r gain_r S_r) <= max_r |gain_r|, so the existing contraction bound on the gains still makes the
        iteration a contraction without dividing any node's input by its own in-degree."""
        aggregate = torch.sparse_coo_tensor(torch.stack([target.long(), source.long()]), torch.ones(len(source)),
                                            (num_graph_nodes, num_graph_nodes)).coalesce()
        vector = torch.full((num_graph_nodes,), num_graph_nodes ** -0.5)
        radius = 0.0
        for _ in range(iterations):
            product = torch.sparse.mm(aggregate, vector[:, None]).squeeze(1)
            norm = float(product.norm())
            if norm == 0.0:
                return 0.0
            vector, previous_radius = product / norm, radius
            radius = norm
            if abs(radius - previous_radius) <= 1e-9 * max(radius, 1.0):
                break
        return radius

    @staticmethod
    def signed_stacked_adjacency(num_graph_nodes: int, relation_types: list[str], edge_source: Tensor, edge_target: Tensor,
                                 edge_relation: Tensor, edge_sign: Tensor, non_propagating_nodes: Tensor | None = None,
                                 cofactor_edges: Tensor | None = None, normalisation: str = "in_degree",
                                 normalisation_exponent: float = 1.0) -> tuple[Tensor, list[str]]:
        """All relations' signed, normalised adjacencies stacked into one sparse [R * N, N] matrix, with the derived
        depletes_substrate relation appended, carrier edges moved to relations of their own when cofactor_edges is
        given, and edges leaving a non-propagating node dropped; returns the matrix and the relation names in stack order."""
        edge_relation = edge_relation.long().clone()
        relation_names = list(relation_types)
        if cofactor_edges is not None and SUBSTRATE_RELATION in relation_types and PRODUCT_RELATION in relation_types:
            cofactor_edges = cofactor_edges.bool()
            for main_relation, carrier_relation in ((SUBSTRATE_RELATION, COSUBSTRATE_RELATION), (PRODUCT_RELATION, COPRODUCT_RELATION)):
                relation_names.append(carrier_relation)
                edge_relation[cofactor_edges & (edge_relation == relation_types.index(main_relation))] = len(relation_names) - 1
        sources, targets, relations, signs = [edge_source.long()], [edge_target.long()], [edge_relation], [edge_sign.float()]
        for forward_relation, depletion_relation in ((SUBSTRATE_RELATION, DEPLETES_SUBSTRATE_RELATION),):  # carriers get none (module docstring)
            if forward_relation not in relation_names:
                continue
            forward_edges = edge_relation == relation_names.index(forward_relation)
            relation_names.append(depletion_relation)
            sources.append(edge_target[forward_edges].long())
            targets.append(edge_source[forward_edges].long())
            relations.append(torch.full((int(forward_edges.sum()),), len(relation_names) - 1, dtype=torch.long))
            signs.append(-torch.ones(int(forward_edges.sum())))
        source, target, relation, sign = torch.cat(sources), torch.cat(targets), torch.cat(relations), torch.cat(signs)
        if non_propagating_nodes is not None:
            keep = ~non_propagating_nodes.bool()[source]
            source, target, relation, sign = source[keep], target[keep], relation[keep], sign[keep]
        for index, name in enumerate(relation_names):
            if name in UNSIGNED_RELATIONS:
                sign = torch.where(relation == index, torch.ones_like(sign), sign)  # the learned gain carries the sign
        if normalisation not in NORMALISATIONS:
            raise ValueError(f"normalisation must be one of {NORMALISATIONS}, not {normalisation!r}")
        row = relation * num_graph_nodes + target
        if normalisation == "spectral":  # one global scale and no per-node divisor; see NORMALISATIONS
            radius = LinearResponseEncoder.spectral_radius_of_unsigned_aggregate(num_graph_nodes, source, target)
            values = sign / max(radius, 1.0)
        else:
            in_degree = torch.zeros(len(relation_names) * num_graph_nodes).index_add(0, row, torch.ones(len(row)))
            relations_feeding_node = (in_degree.view(len(relation_names), num_graph_nodes) > 0).sum(dim=0).clamp_min(1)
            exponent = 1.0 if normalisation == "in_degree" else float(normalisation_exponent)
            values = sign / (in_degree[row] ** exponent * relations_feeding_node[target])
            if exponent != 1.0:
                # row sums of sum_r |S_r| are in_degree^(1 - exponent) / relations_feeding and so can exceed one, which
                # would break the contraction; one global rescale by the largest of them restores it, and unlike the
                # "spectral" constant it divides an already softened matrix, so it stays near one
                largest_row_sum = float(torch.zeros(len(relation_names) * num_graph_nodes).index_add(0, row, values.abs()).max())
                values = values / max(largest_row_sum, 1.0)
        stacked = torch.sparse_coo_tensor(torch.stack([row, source]), values, (len(relation_names) * num_graph_nodes, num_graph_nodes)).coalesce()
        return stacked, relation_names

    def relation_gain(self) -> Tensor:
        """[R, D] gains: positive for signed relations, either sign for unsigned ones, all below the contraction bound."""
        positive = torch.sigmoid(self.gain_logit)
        either_sign = torch.tanh(self.gain_logit)
        return self.maximum_gain * torch.where(self.relation_is_unsigned[:, None], either_sign, positive)

    def sustained_input(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor) -> Tensor:
        """u: [B, N, C], sign * magnitude * input_weight at each perturbed node (index -1 is padding)."""
        batch_size = perturbation_node_index.shape[0]
        valid = (perturbation_node_index >= 0).to(perturbation_sign_and_magnitude.dtype)
        signed_magnitude = perturbation_sign_and_magnitude[:, :, 0] * perturbation_sign_and_magnitude[:, :, 1] * valid
        injected = signed_magnitude[:, :, None] * self.input_weight[None, None, :]
        field = torch.zeros(batch_size, self.num_graph_nodes, self.propagation_channels, device=injected.device, dtype=injected.dtype)
        index = perturbation_node_index.clamp_min(0)[:, :, None].expand(-1, -1, self.propagation_channels)
        return field.scatter_add(1, index, injected)

    def response(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, num_steps: int | None = None) -> Tensor:
        """The propagated response h after num_steps (default num_propagation_steps), [B, N, C], before expansion and gating."""
        sustained = self.sustained_input(perturbation_node_index, perturbation_sign_and_magnitude)
        batch_size = sustained.shape[0]
        gain = self.relation_gain()  # [R, C]
        num_relations = gain.shape[0]
        node_major_input = sustained.permute(1, 0, 2).contiguous()  # [N, B, C]
        state = node_major_input  # h(0) = u, so after k steps the response reaches k edges from the perturbed nodes
        for _ in range(self.num_propagation_steps if num_steps is None else num_steps):
            aggregated = torch.sparse.mm(self.stacked_adjacency, state.reshape(self.num_graph_nodes, batch_size * self.propagation_channels))
            relation_messages = (aggregated.view(num_relations, self.num_graph_nodes, batch_size, self.propagation_channels) * gain[:, None, None, :]).sum(dim=0)
            state = (1.0 - self.damping) * state + self.damping * (relation_messages + node_major_input)
        return state.permute(1, 0, 2)

    def forward(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:  # noqa: ARG002
        """The gated response field [B, N, D]; relation_adjacencies is accepted for the message passing signature and
        ignored, since the signed adjacency is built once from the edges given at construction."""
        gate = torch.sigmoid(self.output_gate(self.node_features))  # [N, D]
        response = self.response(perturbation_node_index, perturbation_sign_and_magnitude)
        if self.response_scale == "signed_log":
            response = torch.sign(response) * torch.log1p(response.abs() / torch.exp(self.log_response_scale))
        return (response @ self.channel_expansion) * gate[None, :, :]

    def perturbation_difference_field(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:  # noqa: ARG002
        """The response is linear in the input, so the unperturbed field is zero and the difference field is the field."""
        return self.forward(perturbation_node_index, perturbation_sign_and_magnitude)
