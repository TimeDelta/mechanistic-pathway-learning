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
assumption A9 in another form.

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

The system is linear in the input, so the unperturbed response is zero and the field is a difference field by
construction: no base state, and no node identity, enters it. Node kinds enter only through an output gate,
field = (h W_expand) * sigmoid(node_features W + b), which lets a pooled readout tell where the response landed
without adding anything where it did not.
"""
from __future__ import annotations

import torch
from torch import Tensor, nn

SUBSTRATE_RELATION = "substrate_of"
DEPLETES_SUBSTRATE_RELATION = "depletes_substrate"
UNSIGNED_RELATIONS = ("binds",)


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
        num_propagation_steps: int = 8,
        propagation_channels: int = 4,
        damping: float = 0.5,
        contraction: float = 0.9,
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
                                                                          non_propagating_nodes)
        self.relation_names = relation_names
        self.register_buffer("stacked_adjacency", stacked_adjacency, persistent=False)
        self.register_buffer("node_features", node_features.to(torch.float32), persistent=False)
        self.register_buffer("relation_is_unsigned", torch.tensor([name in UNSIGNED_RELATIONS for name in relation_names]), persistent=False)
        self.maximum_gain = contraction
        self.propagation_channels = propagation_channels
        self.gain_logit = nn.Parameter(torch.randn(len(relation_names), propagation_channels))  # spread so channels start at different time scales
        self.input_weight = nn.Parameter(torch.randn(propagation_channels) / propagation_channels**0.5)
        self.channel_expansion = nn.Parameter(torch.randn(propagation_channels, node_state_dim) / propagation_channels**0.5)
        self.output_gate = nn.Linear(node_features.shape[1], node_state_dim)

    @staticmethod
    def signed_stacked_adjacency(num_graph_nodes: int, relation_types: list[str], edge_source: Tensor, edge_target: Tensor,
                                 edge_relation: Tensor, edge_sign: Tensor, non_propagating_nodes: Tensor | None = None) -> tuple[Tensor, list[str]]:
        """All relations' signed, normalised adjacencies stacked into one sparse [R * N, N] matrix, with the derived
        depletes_substrate relation appended and edges leaving a non-propagating node dropped; returns the matrix and
        the relation names in stack order."""
        sources, targets, relations, signs = [edge_source.long()], [edge_target.long()], [edge_relation.long()], [edge_sign.float()]
        relation_names = list(relation_types)
        if SUBSTRATE_RELATION in relation_types:
            substrate_edges = edge_relation == relation_types.index(SUBSTRATE_RELATION)
            relation_names.append(DEPLETES_SUBSTRATE_RELATION)
            sources.append(edge_target[substrate_edges].long())
            targets.append(edge_source[substrate_edges].long())
            relations.append(torch.full((int(substrate_edges.sum()),), len(relation_names) - 1, dtype=torch.long))
            signs.append(-torch.ones(int(substrate_edges.sum())))
        source, target, relation, sign = torch.cat(sources), torch.cat(targets), torch.cat(relations), torch.cat(signs)
        if non_propagating_nodes is not None:
            keep = ~non_propagating_nodes.bool()[source]
            source, target, relation, sign = source[keep], target[keep], relation[keep], sign[keep]
        for index, name in enumerate(relation_names):
            if name in UNSIGNED_RELATIONS:
                sign = torch.where(relation == index, torch.ones_like(sign), sign)  # the learned gain carries the sign
        row = relation * num_graph_nodes + target
        in_degree = torch.zeros(len(relation_names) * num_graph_nodes).index_add(0, row, torch.ones(len(row)))
        relations_feeding_node = (in_degree.view(len(relation_names), num_graph_nodes) > 0).sum(dim=0).clamp_min(1)
        values = sign / (in_degree[row] * relations_feeding_node[target])
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
        return (self.response(perturbation_node_index, perturbation_sign_and_magnitude) @ self.channel_expansion) * gate[None, :, :]

    def perturbation_difference_field(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:  # noqa: ARG002
        """The response is linear in the input, so the unperturbed field is zero and the difference field is the field."""
        return self.forward(perturbation_node_index, perturbation_sign_and_magnitude)
