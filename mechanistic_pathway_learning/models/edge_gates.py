"""Learnable edge-presence gates pulled toward an evidence prior (decision 9.57, the owner, 16 October 2026).

Each gated edge e carries a gate g_e = sigmoid(theta_e) in (0, 1) that multiplies the edge's weight wherever an
encoder reads it, initialised at the edge's prior probability p_e (graph/edge_evidence.py), and the training loss adds

    penalty = sum_e n_e * KL(Bernoulli(p_e) || Bernoulli(g_e)) / W

where n_e is the edge's pseudo-count (its tier's unit strength times its tempered record count, times the global
strength) and W is the total loss weight of the labelled pairs the fit trains on, so that one pseudo-count weighs as
much as one unit of label weight. This is soft_constraint_losses.anchored_link_prior_penalty applied to edges. A
well-referenced curated edge has a large n_e and the data must pay to close it; a text-mined edge has a small one and
moves easily. A gate never reaches 0 or 1, so no edge is deleted and the sparsity pattern of every adjacency is fixed.

What a gate is not: a transmission magnitude. The rulings of 12 and 14 October (the occupancy placeholder, the
threshold grid, 9.20) stand: the gate says how sure the sources are that the edge is there, read as the expected
presence of an uncertain edge, and the edge's own weight (a conductance share, a mechanism magnitude) is untouched by
the prior. In the linear-response encoder the weight enters as a share at its target, so gates of one tier at one node
cancel and only relative certainty among a node's inputs matters; in the message-passing encoder the divisor is a
count, so a gate scales its edge's message.

Ungated edges (graph/edge_evidence.py: Human-GEM's metabolic edges, the definitional layers) have gate 1 always.

Hard gates (decision 9.62 (e), the owner, 16 October 2026: "Yes the prior probability should be the cost per edge
absolutely. That's the only valid way. It acting like KLD to the Bernoulli proper is good."). HardConcreteEdgeGates
asks whether an edge is there, a yes-or-no fact, and lets a gate be exactly 0 or exactly 1:

  - the prior that edge e is real is pi_e = sigmoid(n_e logit p_e) (binary_edge_prior: n = 1 gives p, more independent
    evidence at the tier's strength moves it further from one half; the owner: "Sigmoid option for the pseudo count
    point you raised. It's closer to the physiology");
  - the stochastic pass samples hard concrete gates (Louizos, Welling and Kingma 2018, arXiv 1712.01312; gamma -0.1,
    zeta 1.1, beta 2/3, their values), one sample per batch as they did ("We used a single sample of the gate z for each
    minibatch of datapoints during the optimization"), started at P(open) = pi_e, and the penalty is
    KL(Bernoulli(q_e) || Bernoulli(pi_e)) with q_e = P(z_e != 0), their eq. 12: a cost of log((1 - pi)/pi) nats per
    open edge less the gate's entropy;
  - the read-out is binary: an edge is open when q_e >= 1/2;
  - the binary pass (the owner: "Why not include that open probability threshold as a second training pass over just
    that parameter or something? Can be given a proxy gradient that resembles the cliff discontinuity") trains the gates
    alone with the read-out in the forward pass. Its proxy gradient is the cliff itself: the first-order change in the
    loss when the edge flips, with every target it feeds renormalised exactly (flip_effect_estimates), plus the exact
    change in the prior term, logit(pi_e) per closure. A local slope would not do: both encoders read a gate as a share
    of its target's gated in-degree, so that a closed gate is a deleted edge, and a gate that is its target's only input
    then has a slope of zero at every value above 0 (2,105 gated edges on the v5 graph in linear response, 8,325 within
    their relation in message passing), while the deletion it stands for is the largest change there is.
"""
from __future__ import annotations

import numpy as np
import torch
from torch import Tensor, nn

from mechanistic_pathway_learning.models.soft_constraint_losses import anchored_link_prior_penalty

PROBABILITY_FLOOR = 1e-4
# the value gradient of a gated sparse product is gathered in chunks of at most this many floats (about 128 MB)
GRADIENT_CHUNK_ELEMENTS = 1 << 25
# hard concrete (Louizos et al. 2018, section 4: "For all of our experiments we set gamma = -0.1, zeta = 1.1" and
# "beta = 2/3 for the concrete distributions")
HARD_CONCRETE_GAMMA, HARD_CONCRETE_ZETA, HARD_CONCRETE_BETA = -0.1, 1.1, 2.0 / 3.0
EDGE_GATE_KINDS = ("soft", "hard_concrete")
EDGE_GATE_PHASES = ("stochastic", "binary")


def binary_edge_prior(prior_probability, pseudo_count):
    """pi = sigmoid(n logit p): the prior that an edge is real as one yes-or-no fact (decision 9.62, sub-choice (iii)).
    A Beta pseudo-count acts on one shared binary fact only through its mean, so the evidence strength enters here,
    as the number of independent observations at the tier's strength: n = 1 gives p."""
    probability = np.clip(np.asarray(prior_probability, dtype=float), PROBABILITY_FLOOR, 1.0 - PROBABILITY_FLOOR)
    count = np.asarray(pseudo_count, dtype=float)
    log_odds = count * np.log(probability / (1.0 - probability))
    return np.clip(1.0 / (1.0 + np.exp(-log_odds)), PROBABILITY_FLOOR, 1.0 - PROBABILITY_FLOOR)


class GateLayout(nn.Module):
    """Which graph edges are gated, their evidence priors, and which of them the encoder reads; the buffers both kinds
    of gate share, under the names the soft gates have always saved."""

    def __init__(self, gated: np.ndarray, prior_probability: np.ndarray, pseudo_count: np.ndarray):
        super().__init__()
        gated = np.asarray(gated, dtype=bool)
        positions = np.flatnonzero(gated)
        probability = np.clip(np.asarray(prior_probability, dtype=float)[positions], PROBABILITY_FLOOR, 1.0 - PROBABILITY_FLOOR)
        count = np.asarray(pseudo_count, dtype=float)[positions]
        if np.any(count < 0) or not np.all(np.isfinite(probability)):
            raise ValueError("a gated edge needs a finite prior probability and a pseudo-count >= 0")
        gate_position = np.full(len(gated), -1, dtype=np.int64)
        gate_position[positions] = np.arange(len(positions))
        self.num_edges = len(gated)
        self.register_buffer("gate_position", torch.as_tensor(gate_position))
        self.register_buffer("gated_edges", torch.as_tensor(positions, dtype=torch.long))
        self.register_buffer("prior_probability", torch.as_tensor(probability, dtype=torch.float32))
        self.register_buffer("pseudo_count", torch.as_tensor(count, dtype=torch.float32))
        # which gated edges the encoder reads (set_read_edges); a gate it does not read has no data gradient and stays at
        # its prior, so the summary reports the read ones
        self.register_buffer("read_by_encoder", torch.ones(len(positions), dtype=torch.bool))

    def set_read_edges(self, edge_positions: Tensor) -> None:
        """Mark the gated edges among edge_positions (graph edge indices) as read by the encoder, the rest as not."""
        read = torch.zeros(self.num_edges, dtype=torch.bool, device=self.gated_edges.device)
        read[edge_positions.to(read.device).long()] = True
        self.read_by_encoder = read[self.gated_edges]

    def binary_state_version(self) -> int | None:
        """None while the gates are not exactly 0 or 1; soft gates never are."""
        return None


class EdgeGates(GateLayout):
    """One logit per gated edge; edge_values() is the [num_edges] multiplier (1 on ungated edges)."""

    def __init__(self, gated: np.ndarray, prior_probability: np.ndarray, pseudo_count: np.ndarray):
        super().__init__(gated, prior_probability, pseudo_count)
        self.logit = nn.Parameter(torch.logit(self.prior_probability.clone()))

    def gate(self) -> Tensor:
        """[num_gated] the gates themselves."""
        return torch.sigmoid(self.logit)

    def edge_values(self) -> Tensor:
        """[num_edges]: the gate of a gated edge, 1 for every other edge."""
        values = torch.ones(self.num_edges, dtype=self.logit.dtype, device=self.logit.device)
        return values.index_copy(0, self.gated_edges, self.gate())

    def kl_penalty(self) -> Tensor:
        """sum_e n_e KL(Bernoulli(p_e) || Bernoulli(g_e)); 0 at initialisation."""
        return anchored_link_prior_penalty(self.gate(), self.prior_probability, self.pseudo_count)  # (model, prior, count): KL(prior || gate)

    def values_per_edge(self) -> np.ndarray:
        """[num_edges] the gate of each gated edge, NaN for an ungated one (what a run saves beside its predictions)."""
        values = np.full(self.num_edges, np.nan, dtype=np.float32)
        with torch.no_grad():
            values[self.gated_edges.cpu().numpy()] = self.gate().cpu().numpy()
        return values

    def summary(self) -> dict:
        """Over the gates the encoder reads: how far the data moved them from their priors."""
        with torch.no_grad():
            read = self.read_by_encoder
            gate, prior = self.gate()[read], self.prior_probability[read]
            moved = (gate - prior).abs()
            some = bool(gate.numel())
            return {"gated_edges": int(self.gated_edges.numel()), "gated_edges_the_encoder_reads": int(read.sum()),
                    "mean_prior": float(prior.mean()) if some else None, "mean_gate": float(gate.mean()) if some else None,
                    "mean_absolute_move": float(moved.mean()) if some else None, "largest_move": float(moved.max()) if some else None,
                    "moved_more_than_0_1": int((moved > 0.1).sum()), "closed_below_0_5": int((gate < 0.5).sum()),
                    "prior_below_0_5": int((prior < 0.5).sum()), "kl_penalty": float(self.kl_penalty())}


class HardConcreteEdgeGates(GateLayout):
    """Hard concrete gates over the gated edges (decision 9.62 (e); the module docstring gives the rules).

    phase "stochastic": in training mode edge_values() is the current sample, drawn by resample() once per batch; in
    evaluation mode, and in phase "binary" always, it is the read-out, exactly 0 or 1. kl_penalty() is the KL divergence
    of each read gate's open probability to its prior in the stochastic phase and, in the binary phase, -log
    Bernoulli(z; pi) of the read-out, the same per-edge cost with no entropy term."""

    def __init__(self, gated: np.ndarray, prior_probability: np.ndarray, pseudo_count: np.ndarray,
                 beta: float = HARD_CONCRETE_BETA, gamma: float = HARD_CONCRETE_GAMMA, zeta: float = HARD_CONCRETE_ZETA):
        super().__init__(gated, prior_probability, pseudo_count)
        if not (gamma < 0.0 and zeta > 1.0 and beta > 0.0):
            raise ValueError("hard concrete needs gamma < 0, zeta > 1 and beta > 0")
        self.beta, self.gamma, self.zeta = float(beta), float(gamma), float(zeta)
        prior_open = binary_edge_prior(self.prior_probability.numpy(), self.pseudo_count.numpy())
        self.register_buffer("prior_open", torch.as_tensor(prior_open, dtype=torch.float32))
        # started at P(open) = prior: log alpha = logit(prior) + beta log(-gamma / zeta), from eq. 12
        self.log_alpha = nn.Parameter(torch.logit(self.prior_open.clone()) + self.beta * float(np.log(-self.gamma / self.zeta)))
        self.phase = "stochastic"
        self.sample: Tensor | None = None

    def set_phase(self, phase: str) -> None:
        if phase not in EDGE_GATE_PHASES:
            raise ValueError(f"phase must be one of {EDGE_GATE_PHASES}")
        self.phase = phase
        self.sample = None

    def open_log_odds(self) -> Tensor:
        """[num_gated] m = logit P(z != 0) = log alpha - beta log(-gamma / zeta) (eq. 12); open at the read-out when m >= 0."""
        return self.log_alpha - self.beta * float(np.log(-self.gamma / self.zeta))

    def open_probability(self) -> Tensor:
        return torch.sigmoid(self.open_log_odds())

    def read_out(self) -> Tensor:
        """[num_gated] 1.0 where the gate is open (P(open) >= 1/2), 0.0 where it is closed, with no gradient."""
        return (self.open_log_odds().detach() >= 0.0).to(self.log_alpha.dtype)

    def resample(self, generator: torch.Generator | None = None) -> Tensor:
        """Draw one hard concrete sample per gate (call once per batch); returns it and keeps it for edge_values()."""
        uniform = torch.rand(self.log_alpha.shape, generator=generator, device=self.log_alpha.device if generator is None else "cpu").to(self.log_alpha.device)
        uniform = uniform.clamp(1e-6, 1.0 - 1e-6)
        stretched = torch.sigmoid((torch.log(uniform) - torch.log1p(-uniform) + self.log_alpha) / self.beta) * (self.zeta - self.gamma) + self.gamma
        self.sample = stretched.clamp(0.0, 1.0)
        return self.sample

    def gate(self) -> Tensor:
        """[num_gated] what the encoder reads now: the sample while training the stochastic phase, else the read-out."""
        if self.phase == "stochastic" and self.training:
            if self.sample is None:
                self.resample()
            return self.sample
        return self.read_out()

    def edge_values(self) -> Tensor:
        """[num_edges]: the gate of a gated edge, 1 for every other edge."""
        values = torch.ones(self.num_edges, dtype=self.log_alpha.dtype, device=self.log_alpha.device)
        return values.index_copy(0, self.gated_edges, self.gate())

    def open_edges(self) -> Tensor:
        """[num_edges] bool: the read-out, with every ungated edge open."""
        open_edge = torch.ones(self.num_edges, dtype=torch.bool, device=self.log_alpha.device)
        return open_edge.index_copy(0, self.gated_edges, self.read_out() > 0.5)

    def binary_state_version(self) -> int | None:
        """A key that changes whenever the read-out does, while the gates are binary (the binary phase, or evaluation);
        None while a stochastic sample is what the encoder reads."""
        if self.phase == "stochastic" and self.training:
            return None
        state = (self.read_out() > 0.5).cpu().numpy()
        return hash(np.packbits(state).tobytes())

    def kl_penalty(self) -> Tensor:
        """Over the gates the encoder reads: the stochastic phase's sum of KL(Bernoulli(q) || Bernoulli(pi)), q = P(open);
        the binary phase's sum of -log Bernoulli(z; pi) at the read-out z, measured from the read-out at the prior so
        that it is 0 when every gate sits on the side of one half its prior is on."""
        read = self.read_by_encoder
        prior = self.prior_open[read]
        if self.phase == "binary":
            state = self.read_out()[read]
            at_prior = (prior >= 0.5).to(prior.dtype)
            cost = lambda z: -(z * torch.log(prior) + (1.0 - z) * torch.log1p(-prior))  # noqa: E731
            return (cost(state) - cost(at_prior)).sum()
        q = self.open_probability()[read].clamp(PROBABILITY_FLOOR, 1.0 - PROBABILITY_FLOOR)
        return (q * torch.log(q / prior) + (1.0 - q) * torch.log((1.0 - q) / (1.0 - prior))).sum()

    def flip_prior_change(self) -> Tensor:
        """[num_edges] the exact change in the prior term if each gated edge flipped from its read-out: +logit(pi) to close
        an open edge, -logit(pi) to open a closed one (0 for an ungated edge)."""
        state = self.read_out()
        logit_prior = torch.log(self.prior_open / (1.0 - self.prior_open))
        change = torch.where(state > 0.5, logit_prior, -logit_prior)
        values = torch.zeros(self.num_edges, dtype=change.dtype, device=change.device)
        return values.index_copy(0, self.gated_edges, change)

    def apply_flip_gradient(self, flip_objective_change: Tensor, surrogate_width: float) -> None:
        """The binary pass's proxy gradient: with D_e the estimated change in the objective if edge e flipped (data and
        prior, already divided by the penalty denominator), the finite difference along z is J(1) - J(0) = -D_e for an
        open edge and +D_e for a closed one, and the chain rule through a steep sigmoid at the threshold gives
        dJ/dm = (J(1) - J(0)) sigmoid'(m / T) / T. It is added to log_alpha.grad (m and log alpha differ by a constant)."""
        change = flip_objective_change[self.gated_edges].to(self.log_alpha.dtype)
        state = self.read_out()
        finite_difference = torch.where(state > 0.5, -change, change)
        slope = torch.sigmoid(self.open_log_odds().detach() / surrogate_width)
        surrogate = slope * (1.0 - slope) / surrogate_width
        gradient = finite_difference * surrogate
        self.log_alpha.grad = gradient if self.log_alpha.grad is None else self.log_alpha.grad + gradient

    def values_per_edge(self) -> np.ndarray:
        """[num_edges] each gated edge's open probability, NaN for an ungated one."""
        values = np.full(self.num_edges, np.nan, dtype=np.float32)
        with torch.no_grad():
            values[self.gated_edges.cpu().numpy()] = self.open_probability().cpu().numpy()
        return values

    def summary(self) -> dict:
        with torch.no_grad():
            read = self.read_by_encoder
            q, prior, state = self.open_probability()[read], self.prior_open[read], self.read_out()[read]
            some = bool(q.numel())
            closed_against_prior = int(((state < 0.5) & (prior >= 0.5)).sum())
            opened_against_prior = int(((state > 0.5) & (prior < 0.5)).sum())
            return {"kind": "hard_concrete", "phase": self.phase, "gated_edges": int(self.gated_edges.numel()), "gated_edges_the_encoder_reads": int(read.sum()),
                    "mean_prior": float(prior.mean()) if some else None, "mean_open_probability": float(q.mean()) if some else None,
                    "open_at_read_out": int((state > 0.5).sum()), "closed_at_read_out": int((state < 0.5).sum()),
                    "closed_against_prior": closed_against_prior, "opened_against_prior": opened_against_prior,
                    "prior_below_0_5": int((prior < 0.5).sum()), "penalty": float(self.kl_penalty())}


def flip_effect_estimates(group: Tensor, origin: Tensor, numerator: Tensor, denominator_weight: Tensor, value_gradient: Tensor,
                          edge_state: Tensor, num_edges: int) -> Tensor:
    """[num_edges] the first-order change in the loss if each edge flipped from edge_state, its target renormalised.

    The adjacency entries k (one or more per graph edge, origin[k] the edge) take values v_k = a_k z_k / D_g, with a_k
    the numerator (sign times weight), z_k the state of the entry's edge, g = group[k] the entry's normalisation group
    (its target in linear response, its relation and target in message passing) and D_g = sum over the group of c_j z_j
    (c the denominator weight: the edge weight in linear response, 1 in message passing). value_gradient[k] = dL/dv_k.

    Closing an open edge e sets its entries to 0 and scales every other entry of each group g it reaches by
    D_g / (D_g - W_eg), W_eg the weight e brings to g; opening a closed edge sets its entries to a_k / (D_g + W_eg) and
    scales the others by D_g / (D_g + W_eg). The value changes are exact; the loss is extrapolated to first order in them,
    which is where the estimate can be wrong (and why the binary pass ends with exact checks). A group that loses its
    only input loses everything, which the slope of v in z (zero there) would not have shown."""
    device = value_gradient.device
    numerator, denominator_weight = numerator.to(device), denominator_weight.to(device)
    state = edge_state.to(device=device, dtype=numerator.dtype)[origin]
    num_groups = int(group.max()) + 1 if group.numel() else 0
    denominator = torch.zeros(num_groups, dtype=numerator.dtype, device=device).index_add(0, group, denominator_weight * state)
    safe = denominator[group].clamp_min(torch.finfo(numerator.dtype).tiny)
    value = numerator * state / safe
    gradient_times_value = value_gradient * value
    group_sum = torch.zeros(num_groups, dtype=numerator.dtype, device=device).index_add(0, group, gradient_times_value)
    # aggregate each edge's entries per group it reaches
    pair_key = origin.to(device) * num_groups + group
    pairs, pair_of_entry = torch.unique(pair_key, return_inverse=True)
    pair_edge, pair_group = torch.div(pairs, num_groups, rounding_mode="floor"), pairs % num_groups
    weight = torch.zeros(len(pairs), dtype=numerator.dtype, device=device).index_add(0, pair_of_entry, denominator_weight)
    own_gradient_times_value = torch.zeros(len(pairs), dtype=numerator.dtype, device=device).index_add(0, pair_of_entry, gradient_times_value)
    own_gradient_times_numerator = torch.zeros(len(pairs), dtype=numerator.dtype, device=device).index_add(0, pair_of_entry, value_gradient * numerator)
    group_denominator, group_total = denominator[pair_group], group_sum[pair_group]
    is_open = edge_state.to(device=device)[pair_edge] > 0.5
    remaining = group_denominator - weight
    closing = -own_gradient_times_value + torch.where(remaining > 0, weight / remaining.clamp_min(torch.finfo(numerator.dtype).tiny), torch.zeros_like(weight)) \
        * (group_total - own_gradient_times_value)
    enlarged = group_denominator + weight
    opening = own_gradient_times_numerator / enlarged.clamp_min(torch.finfo(numerator.dtype).tiny) - weight / enlarged.clamp_min(torch.finfo(numerator.dtype).tiny) * group_total
    per_pair = torch.where(is_open, closing, opening)
    return torch.zeros(num_edges, dtype=numerator.dtype, device=device).index_add(0, pair_edge, per_pair)


def relation_adjacency_parts(edge_index: Tensor, edge_relation_type: Tensor, num_graph_nodes: int, num_relation_types: int,
                             edge_sign: Tensor | None = None, edge_weight: Tensor | None = None) -> list[dict | None]:
    """Per relation, what the message-passing adjacency is built from, with the graph edge behind each value: the
    coalesced indices, each edge's slot among them, its value before gating (sign x weight / in-degree count, as
    RelationalMessagePassingEncoder.build_relation_adjacencies computes it), its numerator (sign x weight), its
    destination and its edge index. The numerator and destination let a gated encoder divide by the gated in-degree,
    the sum of the gates of the relation's edges into the destination (decision 9.62: a closed gate is a deleted edge)."""
    parts: list[dict | None] = []
    for relation_index in range(num_relation_types):
        edges_of_relation = torch.nonzero(edge_relation_type == relation_index).flatten()
        if not len(edges_of_relation):
            parts.append(None)
            continue
        source, destination = edge_index[0][edges_of_relation], edge_index[1][edges_of_relation]
        in_degree = torch.zeros(num_graph_nodes, dtype=torch.float32).index_add(0, destination, torch.ones(len(destination)))
        numerator = torch.ones(len(destination), dtype=torch.float32)
        if edge_sign is not None:
            numerator = numerator * edge_sign[edges_of_relation].to(numerator.dtype)
        if edge_weight is not None:
            numerator = numerator * edge_weight[edges_of_relation].to(numerator.dtype)
        values = numerator / in_degree[destination]
        key = destination * num_graph_nodes + source
        unique_keys, slot = torch.unique(key, return_inverse=True)  # sorted, which is the coalesced order (row-major)
        indices = torch.stack([torch.div(unique_keys, num_graph_nodes, rounding_mode="floor"), unique_keys % num_graph_nodes])
        parts.append({"indices": indices, "slot": slot, "values": values, "numerator": numerator, "destination": destination,
                      "edges": edges_of_relation, "num_slots": len(unique_keys)})
    return parts


def gated_relation_adjacencies(parts: list[dict | None], edge_values: Tensor, num_graph_nodes: int) -> list[Tensor | None]:
    """The message-passing adjacencies with each edge's sign x weight x gate over the gated in-degree of its relation at
    its destination (decision 9.62: with gates of 0 and 1 this is the adjacency of the graph without the closed edges);
    duplicates of one entry sum, as coalescing sums them. Until 16 October 2026 the divisor was the ungated count, so a
    gate scaled its edge's message instead of taking a share."""
    adjacencies: list[Tensor | None] = []
    for part in parts:
        if part is None:
            adjacencies.append(None)
            continue
        device = edge_values.device
        gate = edge_values[part["edges"].to(device)]
        destination = part["destination"].to(device)
        gated_in_degree = torch.zeros(num_graph_nodes, dtype=gate.dtype, device=device).index_add(0, destination, gate)
        gated = part["numerator"].to(device) * gate / gated_in_degree[destination].clamp_min(torch.finfo(gate.dtype).tiny)
        values = torch.zeros(part["num_slots"], dtype=gated.dtype, device=device).index_add(0, part["slot"].to(device), gated)
        adjacencies.append(torch.sparse_coo_tensor(part["indices"].to(device), values, (num_graph_nodes, num_graph_nodes), is_coalesced=True))
    return adjacencies


class SparseValuesProduct(torch.autograd.Function):
    """sparse(indices, values) @ dense, differentiable in the values without a dense intermediate.

    Why it exists. PyTorch's backward of torch.sparse.mm with respect to a sparse COO matrix forms the dense product
    grad @ dense.T, [rows, columns], and then masks it to the matrix's pattern. With the adjacency's values fixed (every
    run without edge gates) no gradient is asked of them and nothing dense is formed; with gated values it is, and on
    graph_full_neuronal_split_binders_drugs_v4 one linear-response step asked for 15.8 GB, a [K, N] float matrix per
    propagation round. Here the gradient of each stored value is gathered directly,
    d loss / d value_e = sum_m grad[row_e, m] * dense[column_e, m], in chunks of entries, which is O(nnz * m) arithmetic
    and memory bounded by GRADIENT_CHUNK_ELEMENTS. The gradient with respect to the dense input is the transposed matrix
    times grad, over a transposed layout computed once by sparse_layout."""

    @staticmethod
    def forward(ctx, values: Tensor, dense: Tensor, indices: Tensor, transposed_indices: Tensor, transposed_order: Tensor, num_rows: int) -> Tensor:
        matrix = torch.sparse_coo_tensor(indices, values, (num_rows, dense.shape[0]), is_coalesced=True)
        ctx.save_for_backward(values, dense, indices, transposed_indices, transposed_order)
        ctx.num_rows = num_rows
        return torch.sparse.mm(matrix, dense)

    @staticmethod
    def backward(ctx, grad_output: Tensor):
        values, dense, indices, transposed_indices, transposed_order = ctx.saved_tensors
        grad_values = grad_dense = None
        grad_output = grad_output.contiguous()
        if ctx.needs_input_grad[1]:
            transposed = torch.sparse_coo_tensor(transposed_indices, values[transposed_order], (dense.shape[0], ctx.num_rows), is_coalesced=True)
            grad_dense = torch.sparse.mm(transposed, grad_output)
        if ctx.needs_input_grad[0]:
            grad_values = torch.empty_like(values)
            rows, columns = indices[0], indices[1]
            chunk = max(1, GRADIENT_CHUNK_ELEMENTS // max(1, dense.shape[1]))
            for start in range(0, values.numel(), chunk):
                end = start + chunk
                grad_values[start:end] = (grad_output[rows[start:end]] * dense[columns[start:end]]).sum(dim=1)
        return grad_values, grad_dense, None, None, None, None


def sparse_layout(indices: Tensor, num_rows: int) -> tuple[Tensor, Tensor]:
    """(transposed indices, order) of a coalesced pattern: the transposed pattern in coalesced (row-major) order and
    the permutation that carries the values there, for SparseValuesProduct."""
    order = torch.argsort(indices[1] * num_rows + indices[0])
    return torch.stack([indices[1][order], indices[0][order]]), order


def sparse_values_product(indices: Tensor, values: Tensor, dense: Tensor, num_rows: int, transposed_indices: Tensor, transposed_order: Tensor) -> Tensor:
    """[num_rows, m] = sparse(indices, values) @ dense, with the value gradient gathered per entry (SparseValuesProduct)."""
    return SparseValuesProduct.apply(values, dense, indices, transposed_indices, transposed_order, num_rows)


def penalty_denominator(pair_weights: np.ndarray, fit_rows: np.ndarray) -> float:
    """W: the total loss weight of the pairs the fit trains on, so a pseudo-count weighs as one unit of label weight."""
    total = float(np.asarray(pair_weights, dtype=float)[np.asarray(fit_rows, dtype=int)].sum())
    return max(total, 1.0)
