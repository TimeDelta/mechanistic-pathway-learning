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
"""
from __future__ import annotations

import numpy as np
import torch
from torch import Tensor, nn

from mechanistic_pathway_learning.models.soft_constraint_losses import anchored_link_prior_penalty

PROBABILITY_FLOOR = 1e-4
# the value gradient of a gated sparse product is gathered in chunks of at most this many floats (about 128 MB)
GRADIENT_CHUNK_ELEMENTS = 1 << 25


class EdgeGates(nn.Module):
    """One logit per gated edge; edge_values() is the [num_edges] multiplier (1 on ungated edges)."""

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
        self.logit = nn.Parameter(torch.logit(torch.as_tensor(probability, dtype=torch.float32)))

    def set_read_edges(self, edge_positions: Tensor) -> None:
        """Mark the gated edges among edge_positions (graph edge indices) as read by the encoder, the rest as not."""
        read = torch.zeros(self.num_edges, dtype=torch.bool, device=self.gated_edges.device)
        read[edge_positions.to(read.device).long()] = True
        self.read_by_encoder = read[self.gated_edges]

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


def relation_adjacency_parts(edge_index: Tensor, edge_relation_type: Tensor, num_graph_nodes: int, num_relation_types: int,
                             edge_sign: Tensor | None = None, edge_weight: Tensor | None = None) -> list[dict | None]:
    """Per relation, what the message-passing adjacency is built from, with the graph edge behind each value: the
    coalesced indices, each edge's slot among them, its value before gating (sign x weight / in-degree count, as
    RelationalMessagePassingEncoder.build_relation_adjacencies computes it) and its edge index."""
    parts: list[dict | None] = []
    for relation_index in range(num_relation_types):
        edges_of_relation = torch.nonzero(edge_relation_type == relation_index).flatten()
        if not len(edges_of_relation):
            parts.append(None)
            continue
        source, destination = edge_index[0][edges_of_relation], edge_index[1][edges_of_relation]
        in_degree = torch.zeros(num_graph_nodes, dtype=torch.float32).index_add(0, destination, torch.ones(len(destination)))
        values = 1.0 / in_degree[destination]
        if edge_sign is not None:
            values = values * edge_sign[edges_of_relation].to(values.dtype)
        if edge_weight is not None:
            values = values * edge_weight[edges_of_relation].to(values.dtype)
        key = destination * num_graph_nodes + source
        unique_keys, slot = torch.unique(key, return_inverse=True)  # sorted, which is the coalesced order (row-major)
        indices = torch.stack([torch.div(unique_keys, num_graph_nodes, rounding_mode="floor"), unique_keys % num_graph_nodes])
        parts.append({"indices": indices, "slot": slot, "values": values, "edges": edges_of_relation, "num_slots": len(unique_keys)})
    return parts


def gated_relation_adjacencies(parts: list[dict | None], edge_values: Tensor, num_graph_nodes: int) -> list[Tensor | None]:
    """The message-passing adjacencies with each edge's value multiplied by its gate; duplicates of one entry sum, as
    coalescing sums them."""
    adjacencies: list[Tensor | None] = []
    for part in parts:
        if part is None:
            adjacencies.append(None)
            continue
        device = edge_values.device
        gated = part["values"].to(device) * edge_values[part["edges"].to(device)]
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
