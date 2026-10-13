"""How the node descriptors enter an encoder, so that they describe the graph a perturbation spreads through without
standing in for the graph (user, 8 October 2026: keep the descriptors, but stop them overwhelming the graph).

The node features are the structural columns followed by num_descriptor_columns descriptor columns
(experiments/run_main_model.py, node_feature_matrix). Three treatments:

- plain: one linear map of all the columns, the form every run before 8 October 2026 used. With it the perturbed node
  reads its own descriptors: in the message passing encoder through its base state, which the ReLU makes interact
  with the injection, and in the linear-response encoder through its output gate, at the node where the response is
  largest. A readout can then learn "which protein is this" and skip the propagation, the shortcut the descriptors-only
  control measures (Gillis and Pavlidis, PLoS ONE 6, e17258, 2011: "high-quality gene function predictions can be made
  using data that possesses no information on which gene interacts with which").
- seed_masked: the perturbed nodes of each perturbation see their structural columns, and the expression columns of
  SEED_READS_THESE_DESCRIPTOR_PREFIXES, only (the other descriptor columns set to zero, the mean of their node type,
  since each block is standardised within its type); every other node keeps its descriptors. What a perturbation's own protein is reaches the readout only through the nodes around it. On a graph
  split into gene and protein nodes (docs/gene_protein_split.md) the mask also covers the nodes joined to a seed by an
  encodes edge: a knockout seeds the gene and a drug the protein, and without it the knocked-out gene's protein
  descriptors, or the drug target's gene expression, would sit one hop from the seed, where the merged graph hid them.
- zero_init_slow: the descriptor columns get their own linear map, initialised at zero and trained at their own
  learning rate (--descriptor-learning-rate). Adam moves a parameter by about one learning rate per step, so the
  descriptors can change a node's state by at most about that rate per step, and the runs stop early (best epochs 2 to
  11 on the slice with descriptors). A penalty on the descriptor weights was not chosen: under Adam the step does
  not scale with the gradient (Loshchilov and Hutter, ICLR 2019: L2 regularization and weight decay "are equivalent for
  standard stochastic gradient descent ... but as we demonstrate this is not the case for adaptive gradient algorithms,
  such as Adam"), and a penalty acts through the optimum, which early stopping does not reach.

With plain the encoders build exactly the parameters, buffers and random draws they built before, so runs and
checkpoints from before the treatments existed are unchanged.
"""
from __future__ import annotations

import torch
from torch import Tensor

DESCRIPTOR_TREATMENTS = ("plain", "seed_masked", "zero_init_slow")


def check_descriptor_treatment(descriptor_treatment: str, num_descriptor_columns: int, num_feature_columns: int) -> None:
    if descriptor_treatment not in DESCRIPTOR_TREATMENTS:
        raise ValueError(f"descriptor_treatment must be one of {DESCRIPTOR_TREATMENTS}, not {descriptor_treatment!r}")
    if descriptor_treatment != "plain" and not 0 < num_descriptor_columns < num_feature_columns:
        raise ValueError(f"descriptor treatment {descriptor_treatment!r} needs descriptor columns after at least one structural column "
                         f"(got {num_descriptor_columns} descriptor columns of {num_feature_columns})")


def features_without_descriptors(node_features: Tensor, num_descriptor_columns: int,
                                 descriptor_columns_kept: Tensor | None = None) -> Tensor:
    """The node features with the descriptor columns (the last num_descriptor_columns) set to zero.

    descriptor_columns_kept, a boolean over those columns, names the ones a perturbed node still reads
    (graph/node_descriptors.py, SEED_READS_THESE_DESCRIPTOR_PREFIXES: the regional and whole-tissue expression of the
    perturbed gene, by the owner's decision of 13 October 2026, because whether the gene is expressed in a region
    decides whether its loss matters there). None masks every descriptor column, which is what every run before that
    date did."""
    masked = node_features.clone()
    first = node_features.shape[1] - num_descriptor_columns
    masked[:, first:] = 0.0
    if descriptor_columns_kept is not None:
        if descriptor_columns_kept.shape[0] != num_descriptor_columns:
            raise ValueError(f"descriptor_columns_kept must have one entry per descriptor column "
                             f"({descriptor_columns_kept.shape[0]} for {num_descriptor_columns})")
        kept = descriptor_columns_kept.bool()
        masked[:, first:][:, kept] = node_features[:, first:][:, kept]
    return masked


def seed_node_mask(perturbation_node_index: Tensor, num_graph_nodes: int) -> Tensor:
    """[batch_size, num_graph_nodes] boolean, True at each perturbation's perturbed nodes (index -1 is padding).

    The mask is built by adding ones and testing for a positive count, so a padding entry that clamps onto node 0 never
    overwrites a real seed at node 0."""
    valid = (perturbation_node_index >= 0).to(torch.float32)
    counts = torch.zeros(perturbation_node_index.shape[0], num_graph_nodes, device=perturbation_node_index.device)
    return counts.scatter_add(1, perturbation_node_index.clamp_min(0), valid) > 0


def seed_unit_mask(perturbation_node_index: Tensor, num_graph_nodes: int, partner_index: Tensor | None = None) -> Tensor:
    """seed_node_mask, widened to the partners of each seed: partner_index is [2, num_pairs] of (node, partner) pairs,
    both directions given (the encodes edges of a split graph). None, or no pairs, gives seed_node_mask unchanged."""
    seeds = seed_node_mask(perturbation_node_index, num_graph_nodes)
    if partner_index is None or partner_index.numel() == 0:
        return seeds
    reached = torch.zeros(seeds.shape, device=seeds.device).index_add_(1, partner_index[1], seeds[:, partner_index[0]].to(torch.float32))
    return seeds | (reached > 0)
