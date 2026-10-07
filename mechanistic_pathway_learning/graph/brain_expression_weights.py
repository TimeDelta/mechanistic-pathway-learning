"""Per-node cell-class weights for the cell-class propagation channels of the linear-response encoder.

Cell types enter propagation without copying the graph (the user's direction of 7 October 2026: "separate edge types
for different cell types"). Every cell class gets its own propagation channel over the same edges, and in class k
every node's response, and the perturbation's input at the perturbed node, is multiplied by the node's weight w_k at
each step. A response arriving at node t is so scaled by w_k(t), and along a path by the product of the weights of
the nodes it passes through: a reaction a cell class does not express neither holds nor passes the response in that
class's channel. This is the same as giving each class its own copy of every relation with edge weights set by the
destination's expression, but stored as one weight per node and class, so the graph keeps one adjacency and a class
costs one more column in the sparse product (LinearResponseEncoder, cell_class_weights).

The weight of a gene in class k is its expression relative to its highest class, with a half-saturation at the Human
Protein Atlas detection threshold of 1 nCPM:

    w_k(g) = x_k(g) / (max_j x_j(g) + 1)

so a gene's own top class gets nearly 1 (TH in dopaminergic neurons: 388 / 389), a class at a fiftieth of the top gets
about 0.02 (MBP in dopaminergic neurons, where its 112 nCPM is likely ambient RNA against 5,606 in oligodendrocytes), and
a gene at a few nCPM everywhere stays below 1 everywhere (2 nCPM in every class gives 2 / 3). The relative form follows
the specificity of a gene across classes rather than its absolute level, because pseudobulk single-nucleus profiles
carry ambient RNA that gives most genes a few nCPM in every class.

Reactions take their class expression from the Human-GEM gene rule (minimum over "and", maximum over "or", as for the
descriptors: brain_expression_descriptors.evaluate_gene_rule) on the raw nCPM, then the same relative weight. Nodes
that carry no expression get weight 1 in every class, so a class channel does not block what nothing is known about:
metabolites (shared by every cell that holds them), the Reactome protein entities, the membrane potential, genes absent
from the expression tables and reactions whose rule names no known gene (transport, exchange and spontaneous
reactions). An all_cells class of weight 1 everywhere can be added, which gives the original propagation as one of the
channels, so a model with cell classes contains the model without them.

Extracellular metabolites (Human-GEM compartment "e") can be marked as one pool shared by every class
(extracellular_pool_nodes): a transmitter released by one class then reaches the receptors of every other class.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.brain_expression_descriptors import evaluate_gene_rule, gene_node_ensembl_ids

DETECTION_THRESHOLD_NCPM = 1.0
ALL_CELLS_CLASS = "all_cells"
EXTRACELLULAR_COMPARTMENT = "e"


def relative_expression_weights(class_expression: pd.DataFrame, half_saturation: float = DETECTION_THRESHOLD_NCPM) -> pd.DataFrame:
    """x_k / (max_j x_j + half_saturation) per row; rows with no value in any class stay missing."""
    values = class_expression.astype(float).clip(lower=0.0)
    return values.div(values.max(axis=1) + half_saturation, axis=0)


def node_cell_class_weights(nodes: pd.DataFrame, class_expression: pd.DataFrame, symbol_to_ensembl: dict[str, str] | None = None,
                            half_saturation: float = DETECTION_THRESHOLD_NCPM, include_all_cells_class: bool = True) -> pd.DataFrame:
    """node_id-indexed weights in [0, 1], one column per class of class_expression (an Ensembl-indexed nCPM table), plus
    ALL_CELLS_CLASS when include_all_cells_class; see the module docstring for which nodes get what."""
    classes = list(class_expression.columns)
    weights = pd.DataFrame(1.0, index=pd.Index(nodes.node_id, name="node_id"), columns=classes)

    is_gene = (nodes.node_type == "gene").to_numpy()
    gene_ensembl = gene_node_ensembl_ids(nodes, symbol_to_ensembl)[is_gene]
    gene_weights = relative_expression_weights(class_expression, half_saturation).reindex(gene_ensembl.to_numpy())
    weights.loc[is_gene, classes] = gene_weights.fillna(1.0).to_numpy()

    is_reaction = (nodes.node_type == "reaction").to_numpy()
    rules = nodes.gene_reaction_rule[is_reaction].fillna("").to_numpy()
    reaction_expression = pd.DataFrame({class_name: [evaluate_gene_rule(rule, class_expression[class_name].dropna().to_dict()) for rule in rules]
                                        for class_name in classes})
    weights.loc[is_reaction, classes] = relative_expression_weights(reaction_expression, half_saturation).fillna(1.0).to_numpy()
    if include_all_cells_class:
        weights[ALL_CELLS_CLASS] = 1.0
    return weights


def extracellular_pool_nodes(nodes: pd.DataFrame) -> np.ndarray:
    """Boolean per node: an extracellular metabolite (Human-GEM compartment "e")."""
    return ((nodes.node_type == "metabolite") & (nodes.compartment.fillna("") == EXTRACELLULAR_COMPARTMENT)).to_numpy()
