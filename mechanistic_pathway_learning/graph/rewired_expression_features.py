"""Reaction expression features that follow a rewired graph (H2's rewiring control; the user's decision of 8 October 2026).

Two node inputs of the confirmatory models carry gene expression onto reactions through the Human-GEM gene rules: the
reaction_brain descriptor block (brain_expression_descriptors.brain_expression_blocks) and the reactions' cell-class
weights (brain_expression_weights.node_cell_class_weights). A model trained on a rewired graph that keeps them still
learns which genes catalyse each reaction through its features, so the real-minus-rewired difference would leave out
part of what the wiring carries. Here a reaction's rule is rewritten onto its rewired catalysts, and both inputs are
recomputed from the rewritten rules with the functions that built them, so a reaction's expression is that of the genes
the rewired graph says catalyse it.

Rewriting. A degree-preserving rewiring of catalyzed_by (gene or protein entity -> reaction) keeps every reaction's
number of catalysts. Per reaction the catalysts it keeps map to themselves and the others are paired in node order,
old with new (any pairing is one more random choice on top of the rewiring). Each gene token of the rule is replaced by
the Ensembl id of its paired new catalyst, or by a token that matches no gene when that catalyst has no Ensembl id of its
own (a Reactome protein entity, a gene named by symbol only), which evaluate_gene_rule skips as it skips any gene
without expression. The and/or structure of the rule is kept. A rule gene with no catalyzed_by edge to its reaction
(8 mentions on graph_full_neuronal) stays as it is. Intrinsic node properties (protein embeddings, metabolite and
reaction chemistry, a gene's own expression) are not wiring and stay as they are.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.brain_expression_descriptors import RULE_TOKEN, brain_expression_blocks
from mechanistic_pathway_learning.graph.brain_expression_weights import node_cell_class_weights

RULE_OPERATORS = {"(", ")", "and", "or"}
NO_EXPRESSION_TOKEN = "NO_EXPRESSION:{node}"


def substitute_rule_genes(rule: str, replacement: dict[str, str]) -> str:
    """The rule with every gene token found in replacement replaced; operators, parentheses and other genes kept."""
    return " ".join(replacement.get(token, token) if token not in RULE_OPERATORS else token for token in RULE_TOKEN.findall(rule))


def rewired_gene_reaction_rules(nodes: pd.DataFrame, edges_before: np.ndarray, edges_after: np.ndarray, catalysis: np.ndarray) -> tuple[pd.Series, dict]:
    """node_id-indexed rules of the reactions after the rewiring of catalyzed_by, and counts for the record.

    nodes: the graph's nodes table in node-index order (node_id, node_type, ensembl_gene_id, gene_reaction_rule);
    edges_before and edges_after: [2, E] node indices; catalysis: boolean per edge, the catalyzed_by edges."""
    ensembl = nodes.ensembl_gene_id.where(nodes.node_type.isin(("gene", "protein"))).to_numpy()  # protein: the catalyst on a split graph
    node_ids = nodes.node_id.to_numpy()

    def token(node: int) -> str:
        value = ensembl[node]
        return str(value) if isinstance(value, str) and value else NO_EXPRESSION_TOKEN.format(node=node_ids[node])

    catalysts_before, catalysts_after = defaultdict(set), defaultdict(set)
    for (source, target), (new_source, new_target) in zip(edges_before[:, catalysis].T.tolist(), edges_after[:, catalysis].T.tolist()):
        catalysts_before[target].add(source)
        catalysts_after[new_target].add(new_source)
    rules = nodes.gene_reaction_rule.fillna("").to_numpy()
    rewritten = pd.Series(rules, index=pd.Index(node_ids, name="node_id"), dtype=object)
    counts = {"reactions_with_a_rule": 0, "rules_changed": 0, "rule_genes_replaced": 0, "rule_genes_without_an_edge": 0,
              "rule_genes_replaced_by_a_catalyst_without_expression": 0}
    for reaction in np.flatnonzero((nodes.node_type == "reaction").to_numpy() & (rules != "")):
        counts["reactions_with_a_rule"] += 1
        before, after = catalysts_before.get(reaction, set()), catalysts_after.get(reaction, set())
        if len(before) != len(after):
            raise ValueError(f"reaction {node_ids[reaction]} has {len(before)} catalysts before the rewiring and {len(after)} after; not a degree-preserving rewiring")
        kept = before & after
        pairing = {old: old for old in kept} | dict(zip(sorted(before - kept), sorted(after - kept)))
        replacement = {token(old): token(new) for old, new in pairing.items()}
        rule_genes = {gene for gene in RULE_TOKEN.findall(rules[reaction]) if gene not in RULE_OPERATORS}
        counts["rule_genes_without_an_edge"] += len(rule_genes - set(replacement))
        changed = {gene for gene in rule_genes if replacement.get(gene, gene) != gene}
        counts["rule_genes_replaced"] += len(changed)
        counts["rule_genes_replaced_by_a_catalyst_without_expression"] += sum(replacement[gene].startswith("NO_EXPRESSION:") for gene in changed)
        if changed:
            counts["rules_changed"] += 1
            rewritten.iloc[reaction] = substitute_rule_genes(rules[reaction], replacement)
    return rewritten, counts


def rewired_reaction_brain_block(nodes: pd.DataFrame, rules: pd.Series, gene_expression: pd.DataFrame) -> pd.DataFrame:
    """The reaction_brain_* columns recomputed from rules (node_id-indexed) and the Ensembl-indexed expression table the
    descriptors were built from (gene rows are not needed: only the reaction block is returned)."""
    rewired_nodes = nodes.assign(gene_reaction_rule=rules.reindex(nodes.node_id).to_numpy())
    blocks = brain_expression_blocks(rewired_nodes, gene_expression, symbol_to_ensembl=None)
    return blocks[[column for column in blocks.columns if column.startswith("reaction_brain_")]]


def rewired_reaction_cell_class_weights(nodes: pd.DataFrame, rules: pd.Series, class_expression: pd.DataFrame, include_all_cells_class: bool = True) -> pd.DataFrame:
    """Cell-class weights of the reactions (other rows unchanged by a rewiring) recomputed from rules."""
    rewired_nodes = nodes.assign(gene_reaction_rule=rules.reindex(nodes.node_id).to_numpy())
    weights = node_cell_class_weights(rewired_nodes, class_expression, symbol_to_ensembl=None, include_all_cells_class=include_all_cells_class)
    return weights[(nodes.node_type == "reaction").to_numpy()]
