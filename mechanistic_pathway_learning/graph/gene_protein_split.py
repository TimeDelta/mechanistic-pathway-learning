"""Split each merged gene node into a gene node and a protein node (docs/gene_protein_split.md; the user's decision of
8 October 2026 that "The genes shouldn't need descriptors but the proteins should").

In the graphs built so far a gene node stands for the gene and its product. Here every gene node keeps the gene and
gets a protein node joined to it by an encodes edge (gene -> protein):
- regulates_transcription_of runs from the regulator's protein to the target's gene;
- every other relation that touches a gene node (activates, inhibits, binds, catalyzed_by, member_of, cofactor_of,
  oxidatively_modifies, produces_oxidant and the rest) touches its protein instead;
- genes whose product is one reviewed UniProt entry (UniProt lists several primary gene names, as for H3-3A and H3-3B,
  both encoding histone H3.3) share one protein node, with an encodes edge from each gene; edges that become identical
  once rerouted onto the shared node are kept once;
- a gene with several reviewed entries keeps one protein node, because the sources give the entries the same edges
  (OmniPath copies a gene's interactions onto each of its entries), and isoforms are not represented.

Degree is recomputed on the split graph (the user: "the registered degree strata SHOULD change to the new split").
The per-node expression columns of nodes.parquet move to the protein node (a shared protein takes the largest value of
its genes, as a reaction does over isozymes). Descriptor and cell-class weight tables are moved by
move_gene_blocks_to_proteins and cell_class_weights_for_split.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

PROTEIN_NODE_TYPE = "protein"
PROTEIN_NODE_PREFIX = "PROTEIN:"
ENCODES_RELATION = "encodes"
GENE_TARGET_RELATIONS = frozenset({"regulates_transcription_of"})  # the target end of these relations stays at the gene
EXPRESSION_NODE_COLUMNS = ("brain_median_tpm_max", "brain_expression_known", "brain_expressed")
MOVED_DESCRIPTOR_PREFIXES = ("protein_", "gene_brain_")


def shared_entry_gene_sets(uniprot_primary_gene_names: pd.Series) -> list[set[str]]:
    """The gene-symbol sets of the reviewed entries that list several primary gene names ("H3-3A; H3-3B")."""
    sets = []
    for field in uniprot_primary_gene_names.dropna().astype(str):
        symbols = {symbol.strip() for symbol in field.split(";") if symbol.strip()}
        if len(symbols) > 1:
            sets.append(symbols)
    return sets


def protein_node_assignment(nodes: pd.DataFrame, uniprot_primary_gene_names: pd.Series | None = None) -> pd.DataFrame:
    """One row per gene node: gene_node_id, protein_node_id, protein_display_name and shared (the protein node has
    several genes). Genes joined by a shared entry (transitively) share one protein node named after the smallest gene
    node id of the set."""
    genes = nodes[nodes.node_type == "gene"]
    node_of_symbol = dict(zip(genes.gene_symbol.astype(str), genes.node_id))
    parent = {node_id: node_id for node_id in genes.node_id}

    def find(item: str) -> str:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    for symbols in shared_entry_gene_sets(uniprot_primary_gene_names) if uniprot_primary_gene_names is not None else []:
        present = sorted(node_of_symbol[symbol] for symbol in symbols if symbol in node_of_symbol)
        for other in present[1:]:
            parent[find(other)] = find(present[0])
    members: dict[str, list[str]] = {}
    for node_id in genes.node_id:
        members.setdefault(find(node_id), []).append(node_id)
    symbol_of_node = dict(zip(genes.node_id, genes.gene_symbol.astype(str)))
    rows = []
    for group in members.values():
        group = sorted(group)
        protein_node_id = PROTEIN_NODE_PREFIX + symbol_of_node[group[0]]
        display_name = "/".join(symbol_of_node[node_id] for node_id in group)
        for node_id in group:
            rows.append({"gene_node_id": node_id, "protein_node_id": protein_node_id, "protein_display_name": display_name, "shared": len(group) > 1})
    return pd.DataFrame(rows).sort_values("gene_node_id").reset_index(drop=True)


def split_graph(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str], assignment: pd.DataFrame
                ) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The split graph's nodes, edges and relation types (encodes appended last, so earlier relation indices keep their
    meaning), and counts for the record."""
    if ENCODES_RELATION in relation_types:
        raise ValueError("the graph already has an encodes relation; it is split already")
    protein_of_gene = dict(zip(assignment.gene_node_id, assignment.protein_node_id))
    genes = nodes[nodes.node_type == "gene"].set_index("node_id")

    protein_rows = []
    for protein_node_id, group in assignment.groupby("protein_node_id", sort=True):
        members = genes.loc[group.gene_node_id]
        row = {column: np.nan for column in nodes.columns}
        row.update({"node_id": protein_node_id, "node_type": PROTEIN_NODE_TYPE, "display_name": group.protein_display_name.iloc[0],
                    "gene_symbol": members.gene_symbol.iloc[0], "is_currency": False})
        if "ensembl_gene_id" in members.columns:
            with_ensembl = members.ensembl_gene_id.dropna()
            row["ensembl_gene_id"] = with_ensembl.iloc[0] if len(with_ensembl) else None
        if "in_metabolic_layer" in members.columns:
            row["in_metabolic_layer"] = bool((members.in_metabolic_layer == True).any())  # noqa: E712 - NaN counts as False
        if "brain_median_tpm_max" in members.columns:
            values = pd.to_numeric(members.brain_median_tpm_max, errors="coerce")
            row["brain_median_tpm_max"] = values.max() if values.notna().any() else np.nan
        for flag in ("brain_expression_known", "brain_expressed"):
            if flag in members.columns:
                row[flag] = bool((members[flag] == True).any())  # noqa: E712 - NaN counts as False
        protein_rows.append(row)
    split_nodes = nodes.copy()
    is_gene = split_nodes.node_type == "gene"
    for column in EXPRESSION_NODE_COLUMNS:
        if column in split_nodes.columns:
            split_nodes.loc[is_gene, column] = np.nan if column == "brain_median_tpm_max" else False
    split_nodes = pd.concat([split_nodes, pd.DataFrame(protein_rows)], ignore_index=True)

    rerouted = edges.copy()
    keep_target_at_gene = rerouted.relation_type.isin(GENE_TARGET_RELATIONS)
    rerouted["source_id"] = rerouted.source_id.map(lambda node_id: protein_of_gene.get(node_id, node_id))
    rerouted.loc[~keep_target_at_gene, "target_id"] = rerouted.loc[~keep_target_at_gene, "target_id"].map(lambda node_id: protein_of_gene.get(node_id, node_id))
    key_columns = [column for column in ("source_id", "target_id", "relation_type", "sign") if column in rerouted.columns]
    before = len(rerouted)
    rerouted = rerouted.drop_duplicates(subset=key_columns, keep="first")
    encodes = pd.DataFrame({"source_id": assignment.gene_node_id, "target_id": assignment.protein_node_id, "relation_type": ENCODES_RELATION})
    if "sign" in edges.columns:
        encodes["sign"] = 1.0
    if "evidence_source" in edges.columns:
        encodes["evidence_source"] = "gene_protein_split"
    split_edges = pd.concat([rerouted, encodes[[column for column in edges.columns if column in encodes.columns]]], ignore_index=True)

    degree = pd.concat([split_edges.source_id, split_edges.target_id]).value_counts()
    split_nodes["degree"] = split_nodes.node_id.map(degree).fillna(0).astype(int)
    summary = {"gene_nodes": int(is_gene.sum()), "protein_nodes": len(protein_rows),
               "shared_protein_nodes": int(assignment[assignment.shared].protein_node_id.nunique()),
               "genes_on_a_shared_protein_node": int(assignment.shared.sum()),
               "encodes_edges": len(encodes), "edges_merged_onto_a_shared_protein": before - len(rerouted),
               "edges_before": len(edges), "edges_after": len(split_edges),
               "edges_per_relation": split_edges.relation_type.value_counts().to_dict()}
    return split_nodes, split_edges, list(relation_types) + [ENCODES_RELATION], summary


def move_gene_blocks_to_proteins(table: pd.DataFrame, split_nodes: pd.DataFrame, assignment: pd.DataFrame,
                                 protein_descriptors: pd.DataFrame | None = None,
                                 uniprot_primary_gene_names: pd.Series | None = None,
                                 moved_prefixes: tuple[str, ...] = MOVED_DESCRIPTOR_PREFIXES) -> pd.DataFrame:
    """A node_id-indexed descriptor table for the split graph: the moved blocks (the protein block and the gene brain
    expression block) leave the gene nodes and sit on their protein nodes; every other block keeps its rows, and the
    protein nodes get zeros in it.

    A protein node with one gene takes that gene's row. A shared protein node takes the largest value of its genes in
    each column, except that the protein block of a shared node comes from its UniProt entry's row of
    protein_descriptors (indexed by UniProt's primary gene-name field, "H3-3A; H3-3B") when given: the merged graphs
    match that table by single symbol, so these genes had no protein descriptors there (docs/gene_protein_split.md)."""
    moved = [column for column in table.columns if column.startswith(moved_prefixes)]
    protein_ids = split_nodes.node_id[split_nodes.node_type == PROTEIN_NODE_TYPE].tolist()
    result = table.reindex(pd.Index(split_nodes.node_id, name=table.index.name or "node_id")).fillna(0.0)
    gene_values = table.loc[table.index.intersection(assignment.gene_node_id), moved]
    per_protein = gene_values.groupby(assignment.set_index("gene_node_id").protein_node_id.reindex(gene_values.index).to_numpy()).max()
    result.loc[protein_ids, moved] = per_protein.reindex(protein_ids).fillna(0.0).to_numpy()
    result.loc[result.index.intersection(assignment.gene_node_id), moved] = 0.0

    protein_columns = [column for column in moved if column.startswith("protein_") and column != "protein_has_protein_descriptors"]
    if protein_descriptors is not None and uniprot_primary_gene_names is not None and protein_columns:
        entry_key_of_symbol = {}
        for field in uniprot_primary_gene_names.dropna().astype(str):
            symbols = [symbol.strip() for symbol in field.split(";") if symbol.strip()]
            if len(symbols) > 1:
                for symbol in symbols:
                    entry_key_of_symbol.setdefault(symbol, field)
        symbol_of_gene = dict(zip(split_nodes.node_id, split_nodes.gene_symbol))
        source_columns = [column[len("protein_"):] for column in protein_columns]
        for protein_node_id, group in assignment[assignment.shared].groupby("protein_node_id"):
            keys = [entry_key_of_symbol.get(str(symbol_of_gene[gene])) for gene in group.gene_node_id]
            keys = [key for key in keys if key is not None and key in protein_descriptors.index]
            if keys and set(source_columns) <= set(protein_descriptors.columns):
                result.loc[protein_node_id, protein_columns] = protein_descriptors.loc[keys[0], source_columns].to_numpy(dtype=float)
                if "protein_has_protein_descriptors" in result.columns:
                    result.loc[protein_node_id, "protein_has_protein_descriptors"] = 1.0
    return result


def cell_class_weights_for_split(weights: pd.DataFrame, split_nodes: pd.DataFrame, assignment: pd.DataFrame) -> pd.DataFrame:
    """Per-node cell-class weights for the split graph: gene nodes keep theirs (the gene is transcribed where its mRNA
    is), each protein node takes its gene's (the largest over the genes of a shared node), other nodes keep theirs."""
    gene_rows = weights.loc[weights.index.intersection(assignment.gene_node_id)]
    per_protein = gene_rows.groupby(assignment.set_index("gene_node_id").protein_node_id.reindex(gene_rows.index).to_numpy()).max()
    result = weights.reindex(pd.Index(split_nodes.node_id, name=weights.index.name or "node_id"))
    protein_ids = split_nodes.node_id[split_nodes.node_type == PROTEIN_NODE_TYPE].tolist()
    result.loc[protein_ids] = per_protein.reindex(protein_ids).to_numpy()
    return result
