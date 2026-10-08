"""Split each merged gene node into a gene node and its protein nodes (docs/gene_protein_split.md; the user's decisions
of 8 October 2026: "The genes shouldn't need descriptors but the proteins should", then that the gene's descriptors
"would be limited to only the brain expression columns", and that several genes encoding one protein needs "an
equivalent several proteins from one gene").

In the graphs built so far a gene node stands for the gene and its product. Here every gene node keeps the gene and
gets one protein node per reviewed UniProt entry that lists the gene among its primary gene names, joined to it by an
encodes edge (gene -> protein):
- regulates_transcription_of runs from the regulator's protein to the target's gene;
- every other relation that touches a gene node (activates, inhibits, binds, catalyzed_by, member_of, cofactor_of,
  oxidatively_modifies, produces_oxidant and the rest) touches its protein instead;
- genes listed by one entry (H3-3A and H3-3B, both encoding histone H3.3) share its protein node, with an encodes edge
  from each gene; edges that become identical once rerouted onto the shared node are kept once;
- a gene listed by several entries (CDKN2A: p16INK4a and p14ARF) gets one protein node per entry. An edge end at such a
  gene goes to the entries its source names (OmniPath and CollecTRI give the accession of each row, a Reactome entity
  its member accessions); where the source names none of the gene's entries, or does not resolve entries at all
  (Human-GEM gene rules, curated tables), the end goes to every protein node of the gene, since which product is meant
  is not known. Isoforms inside one entry are not represented: no source resolves them;
- a gene with no reviewed entry keeps one protein node, named after the gene, without protein descriptors.

Degree is recomputed on the split graph (the user: "the registered degree strata SHOULD change to the new split"). The
per-node expression columns of nodes.parquet and the gene brain expression descriptor block stay on the gene nodes; the
protein descriptor block sits on the protein nodes, one row per entry (place_descriptors_on_split). Cell-class weights:
cell_class_weights_for_split.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd

PROTEIN_NODE_TYPE = "protein"
PROTEIN_NODE_PREFIX = "PROTEIN:"
ENCODES_RELATION = "encodes"
GENE_TARGET_RELATIONS = frozenset({"regulates_transcription_of"})  # the target end of these relations stays at the gene
EXPRESSION_NODE_COLUMNS = ("brain_median_tpm_max", "brain_expression_known", "brain_expressed")  # stay on the gene nodes
PROTEIN_DESCRIPTOR_PREFIX = "protein_"
PROTEIN_DESCRIPTOR_FLAG = "protein_has_protein_descriptors"
UNIPROT_ENTRY_COLUMN = "Entry"
UNIPROT_GENE_NAMES_COLUMN = "Gene Names (primary)"


def entries_of_gene_symbol(uniprot_entries: pd.DataFrame) -> dict[str, list[str]]:
    """gene symbol -> the reviewed accessions listing it among their primary gene names ("H3-3A; H3-3B"), in table order.
    The whole field of an entry with several names maps to the entry too: the Reactome import named a few gene nodes
    after it ("GENE:HBA1; HBA2"), and those nodes stand for that entry."""
    entries: dict[str, list[str]] = defaultdict(list)
    for entry, field in zip(uniprot_entries[UNIPROT_ENTRY_COLUMN], uniprot_entries[UNIPROT_GENE_NAMES_COLUMN]):
        if isinstance(field, str):
            symbols = [symbol.strip() for symbol in field.split(";") if symbol.strip()]
            for symbol in symbols + ([field.strip()] if len(symbols) > 1 else []):
                if entry not in entries[symbol]:
                    entries[symbol].append(entry)
    return dict(entries)


def interaction_members(symbol_field: str, accession_field: str) -> list[tuple[str, str | None]]:
    """(symbol, accession) of each member of an OmniPath row end; complexes are A_B with accessions COMPLEX:P1_P2. When the
    two lists do not line up the accessions are not used."""
    symbols = [member for member in symbol_field.split("_") if member]
    accessions = [member for member in accession_field.replace("COMPLEX:", "").split("_") if member]
    return list(zip(symbols, accessions)) if len(symbols) == len(accessions) else [(symbol, None) for symbol in symbols]


def entries_named_by_interaction_rows(rows, evidence_source: str, entries_of_symbol: dict[str, list[str]], relation_of_row,
                                      gene_node_prefix: str = "GENE:") -> dict[tuple[str, str, str, str, str], set[str]]:
    """For OmniPath-format rows (source, target, source_genesymbol, target_genesymbol): the accessions each graph edge's
    rows name at an end whose gene has several entries, keyed as split_graph's named_entries. relation_of_row gives the
    relation build_physiology_graph makes of a row; binds edges exist in both directions, as it builds them."""
    named: dict[tuple[str, str, str, str, str], set[str]] = defaultdict(set)
    for row in rows:
        relation = relation_of_row(row)
        for source, source_accession in interaction_members(row["source_genesymbol"], row["source"]):
            for target, target_accession in interaction_members(row["target_genesymbol"], row["target"]):
                if source == target:
                    continue
                directions = [(source, source_accession, target, target_accession)]
                if relation == "binds":
                    directions.append((target, target_accession, source, source_accession))
                for first, first_accession, second, second_accession in directions:
                    key = (gene_node_prefix + first, gene_node_prefix + second, relation, evidence_source)
                    if first_accession and len(entries_of_symbol.get(first, [])) > 1:
                        named[(*key, "source")].add(first_accession)
                    if second_accession and len(entries_of_symbol.get(second, [])) > 1:
                        named[(*key, "target")].add(second_accession)
    return dict(named)


def protein_node_assignment(nodes: pd.DataFrame, uniprot_entries: pd.DataFrame | None = None,
                            drug_target_entries: dict[str, set[str]] | None = None) -> pd.DataFrame:
    """One row per (gene node, protein node) pair: gene_node_id, protein_node_id, protein_display_name, uniprot_entry
    (None for a gene without a reviewed entry), shared (the protein node has several genes), several_proteins (the gene
    has several protein nodes) and drug_target (a drug acting on the gene seeds this protein node).

    A protein node is named after the smallest gene node id among its genes, with ":<accession>" added when one of its
    genes has several entries. drug_target_entries maps a gene symbol to the accessions drug targets name for it
    (ChEMBL target components); a gene with several entries has its drugs seed the named ones, or all of them when no
    drug target names one."""
    genes = nodes[nodes.node_type == "gene"]
    symbol_of_node = dict(zip(genes.node_id, genes.gene_symbol.astype(str)))
    entries = entries_of_gene_symbol(uniprot_entries) if uniprot_entries is not None else {}
    entries_of_gene = {node_id: entries.get(symbol, []) for node_id, symbol in symbol_of_node.items()}
    genes_of_entry: dict[str, list[str]] = defaultdict(list)
    for node_id in sorted(entries_of_gene):
        for entry in entries_of_gene[node_id]:
            genes_of_entry[entry].append(node_id)

    rows = []
    for node_id in sorted(entries_of_gene):
        if not entries_of_gene[node_id]:
            rows.append({"gene_node_id": node_id, "protein_node_id": PROTEIN_NODE_PREFIX + symbol_of_node[node_id],
                         "protein_display_name": symbol_of_node[node_id], "uniprot_entry": None})
    for entry, gene_ids in genes_of_entry.items():
        with_suffix = any(len(entries_of_gene[gene]) > 1 for gene in gene_ids)
        protein_node_id = PROTEIN_NODE_PREFIX + symbol_of_node[gene_ids[0]] + (f":{entry}" if with_suffix else "")
        display_name = "/".join(symbol_of_node[gene] for gene in gene_ids) + (f" ({entry})" if with_suffix else "")
        for gene in gene_ids:
            rows.append({"gene_node_id": gene, "protein_node_id": protein_node_id, "protein_display_name": display_name, "uniprot_entry": entry})
    assignment = pd.DataFrame(rows, columns=["gene_node_id", "protein_node_id", "protein_display_name", "uniprot_entry"])
    if assignment.groupby("protein_node_id").uniprot_entry.nunique(dropna=False).max() > 1:
        raise ValueError("two entries were given one protein node id")
    assignment["shared"] = assignment.groupby("protein_node_id").gene_node_id.transform("size") > 1
    assignment["several_proteins"] = assignment.groupby("gene_node_id").protein_node_id.transform("size") > 1
    drug_target_entries = drug_target_entries or {}
    named_for_gene = [set(drug_target_entries.get(symbol_of_node[gene], set())) & set(entries_of_gene[gene]) for gene in assignment.gene_node_id]
    assignment["drug_target"] = [not several or not named or entry in named
                                 for several, named, entry in zip(assignment.several_proteins, named_for_gene, assignment.uniprot_entry)]
    return assignment.sort_values(["gene_node_id", "protein_node_id"]).reset_index(drop=True)


def split_graph(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str], assignment: pd.DataFrame,
                named_entries: dict[tuple[str, str, str, str, str], set[str]] | None = None
                ) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The split graph's nodes, edges and relation types (encodes appended last, so earlier relation indices keep their
    meaning), and counts for the record.

    named_entries maps an edge end at a gene with several entries, keyed (source_id, target_id, relation_type,
    evidence_source, end) with end "source" or "target", to the accessions its source names (experiments/
    build_gene_protein_split.py reads them from the raw OmniPath, CollecTRI and Reactome files)."""
    if ENCODES_RELATION in relation_types:
        raise ValueError("the graph already has an encodes relation; it is split already")
    named_entries = named_entries or {}
    genes = nodes[nodes.node_type == "gene"].set_index("node_id")
    proteins_of_gene = assignment.groupby("gene_node_id").protein_node_id.agg(list).to_dict()
    entry_of_protein = dict(zip(assignment.protein_node_id, assignment.uniprot_entry))

    protein_rows = []
    for protein_node_id, group in assignment.groupby("protein_node_id", sort=True):
        members = genes.loc[group.gene_node_id]
        row = {column: np.nan for column in nodes.columns}
        row.update({"node_id": protein_node_id, "node_type": PROTEIN_NODE_TYPE, "display_name": group.protein_display_name.iloc[0],
                    "gene_symbol": members.gene_symbol.iloc[0], "is_currency": False, "uniprot_entry": group.uniprot_entry.iloc[0]})
        if "ensembl_gene_id" in members.columns:  # identifies the gene whose expression the rewired descriptors read, not a descriptor
            with_ensembl = members.ensembl_gene_id.dropna()
            row["ensembl_gene_id"] = with_ensembl.iloc[0] if len(with_ensembl) else None
        if "in_metabolic_layer" in members.columns:
            row["in_metabolic_layer"] = bool((members.in_metabolic_layer == True).any())  # noqa: E712 - NaN counts as False
        for flag in ("brain_expression_known", "brain_expressed"):
            if flag in nodes.columns:
                row[flag] = False
        protein_rows.append(row)
    split_nodes = pd.concat([nodes.copy(), pd.DataFrame(protein_rows)], ignore_index=True)

    end_counts = {"single": 0, "named_subset": 0, "named_all": 0, "unresolved_all": 0}

    def proteins_at(node_id: str, key: tuple[str, str, str, str], end: str) -> list[str]:
        proteins = proteins_of_gene.get(node_id)
        if proteins is None:
            return [node_id]
        if len(proteins) == 1:
            end_counts["single"] += 1
            return proteins
        named = named_entries.get((*key, end), set())
        chosen = [protein for protein in proteins if entry_of_protein[protein] in named]
        if not chosen:
            end_counts["unresolved_all"] += 1
            return proteins
        end_counts["named_all" if len(chosen) == len(proteins) else "named_subset"] += 1
        return chosen

    evidence = edges.evidence_source if "evidence_source" in edges.columns else pd.Series("", index=edges.index)
    sources, targets = [], []
    for source_id, target_id, relation, evidence_source in zip(edges.source_id, edges.target_id, edges.relation_type, evidence.astype(str)):
        key = (source_id, target_id, relation, evidence_source)
        sources.append(proteins_at(source_id, key, "source"))
        targets.append([target_id] if relation in GENE_TARGET_RELATIONS else proteins_at(target_id, key, "target"))
    rerouted = edges.copy()
    rerouted["source_id"], rerouted["target_id"] = sources, targets
    rerouted = rerouted.explode("source_id").explode("target_id").reset_index(drop=True)
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
    several = assignment[assignment.several_proteins]
    summary = {"gene_nodes": int((nodes.node_type == "gene").sum()), "protein_nodes": len(protein_rows),
               "protein_nodes_without_a_reviewed_entry": int(assignment[assignment.uniprot_entry.isna()].protein_node_id.nunique()),
               "shared_protein_nodes": int(assignment[assignment.shared].protein_node_id.nunique()),
               "genes_on_a_shared_protein_node": int(assignment[assignment.shared].gene_node_id.nunique()),
               "genes_with_several_protein_nodes": int(several.gene_node_id.nunique()),
               "protein_nodes_of_genes_with_several": int(several.protein_node_id.nunique()),
               "edge_ends_at_genes_with_several_proteins": {key: value for key, value in end_counts.items() if key != "single"},
               "encodes_edges": len(encodes), "edges_before": len(edges), "edges_after_rerouting_before_merging": before,
               "edges_merged_after_rerouting": before - len(rerouted), "edges_after": len(split_edges),
               "edges_per_relation": split_edges.relation_type.value_counts().to_dict()}
    return split_nodes, split_edges, list(relation_types) + [ENCODES_RELATION], summary


def place_descriptors_on_split(table: pd.DataFrame, split_nodes: pd.DataFrame, assignment: pd.DataFrame,
                               protein_descriptors_per_entry: pd.DataFrame | None = None) -> pd.DataFrame:
    """A node_id-indexed descriptor table for the split graph. The protein block (columns protein_*) leaves the gene
    nodes and sits on the protein nodes: each protein node takes its entry's row of protein_descriptors_per_entry
    (indexed by accession, columns without the protein_ prefix), and a protein node without an entry gets zeros and the
    flag at 0. Without the per-entry table a protein node takes the largest value over its genes' rows. Every other
    block keeps its rows (the gene brain expression block stays on the genes, the user's decision), and the protein
    nodes get zeros in it."""
    protein_columns = [column for column in table.columns if column.startswith(PROTEIN_DESCRIPTOR_PREFIX)]
    value_columns = [column for column in protein_columns if column != PROTEIN_DESCRIPTOR_FLAG]
    protein_ids = split_nodes.node_id[split_nodes.node_type == PROTEIN_NODE_TYPE].tolist()
    result = table.reindex(pd.Index(split_nodes.node_id, name=table.index.name or "node_id")).fillna(0.0)
    if not protein_columns:
        return result
    gene_ids = result.index.intersection(assignment.gene_node_id)
    if protein_descriptors_per_entry is None:
        pairs = assignment[assignment.gene_node_id.isin(table.index)]
        per_protein = table.loc[pairs.gene_node_id, protein_columns].set_axis(pairs.protein_node_id.to_numpy()).groupby(level=0).max()
        result.loc[protein_ids, protein_columns] = per_protein.reindex(protein_ids).fillna(0.0).to_numpy()
    else:
        source_columns = [column[len(PROTEIN_DESCRIPTOR_PREFIX):] for column in value_columns]
        entry_of_protein = assignment.drop_duplicates("protein_node_id").set_index("protein_node_id").uniprot_entry.reindex(protein_ids)
        known = entry_of_protein.isin(protein_descriptors_per_entry.index).to_numpy()
        values = np.zeros((len(protein_ids), len(value_columns)))
        values[known] = protein_descriptors_per_entry.loc[entry_of_protein[known], source_columns].to_numpy(dtype=float)
        result.loc[protein_ids, value_columns] = values
        if PROTEIN_DESCRIPTOR_FLAG in protein_columns:
            result.loc[protein_ids, PROTEIN_DESCRIPTOR_FLAG] = known.astype(float)
    result.loc[gene_ids, protein_columns] = 0.0
    return result


def cell_class_weights_for_split(weights: pd.DataFrame, split_nodes: pd.DataFrame, assignment: pd.DataFrame) -> pd.DataFrame:
    """Per-node cell-class weights for the split graph: gene nodes keep theirs (the gene is transcribed where its mRNA
    is), each protein node takes its gene's (the largest over the genes of a shared node), other nodes keep theirs."""
    pairs = assignment[assignment.gene_node_id.isin(weights.index)]
    per_protein = weights.loc[pairs.gene_node_id].set_axis(pairs.protein_node_id.to_numpy()).groupby(level=0).max()
    result = weights.reindex(pd.Index(split_nodes.node_id, name=weights.index.name or "node_id"))
    protein_ids = split_nodes.node_id[split_nodes.node_type == PROTEIN_NODE_TYPE].tolist()
    result.loc[protein_ids] = per_protein.reindex(protein_ids).to_numpy()
    return result
