"""Regulation of transmitter release by auto- and heteroreceptors (graph variant of design section 5.7, on top of the
neuronal layer of graph/reactome_import.py; experiments/build_neuronal_graph_variant.py).

Input: docs/curated_presynaptic_receptors.csv, one row per receptor gene, transmitter population, role (autoreceptor
or heteroreceptor) and site (presynaptic terminal or somatodendritic), with the sign of the receptor's effect on that
population's release and the PubMed record that reports it (curated from abstracts; species and mechanism recorded).

Each receptor gene becomes an activation node RCPT_<gene> ("<gene> activated by <ligand>"): the gene feeds it
(member_of, +1) and so does its endogenous ligand outside the cell (activates, +1; the linearisation of occupancy x
expression around baseline is the sum of the two). The activation node regulates the release reaction of each
population the table names (regulates_release, with the table's sign). Somatodendritic receptors act on firing and
terminal receptors on exocytosis; both end at release here, since the graph has one membrane potential for a generic
neuron, and the site is kept in the evidence. Autoreceptors close a negative loop: release -> transmitter outside ->
D2, 5-HT1A/1B, alpha2A, M2/M4, H3 or mGluR2/3 activation -> less release.

The ligand node is the extracellular copy of the ligand's Human-GEM metabolite (made non-currency, since it is a
signal), or the cytosolic copy when there is no extracellular one (2-arachidonoylglycerol crosses membranes as a
retrograde messenger).
"""
from __future__ import annotations

import pandas as pd

REGULATION_RELATION = "regulates_release"
LIGAND_BASES = {"dopamine": "MAM01736", "serotonin": "MAM02897", "noradrenaline": "MAM02617", "adrenaline": "MAM01290", "acetylcholine": "MAM01260",
                "histamine": "MAM02124", "glutamate": "MAM01974", "GABA": "MAM00970", "glycine": "MAM01986", "adenosine": "MAM01280",
                "2-arachidonoylglycerol": "MAM00635"}
POPULATION_TRANSMITTERS = {"dopaminergic": "MAM01736", "serotonergic": "MAM02897", "noradrenergic": "MAM02617", "cholinergic": "MAM01260",
                           "histaminergic": "MAM02124", "glutamatergic": "MAM01974", "GABAergic": "MAM00970", "glycinergic": "MAM01986"}


def ligand_node(base: str, node_ids: set[str]) -> str:
    for letter in ("e", "c"):
        if f"{base}{letter}" in node_ids:
            return f"{base}{letter}"
    return f"{base}e"


def add_release_regulation(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str], receptor_table: pd.DataFrame,
                           release_reactions_of_transmitter: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The graph with the receptor activation nodes and their regulation of release. receptor_table columns:
    receptor_gene, endogenous_ligand, population, role, site, effect_on_release, pmid.
    release_reactions_of_transmitter: {transmitter base id: release reaction node ids}."""
    node_ids = set(nodes.node_id)
    blank = {column: None for column in nodes.columns}
    new_rows, edge_rows, unplaced = [], [], []
    signs_by_pair: dict[tuple[str, str], set[float]] = {}
    for row in receptor_table.itertuples():
        population_base = POPULATION_TRANSMITTERS.get(row.population)
        if population_base is None or not release_reactions_of_transmitter.get(population_base) or row.endogenous_ligand not in LIGAND_BASES:
            unplaced.append(f"{row.receptor_gene} on {row.population} ({row.endogenous_ligand})")
            continue
        signs_by_pair.setdefault((row.receptor_gene, row.population), set()).add(float(row.effect_on_release))
    conflicting = sorted(f"{gene} on {population}" for (gene, population), signs in signs_by_pair.items() if len(signs) > 1)
    receptors_added = set()
    for row in receptor_table.itertuples():
        pair = (row.receptor_gene, row.population)
        if pair not in signs_by_pair or len(signs_by_pair[pair]) > 1:
            continue
        gene, activation = f"GENE:{row.receptor_gene}", f"RCPT_{row.receptor_gene}"
        ligand = ligand_node(LIGAND_BASES[row.endogenous_ligand], node_ids)
        for node_id, values in ((gene, {"node_type": "gene", "display_name": row.receptor_gene, "gene_symbol": row.receptor_gene, "in_metabolic_layer": False}),
                                (ligand, {"node_type": "metabolite", "display_name": row.endogenous_ligand, "compartment": ligand[-1],
                                          "base_metabolite_id": LIGAND_BASES[row.endogenous_ligand], "is_currency": False, "in_metabolic_layer": True}),
                                (activation, {"node_type": "protein_entity", "display_name": f"{row.receptor_gene} activated by {row.endogenous_ligand}", "compartment": "c", "is_currency": False})):
            if node_id not in node_ids:
                new_rows.append({**blank, "node_id": node_id, "degree": 0, **values})
                node_ids.add(node_id)
        if row.receptor_gene not in receptors_added:
            edge_rows.append({"source_id": gene, "target_id": activation, "relation_type": "member_of", "sign": 1.0, "evidence_source": "receptor expression"})
            edge_rows.append({"source_id": ligand, "target_id": activation, "relation_type": "activates", "sign": 1.0, "evidence_source": "endogenous ligand (curated presynaptic receptors)"})
            receptors_added.add(row.receptor_gene)
        for reaction in release_reactions_of_transmitter[POPULATION_TRANSMITTERS[row.population]]:
            edge_rows.append({"source_id": activation, "target_id": reaction, "relation_type": REGULATION_RELATION, "sign": float(row.effect_on_release),
                              "evidence_source": f"{row.role}, {row.site} (PMID {row.pmid})"})
    new_nodes = pd.concat([nodes, pd.DataFrame(new_rows, columns=nodes.columns)], ignore_index=True)
    ligand_nodes = {ligand_node(base, node_ids) for base in LIGAND_BASES.values()}
    new_nodes.loc[new_nodes.node_id.isin(ligand_nodes), "is_currency"] = False
    new_edges = pd.concat([edges, pd.DataFrame(edge_rows, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict())
    new_edges = new_edges.drop_duplicates(["source_id", "target_id", "relation_type"], ignore_index=True)
    counts = new_edges.source_id.value_counts().add(new_edges.target_id.value_counts(), fill_value=0)
    new_nodes["degree"] = counts.reindex(new_nodes.node_id).fillna(0).astype("int64").to_numpy()
    new_relations = list(relation_types) + ([REGULATION_RELATION] if REGULATION_RELATION not in relation_types else [])
    regulation = new_edges[new_edges.relation_type == REGULATION_RELATION]
    summary = {"receptors": len(receptors_added), "regulation_edges": int(len(regulation)), "inhibitory_edges": int((regulation.sign < 0).sum()),
               "facilitatory_edges": int((regulation.sign > 0).sum()), "conflicting_signs_left_out": conflicting, "rows_without_a_release_reaction": unplaced}
    return new_nodes, new_edges, new_relations, summary
