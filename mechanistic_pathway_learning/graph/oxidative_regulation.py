"""Oxidant sources and redox-sensitive protein targets (graph variant of design section 5.7, on top of the layer of
graph/reactome_import.py; experiments/build_neuronal_graph_variant.py).

Human-GEM and the imported Reactome pathways carry the chemistry of the oxidants (which enzymes make superoxide,
hydrogen peroxide, nitric oxide or peroxynitrite and which remove them) but not what an oxidant does to a protein. Two
curated tables, one PubMed record per row, supply that:

  - docs/curated_oxidant_targets.csv -> oxidatively_modifies edges, oxidant metabolite -> gene, signed by the effect of
    the modification on the protein's activity (columns species, target_gene, effect_sign, pmid).
  - docs/curated_oxidant_sources.csv -> produces_oxidant edges, gene -> oxidant metabolite, +1, for the enzymes whose
    side or main product is an oxidant where Human-GEM's reaction does not carry it (columns source_gene,
    species_produced, pmid).

The oxidant metabolite is the Human-GEM copy in the compartment the table names, else the cytosolic copy; the copies are
made non-currency, so a change in them propagates. Rows whose species has no Human-GEM metabolite, and target genes
absent from the graph, are reported rather than invented.
"""
from __future__ import annotations

import pandas as pd

TARGET_RELATION, SOURCE_RELATION = "oxidatively_modifies", "produces_oxidant"
# curated table species name -> Human-GEM base metabolite id
OXIDANT_BASES = {"hydrogen peroxide": "MAM02041", "superoxide": "MAM02631", "nitric oxide": "MAM02609", "peroxynitrite": "MAM02714",
                 "hydroxyl radical": "MAM02149", "4-hydroxynonenal": "MAM00988", "hypochlorous acid": "MAM02156", "nitrite": "MAM02588",
                 "nitrogen dioxide": "MAM02591"}
COMPARTMENT_WORDS = (("mitochondri", "m"), ("peroxisom", "x"), ("endoplasmic reticulum", "r"), ("nucle", "n"), ("lysosom", "l"),
                     ("extracellular", "e"), ("intermembrane", "i"))


def oxidant_bases(species: str) -> list[str]:
    """Human-GEM base ids for a curated species name, empty when the row names no oxidant Human-GEM has: a species
    without a metabolite (lipid hydroperoxides) or an unresolved one ('reactive oxygen species (unspecified)'). A row
    that names two ('superoxide and hydrogen peroxide') gives both."""
    name = str(species).lower()
    return [base for key, base in OXIDANT_BASES.items() if key in name]


def oxidant_node(base: str, compartment_text: str, node_ids: set[str]) -> str | None:
    """The compartment copy of the oxidant the row names, falling back to the cytosolic copy, or None when the graph
    has no copy of it (nothing in the graph makes or removes that oxidant, so there is nothing to connect)."""
    text = str(compartment_text).lower()
    for word, letter in COMPARTMENT_WORDS:
        if word in text and f"{base}{letter}" in node_ids:
            return f"{base}{letter}"
    return f"{base}c" if f"{base}c" in node_ids else None


def add_oxidative_regulation(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str], target_table: pd.DataFrame,
                             source_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The graph with the curated oxidant target and source edges."""
    node_ids = set(nodes.node_id)
    edge_rows, unresolved_species, absent_genes = [], [], []
    signs_by_pair: dict[tuple[str, str], set[float]] = {}
    for row in target_table.itertuples():
        oxidants = [oxidant_node(base, getattr(row, "compartment", ""), node_ids) for base in oxidant_bases(row.species)]
        oxidants = [oxidant for oxidant in oxidants if oxidant]
        if not oxidants:
            unresolved_species.append(f"{row.species} -> {row.target_gene}")
        elif f"GENE:{row.target_gene}" not in node_ids:
            absent_genes.append(row.target_gene)
        else:
            for oxidant in oxidants:
                signs_by_pair.setdefault((oxidant, row.target_gene), set()).add(float(row.effect_sign))
    conflicting = sorted(f"{oxidant} -> {gene}" for (oxidant, gene), signs in signs_by_pair.items() if len(signs) > 1)
    for (oxidant, gene), signs in signs_by_pair.items():
        if len(signs) == 1:
            pmids = sorted({str(row.pmid) for row in target_table.itertuples() if row.target_gene == gene and oxidant in
                            [node for node in (oxidant_node(base, getattr(row, "compartment", ""), node_ids) for base in oxidant_bases(row.species)) if node]})
            edge_rows.append({"source_id": oxidant, "target_id": f"GENE:{gene}", "relation_type": TARGET_RELATION, "sign": next(iter(signs)),
                              "evidence_source": f"curated oxidative modification (PMID {', '.join(pmids)})"})
    for row in source_table.itertuples():
        oxidants = [oxidant_node(base, row.compartment, node_ids) for base in oxidant_bases(row.species_produced)]
        oxidants = [oxidant for oxidant in oxidants if oxidant]
        if not oxidants:
            unresolved_species.append(f"{row.source_gene} -> {row.species_produced}")
        elif f"GENE:{row.source_gene}" not in node_ids:
            absent_genes.append(row.source_gene)
        else:
            edge_rows += [{"source_id": f"GENE:{row.source_gene}", "target_id": oxidant, "relation_type": SOURCE_RELATION, "sign": 1.0,
                           "evidence_source": f"curated oxidant source (PMID {row.pmid})"} for oxidant in oxidants]
    new_edges = pd.concat([edges, pd.DataFrame(edge_rows, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict())
    new_edges = new_edges.drop_duplicates(["source_id", "target_id", "relation_type"], ignore_index=True)
    new_nodes = nodes.copy()
    new_nodes.loc[new_nodes.base_metabolite_id.isin(set(OXIDANT_BASES.values())), "is_currency"] = False
    counts = new_edges.source_id.value_counts().add(new_edges.target_id.value_counts(), fill_value=0)
    new_nodes["degree"] = counts.reindex(new_nodes.node_id).fillna(0).astype("int64").to_numpy()
    new_relations = list(relation_types) + [relation for relation in (TARGET_RELATION, SOURCE_RELATION) if relation not in relation_types]
    target_edges = new_edges[new_edges.relation_type == TARGET_RELATION]
    summary = {"target_edges": int(len(target_edges)), "inactivating_edges": int((target_edges.sign < 0).sum()), "activating_edges": int((target_edges.sign > 0).sum()),
               "source_edges": int((new_edges.relation_type == SOURCE_RELATION).sum()), "conflicting_signs_left_out": conflicting,
               "species_without_a_human_gem_metabolite": sorted(set(unresolved_species)), "genes_absent_from_the_graph": sorted(set(absent_genes))}
    return new_nodes, new_edges, new_relations, summary
