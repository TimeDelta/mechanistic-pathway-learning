"""BRENDA's human enzyme-regulation records as metabolite-to-enzyme edges (decision 9.5, settled 15 October 2026).

Why. The metabolic layer could receive but almost could not act: of 266,966 edges, 591 ran metabolite to protein and
only **four** of those were `inhibits`, so no metabolite could regulate an enzyme in the graph. BRENDA supplies 4,143
pairs that land on nodes the graph already holds - 1,552 `inhibits`, 482 `activates` and 2,109 `cofactor_of` - over 501
metabolites and 1,364 enzymes (`experiments/map_brenda_regulation_to_graph.py`, `docs/brenda_regulation_mapping.json`).

Licence. "The usage of the BRENDA data is licensed under Creative Commons Attribution License CC BY 4.0", read on the
provider's download page and again at `copy.php`, which the shipped README points to, so the derived table ships here.

Unweighted, and why. The owner asked whether the sources could instead be weighted by a Dawid-Skene fit, estimating a
source's reliability from its overlapping rows and carrying it to the rest. The transfer logic is right - reliability
is a per-source parameter - but the data cannot anchor it. BRENDA meets the existing edge sources on **53 of 4,143
pairs**, 40 of them `cofactor_of` against one source (`UniProt cofactor Mn(2+)`) on one metabolite; and of the other 13,
the existing edge is usually OmniPath `binds`, sign 0, meaning unknown, so BRENDA naming a sign is a specificity gain
that a fit would score as disagreement. Above all, **no edge source ever asserts a negative**: no row says a metabolite
does not inhibit an enzyme, so specificity is unestimable for every source and sensitivity alone rewards volume, which
`evidence_reliability_model` already states for the label sources ("the estimate becomes a decreasing function of how
much the source reports rather than a measurement of anything"). So every pair enters at sign alone and the evidence
quantity travels as recorded provenance instead of as a weight.

One edge per (metabolite, enzyme) pair, from one representative compartment copy. Protein nodes carry no compartment
in this graph while metabolites do, so emitting from every compartment copy of a base would multiply an
uncompartmentalised enzyme's in-degree by the number of compartments and distort item 28's `total_in_degree` divisor.
The representative is the cytosolic copy where the base has one, since BRENDA's assays are cytosolic preparations for
the most part, and otherwise the first compartment present in sorted order. The choice is recorded per pair.
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

DEFAULT_BRENDA_REGULATION_PATH = Path("configs/brenda_metabolite_enzyme_regulation.csv")

# BRENDA's relation names and the sign each carries on the destination's activity. cofactor_of is positive because an
# enzyme without its cofactor does not turn over; it is already a relation of this graph (the manganese extension emits
# it), so nothing new enters the vocabulary.
SIGN_OF_RELATION = {"inhibits": -1.0, "activates": 1.0, "cofactor_of": 1.0}
PREFERRED_COMPARTMENT = "c"


def load_brenda_regulation(regulation_path: Path | None) -> list[dict]:
    """The derived table: base_metabolite_id, uniprot_accession, relation_type, sign, brenda_record_count, ec_numbers,
    brenda_reference_count, enzyme_routes, example_ligand_name."""
    if regulation_path is None or not Path(regulation_path).exists():
        return []
    with open(regulation_path, encoding="utf-8") as regulation_file:
        return list(csv.DictReader(regulation_file))


def representative_metabolite_node(nodes: pd.DataFrame) -> dict[str, str]:
    """base metabolite id -> the one node id its edges are emitted from (cytosol where present, else first in order)."""
    metabolites = nodes[(nodes.node_type == "metabolite") & nodes.base_metabolite_id.notna()]
    representative: dict[str, str] = {}
    for base, group in metabolites.groupby("base_metabolite_id"):
        cytosolic = group[group.compartment == PREFERRED_COMPARTMENT]
        chosen = cytosolic if len(cytosolic) else group.sort_values("node_id")
        representative[str(base)] = str(chosen.node_id.iloc[0])
    return representative


def enzyme_node_of_accession(nodes: pd.DataFrame, gene_of_accession: dict[str, str] | None = None) -> dict[str, str]:
    """UniProt accession -> the node its regulation arrives at.

    Protein nodes where the gene/protein split has already run, gene nodes where it has not, so this module can sit on
    either side of it: the split moves every relation that touches a gene node onto its protein. Protein nodes carry
    their accession in `uniprot_entry`, a column the split writes; before the split the graph has gene nodes only, and
    those carry a gene symbol and no accession (`build_physiology_graph.py` writes no `uniprot_entry`), so an accession
    reaches its gene node through `gene_of_accession`, UniProt's accession -> primary gene name, the map the neuronal
    build already reads for its Reactome layer. A gene node that does carry `uniprot_entry` resolves by it as well.
    """
    entries = nodes["uniprot_entry"] if "uniprot_entry" in nodes.columns else pd.Series([None] * len(nodes), index=nodes.index, dtype=object)
    proteins = nodes[(nodes.node_type == "protein") & entries.notna()]
    if len(proteins):
        return {str(accession): str(node_id) for node_id, accession in zip(proteins.node_id, proteins.uniprot_entry)}
    genes = nodes[nodes.node_type == "gene"]
    node_of_symbol = {str(symbol): str(node_id) for node_id, symbol in zip(genes.node_id, genes.gene_symbol)
                      if isinstance(symbol, str)} if "gene_symbol" in genes.columns else {}
    resolved = {str(accession): node_of_symbol[str(symbol)] for accession, symbol in (gene_of_accession or {}).items()
                if str(symbol) in node_of_symbol}
    gene_entries = entries.reindex(genes.index)
    resolved.update({str(accession): str(node_id) for node_id, accession in zip(genes.node_id, gene_entries) if accession is not None and pd.notna(accession)})
    return resolved


def add_brenda_enzyme_regulation(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str],
                                 regulation_rows: list[dict],
                                 gene_of_accession: dict[str, str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The graph with BRENDA's regulation edges added, and a summary that counts what landed and what did not.

    Nothing is invented: a row whose metabolite or whose accession the graph does not hold is skipped and counted, and
    a pair the graph already carries under the same relation is skipped as a duplicate rather than doubled.
    `gene_of_accession` (UniProt accession -> primary gene name) is what resolves an accession before the gene/protein
    split, when the graph holds gene nodes and no `uniprot_entry` column (enzyme_node_of_accession).
    """
    representative = representative_metabolite_node(nodes)
    enzyme_of_accession = enzyme_node_of_accession(nodes, gene_of_accession)
    existing_pairs = {(str(source), str(target), str(relation))
                      for source, target, relation in zip(edges.source_id, edges.target_id, edges.relation_type)}

    edge_rows: list[dict] = []
    emitted: dict[str, int] = {relation: 0 for relation in SIGN_OF_RELATION}
    skipped = {"metabolite_absent": 0, "enzyme_absent": 0, "relation_unknown": 0, "already_present": 0,
               "duplicate_within_brenda": 0}
    seen_within_brenda: set[tuple[str, str, str]] = set()
    for row in regulation_rows:
        relation = row["relation_type"]
        if relation not in SIGN_OF_RELATION:
            skipped["relation_unknown"] += 1
            continue
        source_node = representative.get(row["base_metabolite_id"])
        target_node = enzyme_of_accession.get(row["uniprot_accession"])
        if source_node is None:
            skipped["metabolite_absent"] += 1
            continue
        if target_node is None:
            skipped["enzyme_absent"] += 1
            continue
        key = (source_node, target_node, relation)
        if key in seen_within_brenda:
            skipped["duplicate_within_brenda"] += 1
            continue
        seen_within_brenda.add(key)
        if key in existing_pairs:
            skipped["already_present"] += 1
            continue
        references = row.get("brenda_reference_count") or "0"
        edge_rows.append({"source_id": source_node, "target_id": target_node, "relation_type": relation,
                          "sign": SIGN_OF_RELATION[relation],
                          "evidence_source": f"BRENDA CC BY 4.0 ({references} references, {row.get('enzyme_routes', 'route unrecorded')})"})
        emitted[relation] += 1

    new_edges = pd.concat([edges, pd.DataFrame(edge_rows, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict()) if edge_rows else edges.copy()
    new_nodes = nodes.copy()
    counts = new_edges.source_id.value_counts().add(new_edges.target_id.value_counts(), fill_value=0)
    new_nodes["degree"] = counts.reindex(new_nodes.node_id).fillna(0).astype("int64").to_numpy()
    new_relations = list(relation_types) + [relation for relation in SIGN_OF_RELATION
                                            if relation not in relation_types and emitted[relation]]
    landed_on = {node_type for node_type in nodes.set_index("node_id").node_type.reindex({row["target_id"] for row in edge_rows}).dropna()}
    summary = {"rows_read": len(regulation_rows), "edges_emitted_by_relation": emitted,
               "edges_emitted": sum(emitted.values()), "skipped": skipped,
               "distinct_metabolites": len({row["source_id"] for row in edge_rows}),
               "distinct_enzymes": len({row["target_id"] for row in edge_rows}),
               "enzyme_nodes": sorted(landed_on) or ["none"],
               "accession_resolution": ("uniprot_entry of the protein nodes (after the gene/protein split)" if "protein" in landed_on
                                        else "UniProt accession -> primary gene name -> gene node (before the split)"),
               "accessions_resolvable": len(enzyme_of_accession),
               "compartment_policy": f"one edge per pair from the {PREFERRED_COMPARTMENT} copy where present, else the first in sorted order",
               "weighting": "none: every pair enters at sign alone, with the reference count and resolution route recorded in evidence_source",
               "licence": "Creative Commons Attribution License CC BY 4.0, stated on the BRENDA download page and at copy.php"}
    return new_nodes, new_edges, new_relations, summary
