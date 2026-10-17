"""Drugs as graph entry nodes (docs/drug_entry_nodes.md; the user's decision of 10 October 2026: "i definitely want
drugs as nodes with an ablation on them").

What is added to a graph:

  - one node per drug perturbation, node type "drug", id DRUG:<perturbation id>;
  - one mechanism edge, relation targets, from the drug node to each node the drug would otherwise be seeded on, with
    the sign the mechanism table assigns (-1, -0.5, 0, 0.5 or 1) in the sign column and the seed's magnitude (a share
    of the drug's one mechanism target among that target's components) in the entry table. targets is the relation
    the design document reserves for this (docs/experiment_design.md section 4.1: "targets (drug to protein, with
    action type and affinity from ChEMBL)"); every graph's relation_types.json has carried it since the first build,
    with no edge, so the mechanism edges add no relation type;
  - for a drug whose plasma carrier is known, one sequestration edge, relation sequesters, sign -1, from the carrier's
    node to the drug node. The existing carriage edges use binds, whose rows all carry sign 0, which reads as moving a
    cargo without acting on it; a carriage edge that is to lower free drug needs a sign.
  - for a drug with a known clearing enzyme, one clearance edge, relation clears, sign -1, from the enzyme's node to the
    drug node (the user, 11 October 2026: "Add the enzyme labels as the carrier twin instead then"). It is the carrier's
    edge with another relation name: more enzyme, less drug. Added only when a table of clearing enzymes is passed
    (configs/drug_clearing_enzymes.csv, experiments/build_drug_clearing_enzymes.py), so a build without one is the build
    it was.

A drug that is itself a graph compound gets no node (the user: "it just has to interact with the endogenous and
exogenous drug nodes consistently"). Its perturbation stays seeded on the compound's own metabolite nodes, and a
carrier's sequestration edge goes to the compound's extracellular node. One molecule is then one node whether it is
made by the body or taken as a drug.

Two properties are kept on purpose, because the measurements of docs/drug_entry_nodes_measured.md depend on them:

  - No existing row changes. Node and edge rows are appended, the degree of an existing node is left as it was, and the
    new relation types go after the old ones. A reader that drops the appended rows has the source graph back, which is
    what makes the ablation arm (seed on the targets, as before) the same inputs as a run on the source graph.
  - A drug node holds no identity. It has a type and a degree and nothing else, so a held-out drug's node is as
    informative as a training drug's and no parameter belongs to one drug.

A neomorphic allele (decision 9.54 (d), evidence/allele_mechanisms.py) enters the same way: one node per curated
allele, node type "allele", id ALLELE:<gene>_<allele>, with one targets edge per cited mechanism of the mutant product
(configs/neomorphic_alleles.csv), to the target's protein nodes on a split graph or its gene node for a transcriptional
effect. The CAG expansion of HTT is such a node: the mutant product is not more of huntingtin's normal activity, so it
cannot be a sign on HTT's own node, and the wild-type copy keeps HTT's node as it was.

A drug node exists only in the perturbation that takes the drug. The graph table holds every drug's node, because a
table has one row per node, but the encoders treat a drug node that is not seeded as absent (state zero, nothing sent),
and they add the mechanism edges without the in-degree divisor other relations use. Otherwise a target's input would be
divided by the number of study drugs that share it, and a gene knockout's response would depend on which drugs the
study happens to contain (relational_message_passing_encoder.py and linear_response_encoder.py, entry edges).
"""
from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

DRUG_NODE_TYPE = "drug"
DRUG_NODE_PREFIX = "DRUG:"
DRUG_MECHANISM_RELATION = "targets"
SEQUESTERS_RELATION = "sequesters"
SEQUESTRATION_SIGN = -1.0
CLEARS_RELATION = "clears"
CLEARANCE_SIGN = -1.0
# relations whose edges end on a drug's node and say "more of the source, less of the drug". The message-passing encoder
# sends along them how far the source has moved from its unperturbed state, so they are silent in a perturbation that
# leaves the carrier or the enzyme at rest (experiments/run_main_model.py, change_carrying_inputs)
CHANGE_CARRYING_RELATIONS: tuple[str, ...] = (SEQUESTERS_RELATION, CLEARS_RELATION)
# relations whose edges leave an entry node: added without an in-degree divisor, and silent unless their source is seeded
ENTRY_RELATIONS: tuple[str, ...] = (DRUG_MECHANISM_RELATION,)
ALLELE_NODE_TYPE = "allele"
ENTRY_NODE_TYPES: tuple[str, ...] = (DRUG_NODE_TYPE, ALLELE_NODE_TYPE)
ALLELE_MECHANISM_EVIDENCE_PREFIX = "neomorphic allele mechanism"
DRUG_ENTRY_FILE = "drug_entry_nodes.parquet"  # one row per mechanism edge; read by experiment_data.load_experiment_data
DRUG_ENTRY_COLUMNS: tuple[str, ...] = ("perturbation_id", "drug_node_id", "target_node_id", "sign", "magnitude")
MECHANISM_EVIDENCE = "drug mechanism (ChEMBL), moved from the perturbation's seeds"
CARRIAGE_EVIDENCE_PREFIX = "plasma carrier of the drug"
CLEARANCE_EVIDENCE_PREFIX = "enzyme that clears the drug"
PLASMA_COMPARTMENT = "e"
CARRIER_GENES: dict[str, tuple[str, ...]] = {"albumin": ("ALB",), "orosomucoid": ("ORM1", "ORM2"), "both": ("ALB", "ORM1", "ORM2")}


def drug_node_id(perturbation_id: str) -> str:
    return f"{DRUG_NODE_PREFIX}{perturbation_id}"


def is_graph_compound(seed_node_ids: list[str], base_metabolite_of_node: dict[str, str], graph_compound_base_ids: list[frozenset[str]]) -> bool:
    """True when every seed is a metabolite node of one listed graph compound (configs/drugs_acting_as_graph_compounds.csv).

    A drug whose mechanism target is a metabolite (a chelator's metal, say) is not one: its seeds are the target's
    nodes, the table does not list them as a drug, and it gets a drug node like any other."""
    if not seed_node_ids:
        return False
    seed_bases = {base_metabolite_of_node.get(node_id) for node_id in seed_node_ids}
    if None in seed_bases:
        return False
    return any(seed_bases <= compound_bases for compound_bases in graph_compound_base_ids)


def add_drug_entry_nodes(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str],
                         drug_seeds: dict[str, tuple[str, list[tuple[str, float, float]]]],
                         graph_compound_base_ids: list[frozenset[str]] | None = None,
                         carriers: pd.DataFrame | None = None,
                         carrier_nodes_of_gene: dict[str, list[str]] | None = None,
                         clearing_enzymes: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame, list[str], pd.DataFrame, dict]:
    """The graph with drug entry nodes, their mechanism edges, the carriers' sequestration edges and, when
    clearing_enzymes is given, the enzymes' clearance edges appended.

    drug_seeds: perturbation id -> (drug name, [(node id, sign, magnitude)]), the seeds the loader would place for that
    drug on this graph (on a gene and protein split graph, the protein nodes of its target).
    graph_compound_base_ids: the base metabolite ids of each drug that is itself a graph compound.
    carriers: one row per drug with a known plasma carrier, columns perturbation_id and carrier ("albumin",
    "orosomucoid" or "both") and, optionally, evidence (a short source label kept in evidence_source).
    carrier_nodes_of_gene: which node carries the binding, GENE:<symbol> by default and the symbol's protein nodes on a
    split graph, as for the plasma carriage edges (plasma_binding.add_plasma_carriage).
    clearing_enzymes: one row per drug and enzyme, columns perturbation_id and enzyme_gene and, optionally, evidence.
    The enzyme's node is found as a carrier's is, through carrier_nodes_of_gene on a split graph. None adds no
    clearance edge and leaves the summary as it was.

    Returns (nodes, edges, relation types, entry table, summary). The entry table has one row per mechanism edge with
    its sign and magnitude, which is what lets the ablation arm rebuild the seeds exactly."""
    node_ids = set(nodes.node_id)
    base_metabolite_of_node = {node_id: str(base) for node_id, base in zip(nodes.node_id, nodes.base_metabolite_id) if isinstance(base, str) and base}
    compartment_of_node = dict(zip(nodes.node_id, nodes.compartment))
    graph_compound_base_ids = list(graph_compound_base_ids or [])
    carrier_of_perturbation: dict[str, tuple[str, str]] = {}
    if carriers is not None:
        evidence_column = carriers["evidence"] if "evidence" in carriers.columns else pd.Series("", index=carriers.index)
        for perturbation_id, carrier, evidence in zip(carriers.perturbation_id, carriers.carrier, evidence_column):
            if carrier not in CARRIER_GENES:
                raise ValueError(f"carrier of {perturbation_id} is {carrier!r}; expected one of {sorted(CARRIER_GENES)}")
            carrier_of_perturbation[str(perturbation_id)] = (str(carrier), str(evidence))

    enzymes_of_perturbation: dict[str, list[tuple[str, str]]] = {}
    if clearing_enzymes is not None:
        evidence_column = clearing_enzymes["evidence"] if "evidence" in clearing_enzymes.columns else pd.Series("", index=clearing_enzymes.index)
        for perturbation_id, enzyme_gene, evidence in zip(clearing_enzymes.perturbation_id, clearing_enzymes.enzyme_gene, evidence_column):
            entry = (str(enzyme_gene), str(evidence))
            if entry not in enzymes_of_perturbation.setdefault(str(perturbation_id), []):
                enzymes_of_perturbation[str(perturbation_id)].append(entry)
    clearance_edges_by_enzyme: dict[str, int] = {}
    clearance_edge_rows: list[dict] = []
    enzymes_without_a_node, drugs_with_clearance = [], []

    blank_node = {column: None for column in nodes.columns}
    node_rows, edge_rows, entry_rows = [], [], []
    attached_to_compound, without_seeds, seeds_absent_from_graph = [], [], []
    carriage_edges_by_carrier: dict[str, int] = {}
    carriers_without_a_node, drugs_with_carriage = [], []
    for perturbation_id in sorted(drug_seeds):
        drug_name, triples = drug_seeds[perturbation_id]
        present = [(node_id, float(sign), float(magnitude)) for node_id, sign, magnitude in triples if node_id in node_ids]
        if len(present) < len(triples):
            seeds_absent_from_graph.append(perturbation_id)
        if not present:
            without_seeds.append(perturbation_id)
            continue
        seed_node_ids = [node_id for node_id, _, _ in present]
        if is_graph_compound(seed_node_ids, base_metabolite_of_node, graph_compound_base_ids):
            attached_to_compound.append(perturbation_id)
            cargo_nodes = [node_id for node_id in seed_node_ids if compartment_of_node.get(node_id) == PLASMA_COMPARTMENT]
            mechanism_degree = 0
        else:
            new_node_id = drug_node_id(perturbation_id)
            if new_node_id in node_ids:
                raise ValueError(f"{new_node_id} is already a node of the graph")
            cargo_nodes = [new_node_id]
            for target_node_id, sign, magnitude in present:
                edge_rows.append({"source_id": new_node_id, "target_id": target_node_id, "relation_type": DRUG_MECHANISM_RELATION,
                                  "sign": sign, "evidence_source": MECHANISM_EVIDENCE})
                entry_rows.append({"perturbation_id": perturbation_id, "drug_node_id": new_node_id, "target_node_id": target_node_id,
                                   "sign": sign, "magnitude": magnitude})
            mechanism_degree = len(present)
        carriage_degree = 0
        if perturbation_id in carrier_of_perturbation and cargo_nodes:
            carrier, evidence = carrier_of_perturbation[perturbation_id]
            evidence_source = f"{CARRIAGE_EVIDENCE_PREFIX} ({carrier}{'; ' + evidence if evidence else ''})"
            for carrier_gene in CARRIER_GENES[carrier]:
                carrier_nodes = [node for node in (carrier_nodes_of_gene or {}).get(carrier_gene, [f"GENE:{carrier_gene}"]) if node in node_ids]
                if not carrier_nodes:
                    carriers_without_a_node.append(carrier_gene)
                for carrier_node in carrier_nodes:
                    for cargo_node in cargo_nodes:
                        edge_rows.append({"source_id": carrier_node, "target_id": cargo_node, "relation_type": SEQUESTERS_RELATION,
                                          "sign": SEQUESTRATION_SIGN, "evidence_source": evidence_source})
                        carriage_edges_by_carrier[carrier_gene] = carriage_edges_by_carrier.get(carrier_gene, 0) + 1
                        carriage_degree += 1
            if carriage_degree:
                drugs_with_carriage.append(perturbation_id)
        clearance_degree = 0
        if perturbation_id in enzymes_of_perturbation and cargo_nodes:
            for enzyme_gene, evidence in enzymes_of_perturbation[perturbation_id]:
                evidence_source = f"{CLEARANCE_EVIDENCE_PREFIX} ({enzyme_gene}{'; ' + evidence if evidence else ''})"
                enzyme_nodes = [node for node in (carrier_nodes_of_gene or {}).get(enzyme_gene, [f"GENE:{enzyme_gene}"]) if node in node_ids]
                if not enzyme_nodes:
                    enzymes_without_a_node.append(enzyme_gene)
                for enzyme_node in enzyme_nodes:
                    for cargo_node in cargo_nodes:
                        clearance_edge_rows.append({"source_id": enzyme_node, "target_id": cargo_node, "relation_type": CLEARS_RELATION,
                                                    "sign": CLEARANCE_SIGN, "evidence_source": evidence_source})
                        clearance_edges_by_enzyme[enzyme_gene] = clearance_edges_by_enzyme.get(enzyme_gene, 0) + 1
                        clearance_degree += 1
            if clearance_degree:
                drugs_with_clearance.append(perturbation_id)
        if perturbation_id not in attached_to_compound:
            # the degree does not count clearance edges, so the node table is the one a build without them writes
            node_rows.append({**blank_node, "node_id": drug_node_id(perturbation_id), "node_type": DRUG_NODE_TYPE, "display_name": drug_name,
                              "degree": mechanism_degree + carriage_degree})
    # after every other appended row, so the graph without clearance edges is this one with its last rows dropped
    edge_rows += clearance_edge_rows

    # column by column, so each keeps the type the source table gave it: concatenating frames would let the appended
    # rows, which are missing in most columns, decide the type of a column
    new_nodes = nodes if not node_rows else pd.DataFrame(
        {column: pd.Series(list(nodes[column]) + [row[column] for row in node_rows], dtype=nodes[column].dtype) for column in nodes.columns})
    new_edges = edges if not edge_rows else pd.concat([edges, pd.DataFrame(edge_rows, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict())
    new_relations = list(relation_types)
    for relation in (DRUG_MECHANISM_RELATION, *CHANGE_CARRYING_RELATIONS):
        if relation not in new_relations and any(row["relation_type"] == relation for row in edge_rows):
            new_relations.append(relation)
    entry_table = pd.DataFrame(entry_rows, columns=list(DRUG_ENTRY_COLUMNS))
    summary = {
        "drug_perturbations": len(drug_seeds),
        "drug_nodes_added": len(node_rows),
        "mechanism_edges_added": int(len(entry_table)),
        "mechanism_edges_by_sign": {str(sign): int(count) for sign, count in entry_table.sign.value_counts().sort_index().items()},
        "drugs_attached_to_a_graph_compound": sorted(attached_to_compound),
        "drugs_without_a_seed_in_this_graph": sorted(without_seeds),
        "drugs_with_a_seed_absent_from_this_graph": sorted(seeds_absent_from_graph),
        "drugs_with_a_known_carrier": int(len(carrier_of_perturbation)),
        "drugs_given_a_sequestration_edge": int(len(drugs_with_carriage)),
        "known_carriers_of_drugs_absent_from_this_graph": sorted(set(carrier_of_perturbation) - set(drug_seeds)),
        "sequestration_edges_added": int(sum(carriage_edges_by_carrier.values())),
        "sequestration_edges_by_carrier_gene": dict(sorted(carriage_edges_by_carrier.items())),
        "carrier_genes_without_a_node": sorted(set(carriers_without_a_node)),
        "relation_types_added": [relation for relation in new_relations if relation not in relation_types],
        "existing_rows_unchanged": True,
    }
    if clearing_enzymes is not None:  # recorded only when a table was passed, so a build without one has the summary it always had
        summary.update({
            "drugs_with_a_known_clearing_enzyme": len(enzymes_of_perturbation),
            "drugs_given_a_clearance_edge": len(drugs_with_clearance),
            "known_enzymes_of_drugs_absent_from_this_graph": sorted(set(enzymes_of_perturbation) - set(drug_seeds)),
            "clearance_edges_added": int(sum(clearance_edges_by_enzyme.values())),
            "clearance_edges_by_enzyme_gene": dict(sorted(clearance_edges_by_enzyme.items())),
            "enzyme_genes_without_a_node": sorted(set(enzymes_without_a_node)),
        })
    return new_nodes, new_edges, new_relations, entry_table, summary


def add_allele_entry_nodes(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str], alleles: Mapping,
                           protein_nodes_of_gene: Mapping[str, list[str]] | None = None) -> tuple[pd.DataFrame, pd.DataFrame, list[str], pd.DataFrame, dict]:
    """The graph with one entry node per curated neomorphic allele (allele_mechanisms.NeomorphicAllele, keyed by gene)
    and its mechanism edges appended. A protein-level mechanism goes to the target gene's protein nodes on a split graph
    (protein_nodes_of_gene, GENE:<symbol> -> its protein nodes) and to its gene node otherwise; a transcriptional one to
    the gene node. Returns (nodes, edges, relation types, entry table rows, summary) as add_drug_entry_nodes does; an
    allele none of whose targets is in the graph gets no node."""
    node_ids = set(nodes.node_id)
    blank_node = {column: None for column in nodes.columns}
    node_rows, edge_rows, entry_rows, missing_targets = [], [], [], []
    for gene in sorted(alleles):
        allele = alleles[gene]
        if allele.node_id in node_ids:
            raise ValueError(f"{allele.node_id} is already a node of the graph")
        degree = 0
        for edge in allele.edges:
            gene_node = f"GENE:{edge.target_symbol}"
            targets = [gene_node] if edge.target_kind == "gene" else list((protein_nodes_of_gene or {}).get(gene_node, [gene_node]))
            targets = [target for target in targets if target in node_ids]
            if not targets:
                missing_targets.append(f"{allele.allele_id} -> {edge.target_symbol}")
            for target in targets:
                edge_rows.append({"source_id": allele.node_id, "target_id": target, "relation_type": DRUG_MECHANISM_RELATION, "sign": edge.sign,
                                  "evidence_source": f"{ALLELE_MECHANISM_EVIDENCE_PREFIX} (PMID {edge.pmid})"})
                entry_rows.append({"perturbation_id": allele.allele_id, "drug_node_id": allele.node_id, "target_node_id": target, "sign": edge.sign,
                                   "magnitude": edge.magnitude})
                degree += 1
        if degree:
            node_rows.append({**blank_node, "node_id": allele.node_id, "node_type": ALLELE_NODE_TYPE, "display_name": allele.label, "degree": degree})
    new_nodes = nodes if not node_rows else pd.DataFrame(
        {column: pd.Series(list(nodes[column]) + [row[column] for row in node_rows], dtype=nodes[column].dtype) for column in nodes.columns})
    new_edges = edges if not edge_rows else pd.concat([edges, pd.DataFrame(edge_rows, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict())
    new_relations = list(relation_types) + ([DRUG_MECHANISM_RELATION] if edge_rows and DRUG_MECHANISM_RELATION not in relation_types else [])
    summary = {"allele_nodes_added": len(node_rows), "allele_mechanism_edges_added": len(edge_rows),
               "alleles": {gene: alleles[gene].allele_id for gene in sorted(alleles)}, "allele_targets_absent_from_this_graph": missing_targets}
    return new_nodes, new_edges, new_relations, pd.DataFrame(entry_rows, columns=list(DRUG_ENTRY_COLUMNS)), summary


def extend_table_for_drug_nodes(table: pd.DataFrame, drug_node_ids: list[str], fill_value: float) -> pd.DataFrame:
    """A per-node table (indexed by node id) with one row per drug node appended, every column set to fill_value.

    Node descriptors take 0: a drug node carries no descriptor, which is what keeps it free of identity. Cell-class
    weights take 1: a drug in plasma reaches every cell class, and a weight of 0 would stop the linear-response encoder
    from holding or passing the drug's response in that class."""
    already_present = [node_id for node_id in drug_node_ids if node_id in table.index]
    if already_present:
        raise ValueError(f"{len(already_present)} drug nodes already have a row, for example {already_present[0]}")
    addition = pd.DataFrame(fill_value, index=pd.Index(drug_node_ids, name=table.index.name), columns=table.columns).astype(table.dtypes.to_dict())
    return pd.concat([table, addition])
