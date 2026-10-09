"""The plasma substrate binders as carriers: binds edges, sign 0, from each extracellular cargo metabolite to the
binder's gene node (docs/plasma_protein_binding.md; the user's decision of 9 October 2026, "All of the plasma substrate
binders should be modeled", and "orosomucoid should also be implemented").

Why the edges exist. Albumin, orosomucoid and the other binders are already gene nodes of the full graph, but not one
of their edges reaches a small molecule they carry, so their physiological job is absent: albumin's carriage of the
non-esterified fatty acids and bilirubin, transthyretin's and thyroxine-binding globulin's carriage of the thyroid
hormones, the steroid-binding globulins' carriage of cortisol and the sex steroids, retinol-binding protein's carriage
of retinol. This module adds that carriage.

Direction and sign follow the graph's own convention for a small molecule that binds a protein: the full graph's 436
small_molecule_protein edges from OmniPath run metabolite -> gene with relation binds, and all 31,118 binds edges carry
sign 0, which is the right reading for a carrier, since it moves a substrate without acting on it. Only the
extracellular copy of a cargo is joined, because plasma is the extracellular compartment.

Two sources, kept apart in evidence_source so each edge says where it came from:
  - UniProt binding-site features (data/raw/uniprot/uniprot_plasma_binders.tsv.gz, fetched by
    experiments/fetch_plasma_binder_annotations.py). Their ligand entries carry ChEBI identifiers, which the same
    ChEBI-to-Human-GEM route the Reactome import uses turns into metabolites. These features are structurally
    evidenced and therefore sparse: they give albumin bilirubin and three metals, corticosteroid-binding globulin
    cortisol, thyroxine-binding globulin and transthyretin thyroxine, and nothing at all for orosomucoid, sex
    hormone-binding globulin, retinol-binding protein or the apolipoproteins.
  - docs/curated_plasma_carriage.csv, one row per (binder, cargo) with its source, PMID, DOI and a quoted sentence,
    for the carriage those features miss. It is deliberately restricted to cargo each binder is a principal carrier
    of: albumin's low-affinity, high-capacity buffering of every steroid is left out rather than making albumin a hub.

The gene node is filtered on the UniProt *primary* gene name, not the synonyms: FBF1, a keratin-binding protein of the
apical junction complex, lists ALB among its gene names and would otherwise be given albumin's cargo.
"""
from __future__ import annotations

import re

import pandas as pd

BINDS_RELATION = "binds"
CARRIAGE_SIGN = 0.0
PLASMA_COMPARTMENT = "e"
UNIPROT_SOURCE = "uniprot_binding_site"
UNIPROT_EVIDENCE = "UniProt binding site (plasma binder carriage)"
CURATED_EVIDENCE = "curated plasma carriage"
METABOLITE_NODE_TYPE = "metabolite"
GENE_NODE_TYPE = "gene"
PRIMARY_GENE_COLUMN = "Gene Names (primary)"
BINDING_SITE_COLUMN = "Binding site"
FEATURE_START = re.compile(r"(?=BINDING\s)")
LIGAND_NAME_PATTERN = re.compile(r'/ligand="([^"]+)"')
LIGAND_CHEBI_PATTERN = re.compile(r'/ligand_id="ChEBI:(CHEBI:\d+)"')


def ligands_of_binding_features(features: str) -> list[tuple[str, str]]:
    """(ligand name, ChEBI id) of each BINDING feature of a UniProt 'Binding site' field.

    The field concatenates features with semicolons, so each feature is taken on its own: a feature carrying a bare
    /ligand="substrate" with no /ligand_id (UniProt writes that for RBP4) must not be paired with the ChEBI identifier
    of the next feature."""
    pairs = []
    for chunk in FEATURE_START.split(features):
        name, chebi = LIGAND_NAME_PATTERN.search(chunk), LIGAND_CHEBI_PATTERN.search(chunk)
        if name and chebi:
            pairs.append((name.group(1), chebi.group(1)))
    return pairs


def uniprot_carriage_rows(uniprot_table: pd.DataFrame, binder_genes) -> pd.DataFrame:
    """(binder_gene, cargo_name, cargo_chebi, source) for every ligand of a binder's UniProt binding-site features.

    Rows are matched on the primary gene name against binder_genes, and a ligand without a ChEBI identifier (UniProt
    writes a bare /ligand="substrate" for some entries) is left out, since nothing can be mapped from it."""
    wanted = set(binder_genes)
    rows = []
    for gene, features in zip(uniprot_table[PRIMARY_GENE_COLUMN], uniprot_table[BINDING_SITE_COLUMN]):
        if gene not in wanted or not isinstance(features, str):
            continue
        for ligand_name, chebi_id in ligands_of_binding_features(features):
            rows.append({"binder_gene": gene, "cargo_name": ligand_name, "cargo_chebi": chebi_id, "source": UNIPROT_SOURCE})
    frame = pd.DataFrame(rows, columns=["binder_gene", "cargo_name", "cargo_chebi", "source"])
    return frame.drop_duplicates(["binder_gene", "cargo_chebi"]).reset_index(drop=True)


def curated_carriage_rows(curated_table: pd.DataFrame) -> pd.DataFrame:
    """The same four columns from docs/curated_plasma_carriage.csv, whose other columns hold the citation."""
    columns = ["binder_gene", "cargo_name", "cargo_chebi", "source"]
    return curated_table.loc[:, columns].drop_duplicates(["binder_gene", "cargo_chebi"]).reset_index(drop=True)


def extracellular_nodes_of_base(nodes: pd.DataFrame) -> dict[str, list[str]]:
    """Human-GEM base metabolite id -> its extracellular metabolite node ids (usually one)."""
    rows = nodes[(nodes.node_type == METABOLITE_NODE_TYPE) & (nodes.compartment == PLASMA_COMPARTMENT) & nodes.base_metabolite_id.notna()]
    mapping: dict[str, list[str]] = {}
    for base, node_id in zip(rows.base_metabolite_id, rows.node_id):
        mapping.setdefault(str(base), []).append(str(node_id))
    return {base: sorted(set(ids)) for base, ids in mapping.items()}


def add_plasma_carriage(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str], carriage_rows: pd.DataFrame,
                        bases_of_chebi, binder_nodes_of_gene: dict[str, list[str]] | None = None) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The graph with the carriage edges added, and a summary of what was and was not connected.

    bases_of_chebi(chebi_id) -> the Human-GEM base metabolite ids of that ChEBI entity (a set, possibly empty); the
    build script supplies the ChEBI route. binder_nodes_of_gene: which node carries the binding, GENE:<symbol> by
    default and the symbol's protein nodes on a gene/protein split graph, where binding is the protein's property and
    the gene node only holds expression. No node is added: every binder is already a node and every cargo that reaches
    plasma already has an extracellular copy, so a row that finds neither is reported, not invented."""
    node_ids = set(nodes.node_id)
    nodes_of_base = extracellular_nodes_of_base(nodes)
    existing = set(zip(edges.source_id, edges.target_id, edges.relation_type))
    edge_rows, connected = [], []
    binders_absent, cargo_unmapped, cargo_without_extracellular_node, already_present = [], [], [], []
    for row in carriage_rows.sort_values(["binder_gene", "cargo_name", "cargo_chebi"]).itertuples(index=False):
        binder_nodes = [node for node in (binder_nodes_of_gene or {}).get(row.binder_gene, [f"GENE:{row.binder_gene}"]) if node in node_ids]
        if not binder_nodes:
            binders_absent.append(row.binder_gene)
            continue
        bases = sorted(bases_of_chebi(row.cargo_chebi) or ())
        if not bases:
            cargo_unmapped.append(f"{row.binder_gene} {row.cargo_name} ({row.cargo_chebi})")
            continue
        cargo_nodes = [node_id for base in bases for node_id in nodes_of_base.get(base, ())]
        if not cargo_nodes:
            cargo_without_extracellular_node.append(f"{row.binder_gene} {row.cargo_name} ({row.cargo_chebi} -> {', '.join(bases)})")
            continue
        evidence = UNIPROT_EVIDENCE if row.source == UNIPROT_SOURCE else CURATED_EVIDENCE
        for binder_node in binder_nodes:
            for cargo_node in cargo_nodes:
                if (cargo_node, binder_node, BINDS_RELATION) in existing:
                    already_present.append(f"{cargo_node} -> {binder_node}")
                    continue
                existing.add((cargo_node, binder_node, BINDS_RELATION))
                edge_rows.append({"source_id": cargo_node, "target_id": binder_node, "relation_type": BINDS_RELATION,
                                  "sign": CARRIAGE_SIGN, "evidence_source": evidence})
                connected.append({"binder_gene": row.binder_gene, "cargo_name": row.cargo_name, "binder_node": binder_node,
                                  "cargo_node": cargo_node, "source": row.source})
    new_edges = (edges if not edge_rows else
                 pd.concat([edges, pd.DataFrame(edge_rows, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict()))
    new_nodes = nodes.copy()
    counts = new_edges.source_id.value_counts().add(new_edges.target_id.value_counts(), fill_value=0)
    new_nodes["degree"] = counts.reindex(new_nodes.node_id).fillna(0).astype("int64").to_numpy()
    connected_frame = pd.DataFrame(connected, columns=["binder_gene", "cargo_name", "binder_node", "cargo_node", "source"])
    summary = {
        "carriage_rows_read": int(len(carriage_rows)),
        "edges_added": len(edge_rows),
        "edges_by_source": connected_frame.source.value_counts().to_dict(),
        "edges_by_binder": connected_frame.binder_gene.value_counts().sort_index().to_dict(),
        "cargo_by_binder": {gene: sorted(set(group.cargo_name)) for gene, group in connected_frame.groupby("binder_gene")},
        "binder_nodes": sorted(set(connected_frame.binder_node)),
        "binder_genes_absent_from_graph": sorted(set(binders_absent)),
        "cargo_not_mapped_to_human_gem": sorted(set(cargo_unmapped)),
        "cargo_without_an_extracellular_node": sorted(set(cargo_without_extracellular_node)),
        "edges_already_in_the_graph": sorted(set(already_present)),
    }
    new_relations = list(relation_types) + ([BINDS_RELATION] if BINDS_RELATION not in relation_types else [])
    return new_nodes, new_edges, new_relations, summary
