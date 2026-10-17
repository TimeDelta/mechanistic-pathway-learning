"""Audit: is every transporter connected all the way down, gene -> transport reaction -> the metabolite in both
compartments, in the slice and the full graph?

Transporters are the reviewed human proteins annotated (any evidence, NOT-qualified rows excluded, propagated through
is_a and part_of) with transmembrane transporter activity, GO:0022857. For each, the audit reports whether the gene is
a node of each graph, how many transport reactions it catalyses (catalyzed_by edges to reactions flagged is_transport)
and which metabolites those reactions move; a transporter that is a node but catalyses no transport reaction is
connected only through signalling or binding edges (full graph) or not at all. A list of transporters that matter for
this project (neurotransmitter, metal, glucose, lactate, amino acid and creatine transporters) is shown in full, with
the metabolites they move, so their substrates can be checked against what the protein carries. Transport reactions
without any catalysing gene are counted too: they move metabolites, but no gene perturbation reaches them.

Usage (the default graphs are the merged full graph and the v5 graph the registered models read):
  python experiments/audit_transporter_connectivity.py --markdown-output docs/transporter_connectivity_audit.md
"""
from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.documentation import provenance_stamp
from mechanistic_pathway_learning.graph.protein_descriptors import read_go_ancestors

TRANSMEMBRANE_TRANSPORTER_ACTIVITY = "GO:0022857"
KEY_TRANSPORTERS = {
    "dopamine, noradrenaline, serotonin": ["SLC6A3", "SLC6A2", "SLC6A4", "SLC18A1", "SLC18A2", "SLC22A3", "SLC29A4"],
    "GABA, glycine, glutamate": ["SLC6A1", "SLC6A11", "SLC6A13", "SLC32A1", "SLC6A5", "SLC6A9", "SLC1A1", "SLC1A2", "SLC1A3", "SLC17A6", "SLC17A7", "SLC17A8"],
    "choline, acetylcholine, histamine": ["SLC5A7", "SLC18A3", "SLC44A1", "SLC22A2"],
    "amino acids at the blood-brain barrier": ["SLC7A5", "SLC3A2", "SLC7A11", "SLC38A1", "SLC38A2", "SLC1A5", "SLC6A19"],
    "glucose, lactate, ketone bodies, creatine": ["SLC2A1", "SLC2A3", "SLC16A1", "SLC16A3", "SLC16A7", "SLC6A8"],
    "manganese, iron, zinc, copper": ["SLC39A14", "SLC39A8", "SLC30A10", "SLC11A2", "SLC40A1", "TFRC", "SLC39A4", "SLC30A8", "ATP7A", "ATP7B", "SLC31A1", "ATP2C1"],
    "potassium, sodium, chloride": ["ATP1A1", "ATP1A2", "ATP1A3", "SLC12A1", "SLC12A2", "SLC12A3", "SLC12A5", "KCNJ1", "CLCNKB"],
    "folate, vitamins, cholesterol": ["SLC46A1", "FOLR1", "SLC19A1", "SLC19A3", "SLC52A2", "SLC5A6", "ABCA1", "ABCG1", "NPC1", "NPC2"],
}


def transporter_genes(uniprot_table: pd.DataFrame, gaf_path: Path, go_ancestors: dict[str, set[str]]) -> set[str]:
    gene_of_accession = dict(zip(uniprot_table["Entry"], uniprot_table["Gene Names (primary)"]))
    genes = set()
    with gzip.open(gaf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[8] != "F" or fields[3].startswith("NOT"):
                continue
            if TRANSMEMBRANE_TRANSPORTER_ACTIVITY in (go_ancestors.get(fields[4]) or {fields[4]}) and isinstance(gene_of_accession.get(fields[1]), str):
                genes.add(gene_of_accession[fields[1]])
    return genes


MEMBERSHIP_DEPTH = 3  # a protein -> its Reactome entities -> the entities those are part of
CATALYSIS_RELATIONS = ("catalyzed_by", "catalyzes_reverse_transport")


def transport_connectivity(nodes: pd.DataFrame, edges: pd.DataFrame) -> tuple[dict[str, dict], int, int]:
    """gene symbol -> {transport_reactions, metabolites moved, other_edges}; plus the number of transport reactions and
    of those without a catalysing gene.

    A gene acts through its protein nodes in a split graph (encodes edges) and through its own node in a merged one,
    and through the Reactome entities those are members of, to MEMBERSHIP_DEPTH levels, since a Reactome transport
    reaction is catalysed by an entity rather than by a gene. Until 17 October 2026 only a gene node's own catalysis
    counted, which on a split graph, where catalysis leaves the protein, found no transporter connected at all."""
    is_transport = dict(zip(nodes.node_id, nodes.is_transport.astype("boolean").fillna(False).astype(bool)))
    display = dict(zip(nodes.node_id, nodes.display_name))
    reaction_metabolites: dict[str, set[str]] = defaultdict(set)
    proteins_of_gene: dict[str, list[str]] = defaultdict(list)
    out_edges: dict[str, list[tuple[str, str]]] = defaultdict(list)
    touching: dict[str, int] = defaultdict(int)
    for source, target, relation in zip(edges.source_id, edges.target_id, edges.relation_type):
        source, target, relation = str(source), str(target), str(relation)
        if relation == "substrate_of" and is_transport.get(target):
            reaction_metabolites[target].add(str(display.get(source, source)))
        elif relation == "product_of" and is_transport.get(source):
            reaction_metabolites[source].add(str(display.get(target, target)))
        if relation == "encodes":
            proteins_of_gene[source].append(target)
            continue
        out_edges[source].append((target, relation))
        if relation not in CATALYSIS_RELATIONS:
            touching[source] += 1
            touching[target] += 1
    catalysts: dict[str, set[str]] = defaultdict(set)
    connectivity: dict[str, dict] = {}
    for gene_node in (str(node) for node in nodes.node_id if str(node).startswith("GENE:")):
        actors = proteins_of_gene.get(gene_node) or [gene_node]
        acting, frontier = set(actors), list(actors)
        for _ in range(MEMBERSHIP_DEPTH):
            frontier = [target for node in frontier for target, relation in out_edges.get(node, []) if relation == "member_of" and target not in acting]
            acting.update(frontier)
        entry = {"transport_reactions": set(), "metabolites": set(), "other_edges": sum(touching.get(actor, 0) for actor in actors)}
        for node in acting:
            for target, relation in out_edges.get(node, []):
                if relation in CATALYSIS_RELATIONS and is_transport.get(target):
                    entry["transport_reactions"].add(target)
                    entry["metabolites"] |= reaction_metabolites.get(target, set())
                    catalysts[target].add(gene_node)
        connectivity[gene_node[5:]] = entry
    transport_reactions = [node for node, flag in is_transport.items() if flag]
    return connectivity, len(transport_reactions), sum(1 for reaction in transport_reactions if not catalysts.get(reaction))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--uniprot-table", type=Path, default=Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz"))
    parser.add_argument("--gaf", type=Path, default=Path("data/raw/gene_ontology/goa_human.gaf.gz"))
    parser.add_argument("--go-obo", type=Path, default=Path("data/raw/gene_ontology/go-basic.obo"))
    parser.add_argument("--graphs", nargs="+", default=["data/processed/graph_full", "data/processed/graph_full_neuronal_split_binders_drugs_v5"])
    parser.add_argument("--markdown-output", type=Path, default=None)
    arguments = parser.parse_args()
    transporters = transporter_genes(pd.read_csv(arguments.uniprot_table, sep="\t"), arguments.gaf, read_go_ancestors(arguments.go_obo))
    stamp = provenance_stamp("experiments/audit_transporter_connectivity.py", [*arguments.graphs, str(arguments.uniprot_table), str(arguments.gaf), str(arguments.go_obo)])
    lines = ["# Transporter connectivity audit (generated by experiments/audit_transporter_connectivity.py)", "", stamp, "",
             f"{len(transporters):,} reviewed human proteins carry transmembrane transporter activity (GO:0022857, propagated). "
             "'Connected down' means the gene catalyses at least one transport reaction, which moves a metabolite between compartment copies: "
             "through its protein node in a split graph, its own node in a merged one, or a Reactome entity it is a member of.", ""]
    per_graph = {}
    for graph in arguments.graphs:
        nodes, edges = pd.read_parquet(Path(graph) / "nodes.parquet"), pd.read_parquet(Path(graph) / "edges.parquet")
        gene_nodes = {node[5:] for node in nodes.node_id if str(node).startswith("GENE:")}
        connectivity, transport_reaction_count, orphan_reactions = transport_connectivity(nodes, edges)
        per_graph[graph] = (gene_nodes, connectivity)
        in_graph = transporters & gene_nodes
        connected = {gene for gene in in_graph if connectivity.get(gene, {}).get("transport_reactions")}
        signalling_only = {gene for gene in in_graph - connected if connectivity.get(gene, {}).get("other_edges")}
        lines += [f"## {graph}", "",
                  f"- transporters that are gene nodes: {len(in_graph):,} of {len(transporters):,}",
                  f"- connected down (catalyse a transport reaction): {len(connected):,}",
                  f"- gene nodes catalysing no transport reaction, linked only by signalling or binding edges: {len(signalling_only):,}",
                  f"- gene nodes catalysing no transport reaction and with no other edge: {len(in_graph - connected - signalling_only):,}",
                  f"- transport reactions: {transport_reaction_count:,}, of which without any catalysing gene: {orphan_reactions:,}", ""]
    lines += ["## Transporters that matter here", "", "| group | gene | " + " | ".join(f"{Path(graph).name}: transport reactions (metabolites moved)" for graph in arguments.graphs) + " |",
              "|---|---|" + "---|" * len(arguments.graphs)]
    for group, genes in KEY_TRANSPORTERS.items():
        for gene in genes:
            cells = []
            for graph in arguments.graphs:
                gene_nodes, connectivity = per_graph[graph]
                entry = connectivity.get(gene)
                if gene not in gene_nodes:
                    cells.append("not a node")
                elif not entry or not entry["transport_reactions"]:
                    cells.append(f"no transport reaction ({entry['other_edges'] if entry else 0} other edges)")
                else:
                    moved = sorted(entry["metabolites"] - {"H+", "Na+", "K+", "Cl-", "H2O", "ATP", "ADP", "Pi"})
                    cells.append(f"{len(entry['transport_reactions'])} ({', '.join(moved[:8])}{', ...' if len(moved) > 8 else ''})")
            lines.append(f"| {group} | {gene} | " + " | ".join(cells) + " |")
    report = "\n".join(lines) + "\n"
    print(report)
    if arguments.markdown_output:
        arguments.markdown_output.write_text(report)


if __name__ == "__main__":
    main()
