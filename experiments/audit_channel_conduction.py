"""Audit: does every ion channel in the graph conduct, that is, catalyse a reaction that moves ions into the
membrane potential or a neuronal ion pool?

The electrical layer (mechanistic_pathway_learning/graph/reactome_import.py) reads a channel's effect from the
reactions it catalyses: a reaction that moves net charge across the plasma membrane feeds VM_c
(changes_membrane_potential), one that moves an ion into a neuronal pool feeds the pool, and VM_c gates the reactions
of voltage-gated pores (voltage_gates). A channel protein that catalyses no such reaction, directly or as a member of
a Reactome entity, can reach the rest of the graph only through signalling or binding edges, whatever its gene does
in a neuron. Decision 9.69 found the CaV1 and CaV3 calcium channels in that state; this audit asks the same of every
channel.

Channels are the reviewed human proteins annotated (any evidence, NOT-qualified rows excluded, propagated through
is_a and part_of) with monoatomic ion channel activity, GO:0005216. Each one that is a node of the graph is put in one
class, the first that holds:

  conducts into VM_c        a reaction it catalyses (or an entity it is a member of, to three levels, catalyses)
                            feeds VM_c;
  moves ions, no VM_c edge  such a reaction moves an ion into a neuronal pool or between compartments, with no
                            membrane-potential edge (an organelle channel, an electroneutral flux);
  catalyses another         it catalyses reactions, none of which moves an inorganic ion;
  signalling or binding     its edges out are signalling, binding or membership of entities that catalyse nothing;
  no edge out               nothing leaves it.

Usage:
  python experiments/audit_channel_conduction.py --graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v5 \\
      --markdown-output docs/channel_conduction_audit.md
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ION_CHANNEL_ACTIVITY = "GO:0005216"  # monoatomic ion channel activity
MEMBERSHIP_DEPTH = 3
INORGANIC_ION_BASES = {"MAM02519": "Na+", "MAM02200": "K+", "MAM01413": "Ca2+", "MAM01442": "Cl-", "MAM02482": "Mg2+", "MAM02039": "H+"}
CLASSES = ("conducts into VM_c", "moves ions, no VM_c edge", "catalyses another reaction", "signalling or binding only", "no edge out")
CATALYSIS_RELATIONS = ("catalyzed_by", "catalyzes_reverse_transport")
# families by gene-symbol prefix, for the table; first match wins
FAMILIES = (
    ("voltage-gated sodium (Nav)", r"^SCN\d+A$"), ("sodium channel beta subunits", r"^SCN\d+B$"),
    ("voltage-gated calcium alpha1 (CaV)", r"^CACNA1"), ("calcium channel auxiliary", r"^CACN[ABG]"),
    ("voltage-gated potassium (Kv, KCNQ, hERG, BK, SK)", r"^KCN[ABCDFGHQSVTMNU]"), ("inward rectifier potassium (Kir)", r"^KCNJ"),
    ("two-pore potassium (K2P)", r"^KCNK"), ("potassium channel auxiliary", r"^KCN[EIR]|^KCNAB|^KCNMB"), ("HCN", r"^HCN"),
    ("GABA-A receptors", r"^GABR"), ("glycine receptors", r"^GLR[AB]"), ("ionotropic glutamate (AMPA, kainate, NMDA, delta)", r"^GRI[ADKN]"),
    ("nicotinic acetylcholine", r"^CHRN"), ("serotonin 5-HT3", r"^HTR3"), ("P2X purinergic", r"^P2RX"), ("acid-sensing (ASIC)", r"^ASIC"),
    ("epithelial sodium (ENaC)", r"^SCNN1"), ("TRP", r"^TRP"), ("chloride channels (CLC, CLCA, bestrophin, anoctamin)", r"^CLC|^BEST|^ANO\d"),
    ("ryanodine and IP3 receptors", r"^RYR|^ITPR"), ("ORAI and STIM", r"^ORAI|^STIM"), ("mechanosensitive (PIEZO, TMC)", r"^PIEZO|^TMC"),
    ("connexins and pannexins", r"^GJ|^PANX"), ("leak and other cation channels", r"^NALCN|^TPCN|^MCOLN|^PKD|^CATSPER|^HVCN|^CALHM|^TMEM"),
)


def family_of(symbol: str) -> str:
    return next((name for name, pattern in FAMILIES if re.match(pattern, symbol)), "other")


def channel_genes(uniprot_table: pd.DataFrame, gaf_path: Path, go_ancestors: dict[str, set[str]], term: str = ION_CHANNEL_ACTIVITY) -> set[str]:
    """Gene symbols of the reviewed human proteins annotated with the term or a descendant (NOT rows excluded)."""
    gene_of_accession = dict(zip(uniprot_table["Entry"], uniprot_table["Gene Names (primary)"]))
    genes = set()
    with gzip.open(gaf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[8] != "F" or fields[3].startswith("NOT"):
                continue
            if term in (go_ancestors.get(fields[4]) or {fields[4]}) and isinstance(gene_of_accession.get(fields[1]), str):
                genes.add(gene_of_accession[fields[1]].split(";")[0].strip())
    return genes


def reaction_ion_movement(edges: pd.DataFrame, node_ids: set[str]) -> dict[str, dict]:
    """Per reaction node: whether it feeds VM_c, whether it moves an inorganic ion (into a neuronal pool, or the same ion
    as substrate and product in two compartments), and whether VM_c gates it."""
    feeds_vm = set(edges.loc[edges.relation_type == "changes_membrane_potential", "source_id"])
    gated = set(edges.loc[edges.relation_type == "voltage_gates", "target_id"])
    to_pool = set(edges.loc[(edges.relation_type == "product_of") & edges.target_id.astype(str).str.startswith("NEURON_"), "source_id"])
    to_pool |= set(edges.loc[edges.relation_type == "changes_ion_pool", "source_id"])

    def ion_base(node: str) -> str | None:
        node = str(node)
        if node.startswith("NEURON_"):
            node = node[len("NEURON_"):]
        base = node[:-1]
        return base if base in INORGANIC_ION_BASES else None

    substrates, products = defaultdict(set), defaultdict(set)
    for source, target, relation in zip(edges.source_id, edges.target_id, edges.relation_type):
        if relation == "substrate_of" and ion_base(source):
            substrates[target].add((ion_base(source), str(source)[-1]))
        elif relation == "product_of" and ion_base(target):
            products[source].add((ion_base(target), str(target)[-1]))
    moves_ion = {reaction for reaction in set(substrates) & set(products)
                 if any(base == other_base and letter != other_letter for base, letter in substrates[reaction] for other_base, other_letter in products[reaction])}
    reactions = set(edges.loc[edges.relation_type.isin(CATALYSIS_RELATIONS), "target_id"])
    return {reaction: {"feeds_vm": reaction in feeds_vm, "moves_ion": reaction in moves_ion or reaction in to_pool, "voltage_gated": reaction in gated}
            for reaction in reactions | feeds_vm}


def classify_channels(nodes: pd.DataFrame, edges: pd.DataFrame, genes: set[str]) -> pd.DataFrame:
    """One row per channel gene that is a node: its protein nodes, the reactions they catalyse (directly or through
    entities they are members of), and its class."""
    node_ids = set(nodes.node_id)
    encodes = edges[edges.relation_type == "encodes"]
    proteins_of_gene = defaultdict(list)
    for source, target in zip(encodes.source_id, encodes.target_id):
        proteins_of_gene[str(source)].append(str(target))
    out_edges = defaultdict(list)
    for source, target, relation in zip(edges.source_id, edges.target_id, edges.relation_type):
        out_edges[str(source)].append((str(target), str(relation)))
    movement = reaction_ion_movement(edges, node_ids)
    rows = []
    for gene in sorted(genes):
        gene_node = f"GENE:{gene}"
        if gene_node not in node_ids:
            continue
        actors = proteins_of_gene.get(gene_node) or [gene_node]  # a merged graph has no protein nodes
        frontier, seen = list(actors), set(actors)
        for _ in range(MEMBERSHIP_DEPTH):
            next_frontier = []
            for node in frontier:
                for target, relation in out_edges.get(node, []):
                    if relation == "member_of" and target not in seen:
                        seen.add(target)
                        next_frontier.append(target)
            frontier = next_frontier
        catalysed = sorted({target for node in seen for target, relation in out_edges.get(node, []) if relation in CATALYSIS_RELATIONS})
        actor_edges = Counter(relation for node in actors for _, relation in out_edges.get(node, []))
        if any(movement.get(reaction, {}).get("feeds_vm") for reaction in catalysed):
            klass = CLASSES[0]
        elif any(movement.get(reaction, {}).get("moves_ion") for reaction in catalysed):
            klass = CLASSES[1]
        elif catalysed:
            klass = CLASSES[2]
        elif actor_edges:
            klass = CLASSES[3]
        else:
            klass = CLASSES[4]
        rows.append({"gene": gene, "family": family_of(gene), "class": klass, "protein_nodes": ";".join(actors),
                     "reactions_catalysed": len(catalysed), "voltage_gated_reactions": sum(movement.get(reaction, {}).get("voltage_gated", False) for reaction in catalysed),
                     "edges_out_of_the_protein": ";".join(f"{relation}:{count}" for relation, count in sorted(actor_edges.items())),
                     "reactions": ";".join(catalysed[:6])})
    return pd.DataFrame(rows, columns=["gene", "family", "class", "protein_nodes", "reactions_catalysed", "voltage_gated_reactions",
                                       "edges_out_of_the_protein", "reactions"])


def render(table: pd.DataFrame, graph_dir: Path, stamp: str, channel_count: int) -> str:
    lines = ["# Ion channel conduction audit (generated by experiments/audit_channel_conduction.py — do not edit by hand)", "", stamp, "",
             f"Graph `{graph_dir}`. {channel_count:,} reviewed human proteins carry monoatomic ion channel activity ({ION_CHANNEL_ACTIVITY}, propagated); "
             f"{len(table):,} of their genes are nodes of the graph. A channel conducts when a reaction it catalyses, directly or as a member of a "
             "Reactome entity, feeds the membrane potential (`changes_membrane_potential` into VM_c); otherwise its gene reaches the rest of the graph "
             "only through what else it touches. The classes are defined in the script's docstring.", "",
             "## By class", "", "| class | channels |", "| --- | --- |"]
    for klass in CLASSES:
        lines.append(f"| {klass} | {int((table['class'] == klass).sum()):,} |")
    lines += ["", "## By family", "", "| family | channels | conduct into VM_c | move ions, no VM_c edge | catalyse another reaction | signalling or binding only | no edge out |",
              "| --- | --- | --- | --- | --- | --- | --- |"]
    for family, group in sorted(table.groupby("family"), key=lambda item: -len(item[1])):
        counts = [int((group["class"] == klass).sum()) for klass in CLASSES]
        lines.append(f"| {family} | {len(group)} | " + " | ".join(str(count) for count in counts) + " |")
    lines += ["", "## Channels that do not conduct into VM_c", "", "| gene | family | class | edges out of the protein | reactions catalysed (first six) |", "| --- | --- | --- | --- | --- |"]
    for _, row in table[table["class"] != CLASSES[0]].sort_values(["family", "gene"]).iterrows():
        lines.append(f"| {row.gene} | {row.family} | {row['class']} | {row.edges_out_of_the_protein or '-'} | {row.reactions or '-'} |")
    return "\n".join(lines) + "\n"


def main() -> None:
    from mechanistic_pathway_learning.documentation import provenance_stamp
    from mechanistic_pathway_learning.graph.protein_descriptors import read_go_ancestors

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders_drugs_v5"))
    parser.add_argument("--uniprot-table", type=Path, default=Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz"))
    parser.add_argument("--go-annotations", type=Path, default=Path("data/raw/gene_ontology/goa_human.gaf.gz"))
    parser.add_argument("--go-ontology", type=Path, default=Path("data/raw/gene_ontology/go-basic.obo"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/channel_conduction_audit.md"))
    parser.add_argument("--table-output", type=Path, default=None, help="CSV, one row per channel gene in the graph")
    arguments = parser.parse_args()
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet")
    genes = channel_genes(pd.read_csv(arguments.uniprot_table, sep="\t", usecols=["Entry", "Gene Names (primary)"]), arguments.go_annotations,
                          read_go_ancestors(arguments.go_ontology))
    table = classify_channels(nodes, edges, genes)
    stamp = provenance_stamp("experiments/audit_channel_conduction.py", [str(arguments.graph_dir), str(arguments.uniprot_table), str(arguments.go_annotations),
                                                                        str(arguments.go_ontology)])
    arguments.markdown_output.write_text(render(table, arguments.graph_dir, stamp, len(genes)))
    if arguments.table_output is not None:
        table.to_csv(arguments.table_output, index=False)
    print(json.dumps({"channels": len(genes), "in_graph": len(table), "by_class": table["class"].value_counts().to_dict()}, indent=1))


if __name__ == "__main__":
    main()
