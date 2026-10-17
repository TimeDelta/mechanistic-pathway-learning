"""Write docs/graph_summary.md: what the physiology graph holds, read from its files, and the evidence tables it is
trained against.

Every number on the page is computed from the files named in its stamp, so the page says what the data on disk are
and not what a session remembers (CLAUDE.md, "What is current is generated, not remembered"). The first graph named
is the one the registered models read; any further graph is summarised by how it differs from the first.

Usage:
  python experiments/summarise_graph.py
  python experiments/summarise_graph.py --graph-dirs data/processed/graph_full_neuronal_split_binders_drugs_v4 \\
      --evidence-dirs data/processed/evidence_full_v4 data/processed/evidence_full_v5 --output docs/graph_summary.md
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.documentation import provenance_stamp
from mechanistic_pathway_learning.graph.edge_evidence import edge_priors, read_edge_tiers
from mechanistic_pathway_learning.graph.tag_currency_metabolites import currency_flags

DEFAULT_GRAPHS = ["data/processed/graph_full_neuronal_split_binders_drugs_v4", "data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes"]
DEFAULT_EVIDENCE = ["data/processed/evidence_full_v4", "data/processed/evidence_full_v5"]
LABEL_GRADES = ("A", "B")


def short_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def source_family(evidence_source: str) -> str:
    """The head of an edge's evidence_source, without the parenthesised detail (counts, PubMed ids, carriers)."""
    text = str(evidence_source)
    return text.split(" (")[0].split(";")[0].strip() or "(none)"


def markdown_table(header: list[str], rows: list[list]) -> list[str]:
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(cell) for cell in row) + " |" for row in rows]
    return lines


def count(value) -> str:
    return f"{int(value):,}"


def node_section(nodes: pd.DataFrame, edges: pd.DataFrame) -> list[str]:
    lines = ["## Nodes", ""]
    degree_in = edges.target_id.value_counts()
    degree_out = edges.source_id.value_counts()
    rows = []
    for node_type, group in sorted(nodes.groupby("node_type"), key=lambda item: -len(item[1])):
        total = group.node_id.map(lambda node: degree_in.get(node, 0) + degree_out.get(node, 0))
        rows.append([node_type, count(len(group)), f"{100 * len(group) / len(nodes):.1f}", f"{np.median(total):.0f}", f"{np.percentile(total, 90):.0f}", count(total.max())])
    lines += markdown_table(["node type", "nodes", "percent", "median degree", "90th percentile", "largest"], rows)
    flags = []
    for column, label in (("is_currency", "currency metabolites (they receive but send nothing in the linear-response encoder)"),
                          ("is_transport", "transport reactions"), ("reversible", "reversible reactions"),
                          ("brain_expressed", "genes and proteins expressed in at least one brain tissue (GTEx)"),
                          ("in_metabolic_layer", "nodes of the metabolic layer")):
        if column == "is_currency" and column in nodes.columns:  # metabolites only, as the loader reads it (the fourth currency defect)
            flags.append(f"{count(int(currency_flags(nodes).sum()))} {label}")
        elif column in nodes.columns:
            flags.append(f"{count((nodes[column] == True).sum())} {label}")  # noqa: E712 - the columns hold NaN
    lines += ["", "Of these: " + "; ".join(flags) + "."]
    if "compartment" in nodes.columns:
        compartments = nodes.compartment.fillna("").astype(str)
        single = compartments[(compartments != "") & ~compartments.str.contains(";")].value_counts()
        lines += ["", f"Compartments of the nodes that have one ({count((compartments != '').sum())} nodes; "
                  f"{count(compartments.str.contains(';').sum())} transport reactions span several): "
                  + ", ".join(f"{name} {count(value)}" for name, value in single.items()) + "."]
    total_degree = degree_in.add(degree_out, fill_value=0).sort_values(ascending=False)
    names = dict(zip(nodes.node_id, nodes.display_name.fillna(nodes.node_id))) if "display_name" in nodes.columns else {}
    types = dict(zip(nodes.node_id, nodes.node_type))
    hubs = [f"{names.get(node, node)} ({types.get(node, '?')}, {count(value)})" for node, value in total_degree.head(10).items()]
    lines += ["", "The ten nodes with the most edges: " + "; ".join(hubs) + "."]
    return lines


def relation_section(edges: pd.DataFrame, relation_types: list[str]) -> list[str]:
    lines = ["## Edges by relation", ""]
    rows = []
    for relation in sorted(relation_types, key=lambda name: -(edges.relation_type == name).sum()):
        group = edges[edges.relation_type == relation]
        if not len(group):
            rows.append([relation, "0", "", "", "", ""])
            continue
        sign = group["sign"].fillna(1.0) if "sign" in group.columns else pd.Series(1.0, index=group.index)
        families = group.evidence_source.map(source_family).value_counts() if "evidence_source" in group.columns else pd.Series(dtype=int)
        sources = ", ".join(f"{name} {count(value)}" for name, value in families.head(3).items()) + (f", {len(families) - 3} more" if len(families) > 3 else "")
        weighted = int((group["weight"].fillna(1.0) != 1.0).sum()) if "weight" in group.columns else 0
        rows.append([relation, count(len(group)), count((sign > 0).sum()), count((sign < 0).sum()), count((sign == 0).sum()) if (sign == 0).any() else "",
                     sources + (f"; {count(weighted)} with a weight other than 1" if weighted else "")])
    lines += markdown_table(["relation", "edges", "positive sign", "negative sign", "sign 0", "sources (most edges first)"], rows)
    if "stoichiometric_coefficient" in edges.columns:
        coefficient = edges.stoichiometric_coefficient.dropna()
        lines += ["", f"Stoichiometric coefficients: {count((coefficient > 1).sum())} reaction edges above 1 and {count((coefficient < 1).sum())} below 1 "
                  f"(of {count(len(coefficient))} that state one)."]
    return lines


def drug_section(graph: Path) -> list[str]:
    path = graph / "drug_entry_summary.json"
    if not path.exists():
        return []
    summary = json.loads(path.read_text())
    signs = ", ".join(f"{sign} {count(value)}" for sign, value in sorted(summary.get("mechanism_edges_by_sign", {}).items(), key=lambda item: float(item[0])))
    carriers = ", ".join(f"{gene} {count(value)}" for gene, value in summary.get("sequestration_edges_by_carrier_gene", {}).items())
    lines = ["## The drug entry layer", "",
             f"{count(summary.get('drug_nodes_added', 0))} drug nodes, each present only in a perturbation seeded on it, with "
             f"{count(summary.get('mechanism_edges_added', 0))} mechanism edges (`targets`, by sign times magnitude: {signs}) and "
             f"{count(summary.get('sequestration_edges_added', 0))} sequestration edges from plasma carriers ({carriers}) for "
             f"{count(summary.get('drugs_given_a_sequestration_edge', 0))} drugs with a known carrier."]
    if summary.get("clearance_edges_added") is not None:
        lines += ["", f"{count(summary['clearance_edges_added'])} clearance edges from clearing enzymes."]
    return lines


def edge_evidence_section(graph: Path, tiers_path: Path) -> list[str]:
    path = graph / "edge_evidence.parquet"
    if not path.exists():
        return []
    evidence = pd.read_parquet(path)
    gated, probability, strength = edge_priors(evidence, read_edge_tiers(tiers_path))
    rows = []
    for (family, kind), group in evidence.assign(gated=gated, probability=probability).groupby(["family", "kind"]):
        rows.append([family, kind, count(len(group)), "yes" if group.gated.any() else "no", f"{group.probability[group.gated].mean():.2f}" if group.gated.any() else ""])
    rows.sort(key=lambda row: -int(row[2].replace(",", "")))
    lines = ["## Evidence per edge and the edge priors (decision 9.57)", "",
             f"{count(gated.sum())} of {count(len(evidence))} edges carry an edge prior under `{tiers_path}`; "
             f"{count(evidence.contested.sum())} of them lie between proteins IntAct curators recorded as not interacting.", ""]
    lines += markdown_table(["family", "evidence kind", "edges", "gated", "mean prior"], rows)
    return lines


def evidence_section(evidence_directory: Path) -> list[str]:
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet")
    induces = records[records.relation == "induces"]
    if "positive_report_count" in induces.columns:
        induces = induces[induces.positive_report_count > 0]
    labels = induces[induces.grade.isin(LABEL_GRADES)]
    pairs = labels.drop_duplicates(["perturbation_id", "symptom"])
    by_type = induces.drop_duplicates("perturbation_id").perturbation_type.value_counts()
    labelled_by_type = labels.drop_duplicates("perturbation_id").perturbation_type.value_counts()
    per_symptom = pairs.symptom.value_counts()
    lines = [f"### `{evidence_directory}`", "",
             f"{count(len(records))} records; {count(induces.perturbation_id.nunique())} perturbations with a positive record of a symptom they induce ("
             + ", ".join(f"{kind} {count(value)}" for kind, value in by_type.items()) + f"), over {count(induces.symptom.nunique())} symptoms. "
             f"Label pairs (grades {', '.join(LABEL_GRADES)}, before any label selection): {count(len(pairs))} over "
             f"{count(labels.perturbation_id.nunique())} perturbations (" + ", ".join(f"{kind} {count(value)}" for kind, value in labelled_by_type.items()) + ").", "",
             "Label pairs per symptom, most first: " + ", ".join(f"{symptom} {count(value)}" for symptom, value in per_symptom.items()) + "."]
    return lines


def variant_section(first: Path, other: Path) -> list[str]:
    nodes_a, nodes_b = pd.read_parquet(first / "nodes.parquet"), pd.read_parquet(other / "nodes.parquet")
    edges_a, edges_b = pd.read_parquet(first / "edges.parquet"), pd.read_parquet(other / "edges.parquet")
    relation_difference = edges_b.relation_type.value_counts().subtract(edges_a.relation_type.value_counts(), fill_value=0)
    changed = ", ".join(f"{relation} {int(value):+,}" for relation, value in relation_difference[relation_difference != 0].items())
    return [f"### `{other}` against the first", "",
            f"{count(len(nodes_b))} nodes ({int(len(nodes_b) - len(nodes_a)):+,}) and {count(len(edges_b))} edges ({int(len(edges_b) - len(edges_a)):+,}); "
            f"the relations that differ: {changed or 'none'}."]


def summary_page(graphs: list[Path], evidence_directories: list[Path], tiers_path: Path, inputs: list[str]) -> str:
    first = graphs[0]
    nodes, edges = pd.read_parquet(first / "nodes.parquet"), pd.read_parquet(first / "edges.parquet")
    relation_types = json.loads((first / "relation_types.json").read_text())
    lines = [provenance_stamp("experiments/summarise_graph.py", inputs), "", "# The physiology graph", "",
             f"The graph the registered models read is `{first}`: **{count(len(nodes))} nodes of {nodes.node_type.nunique()} types and "
             f"{count(len(edges))} edges of {len(relation_types)} relations** (`nodes.parquet` {short_sha256(first / 'nodes.parquet')}, "
             f"`edges.parquet` {short_sha256(first / 'edges.parquet')}). A gene and the protein it encodes are separate nodes joined by an "
             f"`encodes` edge; a drug is a node of its own that exists only in a perturbation seeded on it.", ""]
    lines += node_section(nodes, edges) + [""] + relation_section(edges, relation_types) + [""]
    lines += drug_section(first) + [""] + edge_evidence_section(first, tiers_path) + [""]
    if len(graphs) > 1:
        lines += ["## Other graph directories", ""]
        for other in graphs[1:]:
            lines += variant_section(first, other) + [""]
    lines += ["## The evidence tables the graph is trained against", ""]
    for directory in evidence_directories:
        lines += evidence_section(directory) + [""]
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dirs", type=Path, nargs="+", default=[Path(path) for path in DEFAULT_GRAPHS])
    parser.add_argument("--evidence-dirs", type=Path, nargs="+", default=[Path(path) for path in DEFAULT_EVIDENCE])
    parser.add_argument("--tiers", type=Path, default=Path("configs/edge_prior_tiers.csv"))
    parser.add_argument("--output", type=Path, default=Path("docs/graph_summary.md"))
    arguments = parser.parse_args()
    inputs = [str(path) for path in arguments.graph_dirs] + [str(path) for path in arguments.evidence_dirs] + [str(arguments.tiers)]
    page = summary_page(arguments.graph_dirs, arguments.evidence_dirs, arguments.tiers, inputs)
    arguments.output.write_text(page, encoding="utf-8")
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
