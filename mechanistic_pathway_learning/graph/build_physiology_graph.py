"""Physiology graph construction (design section 4.1).

Layers and sources:
  metabolic      Human-GEM (SBML from github.com/SysBioChalmers/Human-GEM) via cobra:
                 metabolite nodes (compartment-specific), reaction nodes, gene nodes,
                 edges substrate_of (metabolite -> reaction), product_of (reaction -> metabolite),
                 catalyzed_by (gene -> reaction, from the gene-protein-reaction rule; the AND/OR
                 structure of the rule is kept as a reaction attribute, not expanded, in version 1).
                 Reversible reactions get both orientations of their substrate and product edges.
                 A reaction whose metabolites span more than one compartment is tagged is_transport;
                 those reactions are the only edges that connect compartments.
  signaling      OmniPath: binds, activates, inhibits with signs                 (not yet built)
  transcription  CollecTRI: regulates_transcription_of, signed                   (not yet built)
  pharmacology   ChEMBL mechanisms: targets, with action type                    (joined downstream)
  attributes     currency tag per metabolite (tag_currency_metabolites); brain expression (not yet)

Hard structural constraints are enforced here, not in the model: no other edge types,
no symptom-symptom edges, no disease nodes, compartment changes only via transport
reactions, pharmacological edges only with measured activity.

Outputs (output_directory, default data/processed/graph/):
  nodes.parquet          node_id, node_type, display_name, compartment, is_currency, degree,
                         gene_symbol (genes), is_transport and reversible (reactions), base_metabolite_id
  edges.parquet          source_id, target_id, relation_type, sign, evidence_source
  relation_types.json    ordered list used as the encoder's relation index
  graph_summary.json     counts for the Phase 1 quality checks
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict, deque
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.tag_currency_metabolites import tag_currency_metabolites

RELATION_TYPES: list[str] = [
    "substrate_of",
    "product_of",
    "catalyzed_by",
    "binds",
    "activates",
    "inhibits",
    "regulates_transcription_of",
    "targets",
]


def load_human_gem(model_path: Path):
    import pickle

    if model_path.suffix == ".pkl":
        with open(model_path, "rb") as cached:
            return pickle.load(cached)
    from cobra.io import read_sbml_model

    return read_sbml_model(str(model_path))


def gene_symbols_from_table(genes_table_path: Path | None) -> dict[str, str]:
    if genes_table_path is None or not genes_table_path.exists():
        return {}
    symbol_by_ensembl: dict[str, str] = {}
    with open(genes_table_path, encoding="utf-8") as genes_file:
        for row in csv.DictReader(genes_file, delimiter="\t", quotechar='"'):
            symbols = [symbol.strip() for symbol in (row.get("geneSymbols") or "").split(";") if symbol.strip()]
            if symbols:
                symbol_by_ensembl[row["genes"]] = symbols[0]
    return symbol_by_ensembl


def subsystem_by_reaction(model) -> dict[str, str]:
    subsystems: dict[str, str] = {}
    for group in getattr(model, "groups", []):
        for member in group.members:
            subsystems.setdefault(member.id, group.name)
    return subsystems


def build_metabolic_layer(model, symbol_by_ensembl: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    subsystems = subsystem_by_reaction(model)
    node_rows: list[dict] = []
    edge_rows: list[dict] = []
    for metabolite in model.metabolites:
        node_rows.append({
            "node_id": metabolite.id,
            "node_type": "metabolite",
            "display_name": metabolite.name,
            "compartment": metabolite.compartment,
            "base_metabolite_id": metabolite.id[:-1] if metabolite.id[-1].isalpha() else metabolite.id,
            "degree": len(metabolite.reactions),
        })
    for gene in model.genes:
        node_rows.append({
            "node_id": gene.id,
            "node_type": "gene",
            "display_name": symbol_by_ensembl.get(gene.id, gene.id),
            "gene_symbol": symbol_by_ensembl.get(gene.id),
            "degree": len(gene.reactions),
        })
    for reaction in model.reactions:
        compartments = {metabolite.compartment for metabolite in reaction.metabolites}
        node_rows.append({
            "node_id": reaction.id,
            "node_type": "reaction",
            "display_name": reaction.name,
            "compartment": ";".join(sorted(compartments)),
            "is_transport": len(compartments) > 1,
            "reversible": reaction.reversibility,
            "gene_reaction_rule": reaction.gene_reaction_rule,
            "subsystem": subsystems.get(reaction.id, ""),
            "degree": len(reaction.metabolites) + len(reaction.genes),
        })
        for metabolite, coefficient in reaction.metabolites.items():
            is_substrate = coefficient < 0
            edge_rows.append({"source_id": metabolite.id if is_substrate else reaction.id, "target_id": reaction.id if is_substrate else metabolite.id, "relation_type": "substrate_of" if is_substrate else "product_of", "sign": 1.0, "evidence_source": "Human-GEM"})
            if reaction.reversibility:
                edge_rows.append({"source_id": reaction.id if is_substrate else metabolite.id, "target_id": metabolite.id if is_substrate else reaction.id, "relation_type": "product_of" if is_substrate else "substrate_of", "sign": 1.0, "evidence_source": "Human-GEM"})
        for gene in reaction.genes:
            edge_rows.append({"source_id": gene.id, "target_id": reaction.id, "relation_type": "catalyzed_by", "sign": 1.0, "evidence_source": "Human-GEM GPR"})
    nodes = pd.DataFrame(node_rows)
    edges = pd.DataFrame(edge_rows).drop_duplicates()
    return nodes, edges


def tag_currency(nodes: pd.DataFrame, degree_quantile_threshold: float) -> pd.DataFrame:
    metabolites = nodes[nodes.node_type == "metabolite"]
    tagged = tag_currency_metabolites(
        dict(zip(metabolites.node_id, metabolites.display_name)),
        dict(zip(metabolites.node_id, metabolites.degree)),
        degree_quantile_threshold=degree_quantile_threshold,
    )
    nodes = nodes.copy()
    nodes["is_currency"] = nodes.node_id.isin(tagged)
    return nodes


def largest_connected_component_fraction(nodes: pd.DataFrame, edges: pd.DataFrame, exclude_currency: bool) -> float:
    excluded = set(nodes.loc[nodes["is_currency"], "node_id"]) if exclude_currency and "is_currency" in nodes else set()
    adjacency: dict[str, list[str]] = defaultdict(list)
    for source, target in zip(edges.source_id, edges.target_id):
        if source in excluded or target in excluded:
            continue
        adjacency[source].append(target)
        adjacency[target].append(source)
    all_nodes = [node_id for node_id in nodes.node_id if node_id not in excluded]
    seen: set[str] = set()
    largest = 0
    for start in all_nodes:
        if start in seen:
            continue
        queue = deque([start])
        seen.add(start)
        size = 0
        while queue:
            current = queue.popleft()
            size += 1
            for neighbor in adjacency.get(current, ()):
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        largest = max(largest, size)
    return largest / max(1, len(all_nodes))


def summarize(nodes: pd.DataFrame, edges: pd.DataFrame) -> dict:
    metabolites = nodes[nodes.node_type == "metabolite"]
    reactions = nodes[nodes.node_type == "reaction"]
    return {
        "nodes_by_type": {key: int(value) for key, value in nodes.node_type.value_counts().items()},
        "edges_by_relation": {key: int(value) for key, value in edges.relation_type.value_counts().items()},
        "metabolites_by_compartment": {key: int(value) for key, value in metabolites.compartment.value_counts().items()},
        "currency_metabolites": int(metabolites.is_currency.sum()),
        "transport_reactions": int(reactions.is_transport.sum()),
        "reversible_reactions": int(reactions.reversible.sum()),
        "reactions_without_gene": int((reactions.gene_reaction_rule.fillna("") == "").sum()),
        "largest_component_fraction_all": largest_connected_component_fraction(nodes, edges, exclude_currency=False),
        "largest_component_fraction_without_currency": largest_connected_component_fraction(nodes, edges, exclude_currency=True),
    }


def build_physiology_graph(model_path: Path, output_directory: Path, genes_table_path: Path | None = None, degree_quantile_threshold: float = 0.995) -> dict:
    model = load_human_gem(model_path)
    nodes, edges = build_metabolic_layer(model, gene_symbols_from_table(genes_table_path))
    nodes = tag_currency(nodes, degree_quantile_threshold)
    output_directory.mkdir(parents=True, exist_ok=True)
    nodes.to_parquet(output_directory / "nodes.parquet", index=False)
    edges.to_parquet(output_directory / "edges.parquet", index=False)
    (output_directory / "relation_types.json").write_text(json.dumps(RELATION_TYPES))
    summary = summarize(nodes, edges)
    (output_directory / "graph_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model-path", type=Path, default=Path("data/raw/Human-GEM/model/Human-GEM.xml"))
    parser.add_argument("--genes-table", type=Path, default=Path("data/raw/Human-GEM/model/genes.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--degree-quantile-threshold", type=float, default=0.995)
    arguments = parser.parse_args()
    summary = build_physiology_graph(arguments.model_path, arguments.output_dir, arguments.genes_table, arguments.degree_quantile_threshold)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
