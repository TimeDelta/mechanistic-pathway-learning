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
  signaling      OmniPath post-translational interactions (omnipath, pathwayextra, ligrecextra):
                 activates or inhibits by the consensus sign, binds when the sign is unknown
                 (binds edges are added in both directions). Protein complexes written as
                 A_B_C are expanded to their member genes.
  small molecule OmniPath small_molecule_protein interactions joined to Human-GEM metabolites
                 by PubChem identifier (metabolites.tsv): binds edges from a metabolite node
                 (its extracellular instance when present, else every compartment instance)
                 to the receptor or enzyme node, so dopamine reaches DRD1-5 and so on.
  transcription  CollecTRI regulons: regulates_transcription_of with the consensus sign.
  pharmacology   ChEMBL mechanisms: targets, with action type                    (joined downstream)
  attributes     currency tag per metabolite (tag_currency_metabolites); brain expression (not yet)

Gene and protein nodes share one identifier space, GENE:<symbol>: metabolic genes with a
symbol in genes.tsv are remapped from their Ensembl identifier, and signaling or
transcription participants that are not in Human-GEM become gene nodes with
in_metabolic_layer false. Metabolic genes without a symbol keep the Ensembl identifier.

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


def subsystem_by_reaction(model, model_yaml_path: Path | None = None) -> dict[str, str]:
    """Reaction id -> subsystem. SBML groups when the file has them; else the Human-GEM yml, whose reaction
    entries carry a subsystem list (the first entry is kept), parsed line by line to avoid loading the whole document."""
    subsystems: dict[str, str] = {}
    for group in getattr(model, "groups", []):
        for member in group.members:
            subsystems.setdefault(member.id, group.name)
    if subsystems or model_yaml_path is None or not model_yaml_path.exists():
        return subsystems
    current_reaction: str | None = None
    in_subsystem_list = False
    with open(model_yaml_path, encoding="utf-8") as yaml_file:
        for raw_line in yaml_file:
            line = raw_line.strip()
            if line.startswith("- id:"):
                current_reaction = line[len("- id:"):].strip().strip('"').strip("'")
                in_subsystem_list = False
            elif line.startswith("- subsystem:") or line.startswith("subsystem:"):
                remainder = line.split(":", 1)[1].strip()
                in_subsystem_list = remainder == ""
                if remainder and current_reaction is not None:
                    subsystems.setdefault(current_reaction, remainder.strip("[]").split(",")[0].strip().strip('"').strip("'"))
            elif in_subsystem_list and line.startswith("- ") and current_reaction is not None:
                subsystems.setdefault(current_reaction, line[2:].strip().strip('"').strip("'"))
                in_subsystem_list = False
            elif in_subsystem_list and not line.startswith("- "):
                in_subsystem_list = False
    return subsystems


def gene_node_id(symbol: str) -> str:
    return f"GENE:{symbol}"


def build_metabolic_layer(model, symbol_by_ensembl: dict[str, str], model_yaml_path: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    subsystems = subsystem_by_reaction(model, model_yaml_path)
    node_rows: list[dict] = []
    edge_rows: list[dict] = []

    def metabolic_gene_node(gene_identifier: str) -> str:
        symbol = symbol_by_ensembl.get(gene_identifier)
        return gene_node_id(symbol) if symbol else gene_identifier

    for metabolite in model.metabolites:
        node_rows.append({
            "node_id": metabolite.id,
            "node_type": "metabolite",
            "display_name": metabolite.name,
            "compartment": metabolite.compartment,
            "base_metabolite_id": metabolite.id[:-1] if metabolite.id[-1].isalpha() else metabolite.id,
            "degree": len(metabolite.reactions),
        })
    seen_gene_nodes: set[str] = set()
    for gene in model.genes:
        node_identifier = metabolic_gene_node(gene.id)
        if node_identifier in seen_gene_nodes:
            continue
        seen_gene_nodes.add(node_identifier)
        node_rows.append({
            "node_id": node_identifier,
            "node_type": "gene",
            "display_name": symbol_by_ensembl.get(gene.id, gene.id),
            "gene_symbol": symbol_by_ensembl.get(gene.id),
            "ensembl_gene_id": gene.id,
            "in_metabolic_layer": True,
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
            edge_rows.append({"source_id": metabolic_gene_node(gene.id), "target_id": reaction.id, "relation_type": "catalyzed_by", "sign": 1.0, "evidence_source": "Human-GEM GPR"})
    nodes = pd.DataFrame(node_rows)
    edges = pd.DataFrame(edge_rows).drop_duplicates()
    return nodes, edges


def expand_complex(symbol_field: str) -> list[str]:
    """OmniPath writes complexes as A_B_C; return the member symbols (a plain symbol returns itself)."""
    return [member for member in symbol_field.split("_") if member]


def interaction_sign(row: dict) -> tuple[str, float]:
    stimulation = row.get("consensus_stimulation") == "True" or (row.get("is_stimulation") == "True" and row.get("is_inhibition") != "True")
    inhibition = row.get("consensus_inhibition") == "True" or (row.get("is_inhibition") == "True" and row.get("is_stimulation") != "True")
    if stimulation and not inhibition:
        return "activates", 1.0
    if inhibition and not stimulation:
        return "inhibits", -1.0
    return "binds", 0.0


def add_gene_nodes(nodes: pd.DataFrame, symbols: set[str]) -> pd.DataFrame:
    existing = set(nodes.node_id)
    new_rows = [{"node_id": gene_node_id(symbol), "node_type": "gene", "display_name": symbol, "gene_symbol": symbol, "in_metabolic_layer": False, "degree": 0}
                for symbol in sorted(symbols) if gene_node_id(symbol) not in existing]
    return pd.concat([nodes, pd.DataFrame(new_rows)], ignore_index=True) if new_rows else nodes


def build_signaling_layer(nodes: pd.DataFrame, interactions_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    edge_rows: list[dict] = []
    symbols: set[str] = set()
    with open(interactions_path, encoding="utf-8") as interactions_file:
        for row in csv.DictReader(interactions_file, delimiter="\t"):
            relation, sign = interaction_sign(row)
            sources = expand_complex(row["source_genesymbol"])
            targets = expand_complex(row["target_genesymbol"])
            symbols.update(sources + targets)
            for source in sources:
                for target in targets:
                    if source == target:
                        continue
                    edge_rows.append({"source_id": gene_node_id(source), "target_id": gene_node_id(target), "relation_type": relation, "sign": sign, "evidence_source": "OmniPath"})
                    if relation == "binds":
                        edge_rows.append({"source_id": gene_node_id(target), "target_id": gene_node_id(source), "relation_type": relation, "sign": sign, "evidence_source": "OmniPath"})
    return add_gene_nodes(nodes, symbols), pd.DataFrame(edge_rows).drop_duplicates()


def build_transcription_layer(nodes: pd.DataFrame, regulons_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    edge_rows: list[dict] = []
    symbols: set[str] = set()
    with open(regulons_path, encoding="utf-8") as regulons_file:
        for row in csv.DictReader(regulons_file, delimiter="\t"):
            _, sign = interaction_sign(row)
            sources = expand_complex(row["source_genesymbol"])
            targets = expand_complex(row["target_genesymbol"])
            symbols.update(sources + targets)
            for source in sources:
                for target in targets:
                    if source != target:
                        edge_rows.append({"source_id": gene_node_id(source), "target_id": gene_node_id(target), "relation_type": "regulates_transcription_of", "sign": sign if sign != 0.0 else 1.0, "evidence_source": "CollecTRI"})
    return add_gene_nodes(nodes, symbols), pd.DataFrame(edge_rows).drop_duplicates()


def build_small_molecule_layer(nodes: pd.DataFrame, small_molecule_path: Path, metabolites_table_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metabolites_by_pubchem: dict[str, list[str]] = defaultdict(list)
    with open(metabolites_table_path, encoding="utf-8") as metabolites_file:
        for row in csv.DictReader(metabolites_file, delimiter="\t"):
            pubchem_id = (row.get("metPubChemID") or "").strip()
            if pubchem_id:
                metabolites_by_pubchem[pubchem_id].append(row["mets"])
    known_nodes = set(nodes.node_id)
    edge_rows: list[dict] = []
    symbols: set[str] = set()
    with open(small_molecule_path, encoding="utf-8") as small_molecule_file:
        for row in csv.DictReader(small_molecule_file, delimiter="\t"):
            metabolite_instances = [instance for instance in metabolites_by_pubchem.get(row["source"], []) if instance in known_nodes]
            if not metabolite_instances:
                continue
            extracellular = [instance for instance in metabolite_instances if instance.endswith("e")]
            chosen_instances = extracellular or metabolite_instances
            relation, sign = interaction_sign(row)
            for target in expand_complex(row["target_genesymbol"]):
                symbols.add(target)
                for instance in chosen_instances:
                    edge_rows.append({"source_id": instance, "target_id": gene_node_id(target), "relation_type": relation if relation != "binds" else "binds", "sign": sign, "evidence_source": "OmniPath small_molecule_protein"})
    return add_gene_nodes(nodes, symbols), pd.DataFrame(edge_rows).drop_duplicates()


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
    genes = nodes[nodes.node_type == "gene"]
    return {
        "nodes_by_type": {key: int(value) for key, value in nodes.node_type.value_counts().items()},
        "gene_nodes_in_metabolic_layer": int(genes["in_metabolic_layer"].eq(True).sum()) if "in_metabolic_layer" in genes else 0,
        "edges_by_relation": {key: int(value) for key, value in edges.relation_type.value_counts().items()},
        "metabolites_by_compartment": {key: int(value) for key, value in metabolites.compartment.value_counts().items()},
        "currency_metabolites": int(metabolites.is_currency.sum()),
        "transport_reactions": int(reactions.is_transport.sum()),
        "reversible_reactions": int(reactions.reversible.sum()),
        "reactions_without_gene": int((reactions.gene_reaction_rule.fillna("") == "").sum()),
        "largest_component_fraction_all": largest_connected_component_fraction(nodes, edges, exclude_currency=False),
        "largest_component_fraction_without_currency": largest_connected_component_fraction(nodes, edges, exclude_currency=True),
    }


def build_physiology_graph(
    model_path: Path,
    output_directory: Path,
    genes_table_path: Path | None = None,
    degree_quantile_threshold: float = 0.995,
    signaling_interactions_path: Path | None = None,
    transcription_regulons_path: Path | None = None,
    small_molecule_path: Path | None = None,
    metabolites_table_path: Path | None = None,
    model_yaml_path: Path | None = None,
) -> dict:
    model = load_human_gem(model_path)
    nodes, edges = build_metabolic_layer(model, gene_symbols_from_table(genes_table_path), model_yaml_path)
    edge_tables = [edges]
    if signaling_interactions_path is not None and signaling_interactions_path.exists():
        nodes, signaling_edges = build_signaling_layer(nodes, signaling_interactions_path)
        edge_tables.append(signaling_edges)
    if transcription_regulons_path is not None and transcription_regulons_path.exists():
        nodes, transcription_edges = build_transcription_layer(nodes, transcription_regulons_path)
        edge_tables.append(transcription_edges)
    if small_molecule_path is not None and small_molecule_path.exists() and metabolites_table_path is not None:
        nodes, small_molecule_edges = build_small_molecule_layer(nodes, small_molecule_path, metabolites_table_path)
        edge_tables.append(small_molecule_edges)
    edges = pd.concat(edge_tables, ignore_index=True).drop_duplicates()
    edges = edges[edges.source_id.isin(set(nodes.node_id)) & edges.target_id.isin(set(nodes.node_id))]
    degree = pd.concat([edges.source_id, edges.target_id]).value_counts()
    nodes = nodes.copy()
    nodes["degree"] = nodes.node_id.map(degree).fillna(0).astype(int)
    nodes = tag_currency(nodes, degree_quantile_threshold)
    output_directory.mkdir(parents=True, exist_ok=True)
    import os

    for name, table in (("nodes.parquet", nodes), ("edges.parquet", edges)):  # atomic replace so a concurrent reader never sees a partial file
        table.to_parquet(output_directory / (name + ".tmp"), index=False)
        os.replace(output_directory / (name + ".tmp"), output_directory / name)
    (output_directory / "relation_types.json").write_text(json.dumps(RELATION_TYPES))
    summary = summarize(nodes, edges)
    reactions = nodes[nodes.node_type == "reaction"]
    summary["reactions_with_subsystem"] = int((reactions.subsystem.fillna("") != "").sum())
    summary["distinct_subsystems"] = int(reactions.subsystem.fillna("").replace("", pd.NA).nunique())
    (output_directory / "graph_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model-path", type=Path, default=Path("data/raw/Human-GEM/model/Human-GEM.xml"))
    parser.add_argument("--genes-table", type=Path, default=Path("data/raw/Human-GEM/model/genes.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--degree-quantile-threshold", type=float, default=0.995)
    parser.add_argument("--signaling", type=Path, default=Path("data/raw/omnipath/omnipath_interactions.tsv"))
    parser.add_argument("--transcription", type=Path, default=Path("data/raw/omnipath/collectri_interactions.tsv"))
    parser.add_argument("--small-molecule", type=Path, default=Path("data/raw/omnipath/small_molecule_protein.tsv"))
    parser.add_argument("--metabolites-table", type=Path, default=Path("data/raw/Human-GEM/model/metabolites.tsv"))
    parser.add_argument("--model-yaml", type=Path, default=Path("data/raw/Human-GEM/model/Human-GEM.yml"), help="source of reaction subsystems when the SBML has no groups")
    arguments = parser.parse_args()
    summary = build_physiology_graph(arguments.model_path, arguments.output_dir, arguments.genes_table, arguments.degree_quantile_threshold,
                                     arguments.signaling, arguments.transcription, arguments.small_molecule, arguments.metabolites_table, arguments.model_yaml)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
