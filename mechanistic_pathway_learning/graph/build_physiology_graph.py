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
Brain expression (GTEx median TPM, --brain-expression) annotates gene nodes with
brain_median_tpm_max, brain_expression_known and brain_expressed, gives reactions the largest
value over their genes and drops regulon edges of factors known not to be brain expressed.
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


GTEX_BRAIN_TISSUE_PREFIX = "Brain_"


def load_brain_expression(gct_path: Path) -> pd.DataFrame:
    """Per-gene brain expression from a GTEx median-TPM GCT file (design section 4.1, reference [30]).

    Returns one row per gene with ensembl_gene_id (version suffix removed), gene_symbol,
    brain_median_tpm_max (the largest median TPM over the brain tissues), brain_median_tpm_mean and
    brain_tissue_count. GTEx v10 has 13 brain tissues; any column whose name starts with Brain_ counts.
    """
    table = pd.read_csv(gct_path, sep="\t", skiprows=2)
    brain_columns = [column for column in table.columns if column.startswith(GTEX_BRAIN_TISSUE_PREFIX)]
    if not brain_columns:
        raise ValueError(f"{gct_path} has no columns starting with {GTEX_BRAIN_TISSUE_PREFIX}")
    expression = pd.DataFrame({
        "ensembl_gene_id": table["Name"].astype(str).str.split(".").str[0],
        "gene_symbol": table["Description"].astype(str),
        "brain_median_tpm_max": table[brain_columns].max(axis=1).astype(float),
        "brain_median_tpm_mean": table[brain_columns].mean(axis=1).astype(float),
    })
    expression["brain_tissue_count"] = len(brain_columns)
    # a symbol can map to several Ensembl ids (and vice versa); keep the highest brain expression per key
    return expression.sort_values("brain_median_tpm_max", ascending=False)


def annotate_brain_expression(nodes: pd.DataFrame, expression: pd.DataFrame, minimum_tpm: float) -> pd.DataFrame:
    """Add brain_median_tpm_max, brain_expression_known and brain_expressed to gene nodes.

    A gene node is matched by its Ensembl id when it has one (metabolic genes) and by symbol otherwise
    (signaling and transcription genes). brain_expressed is median TPM >= minimum_tpm in at least one
    brain tissue; it is False when the gene is known and below the threshold and also False when the
    gene is not in GTEx, with brain_expression_known separating the two cases. Reaction nodes inherit the
    largest value over the genes that catalyze them, so a reaction's brain weight is available to readers.
    """
    nodes = nodes.copy()
    by_ensembl = expression.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["brain_median_tpm_max"]
    by_symbol = expression.drop_duplicates("gene_symbol").set_index("gene_symbol")["brain_median_tpm_max"]
    gene_mask = nodes.node_type == "gene"
    ensembl_values = nodes.loc[gene_mask, "ensembl_gene_id"].map(by_ensembl) if "ensembl_gene_id" in nodes else pd.Series(index=nodes.index[gene_mask], dtype=float)
    symbol_values = nodes.loc[gene_mask, "gene_symbol"].map(by_symbol)
    brain_tpm = ensembl_values.where(ensembl_values.notna(), symbol_values)
    nodes["brain_median_tpm_max"] = pd.NA
    nodes.loc[gene_mask, "brain_median_tpm_max"] = brain_tpm.values
    nodes["brain_median_tpm_max"] = pd.to_numeric(nodes["brain_median_tpm_max"], errors="coerce")
    nodes["brain_expression_known"] = nodes["brain_median_tpm_max"].notna()
    nodes["brain_expressed"] = nodes["brain_median_tpm_max"].fillna(-1.0) >= minimum_tpm
    return nodes


def propagate_brain_expression_to_reactions(nodes: pd.DataFrame, edges: pd.DataFrame) -> pd.DataFrame:
    """Reaction nodes take the largest brain_median_tpm_max over genes linked by catalyzed_by edges."""
    nodes = nodes.copy()
    catalysis = edges[edges.relation_type == "catalyzed_by"]
    gene_values = nodes.set_index("node_id")["brain_median_tpm_max"]
    gene_side = catalysis.target_id.where(catalysis.target_id.isin(gene_values.index), catalysis.source_id)
    reaction_side = catalysis.source_id.where(catalysis.target_id.isin(gene_values.index), catalysis.target_id)
    per_reaction = pd.DataFrame({"reaction": reaction_side.values, "value": gene_side.map(gene_values).values}).dropna().groupby("reaction")["value"].max()
    reaction_mask = nodes.node_type == "reaction"
    nodes.loc[reaction_mask, "brain_median_tpm_max"] = nodes.loc[reaction_mask, "node_id"].map(per_reaction).values
    nodes["brain_expression_known"] = nodes["brain_median_tpm_max"].notna()
    return nodes


def restrict_transcription_edges_to_brain_expressed(nodes: pd.DataFrame, edges: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Drop regulates_transcription_of edges whose transcription factor is known not to be brain expressed
    (design section 3.3). Factors absent from GTEx keep their edges, because absence from the expression
    table is a mapping gap rather than evidence of absence; the counts of both cases are returned.
    """
    known = nodes.set_index("node_id")["brain_expression_known"]
    expressed = nodes.set_index("node_id")["brain_expressed"]
    is_transcription = edges.relation_type == "regulates_transcription_of"
    factor_known = edges.source_id.map(known).fillna(False).astype(bool)
    factor_expressed = edges.source_id.map(expressed).fillna(False).astype(bool)
    dropped = is_transcription & factor_known & ~factor_expressed
    counts = {
        "transcription_edges_before": int(is_transcription.sum()),
        "transcription_edges_dropped_factor_not_brain_expressed": int(dropped.sum()),
        "transcription_edges_kept_factor_expression_unknown": int((is_transcription & ~factor_known).sum()),
        "transcription_factors_dropped": int(edges.loc[dropped, "source_id"].nunique()),
        "transcription_factors_expression_unknown": int(edges.loc[is_transcription & ~factor_known, "source_id"].nunique()),
    }
    return edges[~dropped].copy(), counts


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
    brain_expression_path: Path | None = None,
    brain_expression_minimum_tpm: float = 1.0,
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
    brain_expression_counts: dict = {}
    if brain_expression_path is not None and brain_expression_path.exists():
        nodes = annotate_brain_expression(nodes, load_brain_expression(brain_expression_path), brain_expression_minimum_tpm)
        edges, brain_expression_counts = restrict_transcription_edges_to_brain_expressed(nodes, edges)
        nodes = propagate_brain_expression_to_reactions(nodes, edges)
        genes = nodes[nodes.node_type == "gene"]
        brain_expression_counts.update({
            "brain_expression_minimum_tpm": brain_expression_minimum_tpm,
            "gene_nodes_with_brain_expression_known": int(genes.brain_expression_known.sum()),
            "gene_nodes_brain_expressed": int(genes.brain_expressed.sum()),
        })
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
    summary["brain_expression"] = brain_expression_counts
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
    parser.add_argument("--brain-expression", type=Path, default=Path("data/raw/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"), help="GTEx median-TPM GCT file; gene nodes get brain expression attributes and transcription factors known not to be brain expressed lose their regulon edges")
    parser.add_argument("--brain-expression-minimum-tpm", type=float, default=1.0, help="median TPM in at least one brain tissue that counts as brain expression evidence")
    arguments = parser.parse_args()
    summary = build_physiology_graph(arguments.model_path, arguments.output_dir, arguments.genes_table, arguments.degree_quantile_threshold,
                                     arguments.signaling, arguments.transcription, arguments.small_molecule, arguments.metabolites_table, arguments.model_yaml,
                                     arguments.brain_expression, arguments.brain_expression_minimum_tpm)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
