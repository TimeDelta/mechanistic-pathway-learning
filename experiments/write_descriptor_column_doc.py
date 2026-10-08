"""Write docs/node_descriptor_columns.md: the name of every column of the final node descriptor table (the full neuronal
graph with brain expression), which node type fills it, and a name for each protein component.

The column lists are read from the table itself, and a column is listed for a node type only when at least one node of
that type has a non-zero value in it, so the document cannot drift from the table. The protein component names come
from experiments/name_protein_descriptor_components.py. The table's SHA-256 is recorded; when it changes, rerun both
scripts, raise DOCUMENT_VERSION and add a line to CHANGES.

    python experiments/write_descriptor_column_doc.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import date
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.node_descriptors import descriptor_blocks

DOCUMENT_VERSION = 3
CHANGES = [
    (1, "2026-10-08", "First version: the 137 descriptor columns of the brain-expression tables and the 83 of the graph directories' own tables."),
    (2, "2026-10-08", "Only the final table (full neuronal graph with brain expression, used by the confirmatory configurations); "
                      "a name for each of the 64 protein components from the annotations it correlates with most."),
    (3, "2026-10-08", "The final table is the split one (the user's decision): the full neuronal graph with genes and proteins as "
                      "separate nodes, so the protein components sit on the protein nodes."),
]
# The final descriptor table, the graph directory whose nodes it describes, and where it is used.
FINAL_TABLE = "data/processed/node_descriptors/full_neuronal_split_descriptors_brain_expression.parquet"
FINAL_GRAPH_DIR = "data/processed/graph_full_neuronal_split"
FINAL_TABLE_USE = "the full graph with genes and proteins split (docs/gene_protein_split.md)"
COMPONENT_NAMES = Path("data/processed/node_descriptors/protein_rrr_component_names.json")  # experiments/name_protein_descriptor_components.py
EC_CLASS_NAMES = {1: "oxidoreductases", 2: "transferases", 3: "hydrolases", 4: "lyases", 5: "isomerases", 6: "ligases", 7: "translocases"}
METABOLITE_MEANINGS = {
    "metabolite_log_molecular_weight": "log molecular weight (RDKit MolWt), standardised",
    "metabolite_logp": "Crippen logP, standardised",
    "metabolite_net_charge": "net charge in Human-GEM (SBML fbc:charge, near pH 7.3), standardised; known for every metabolite",
    "metabolite_log_polar_surface_area": "log(1 + topological polar surface area), standardised",
    "metabolite_log_hydrogen_bond_donors": "log(1 + hydrogen-bond donors), standardised",
    "metabolite_log_hydrogen_bond_acceptors": "log(1 + hydrogen-bond acceptors), standardised",
    "metabolite_log_rotatable_bonds": "log(1 + rotatable bonds), standardised",
    "metabolite_log_rings": "log(1 + rings), standardised",
    "metabolite_has_structure": "1 when a complete SMILES gave the properties above; 0 leaves them at the mean",
    "metabolite_partial_structure": "1 when the SMILES has R-group atoms (*), e.g. the acyl chain of an acyl-CoA",
}
BRAIN_MEANINGS = {
    "gtex_brain_max": "largest GTEx v10 median TPM over the 13 brain tissues",
    "gtex_other_max": "largest GTEx v10 median TPM over every other tissue",
    "has_expression": "1 when the gene (for a reaction: a gene in its rule) is in the Human Protein Atlas region table",
}


def meaning(column: str, component_names: dict | None = None) -> str:
    if column in METABOLITE_MEANINGS:
        return METABOLITE_MEANINGS[column]
    if column.startswith("reaction_ec_class_"):
        number = int(column.rsplit("_", 1)[1])
        return f"1 when an annotated EC number of the reaction starts with {number} ({EC_CLASS_NAMES[number]})"
    if column == "reaction_has_ec":
        return "1 when the reaction has an annotated EC number"
    if column.startswith("protein_rrr_"):
        if component_names and column in component_names:
            return f"{component_names[column]['name']} (see Protein components)"
        return (f"component {column.rsplit('_', 1)[1]} of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto "
                "EC, GO function, GO component, UniProt location and Pfam annotations; standardised")
    if column == "protein_has_protein_descriptors":
        return "1 when the node's protein has an ESM-2 embedding; 0 leaves the components at 0"
    for prefix, owner in (("gene_brain_", "gene"), ("reaction_brain_", "reaction")):
        if column.startswith(prefix):
            name = column[len(prefix):]
            via_rule = "" if owner == "gene" else "; from the gene rule (minimum over and, maximum over or)"
            if name in BRAIN_MEANINGS:
                return BRAIN_MEANINGS[name] + ("" if name == "has_expression" else ", log1p, standardised") + via_rule
            if name.startswith("region_"):
                return f"Human Protein Atlas consensus nTPM in the {name[len('region_'):].replace('_', ' ')}, log1p, standardised{via_rule}"
            if name == "class_dopaminergic_neuron":
                return ("dopaminergic cluster 395 of Siletti et al. 2023 (CELLxGENE counts put on the HPA nCPM scale), log1p, "
                        f"standardised{via_rule}")
            if name.startswith("class_"):
                return (f"Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class "
                        f"{name[len('class_'):].replace('_', ' ')}, log1p, standardised{via_rule}")
    return ""


def read_descriptor_table(path: str) -> pd.DataFrame:
    table = pd.read_parquet(path)
    return table.set_index("node_id") if "node_id" in table.columns else table


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def columns_filled_per_node_type(table: pd.DataFrame, nodes: pd.DataFrame) -> dict[str, list[str]]:
    """Per node type, the descriptor columns with a non-zero value on at least one node of that type, in table order."""
    values = table.reindex(nodes.node_id)
    if values.isna().any(axis=None):
        raise ValueError("a graph node has no row in the descriptor table")
    filled = {}
    for node_type in sorted(nodes.node_type.unique()):
        is_type = (nodes.node_type == node_type).to_numpy()
        nonzero = (values.to_numpy()[is_type] != 0).any(axis=0)
        filled[node_type] = [column for column, used in zip(values.columns, nonzero) if used]
    return filled


def structural_feature_names(nodes: pd.DataFrame) -> list[str]:
    """The names of ExperimentData.structural_node_features() for these nodes, in its order (the code builds them by
    position: one-hot node type, multi-hot compartment, then six columns)."""
    node_types = sorted(nodes.node_type.unique())
    compartments = sorted({part for value in nodes.compartment.fillna("") for part in str(value).split(";") if part})
    return ([f"type_{name}" for name in node_types] + [f"compartment_{name}" for name in compartments]
            + ["log1p_degree", "is_currency", "is_transport", "is_reversible", "log1p_gtex_brain_median_tpm_max", "brain_expressed"])


def code_list(columns: list[str], component_names: dict | None = None) -> str:
    return "\n".join(f"- `{column}`" + (f": {meaning(column, component_names)}" if meaning(column, component_names) else "") for column in columns)


def pole_text(entries: list[dict]) -> str:
    return "; ".join(f"{entry['name']} ({entry['block']}, {entry['correlation']:+.2f})" for entry in entries)


def component_section(component_names: dict, names_file: dict) -> list[str]:
    lines = ["## Protein components", "",
             "`protein_rrr_1` to `protein_rrr_64` are the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, "
             "GO function, GO component, UniProt location and Pfam annotations (mechanistic_pathway_learning/graph/protein_descriptors.py, "
             "docs/protein_descriptor_report.md), averaged over a gene's reviewed entries and standardised. Each is a direction in "
             "annotation space and mixes many annotations; component 1 holds the most predicted annotation variance. The name gives the "
             "annotation with the largest positive correlation (high) and the most negative correlation (low) with the component over the "
             f"{names_file['proteins_correlated']:,} proteins of the fit, and each pole lists up to four annotations (one of any group whose "
             "proteins nearly coincide, as with a GO term and its parent). A high value means the protein looks like the positive pole, a "
             "low one like the negative pole. The correlations show how well a name fits: below about 0.3 the name is the strongest "
             "of weak associations. Generated by experiments/name_protein_descriptor_components.py; the column identifiers in the table "
             "are unchanged.", "",
             "| column | name | largest abs. r | positive pole | negative pole |", "|---|---|---|---|---|"]
    for column, entry in component_names.items():
        lines.append(f"| `{column}` | {entry['name']} | {entry['largest_absolute_correlation']:.2f} | {pole_text(entry['positive'])} | {pole_text(entry['negative'])} |")
    return lines + [""]


def placement_section(table_path: str, graph_dir: str, use: str) -> list[str]:
    table = read_descriptor_table(table_path)
    nodes = pd.read_parquet(Path(graph_dir) / "nodes.parquet", columns=["node_id", "node_type"])
    filled = columns_filled_per_node_type(table, nodes)
    blocks = descriptor_blocks(table.columns)
    node_counts = nodes.node_type.value_counts()
    lines = ["## Placement", "", "Cells give how many of a block's columns are non-zero on at least one node of the type (blank: none).", "",
             "| node type | nodes | " + " | ".join(block for block in blocks if blocks[block]) + " |",
             "|---|---|" + "---|" * sum(1 for block in blocks if blocks[block])]
    for node_type, columns in filled.items():
        cells = []
        for block, block_columns in blocks.items():
            if not block_columns:
                continue
            used = [column for column in block_columns if column in columns]
            cells.append("" if not used else ("all" if len(used) == len(block_columns) else f"{len(used)} of {len(block_columns)}"))
        lines.append(f"| {node_type} | {node_counts[node_type]} | " + " | ".join(cells) + " |")
    return lines + [""]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--output", type=Path, default=Path("docs/node_descriptor_columns.md"))
    arguments = parser.parse_args()

    table = read_descriptor_table(FINAL_TABLE)
    nodes = pd.read_parquet(Path(FINAL_GRAPH_DIR) / "nodes.parquet", columns=["node_id", "node_type", "compartment"])
    filled = columns_filled_per_node_type(table, nodes)
    unused = [column for column in table.columns if not any(column in columns for columns in filled.values())]
    names_file = json.loads(COMPONENT_NAMES.read_text())
    component_names = names_file["components"]
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()

    lines = [
        "# Node descriptor columns", "",
        f"Version {DOCUMENT_VERSION}, generated {date.today().isoformat()} by experiments/write_descriptor_column_doc.py from the code at commit {commit}. "
        "Rerun the script after any change to the table; it reads the column names and their placement from the table.", "",
        f"This document covers the final descriptor table only: `{FINAL_TABLE}` (graph `{FINAL_GRAPH_DIR}`, the full graph with the "
        f"neuronal variant and with genes and proteins as separate nodes), used by {FINAL_TABLE_USE}. It holds every block: protein "
        "components on the protein nodes, brain expression by GTEx tissue, Human Protein Atlas region and cell class (with the "
        "dopaminergic neuron class) on the gene and reaction nodes, metabolite properties and reaction EC classes. "
        f"{table.shape[1]} columns; SHA-256 `{file_sha256(Path(FINAL_TABLE))[:16]}`.", "",
        "Every node gets the same descriptor columns. Each node type fills its own block and is 0 in the others, so the "
        "encoder's input layer acts as one linear map per node type and nothing is learned per node "
        "(mechanistic_pathway_learning/graph/node_descriptors.py). A column is listed under a node type below when at least one "
        "node of that type has a non-zero value in it.", "",
        "## Columns per node type", "",
    ]
    for node_type, columns in filled.items():
        lines += [f"### {node_type} ({len(columns)} columns)", "", code_list(columns, component_names) if columns else "No descriptor column; all 0.", ""]
    if unused:
        lines += ["### Columns no node type fills", "", code_list(unused, component_names), ""]
    lines += component_section(component_names, names_file)
    lines += ["## Structural features before the descriptors", "",
              "The encoders read the descriptors after the structural features of ExperimentData.structural_node_features() "
              "(experiments/run_main_model.py, node_feature_matrix). The code builds these by position, without names; the names "
              "below are given here for reading, in order:", "",
              "\n".join(f"{position + 1}. `{name}`" for position, name in enumerate(structural_feature_names(nodes))), "",
              "## Blocks", "",
              "experiments/run_main_model.py --drop-descriptor-blocks leaves out whole blocks, by column prefix "
              "(DESCRIPTOR_BLOCK_PREFIXES in node_descriptors.py):", ""]
    for block, columns in descriptor_blocks(table.columns).items():
        lines.append(f"- `{block}`: {len(columns)} columns")
    lines += [""] + placement_section(FINAL_TABLE, FINAL_GRAPH_DIR, FINAL_TABLE_USE)
    lines += ["## Changes", "", "| version | date | change |", "|---|---|---|"] + [f"| {version} | {day} | {change} |" for version, day, change in CHANGES]
    arguments.output.write_text("\n".join(lines) + "\n")
    print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
