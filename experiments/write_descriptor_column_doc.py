"""Write docs/node_descriptor_columns.md: the name of every node descriptor column, which node type fills it, and where
each descriptor table puts it.

The column lists are read from the descriptor tables themselves, and a column is listed for a node type only when at
least one node of that type has a non-zero value in it, so the document cannot drift from the tables. Each table's
SHA-256 is recorded; when a table changes, rerun this script, raise DOCUMENT_VERSION and add a line to CHANGES.

    python experiments/write_descriptor_column_doc.py
"""
from __future__ import annotations

import argparse
import hashlib
import subprocess
from datetime import date
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.node_descriptors import descriptor_blocks

DOCUMENT_VERSION = 1
CHANGES = [
    (1, "2026-10-08", "First version: the 137 descriptor columns of the brain-expression tables and the 83 of the graph directories' own tables."),
]
# (descriptor table, graph directory whose nodes it describes, where it is used)
TABLES = [
    ("data/processed/node_descriptors/full_neuronal_descriptors_brain_expression.parquet", "data/processed/graph_full_neuronal",
     "the confirmatory configurations (FULL_GRAPH_NODE_PROPERTIES in experiments/run_main_model_batch.py)"),
    ("data/processed/node_descriptors/full_neuronal_descriptors_brain_expression_resolved.parquet", "data/processed/graph_full_neuronal",
     "merged full graph with the protein block resolved per entry (twin of the full split)"),
    ("data/processed/node_descriptors/full_neuronal_split_descriptors_brain_expression.parquet", "data/processed/graph_full_neuronal_split",
     "full graph with genes and proteins split"),
    ("data/processed/node_descriptors/slice_descriptors_brain_expression.parquet", "data/processed/graph", "slice"),
    ("data/processed/node_descriptors/slice_descriptors_brain_expression_resolved.parquet", "data/processed/graph",
     "slice with the protein block resolved per entry (twin of the split slice)"),
    ("data/processed/node_descriptors/slice_split_descriptors_brain_expression.parquet", "data/processed/graph_split", "slice with genes and proteins split"),
    ("data/processed/graph/node_descriptors.parquet", "data/processed/graph", "slice configurations without brain expression"),
    ("data/processed/graph_full_neuronal/node_descriptors.parquet", "data/processed/graph_full_neuronal", "full graph without brain expression"),
    ("data/processed/graph_split/node_descriptors.parquet", "data/processed/graph_split", "split slice without brain expression"),
    ("data/processed/graph_full_neuronal_split/node_descriptors.parquet", "data/processed/graph_full_neuronal_split", "full split without brain expression"),
]
REFERENCE_TABLE = TABLES[0][0]
SPLIT_REFERENCE_TABLE = TABLES[2][0]
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


def meaning(column: str) -> str:
    if column in METABOLITE_MEANINGS:
        return METABOLITE_MEANINGS[column]
    if column.startswith("reaction_ec_class_"):
        number = int(column.rsplit("_", 1)[1])
        return f"1 when an annotated EC number of the reaction starts with {number} ({EC_CLASS_NAMES[number]})"
    if column == "reaction_has_ec":
        return "1 when the reaction has an annotated EC number"
    if column.startswith("protein_rrr_"):
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


def code_list(columns: list[str]) -> str:
    return "\n".join(f"- `{column}`" + (f": {meaning(column)}" if meaning(column) else "") for column in columns)


def table_section(table_path: str, graph_dir: str, use: str) -> list[str]:
    table = read_descriptor_table(table_path)
    nodes = pd.read_parquet(Path(graph_dir) / "nodes.parquet", columns=["node_id", "node_type"])
    filled = columns_filled_per_node_type(table, nodes)
    blocks = descriptor_blocks(table.columns)
    node_counts = nodes.node_type.value_counts()
    lines = [f"### `{table_path}`", "", f"Used by: {use}. Graph: `{graph_dir}`. {table.shape[1]} columns. SHA-256 `{file_sha256(Path(table_path))[:16]}`.", "",
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

    reference_table = read_descriptor_table(REFERENCE_TABLE)
    reference_nodes = pd.read_parquet(Path(TABLES[0][1]) / "nodes.parquet", columns=["node_id", "node_type", "compartment"])
    split_table = read_descriptor_table(SPLIT_REFERENCE_TABLE)
    split_nodes = pd.read_parquet(Path(TABLES[2][1]) / "nodes.parquet", columns=["node_id", "node_type"])
    filled = columns_filled_per_node_type(reference_table, reference_nodes)
    filled_split = columns_filled_per_node_type(split_table, split_nodes)
    unused = [column for column in reference_table.columns if not any(column in columns for columns in filled.values())]
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()

    lines = [
        "# Node descriptor columns", "",
        f"Version {DOCUMENT_VERSION}, generated {date.today().isoformat()} by experiments/write_descriptor_column_doc.py from the code at commit {commit}. "
        "Rerun the script after any descriptor table changes; it reads the column names and their placement from the tables.", "",
        "Every node gets the same descriptor columns. Each node type fills its own block and is 0 in the others, so the "
        "encoder's input layer acts as one linear map per node type and nothing is learned per node "
        "(mechanistic_pathway_learning/graph/node_descriptors.py). A column is listed under a node type below when at least one "
        "node of that type has a non-zero value in it.", "",
        f"## Columns per node type in the confirmatory table (`{REFERENCE_TABLE}`)", "",
    ]
    for node_type, columns in filled.items():
        lines += [f"### {node_type} ({len(columns)} columns)", "", code_list(columns) if columns else "No descriptor column; all 0.", ""]
    if unused:
        lines += ["### Columns no node type fills in this table", "", code_list(unused), ""]
    moved = {node_type: columns for node_type, columns in filled_split.items() if columns != filled.get(node_type)}
    lines += [f"## Where the split graphs differ (`{SPLIT_REFERENCE_TABLE}`)", "",
              "With genes and proteins split (docs/gene_protein_split.md), the protein block sits on the protein nodes and the "
              "gene nodes keep the brain-expression block. Node types whose columns differ from the table above:", ""]
    for node_type, columns in moved.items():
        filled_blocks = [f"`{block}`" for block, block_columns in descriptor_blocks(columns).items() if block_columns]
        lines.append(f"- {node_type}: {len(columns)} columns" + (f", block {' and '.join(filled_blocks)}" if columns else ""))
    lines += ["", "## Structural features before the descriptors", "",
              "The encoders read the descriptors after the structural features of ExperimentData.structural_node_features() "
              "(experiments/run_main_model.py, node_feature_matrix). The code builds these by position, without names; the names "
              f"below are given here for reading, in order, for `{TABLES[0][1]}`:", "",
              "\n".join(f"{position + 1}. `{name}`" for position, name in enumerate(structural_feature_names(reference_nodes))), "",
              "## Blocks", "",
              "experiments/run_main_model.py --drop-descriptor-blocks leaves out whole blocks, by column prefix "
              "(DESCRIPTOR_BLOCK_PREFIXES in node_descriptors.py):", ""]
    for block, columns in descriptor_blocks(reference_table.columns).items():
        lines.append(f"- `{block}`: {len(columns)} columns")
    lines += ["", "## Placement in every descriptor table", "",
              "Cells give how many of a block's columns are non-zero on at least one node of the type (blank: none).", ""]
    for table_path, graph_dir, use in TABLES:
        if Path(table_path).exists() and (Path(graph_dir) / "nodes.parquet").exists():
            lines += table_section(table_path, graph_dir, use)
        else:
            lines += [f"### `{table_path}`", "", "Not present when this version was generated.", ""]
    lines += ["## Changes", "", "| version | date | change |", "|---|---|---|"] + [f"| {version} | {day} | {change} |" for version, day, change in CHANGES]
    arguments.output.write_text("\n".join(lines) + "\n")
    print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
