"""Write a split copy of a graph: each gene node gets a protein node, joined by an encodes edge
(mechanistic_pathway_learning/graph/gene_protein_split.py; docs/gene_protein_split.md).

Reads a graph directory and writes a new one. It refuses to write into a directory that already holds a graph unless
--overwrite is given, so a graph a running job reads is never regenerated. It writes:
- nodes.parquet, edges.parquet and relation_types.json, with encodes appended last and degree recomputed;
- gene_to_protein.parquet, which tells experiment_data.load_experiment_data to seed a drug's target proteins and a
  knockout's gene;
- node_descriptors.parquet, if the source graph has one, with the protein and gene brain blocks moved to the protein
  nodes;
- split_summary.json.
Each --descriptor-table and --cell-class-weights pair (source:destination) is moved the same way into a new file.

Usage:
  python experiments/build_gene_protein_split.py --graph-dir data/processed/graph --output-dir data/processed/graph_split \\
      --descriptor-table data/processed/node_descriptors/slice_descriptors_brain_expression.parquet:data/processed/node_descriptors/slice_split_descriptors_brain_expression.parquet \\
      --cell-class-weights data/processed/cell_class_weights/slice_cell_class_weights.parquet:data/processed/cell_class_weights/slice_split_cell_class_weights.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import GENE_TO_PROTEIN_FILE
from mechanistic_pathway_learning.graph.gene_protein_split import (
    PROTEIN_NODE_TYPE, cell_class_weights_for_split, move_gene_blocks_to_proteins, protein_node_assignment, split_graph)


def source_and_destination(pair: str) -> tuple[Path, Path]:
    source, separator, destination = pair.partition(":")
    if not separator or not destination:
        raise argparse.ArgumentTypeError(f"expected source:destination, got {pair!r}")
    return Path(source), Path(destination)


def refuse_existing(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise SystemExit(f"{path} exists; a split is written to a new place (pass --overwrite only when no job reads it)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--uniprot-table", type=Path, default=Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz"))
    parser.add_argument("--protein-descriptors", type=Path, default=Path("data/processed/node_descriptors/protein_descriptors.parquet"),
                        help="per-entry protein descriptors (experiments/build_protein_descriptors.py), for the protein nodes shared by several genes")
    parser.add_argument("--descriptor-table", type=source_and_destination, action="append", default=[])
    parser.add_argument("--cell-class-weights", type=source_and_destination, action="append", default=[])
    parser.add_argument("--overwrite", action="store_true")
    arguments = parser.parse_args()

    refuse_existing(arguments.output_dir / "nodes.parquet", arguments.overwrite)
    for _, destination in arguments.descriptor_table + arguments.cell_class_weights:
        refuse_existing(destination, arguments.overwrite)

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet")
    relation_types = json.loads((arguments.graph_dir / "relation_types.json").read_text())
    primary_gene_names = pd.read_csv(arguments.uniprot_table, sep="\t", usecols=["Gene Names (primary)"])["Gene Names (primary)"]
    assignment = protein_node_assignment(nodes, primary_gene_names)
    split_nodes, split_edges, split_relations, summary = split_graph(nodes, edges, relation_types, assignment)

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    split_nodes.to_parquet(arguments.output_dir / "nodes.parquet", index=False)
    split_edges.to_parquet(arguments.output_dir / "edges.parquet", index=False)
    (arguments.output_dir / "relation_types.json").write_text(json.dumps(split_relations, indent=1) + "\n")
    assignment.to_parquet(arguments.output_dir / GENE_TO_PROTEIN_FILE, index=False)

    protein_descriptors = pd.read_parquet(arguments.protein_descriptors) if arguments.protein_descriptors.exists() else None
    tables = list(arguments.descriptor_table)
    if (arguments.graph_dir / "node_descriptors.parquet").exists():
        tables.insert(0, (arguments.graph_dir / "node_descriptors.parquet", arguments.output_dir / "node_descriptors.parquet"))
    summary["descriptor_tables"] = {}
    for source, destination in tables:
        table = pd.read_parquet(source)
        moved = move_gene_blocks_to_proteins(table, split_nodes, assignment, protein_descriptors, primary_gene_names)
        destination.parent.mkdir(parents=True, exist_ok=True)
        moved.to_parquet(destination)
        is_protein = (split_nodes.set_index("node_id").node_type.reindex(moved.index) == PROTEIN_NODE_TYPE).to_numpy()
        gene_rows = (split_nodes.set_index("node_id").node_type.reindex(moved.index) == "gene").to_numpy()
        flags = [column for column in ("protein_has_protein_descriptors", "gene_brain_has_expression") if column in moved.columns]
        summary["descriptor_tables"][str(destination)] = {
            "source": str(source), "rows": len(moved), "columns": len(moved.columns),
            "nonzero_gene_rows": int((moved[gene_rows] != 0).any(axis=1).sum()),
            **{f"protein_nodes_with_{flag}": int(moved.loc[is_protein, flag].sum()) for flag in flags},
            **{f"gene_nodes_with_{flag}_before": int(table.loc[table.index.intersection(assignment.gene_node_id), flag].sum()) for flag in flags}}
    summary["cell_class_weights"] = {}
    for source, destination in arguments.cell_class_weights:
        weights = cell_class_weights_for_split(pd.read_parquet(source), split_nodes, assignment)
        destination.parent.mkdir(parents=True, exist_ok=True)
        weights.to_parquet(destination)
        summary["cell_class_weights"][str(destination)] = {"source": str(source), "rows": len(weights), "rows_with_a_missing_value": int(weights.isna().any(axis=1).sum())}

    is_gene = split_nodes.node_type == "gene"
    is_protein = split_nodes.node_type == PROTEIN_NODE_TYPE
    summary.update({
        "graph_dir": str(arguments.graph_dir), "output_dir": str(arguments.output_dir),
        "nodes_before": len(nodes), "nodes_after": len(split_nodes),
        "gene_degree_quantiles": {str(q): float(np.quantile(split_nodes.degree[is_gene], q)) for q in (0.25, 0.5, 0.75, 0.95)},
        "protein_degree_quantiles": {str(q): float(np.quantile(split_nodes.degree[is_protein], q)) for q in (0.25, 0.5, 0.75, 0.95)},
        "gene_nodes_of_degree_one": int((split_nodes.degree[is_gene] == 1).sum()),
        "shared_protein_nodes": sorted(assignment[assignment.shared].protein_display_name.unique().tolist()),
    })
    (arguments.output_dir / "split_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key not in ("shared_protein_nodes", "edges_per_relation")}, indent=1, default=str))


if __name__ == "__main__":
    main()
