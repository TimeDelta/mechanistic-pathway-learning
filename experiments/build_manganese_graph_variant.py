"""Build the manganese graph variant (mechanistic_pathway_learning/graph/manganese_extension.py) of a graph directory:
<graph-dir>_manganese with nodes, edges and relation types extended and a summary of what was connected.

Input: data/raw/uniprot/uniprot_human_manganese_cofactor.tsv.gz, the reviewed human proteins whose UniProt cofactor
annotation names Mn(2+) (pinned beside it). Node descriptors of the variant, when needed, are built with
experiments/build_node_descriptors.py --graph-dir <graph-dir>_manganese.

Usage:
  python experiments/build_manganese_graph_variant.py --graph-dirs data/processed/graph data/processed/graph_full
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.manganese_extension import add_manganese, manganese_enzymes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dirs", type=Path, nargs="+", default=[Path("data/processed/graph"), Path("data/processed/graph_full")])
    parser.add_argument("--cofactor-table", type=Path, default=Path("data/raw/uniprot/uniprot_human_manganese_cofactor.tsv.gz"))
    arguments = parser.parse_args()
    required, alternative = manganese_enzymes(pd.read_csv(arguments.cofactor_table, sep="\t"))
    for graph_directory in arguments.graph_dirs:
        output_directory = graph_directory.parent / f"{graph_directory.name}_manganese"
        output_directory.mkdir(parents=True, exist_ok=True)
        nodes, edges = pd.read_parquet(graph_directory / "nodes.parquet"), pd.read_parquet(graph_directory / "edges.parquet")
        relation_types = json.loads((graph_directory / "relation_types.json").read_text())
        new_nodes, new_edges, new_relations, summary = add_manganese(nodes, edges, relation_types, required)
        new_nodes.to_parquet(output_directory / "nodes.parquet")
        new_edges.to_parquet(output_directory / "edges.parquet")
        (output_directory / "relation_types.json").write_text(json.dumps(new_relations, indent=1) + "\n")
        summary.update(source_graph=str(graph_directory), proteins_requiring_manganese=len(required), proteins_accepting_manganese_among_other_metals=len(alternative),
                       nodes_added=len(new_nodes) - len(nodes), edges_added=len(new_edges) - len(edges))
        (output_directory / "manganese_extension_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
        print(output_directory, json.dumps({key: value for key, value in summary.items() if key != "cofactor_enzymes_absent_from_graph"}, indent=1))


if __name__ == "__main__":
    main()
