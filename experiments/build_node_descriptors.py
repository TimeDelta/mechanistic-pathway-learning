"""Assemble the fixed node descriptors of one graph (mechanistic_pathway_learning/graph/node_descriptors.py) into
<graph-dir>/node_descriptors.parquet, indexed by node_id, for --node-descriptors of experiments/run_main_model.py.

Blocks: metabolite physicochemical properties (RDKit on the Human-GEM SMILES, SBML charge), reaction EC classes (SBML),
the protein descriptors of experiments/build_protein_descriptors.py (reduced-rank regression of ESM-2 embeddings
onto function annotations; --protein-descriptors, skipped with --no-proteins), and the annotation vectors of the
Reactome protein entities of experiments/build_complex_descriptors.py (--complex-descriptors, off by default; they
carry the protein descriptors' columns, so they share that block). A summary goes to
<graph-dir>/node_descriptors_summary.json.

Usage:
  python experiments/build_node_descriptors.py --graph-dir data/processed/graph
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.node_descriptors import assemble_node_descriptor_table, metabolite_descriptors, reaction_enzyme_classes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--human-gem-metabolites", type=Path, default=Path("data/raw/Human-GEM/model/metabolites.tsv"))
    parser.add_argument("--human-gem-sbml", type=Path, default=Path("data/raw/Human-GEM/model/Human-GEM.xml"))
    parser.add_argument("--protein-descriptors", type=Path, default=Path("data/processed/node_descriptors/protein_descriptors.parquet"))
    parser.add_argument("--no-proteins", action="store_true")
    parser.add_argument("--complex-descriptors", type=Path, default=None,
                        help="annotation vectors of the Reactome protein entities (data/processed/node_descriptors/complex_descriptors.parquet)")
    arguments = parser.parse_args()
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    sbml_text = arguments.human_gem_sbml.read_text()
    metabolite_table = metabolite_descriptors(pd.read_csv(arguments.human_gem_metabolites, sep="\t"), sbml_text)
    reaction_table = reaction_enzyme_classes(sbml_text)
    protein_table = None if arguments.no_proteins else pd.read_parquet(arguments.protein_descriptors)
    complex_table = None if arguments.complex_descriptors is None else pd.read_parquet(arguments.complex_descriptors)
    table = assemble_node_descriptor_table(nodes, metabolite_table, reaction_table, protein_table, complex_table)
    entity_flag = next((column for column in ("protein_has_complex_descriptors", "complex_has_complex_descriptors") if column in table.columns), None)
    table.to_parquet(arguments.graph_dir / "node_descriptors.parquet")
    by_type = nodes.node_type.to_numpy()
    summary = {
        "nodes": len(table), "columns": list(table.columns),
        "metabolite_nodes_with_structure": int(table.loc[by_type == "metabolite", "metabolite_has_structure"].sum()),
        "metabolite_nodes_with_partial_structure": int(table.loc[by_type == "metabolite", "metabolite_partial_structure"].sum()),
        "metabolite_nodes": int((by_type == "metabolite").sum()),
        "reaction_nodes_with_ec": int(table.loc[by_type == "reaction", "reaction_has_ec"].sum()), "reaction_nodes": int((by_type == "reaction").sum()),
        "gene_nodes_with_protein_descriptors": None if protein_table is None else int(table.loc[by_type == "gene", "protein_has_protein_descriptors"].sum()),
        "gene_nodes": int((by_type == "gene").sum()),
        "protein_descriptors": None if protein_table is None else str(arguments.protein_descriptors),
        "complex_descriptors": None if complex_table is None else str(arguments.complex_descriptors),
        "protein_entity_nodes_with_descriptors": None if entity_flag is None else int(table.loc[by_type == "protein_entity", entity_flag].sum()),
        "protein_entity_nodes": int((by_type == "protein_entity").sum()),
    }
    (arguments.graph_dir / "node_descriptors_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "columns"}, indent=2), f"\n{len(table.columns)} descriptor columns")


if __name__ == "__main__":
    main()
