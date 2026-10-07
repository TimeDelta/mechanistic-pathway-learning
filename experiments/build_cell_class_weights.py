"""Per-node cell-class weights for --cell-class-weights of experiments/run_main_model.py (the cell-class propagation
channels: mechanistic_pathway_learning/graph/brain_expression_weights.py).

Classes: the 10 Human Protein Atlas single-nucleus classes (brain_expression_descriptors.hpa_cell_class_table), the
dopaminergic class (experiments/build_dopaminergic_expression.py) and an all_cells class of weight 1. Writes a parquet
indexed by node_id with one column per class, outside the graph directory, and a JSON summary with the weight
distribution per class and per node type and the weights of a few anchor genes.

Usage:
  python experiments/build_cell_class_weights.py --graph-dir data/processed/graph_full_neuronal \
      --output data/processed/cell_class_weights/full_neuronal_cell_class_weights.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.brain_expression_descriptors import (
    add_dopaminergic_class, gene_node_ensembl_ids, hpa_cell_class_table, read_zipped_table, unambiguous_symbol_to_ensembl)
from mechanistic_pathway_learning.graph.brain_expression_weights import (
    DETECTION_THRESHOLD_NCPM, extracellular_pool_nodes, node_cell_class_weights)

ANCHOR_GENES = ["TH", "SLC6A3", "DDC", "GAD1", "SLC17A7", "AQP4", "MBP", "PAH", "OTC", "GCH1", "MAOA", "DRD2"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--hpa-region", type=Path, default=Path("data/raw/hpa/rna_brain_region_hpa/rna_brain_region_hpa.tsv.zip"),
                        help="only its symbol-to-Ensembl table is used, to place gene nodes named by symbol")
    parser.add_argument("--hpa-cluster-type", type=Path, default=Path("data/raw/hpa/rna_single_nuclei_cluster_type/rna_single_nuclei_cluster_type.tsv.zip"))
    parser.add_argument("--dopaminergic-expression", type=Path, default=Path("data/processed/brain_expression/dopaminergic_siletti_cluster395.parquet"))
    parser.add_argument("--half-saturation", type=float, default=DETECTION_THRESHOLD_NCPM)
    parser.add_argument("--no-all-cells-class", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("data/processed/cell_class_weights/slice_cell_class_weights.parquet"))
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    region_table = read_zipped_table(arguments.hpa_region)
    symbol_to_ensembl = unambiguous_symbol_to_ensembl(region_table["Gene name"], region_table["Gene"])
    class_expression = add_dopaminergic_class(hpa_cell_class_table(read_zipped_table(arguments.hpa_cluster_type)),
                                              pd.read_parquet(arguments.dopaminergic_expression)["ncpm_aligned"])
    class_expression.columns = [column.removeprefix("class_") for column in class_expression.columns]
    weights = node_cell_class_weights(nodes, class_expression, symbol_to_ensembl, arguments.half_saturation,
                                      include_all_cells_class=not arguments.no_all_cells_class)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    weights.to_parquet(arguments.output)

    node_type = nodes.set_index("node_id").node_type.reindex(weights.index)
    is_gene = (node_type == "gene").to_numpy()
    gene_symbols = nodes.set_index("node_id").gene_symbol.reindex(weights.index)
    gene_ensembl = gene_node_ensembl_ids(nodes, symbol_to_ensembl)
    anchors = {}
    for symbol in ANCHOR_GENES:
        rows = weights[is_gene & (gene_symbols == symbol).to_numpy()]
        anchors[symbol] = None if rows.empty else {column: round(float(value), 3) for column, value in rows.iloc[0].items()}
    summary = {
        "output": str(arguments.output), "graph_dir": str(arguments.graph_dir), "classes": list(weights.columns),
        "half_saturation_ncpm": arguments.half_saturation,
        "gene_nodes": int(is_gene.sum()), "gene_nodes_with_expression": int(gene_ensembl[is_gene].isin(class_expression.index).sum()),
        "extracellular_pool_nodes": int(extracellular_pool_nodes(nodes).sum()),
        "mean_weight_by_node_type": {str(name): {column: round(float(value), 3) for column, value in group.mean().items()}
                                     for name, group in weights.groupby(node_type.to_numpy())},
        "fraction_below_0.1_by_node_type": {str(name): {column: round(float(value), 3) for column, value in (group < 0.1).mean().items()}
                                            for name, group in weights.groupby(node_type.to_numpy())},
        "anchor_genes": anchors,
    }
    arguments.output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: summary[key] for key in ("output", "classes", "gene_nodes", "gene_nodes_with_expression", "extracellular_pool_nodes")}, indent=2))
    print(json.dumps(summary["fraction_below_0.1_by_node_type"], indent=1))
    print(json.dumps({symbol: summary["anchor_genes"][symbol] for symbol in ("TH", "AQP4", "PAH")}, indent=1))


if __name__ == "__main__":
    main()
