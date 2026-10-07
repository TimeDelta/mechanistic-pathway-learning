"""Append brain region and brain cell-class expression (mechanistic_pathway_learning/graph/brain_expression_descriptors.py)
to a graph's node descriptor table, for --node-descriptors of experiments/run_main_model.py.

Writes a new table outside the graph directory, so the graph that running jobs read is not touched, and a JSON summary
beside it with the coverage of the gene and reaction blocks and the values of a few anchor genes.

Usage:
  python experiments/build_brain_expression_descriptors.py --graph-dir data/processed/graph \
      --output data/processed/node_descriptors/slice_descriptors_brain_expression.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.brain_expression_descriptors import (
    brain_expression_blocks, gene_node_ensembl_ids, gtex_brain_and_other_maximum, hpa_cell_class_table, hpa_region_table,
    read_zipped_table, unambiguous_symbol_to_ensembl)

ANCHOR_GENES = ["PAH", "OTC", "CPS1", "GCH1", "TH", "DDC", "DBH", "MAOA", "COMT", "SLC6A3", "SLC18A2", "GAD1", "ALDH5A1",
                "ATP1A3", "SLC12A5"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--base-descriptors", type=Path, default=None,
                        help="descriptor table to extend (default <graph-dir>/node_descriptors.parquet); --no-base-descriptors gives the brain blocks alone")
    parser.add_argument("--no-base-descriptors", action="store_true")
    parser.add_argument("--gtex", type=Path, default=Path("data/raw/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"))
    parser.add_argument("--hpa-region", type=Path, default=Path("data/raw/hpa/rna_brain_region_hpa/rna_brain_region_hpa.tsv.zip"))
    parser.add_argument("--hpa-cluster-type", type=Path, default=Path("data/raw/hpa/rna_single_nuclei_cluster_type/rna_single_nuclei_cluster_type.tsv.zip"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/node_descriptors/slice_descriptors_brain_expression.parquet"))
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    region_table = read_zipped_table(arguments.hpa_region)
    gene_expression = pd.concat([
        gtex_brain_and_other_maximum(arguments.gtex),
        hpa_region_table(region_table),
        hpa_cell_class_table(read_zipped_table(arguments.hpa_cluster_type)),
    ], axis=1)
    symbol_to_ensembl = unambiguous_symbol_to_ensembl(region_table["Gene name"], region_table["Gene"])
    blocks = brain_expression_blocks(nodes, gene_expression, symbol_to_ensembl)
    gene_ensembl = gene_node_ensembl_ids(nodes, symbol_to_ensembl)[(nodes.node_type == "gene").to_numpy()]
    if arguments.no_base_descriptors:
        table = blocks
    else:
        base = pd.read_parquet(arguments.base_descriptors or arguments.graph_dir / "node_descriptors.parquet")
        table = pd.concat([base.loc[blocks.index], blocks], axis=1)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(arguments.output)

    is_gene = (nodes.node_type == "gene").to_numpy()
    is_reaction = (nodes.node_type == "reaction").to_numpy()
    anchors = {}
    node_symbol_to_ensembl = dict(zip(nodes.gene_symbol[is_gene], gene_node_ensembl_ids(nodes, symbol_to_ensembl)[is_gene]))
    for symbol in ANCHOR_GENES:
        ensembl = node_symbol_to_ensembl.get(symbol)
        if ensembl is None or ensembl not in gene_expression.index:
            anchors[symbol] = None
            continue
        raw = gene_expression.loc[ensembl]
        anchors[symbol] = {column: round(float(raw[column]), 1) for column in ["gtex_brain_max", "gtex_other_max"]}
        region_columns = [column for column in raw.index if column.startswith("region_")]
        class_columns = [column for column in raw.index if column.startswith("class_")]
        anchors[symbol]["top_region"] = str(raw[region_columns].astype(float).idxmax()) if raw[region_columns].notna().any() else None
        anchors[symbol]["top_cell_class"] = str(raw[class_columns].astype(float).idxmax()) if raw[class_columns].notna().any() else None
    summary = {
        "output": str(arguments.output), "graph_dir": str(arguments.graph_dir), "nodes": len(table), "columns": len(table.columns),
        "brain_columns": [column for column in blocks.columns],
        "gene_nodes": int(is_gene.sum()), "gene_nodes_without_own_ensembl_id": int(nodes.ensembl_gene_id[is_gene].isna().sum()),
        "gene_nodes_resolved_by_symbol": int((nodes.ensembl_gene_id[is_gene].isna() & gene_ensembl.notna().to_numpy()).sum()),
        "gene_nodes_with_expression": int(blocks.loc[is_gene, "gene_brain_has_expression"].sum()),
        "reaction_nodes": int(is_reaction.sum()),
        "reaction_nodes_with_rule": int((nodes.gene_reaction_rule[is_reaction].fillna("") != "").sum()),
        "reaction_nodes_with_expression": int(blocks.loc[is_reaction, "reaction_brain_has_expression"].sum()),
        "anchor_genes": anchors,
        "finite": bool(np.isfinite(table.to_numpy(dtype=float)).all()),
    }
    arguments.output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "brain_columns"}, indent=2))


if __name__ == "__main__":
    main()
