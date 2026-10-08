"""Write the Ensembl-indexed expression tables the brain descriptors and the cell-class weights were built from, so that a
rewired run can recompute its reactions' expression from rewired gene rules (mechanistic_pathway_learning/graph/
rewired_expression_features.py; run_main_model.py --recompute-reaction-expression-on-rewired-graph).

The tables are assembled exactly as experiments/build_brain_expression_descriptors.py and experiments/build_cell_class_weights.py
assemble them (GTEx brain and other maxima, HPA brain regions, HPA single-nucleus classes and the dopaminergic class).
Before writing, the script recomputes the reaction rows of the descriptor table and of the weights from the graph's real
rules and refuses to write unless both equal the stored files, so the recomputation in a rewired run starts from the
inputs the real run read. Writes new files only; no file a running job reads is touched.

Usage:
  python experiments/write_expression_tables.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.brain_expression_descriptors import (
    add_dopaminergic_class, gtex_brain_and_other_maximum, hpa_cell_class_table, hpa_region_table, read_zipped_table)
from mechanistic_pathway_learning.graph.brain_expression_weights import ALL_CELLS_CLASS
from mechanistic_pathway_learning.graph.rewired_expression_features import rewired_reaction_brain_block, rewired_reaction_cell_class_weights


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--gtex", type=Path, default=Path("data/raw/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"))
    parser.add_argument("--hpa-region", type=Path, default=Path("data/raw/hpa/rna_brain_region_hpa/rna_brain_region_hpa.tsv.zip"))
    parser.add_argument("--hpa-cluster-type", type=Path, default=Path("data/raw/hpa/rna_single_nuclei_cluster_type/rna_single_nuclei_cluster_type.tsv.zip"))
    parser.add_argument("--dopaminergic-expression", type=Path, default=Path("data/processed/brain_expression/dopaminergic_siletti_cluster395.parquet"))
    parser.add_argument("--descriptors", type=Path, default=Path("data/processed/node_descriptors/full_neuronal_descriptors_brain_expression.parquet"))
    parser.add_argument("--cell-class-weights", type=Path, default=Path("data/processed/cell_class_weights/full_neuronal_cell_class_weights.parquet"))
    parser.add_argument("--gene-expression-output", type=Path, default=Path("data/processed/brain_expression/gene_expression_for_descriptors.parquet"))
    parser.add_argument("--class-expression-output", type=Path, default=Path("data/processed/brain_expression/class_expression_for_weights.parquet"))
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    cell_class_table = add_dopaminergic_class(hpa_cell_class_table(read_zipped_table(arguments.hpa_cluster_type)),
                                              pd.read_parquet(arguments.dopaminergic_expression)["ncpm_aligned"])
    gene_expression = pd.concat([gtex_brain_and_other_maximum(arguments.gtex), hpa_region_table(read_zipped_table(arguments.hpa_region)), cell_class_table], axis=1)
    class_expression = cell_class_table.copy()
    class_expression.columns = [column.removeprefix("class_") for column in class_expression.columns]

    rules = pd.Series(nodes.gene_reaction_rule.fillna("").to_numpy(), index=pd.Index(nodes.node_id, name="node_id"))
    is_reaction = (nodes.node_type == "reaction").to_numpy()
    stored_descriptors = pd.read_parquet(arguments.descriptors)
    recomputed_block = rewired_reaction_brain_block(nodes, rules, gene_expression).loc[nodes.node_id[is_reaction]]
    stored_block = stored_descriptors.loc[nodes.node_id[is_reaction], recomputed_block.columns]
    descriptor_difference = float(np.abs(recomputed_block.to_numpy(dtype=float) - stored_block.to_numpy(dtype=float)).max())
    stored_weights = pd.read_parquet(arguments.cell_class_weights)
    recomputed_weights = rewired_reaction_cell_class_weights(nodes, rules, class_expression, include_all_cells_class=ALL_CELLS_CLASS in stored_weights.columns)
    weight_difference = float(np.abs(recomputed_weights[stored_weights.columns].to_numpy(dtype=float)
                                     - stored_weights.loc[recomputed_weights.index].to_numpy(dtype=float)).max())
    summary = {"graph_dir": str(arguments.graph_dir), "reactions": int(is_reaction.sum()),
               "largest_difference_reaction_brain_block": descriptor_difference, "largest_difference_reaction_cell_class_weights": weight_difference,
               "gene_expression_rows": len(gene_expression), "class_expression_rows": len(class_expression)}
    print(json.dumps(summary, indent=1))
    if descriptor_difference > 1e-5 or weight_difference > 1e-6:
        raise SystemExit("the tables do not reproduce the stored reaction rows; not written")
    for table, path in ((gene_expression, arguments.gene_expression_output), (class_expression, arguments.class_expression_output)):
        if path.exists():
            raise SystemExit(f"{path} exists; written once (rewired runs record its hash)")
        path.parent.mkdir(parents=True, exist_ok=True)
        table.to_parquet(path)
    arguments.gene_expression_output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=1) + "\n")


if __name__ == "__main__":
    main()
