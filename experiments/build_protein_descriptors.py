"""Fit the protein descriptors (mechanistic_pathway_learning/graph/protein_descriptors.py) from the ESM-2 embeddings of
experiments/compute_protein_embeddings.py and write them per gene symbol, with a report.

Steps:
  1. Embeddings of all reviewed human proteins (batches of data/processed/protein_embeddings).
  2. Annotation targets (EC to level 2, GO function and component without phenotype-derived evidence, UniProt
     location, Pfam), standardised within the proteins annotated in each block and weighted per block.
  3. The fit uses the proteins with any GO annotation; the ridge penalty is chosen by five-fold held-out R^2 at full
     rank, then the rank as the smallest whose held-out R^2 reaches 95 percent of the best on the grid.
  4. Descriptors of every protein; a gene with several reviewed entries gets their mean.
  5. ICA rotation of the same subspace (seed 0) for reading, its stability against seeds 1 to 4, and for each
     independent component the annotations it correlates with most.
Outputs: data/processed/node_descriptors/protein_descriptors.parquet (rrr_1 ... rrr_k),
protein_descriptors_ica.parquet (ica_1 ... ica_k), protein_descriptor_fit.json and docs/protein_descriptor_report.md.

Usage:
  OMP_NUM_THREADS=2 python experiments/build_protein_descriptors.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.protein_descriptors import (
    RANK_GRID,
    RELATIVE_RIDGE_GRID,
    annotation_target_matrix,
    choose_rank,
    cross_validated_r_squared,
    fit_reduced_rank_regression,
    go_annotations,
    ica_rotation,
    matched_component_stability,
    read_go_ancestors,
    weighted_standardised_targets,
)

ICA_STABILITY_SEEDS = (1, 2, 3, 4)
TOP_ANNOTATIONS_PER_COMPONENT = 5


def load_embeddings(embedding_directory: Path) -> pd.DataFrame:
    frames = []
    for path in sorted(embedding_directory.glob("batch*.npz")):
        if ".part" in path.name:
            continue
        batch = np.load(path)
        frames.append(pd.DataFrame(batch["embeddings"], index=batch["accessions"]))
    return pd.concat(frames)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--uniprot-table", type=Path, default=Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz"))
    parser.add_argument("--gaf", type=Path, default=Path("data/raw/gene_ontology/goa_human.gaf.gz"))
    parser.add_argument("--go-obo", type=Path, default=Path("data/raw/gene_ontology/go-basic.obo"))
    parser.add_argument("--embedding-dir", type=Path, default=Path("data/processed/protein_embeddings"))
    parser.add_argument("--slice-nodes", type=Path, default=Path("data/processed/graph/nodes.parquet"))
    parser.add_argument("--full-nodes", type=Path, default=Path("data/processed/graph_full/nodes.parquet"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/node_descriptors"))
    parser.add_argument("--report", type=Path, default=Path("docs/protein_descriptor_report.md"))
    arguments = parser.parse_args()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    uniprot = pd.read_csv(arguments.uniprot_table, sep="\t")
    embeddings = load_embeddings(arguments.embedding_dir)
    uniprot = uniprot[uniprot.Entry.isin(embeddings.index)].reset_index(drop=True)
    targets, block_columns, observed = annotation_target_matrix(uniprot, go_annotations(arguments.gaf, read_go_ancestors(arguments.go_obo)))
    target_values = weighted_standardised_targets(targets, block_columns, observed)
    fit_rows = (observed["go_function"] | observed["go_component"]).to_numpy()
    fit_embeddings = embeddings.loc[uniprot.Entry[fit_rows]].to_numpy(dtype=np.float64)
    fit_targets = target_values[fit_rows]

    full_rank = min(fit_embeddings.shape[1], fit_targets.shape[1])
    r_squared_by_alpha = {alpha: cross_validated_r_squared(fit_embeddings, fit_targets, [full_rank], alpha)[full_rank] for alpha in RELATIVE_RIDGE_GRID}
    relative_alpha = max(r_squared_by_alpha, key=r_squared_by_alpha.get)
    r_squared_by_rank = cross_validated_r_squared(fit_embeddings, fit_targets, [rank for rank in RANK_GRID if rank <= full_rank], relative_alpha)
    rank = choose_rank(r_squared_by_rank)

    model = fit_reduced_rank_regression(fit_embeddings, fit_targets, rank, relative_alpha)
    descriptors = pd.DataFrame(model.transform(embeddings.loc[uniprot.Entry].to_numpy(dtype=np.float64)), index=uniprot.Entry,
                               columns=[f"rrr_{index + 1}" for index in range(rank)])
    ica = ica_rotation(descriptors.to_numpy(), seed=0)
    independent = pd.DataFrame(ica.transform(descriptors.to_numpy()), index=uniprot.Entry, columns=[f"ica_{index + 1}" for index in range(rank)])
    stability = {seed: matched_component_stability(independent.to_numpy(), ica_rotation(descriptors.to_numpy(), seed).transform(descriptors.to_numpy()))
                 for seed in ICA_STABILITY_SEEDS}

    gene_of = uniprot.set_index("Entry")["Gene Names (primary)"]
    per_gene = descriptors.groupby(gene_of.reindex(descriptors.index).to_numpy()).mean()
    per_gene_ica = independent.groupby(gene_of.reindex(independent.index).to_numpy()).mean()
    per_gene.index.name = per_gene_ica.index.name = "gene_symbol"
    per_gene.to_parquet(arguments.output_dir / "protein_descriptors.parquet")
    per_gene_ica.to_parquet(arguments.output_dir / "protein_descriptors_ica.parquet")

    annotated = fit_rows  # proteins with a GO annotation; unannotated blocks sit at the mean
    target_frame = pd.DataFrame(target_values, index=uniprot.Entry, columns=targets.columns)
    top_annotations = {}
    for column in independent.columns:
        correlation = target_frame[annotated].corrwith(independent[column][annotated])
        strongest = correlation.abs().sort_values(ascending=False).index[:TOP_ANNOTATIONS_PER_COMPONENT]
        top_annotations[column] = [f"{name} ({correlation[name]:+.2f})" for name in strongest]
    coverage = {}
    for label, nodes_path in (("slice", arguments.slice_nodes), ("full data", arguments.full_nodes)):
        genes = set(pd.read_parquet(nodes_path).query("node_type == 'gene'").gene_symbol.dropna())
        coverage[label] = (len(genes & set(per_gene.index)), len(genes))
    fit = {"proteins_embedded": len(embeddings), "proteins_in_fit": int(fit_rows.sum()), "targets": int(targets.shape[1]),
           "targets_per_block": {block: len(columns) for block, columns in block_columns.items()},
           "held_out_r_squared_by_relative_alpha_full_rank": r_squared_by_alpha, "relative_alpha": relative_alpha,
           "held_out_r_squared_by_rank": r_squared_by_rank, "rank": rank, "ica_stability_by_seed": stability,
           "gene_coverage": {label: {"covered": covered, "genes": total} for label, (covered, total) in coverage.items()}}
    (arguments.output_dir / "protein_descriptor_fit.json").write_text(json.dumps(fit, indent=2) + "\n")

    lines = [
        "# Protein descriptors (generated by experiments/build_protein_descriptors.py)",
        "",
        f"ESM-2 (35M) embeddings of {len(embeddings):,} reviewed human proteins, reduced by reduced-rank regression onto {targets.shape[1]:,} "
        f"function annotations ({', '.join(f'{block} {len(columns)}' for block, columns in block_columns.items())}), fitted on the "
        f"{int(fit_rows.sum()):,} proteins with a GO annotation. No symptom label enters; the directions are the same for every split.",
        "",
        "Held-out R^2 of the weighted targets (five folds over proteins):",
        "",
        "| ridge alpha / proteins | " + " | ".join(f"{alpha:g}" for alpha in r_squared_by_alpha) + " |",
        "|---|" + "---|" * len(r_squared_by_alpha),
        f"| full rank ({full_rank}) | " + " | ".join(f"{value:.3f}" for value in r_squared_by_alpha.values()) + " |",
        "",
        "| rank | " + " | ".join(str(value) for value in r_squared_by_rank) + " |",
        "|---|" + "---|" * len(r_squared_by_rank),
        f"| held-out R^2 (alpha {relative_alpha:g}) | " + " | ".join(f"{value:.3f}" for value in r_squared_by_rank.values()) + " |",
        "",
        ("Note: the chosen ridge penalty is at the edge of the grid. " if relative_alpha in (min(RELATIVE_RIDGE_GRID), max(RELATIVE_RIDGE_GRID)) else "")
        + f"Chosen rank {rank}: the smallest reaching 95 percent of the best held-out R^2. Gene coverage: "
        + "; ".join(f"{label} {covered:,} of {total:,} gene nodes" for label, (covered, total) in coverage.items())
        + " (the rest get zeros and has_protein_descriptors 0).",
        "",
        f"ICA rotation of the same {rank}-dimensional subspace (for reading; the model's input layer sees the same subspace either way). "
        "Stability against other seeds, mean absolute correlation of matched components: "
        + ", ".join(f"seed {seed} {value:.3f}" for seed, value in stability.items()) + ".",
        "",
        "| component | annotations it correlates with most (Pearson r over the proteins of the fit) |",
        "|---|---|",
    ] + [f"| {component} | {'; '.join(annotations)} |" for component, annotations in top_annotations.items()]
    arguments.report.write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:20]))


if __name__ == "__main__":
    main()
