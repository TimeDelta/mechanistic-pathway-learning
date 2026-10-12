"""Options for representing tissue type in the graph, with the coverage each one reaches.

The owner, 12 October 2026: "The tissue type of a particular mechanism matters and if that isn't represented in the
graph, I can't calculate when the tissue type is fully saturated regardless of measuring per-class spread per
propagation so I would be relegated to only the steady state measurement. Can you please look into options for tissue
type representation?"

What a saturation reading needs is a map from every node to the tissues it belongs to, so that "how much of tissue T
has moved" is a quantity. The graph carries nothing of the kind today: one scalar per gene, brain_median_tpm_max, the
maximum median expression over GTEx's brain tissues, and a boolean derived from it. This script measures the options
that the data already on disk would support, and the one number that decides whether any of them answers the
question: how many tissues a gene is expressed in. A gene expressed everywhere cannot tell one tissue from another,
so if the median gene is expressed in most tissues then a tissue partition of the nodes does not exist and only a
weighted reading does.

Sources measured, all already downloaded:
  GTEx v10 gene median TPM           68 tissues (14 of them brain regions), by Ensembl gene id
  HPA rna_brain_region_hpa           13 brain regions, by Ensembl gene id and gene name
  HPA rna_single_nuclei_cluster_type 34 brain cell types, by Ensembl gene id and gene name
  the graph's own columns            brain_median_tpm_max, brain_expressed, compartment, and the 11 brain cell
                                     class weights of the node descriptor table where that is built

Coverage is measured over the node kinds separately, because they differ completely: a gene node can be matched to a
tissue table directly, a reaction node only through its gene reaction rule, and a metabolite node only through the
reactions it takes part in.

Counting only: no label, table, graph or model input is changed.

Usage:
  python experiments/label_review/tissue_representation_options.py
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

GTEX_MEDIAN_TPM = Path("data/raw/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")
HPA_BRAIN_REGION = Path("data/raw/hpa/rna_brain_region_hpa/rna_brain_region_hpa.tsv.zip")
HPA_SINGLE_NUCLEI = Path("data/raw/hpa/rna_single_nuclei_cluster_type/rna_single_nuclei_cluster_type.tsv.zip")
EXPRESSED_TPM = 1.0  # the usual cut-off for "expressed" in a bulk median TPM table


def mechanism_target_genes(evidence_directory: Path, nodes: pd.DataFrame) -> set[str]:
    """Every gene symbol that a drug perturbation of the evidence table seeds, that is, every mechanism target."""
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet",
                              columns=["perturbation_type", "perturbation_nodes"])
    symbol_of_node = dict(zip(nodes.node_id, nodes.gene_symbol))
    genes = set()
    for seeded_nodes in records.perturbation_nodes[records.perturbation_type == "drug"]:
        for node_id, _sign, _magnitude in json.loads(seeded_nodes):
            symbol = symbol_of_node.get(node_id)
            if isinstance(symbol, str) and symbol:
                genes.add(symbol)
    return genes


def gtex_median_tpm(path: Path) -> pd.DataFrame:
    """Gene median TPM per tissue, indexed by the unversioned Ensembl gene id, with the gene symbol kept."""
    table = pd.read_csv(path, sep="\t", skiprows=2, low_memory=False)
    table["ensembl_gene_id"] = table.Name.str.split(".").str[0]
    table = table.rename(columns={"Description": "gene_symbol"})
    return table.set_index("ensembl_gene_id")


def hpa_table(path: Path, group_column: str, value_column: str) -> pd.DataFrame:
    """A long HPA table read from its zip, as gene by group."""
    with zipfile.ZipFile(path) as archive:
        member = archive.namelist()[0]
        with archive.open(member) as handle:
            reader = csv.DictReader(io.TextIOWrapper(handle, encoding="utf-8"), delimiter="\t")
            rows = [(row["Gene"], row["Gene name"], row[group_column], row[value_column]) for row in reader]
    long_table = pd.DataFrame(rows, columns=["ensembl_gene_id", "gene_symbol", "group", "value"])
    long_table["value"] = pd.to_numeric(long_table.value, errors="coerce")
    return long_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"),
                        help="read only for the list of genes the drug perturbations seed")
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet", columns=["source_id", "target_id", "relation_type"])
    print(f"Graph: {arguments.graph_dir}, {len(nodes)} nodes")
    print("  nodes by type: " + ", ".join(f"{kind} {count}" for kind, count in nodes.node_type.value_counts().items()))
    print(f"  with a gene symbol: {int(nodes.gene_symbol.notna().sum())}; "
          f"with an Ensembl gene id: {int(nodes.ensembl_gene_id.notna().sum())}")
    print(f"  with brain_median_tpm_max, the one tissue number the graph has: {int(nodes.brain_median_tpm_max.notna().sum())}")

    gtex = gtex_median_tpm(GTEX_MEDIAN_TPM)
    tissue_columns = [column for column in gtex.columns if column not in ("Name", "gene_symbol")]
    brain_columns = [column for column in tissue_columns if column.startswith("Brain_")]
    print(f"\nGTEx v10: {len(gtex)} genes, {len(tissue_columns)} tissues, {len(brain_columns)} of them brain regions")

    # Matching the graph's gene nodes to GTEx, by Ensembl id first and by symbol for the rest.
    gene_nodes = nodes[nodes.gene_symbol.notna()]
    by_ensembl = gene_nodes.ensembl_gene_id.isin(gtex.index)
    symbol_to_ensembl = gtex.reset_index().groupby("gene_symbol").ensembl_gene_id.first()
    by_symbol = gene_nodes.gene_symbol.isin(symbol_to_ensembl.index)
    matched = by_ensembl | by_symbol
    print(f"  graph nodes with a gene symbol: {len(gene_nodes)}")
    print(f"    matched to GTEx by Ensembl id: {int(by_ensembl.sum())}")
    print(f"    matched by gene symbol (either route): {int(by_symbol.sum())}")
    print(f"    matched by either: {int(matched.sum())} ({matched.mean():.1%})")

    # The number that decides whether a tissue partition exists at all.
    resolved_ensembl = np.where(by_ensembl, gene_nodes.ensembl_gene_id,
                                gene_nodes.gene_symbol.map(symbol_to_ensembl))
    present = pd.Series(resolved_ensembl, index=gene_nodes.index).dropna()
    present = present[present.isin(gtex.index)]
    expression = gtex.loc[present.values, tissue_columns].apply(pd.to_numeric, errors="coerce")
    tissues_expressed = (expression >= EXPRESSED_TPM).sum(axis=1)
    print(f"\n  how many of the {len(tissue_columns)} tissues each matched gene reaches at or above {EXPRESSED_TPM} TPM:")
    print(f"    genes measured: {len(tissues_expressed)}; median {int(tissues_expressed.median())}, "
          f"quartiles {[int(value) for value in np.percentile(tissues_expressed, [25, 75])]}")
    print(f"    in every tissue: {int((tissues_expressed == len(tissue_columns)).sum())}; "
          f"in one tissue only: {int((tissues_expressed == 1).sum())}; "
          f"in none: {int((tissues_expressed == 0).sum())}")
    print(f"    in more than half the tissues: {int((tissues_expressed > len(tissue_columns) / 2).sum())} "
          f"({(tissues_expressed > len(tissue_columns) / 2).mean():.1%} of the matched genes)")
    largest_share = (expression.max(axis=1) / expression.sum(axis=1).replace(0, np.nan))
    print(f"    share of a gene's total expression taken by its largest tissue: median {largest_share.median():.3f} "
          f"against {1 / len(tissue_columns):.3f} for an even spread")

    for name, path, group_column, value_column in (
            ("HPA brain regions", HPA_BRAIN_REGION, "Brain region", "nTPM"),
            ("HPA single-nuclei cell types", HPA_SINGLE_NUCLEI, "Cluster type", "nCPM")):
        long_table = hpa_table(path, group_column, value_column)
        groups = sorted(long_table.group.unique())
        genes = set(long_table.ensembl_gene_id) | set(long_table.gene_symbol)
        covered = gene_nodes.ensembl_gene_id.isin(genes) | gene_nodes.gene_symbol.isin(genes)
        print(f"\n{name}: {len(groups)} groups, {long_table.ensembl_gene_id.nunique()} genes")
        print(f"  graph gene nodes covered: {int(covered.sum())} ({covered.mean():.1%})")

    # The metabolic layer: a reaction can inherit a tissue from its gene reaction rule, and a metabolite from the
    # reactions it takes part in. These are the only routes, because neither kind of node has a gene of its own.
    reaction_nodes = nodes[nodes.node_type == "reaction"] if "reaction" in set(nodes.node_type) else nodes[nodes.gene_reaction_rule.notna()]
    with_rule = reaction_nodes.gene_reaction_rule.notna() & (reaction_nodes.gene_reaction_rule.astype(str).str.strip() != "")
    print("\nThe metabolic layer, the route for nodes that have no gene of their own:")
    print(f"  reaction nodes: {len(reaction_nodes)}; with a gene reaction rule: {int(with_rule.sum())} ({with_rule.mean():.1%})")

    reaction_ids = set(reaction_nodes.node_id[with_rule])
    metabolite_nodes = nodes[nodes.base_metabolite_id.notna()]
    touching = set(edges.source_id[edges.target_id.isin(reaction_ids)]) | set(edges.target_id[edges.source_id.isin(reaction_ids)])
    reachable = metabolite_nodes.node_id.isin(touching)
    print(f"  metabolite nodes: {len(metabolite_nodes)}; taking part in at least one reaction that has a gene rule: "
          f"{int(reachable.sum())} ({reachable.mean():.1%})")

    # The part that bears on the owner's sentence "the tissue type of a particular mechanism matters": not every node
    # needs a tissue for that, only the drugs' mechanism targets, so this is their tissue specificity measured on its
    # own. The graph's current answer for a target is a two-way choice, brain or plasma, from brain_expressed.
    if arguments.evidence_dir.exists():
        target_genes = sorted(mechanism_target_genes(arguments.evidence_dir, nodes))
        in_gtex = [gene for gene in target_genes if gene in symbol_to_ensembl.index]
        target_expression = gtex.loc[symbol_to_ensembl.loc[in_gtex].values, tissue_columns].apply(pd.to_numeric, errors="coerce")
        target_expression.index = in_gtex
        totals = target_expression.sum(axis=1).replace(0, np.nan)
        largest = target_expression.max(axis=1) / totals
        dominant_tissue = target_expression.idxmax(axis=1)
        brain_share = target_expression[brain_columns].sum(axis=1) / totals
        print(f"\nThe drugs' mechanism target genes: {len(target_genes)} named, {len(in_gtex)} found in GTEx")
        print(f"  share of total expression in the single largest tissue: median {largest.median():.3f}, "
              f"above 0.25 for {int((largest > 0.25).sum())} genes, above 0.5 for {int((largest > 0.5).sum())}")
        print(f"  share of total expression in the {len(brain_columns)} brain regions together: median {brain_share.median():.3f}, "
              f"above 0.5 for {int((brain_share > 0.5).sum())} genes, below 0.1 for {int((brain_share < 0.1).sum())}")
        print("  the ten targets most concentrated in one tissue:")
        for gene in largest.sort_values(ascending=False).head(10).index:
            print(f"    {gene}: {largest[gene]:.2f} of its expression in {dominant_tissue[gene]}, "
                  f"brain share {brain_share[gene]:.2f}")
        print("  the ten least concentrated, where a tissue label would say least:")
        for gene in largest.sort_values().head(10).index:
            print(f"    {gene}: {largest[gene]:.2f} of its expression in {dominant_tissue[gene]}, "
                  f"brain share {brain_share[gene]:.2f}")

        # Whether a tissue share would change the site a target's concentration is read at, which is the only place
        # a tissue label touches a number today. brain_expressed is a threshold on the maximum brain median TPM, and
        # a gene can clear that threshold while almost all of its expression sits in another tissue.
        brain_flag_of_gene = dict(zip(nodes.gene_symbol[nodes.node_type == "gene"],
                                      nodes.brain_expressed[nodes.node_type == "gene"]))
        flagged_brain = pd.Series({gene: bool(brain_flag_of_gene.get(gene) is True) for gene in in_gtex})
        would_differ_towards_plasma = flagged_brain & (brain_share < 0.25)
        would_differ_towards_brain = (~flagged_brain) & (brain_share > 0.5)
        print("\n  against the graph's brain_expressed boolean, which is what picks the site today:")
        print(f"    targets flagged brain-expressed: {int(flagged_brain.sum())} of {len(in_gtex)}")
        print(f"    flagged brain yet under a quarter of their expression in brain: {int(would_differ_towards_plasma.sum())} "
              f"({sorted(brain_share.index[would_differ_towards_plasma])})")
        print(f"    not flagged yet over half their expression in brain: {int(would_differ_towards_brain.sum())} "
              f"({sorted(brain_share.index[would_differ_towards_brain])})")
    else:
        print(f"\n(no evidence table at {arguments.evidence_dir}, so the mechanism targets were not measured)")

    compartments = nodes.compartment.dropna()
    print(f"\nThe compartment column, for comparison: {int(compartments.notna().sum())} nodes carry one, "
          f"{compartments.nunique()} distinct strings")
    print(f"  nodes with a gene symbol that carry a compartment: "
          f"{int(nodes.compartment.notna()[nodes.gene_symbol.notna()].sum())}")


if __name__ == "__main__":
    main()
