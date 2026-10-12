"""How much of the study is already enzymes, and what a curated enzyme table would add.

The owner, 12 October 2026: "I can't realistically claim physiological mechanistic pathways without enzymes and
without carriers. I could probably live without carriers but a lot of very important meds act on enzymes themselves
... Even if we have to curate our own enzymes data table, I think it would be worth it but I'm willing to hear a
dissenting opinion if you have one."

Three different things were being called "enzymes" in that exchange, and they have different answers, so this script
counts each separately:

  1. enzymes as nodes and edges of the physiology graph. Every metabolic reaction in Human-GEM is enzyme-catalysed
     and its catalysing genes are nodes, so this is a question about how much of the graph is enzymatic.
  2. enzymes as drug mechanism targets. Monoamine oxidase inhibitors, cholinesterase inhibitors and COMT inhibitors
     act on enzymes, and every drug perturbation trains the mechanism edges, so these already carry gradient.
  3. enzymes as labelled gene perturbations. A labelled knockout or gain of function of an enzyme gene trains the
     graph around that enzyme directly. This is the count that decides whether the study learns enzyme physiology.

A fourth thing is separate from all three: the `clears` relation, an enzyme that metabolises a drug. That is drug
metabolism rather than pathway physiology, and it is the only one of the four that a curated table of drug-clearing
enzymes would add to.

An enzyme is a gene whose reviewed UniProt entry carries an EC number. Beside it, membership of a Human-GEM gene
reaction rule is reported, which is the narrower question of whether the graph's own metabolic layer uses the gene.

Counting only: no label, table, graph or model input is changed.

Usage:
  python experiments/label_review/enzymes_in_the_graph.py
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

UNIPROT_REVIEWED = Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz")
RULE_TOKEN = re.compile(r"[A-Za-z0-9_.-]+")


def enzyme_gene_symbols(uniprot_path: Path) -> set[str]:
    """Primary gene symbols of reviewed human entries that carry an EC number."""
    table = pd.read_csv(uniprot_path, sep="\t", compression="gzip",
                        usecols=["Gene Names (primary)", "EC number"], low_memory=False)
    with_ec = table[table["EC number"].notna() & (table["EC number"].astype(str).str.strip() != "")]
    symbols = set()
    for names in with_ec["Gene Names (primary)"].dropna():
        for symbol in str(names).split(";"):
            if symbol.strip():
                symbols.add(symbol.strip())
    return symbols


def genes_named_in_a_reaction_rule(nodes: pd.DataFrame) -> set[str]:
    """Every identifier that appears in a Human-GEM gene reaction rule, which is how the metabolic layer names its
    catalysing genes (Ensembl gene ids, with a few symbols)."""
    named = set()
    for rule in nodes.gene_reaction_rule.dropna():
        for token in RULE_TOKEN.findall(str(rule)):
            if token.lower() not in ("and", "or"):
                named.add(token)
    return named


def mechanism_and_labelled_genes(evidence_directory: Path, nodes: pd.DataFrame) -> tuple[dict[str, set[str]], set[str]]:
    """({drug label -> its mechanism target gene symbols}, the gene symbols that have a labelled gene perturbation)."""
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet",
                              columns=["perturbation_type", "perturbation_label", "perturbation_nodes"])
    symbol_of_node = dict(zip(nodes.node_id, nodes.gene_symbol))
    mechanism_genes_of_drug: dict[str, set[str]] = {}
    labelled_gene_symbols: set[str] = set()
    for row in records.drop_duplicates(subset=["perturbation_type", "perturbation_label"]).itertuples(index=False):
        symbols = {symbol_of_node.get(node_id) for node_id, _sign, _magnitude in json.loads(row.perturbation_nodes)}
        symbols = {symbol for symbol in symbols if isinstance(symbol, str) and symbol}
        if row.perturbation_type == "drug":
            mechanism_genes_of_drug[row.perturbation_label] = symbols
        else:
            labelled_gene_symbols |= symbols
    return mechanism_genes_of_drug, labelled_gene_symbols


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--clearing-enzymes", type=Path, default=Path("configs/drug_clearing_enzymes.csv"))
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    enzymes = enzyme_gene_symbols(UNIPROT_REVIEWED)
    in_a_rule = genes_named_in_a_reaction_rule(nodes)
    print(f"Reviewed human UniProt entries with an EC number: {len(enzymes)} gene symbols")

    print("\n1. Enzymes as nodes and edges of the graph")
    gene_nodes = nodes[nodes.node_type == "gene"]
    protein_nodes = nodes[nodes.node_type == "protein"]
    print(f"  gene nodes: {len(gene_nodes)}; of them enzymes by EC number: "
          f"{int(gene_nodes.gene_symbol.isin(enzymes).sum())} "
          f"({gene_nodes.gene_symbol.isin(enzymes).mean():.1%})")
    print(f"  protein nodes: {len(protein_nodes)}; of them enzymes by EC number: "
          f"{int(protein_nodes.gene_symbol.isin(enzymes).sum())}")
    print(f"  gene nodes named in a Human-GEM gene reaction rule, so catalysing a reaction the graph holds: "
          f"{int(gene_nodes.ensembl_gene_id.isin(in_a_rule).sum() + gene_nodes.gene_symbol.isin(in_a_rule).sum())}")
    reactions = nodes[nodes.node_type == "reaction"]
    with_rule = reactions.gene_reaction_rule.notna() & (reactions.gene_reaction_rule.astype(str).str.strip() != "")
    print(f"  reaction nodes: {len(reactions)}; with a gene reaction rule naming their enzymes: {int(with_rule.sum())}")

    mechanism_genes_of_drug, labelled_gene_symbols = mechanism_and_labelled_genes(arguments.evidence_dir, nodes)
    all_mechanism_genes = set().union(*mechanism_genes_of_drug.values()) if mechanism_genes_of_drug else set()
    enzyme_mechanism_genes = sorted(all_mechanism_genes & enzymes)
    drugs_on_an_enzyme = sorted(label for label, symbols in mechanism_genes_of_drug.items() if symbols & enzymes)
    print("\n2. Enzymes as drug mechanism targets, which every drug perturbation already trains")
    print(f"  mechanism target genes: {len(all_mechanism_genes)}; of them enzymes: {len(enzyme_mechanism_genes)}")
    print(f"    {enzyme_mechanism_genes}")
    print(f"  drugs whose mechanism includes an enzyme: {len(drugs_on_an_enzyme)} of {len(mechanism_genes_of_drug)}")
    print(f"    {drugs_on_an_enzyme}")

    enzyme_labelled_genes = sorted(labelled_gene_symbols & enzymes)
    print("\n3. Enzymes as labelled gene perturbations, which train the graph around the enzyme itself")
    print(f"  gene perturbations: {len(labelled_gene_symbols)} gene symbols; of them enzymes: "
          f"{len(enzyme_labelled_genes)} ({len(enzyme_labelled_genes) / max(len(labelled_gene_symbols), 1):.1%})")
    ensembl_of_symbol = dict(zip(nodes.gene_symbol, nodes.ensembl_gene_id))
    in_a_rule_by_symbol = {symbol for symbol in enzyme_labelled_genes
                           if symbol in in_a_rule or ensembl_of_symbol.get(symbol) in in_a_rule}
    print(f"  of those, catalysing a reaction the graph holds (named in a Human-GEM reaction rule): "
          f"{len(in_a_rule_by_symbol)}")
    print(f"  the first thirty: {enzyme_labelled_genes[:30]}")

    print("\n4. The clears relation, drug metabolism, which is the only one a curated table would add to")
    if arguments.clearing_enzymes.exists():
        clearing = pd.read_csv(arguments.clearing_enzymes)
        with_an_edge = clearing[clearing.gives_an_edge.astype(str).str.lower() == "true"]
        clearing_enzymes = sorted(with_an_edge.enzyme_gene.dropna().unique())
        print(f"  {arguments.clearing_enzymes}: {len(clearing)} rows read, {len(with_an_edge)} give an edge, "
              f"{len(clearing_enzymes)} enzymes, {with_an_edge.drug.nunique()} drugs")
        print(f"  the enzymes: {clearing_enzymes}")
        print(f"  how many are already gene nodes of the graph: "
              f"{len([symbol for symbol in clearing_enzymes if symbol in set(gene_nodes.gene_symbol)])}")
        print(f"  how many carry an EC number, so are enzymes by the definition used above: "
              f"{len([symbol for symbol in clearing_enzymes if symbol in enzymes])}")
        print(f"  how many have a labelled gene perturbation of their own, which is what would train them: "
              f"{len(set(clearing_enzymes) & labelled_gene_symbols)} "
              f"({sorted(set(clearing_enzymes) & labelled_gene_symbols)})")
    else:
        print(f"  ({arguments.clearing_enzymes} not found)")


if __name__ == "__main__":
    main()
