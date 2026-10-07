"""Coverage of the laboratory-abnormality labels (mechanistic_pathway_learning/evidence/laboratory_abnormality_labels.py):
how many gene-term annotations of the slice and full-data genes name a metabolite change that reaches a metabolite of
the graph, by mapping route, direction and fluid. Decides whether the labels can supervise the perturbation field.

Usage:
  python experiments/check_laboratory_label_coverage.py --markdown-output docs/laboratory_label_coverage.md
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.laboratory_abnormality_labels import (
    generic_metabolite_ids,
    human_gem_metabolites_by_chebi,
    map_definition_to_metabolites,
    parse_hpo_chemical_definitions,
    read_chebi_relations,
)

DIRECTION_NAMES = {1: "increased", -1: "decreased", 0: "abnormal, no direction"}


def gene_symbols_of_evidence(evidence_path: Path) -> set[str]:
    records = pd.read_parquet(evidence_path)
    genes = records[records.perturbation_type == "gene"].perturbation_id
    return {perturbation_id.split(":", 1)[-1] for perturbation_id in genes.unique()}


def coverage_rows(label_table: pd.DataFrame, gene_symbols: set[str], graph_base_metabolites: set[str]) -> dict:
    rows = label_table[label_table.gene_symbol.isin(gene_symbols)]
    in_graph = rows[rows.metabolites_in_graph.map(len) > 0]
    return {
        "genes": len(gene_symbols),
        "gene-term pairs with a chemical definition": len(rows),
        "genes with such a pair": rows.gene_symbol.nunique(),
        "pairs mapped to a Human-GEM metabolite": int(rows.route.isin(["equivalent", "specific_form"]).sum()),
        "pairs left out as class terms": int((rows.route == "class_term").sum()),
        "  of which through an equivalent form": int((rows.route == "equivalent").sum()),
        "  of which through a specific form": int((rows.route == "specific_form").sum()),
        "pairs reaching a metabolite of the graph": len(in_graph),
        "genes with at least one such pair": in_graph.gene_symbol.nunique(),
        "distinct graph metabolites labelled": len(set().union(*in_graph.metabolites_in_graph)) if len(in_graph) else 0,
        "signed pairs reaching the graph (increased / decreased)": f"{int((in_graph.direction == 1).sum())} / {int((in_graph.direction == -1).sum())}",
        "fluids of pairs reaching the graph": ", ".join(f"{fluid} {count}" for fluid, count in Counter(in_graph.fluid).most_common()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hpo-owl", type=Path, default=Path("data/raw/hpo/hp-base.owl"))
    parser.add_argument("--hpo-annotations", type=Path, default=Path("data/raw/hpo/genes_to_phenotype.txt"))
    parser.add_argument("--chebi-obo", type=Path, default=Path("data/raw/chebi/chebi_core.obo.gz"))
    parser.add_argument("--human-gem-metabolites", type=Path, default=Path("data/raw/Human-GEM/model/metabolites.tsv"))
    parser.add_argument("--slice-graph", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--slice-evidence", type=Path, default=Path("data/processed/evidence/evidence_records.parquet"))
    parser.add_argument("--full-graph", type=Path, default=Path("data/processed/graph_full"))
    parser.add_argument("--full-evidence", type=Path, default=Path("data/processed/evidence_full/evidence_records.parquet"))
    parser.add_argument("--markdown-output", type=Path, default=None)
    arguments = parser.parse_args()

    definitions = parse_hpo_chemical_definitions(arguments.hpo_owl.read_text())
    names, neighbours, children = read_chebi_relations(arguments.chebi_obo)
    metabolites_by_chebi = human_gem_metabolites_by_chebi(pd.read_csv(arguments.human_gem_metabolites, sep="\t"))
    excluded_metabolites = frozenset(generic_metabolite_ids(pd.read_parquet(arguments.full_graph / "nodes.parquet")))
    mapping_by_term = {}
    for definition in definitions:
        metabolites, route = map_definition_to_metabolites(definition, names, neighbours, children, metabolites_by_chebi, excluded_metabolites)
        mapping_by_term[definition.hpo_id] = (definition, metabolites, route)

    annotations = pd.read_csv(arguments.hpo_annotations, sep="\t").drop_duplicates(["gene_symbol", "hpo_id"])
    annotations = annotations[annotations.hpo_id.isin(mapping_by_term)]
    sections = []
    for label, graph_directory, evidence_path in (("slice", arguments.slice_graph, arguments.slice_evidence), ("full data", arguments.full_graph, arguments.full_evidence)):
        nodes = pd.read_parquet(graph_directory / "nodes.parquet")
        graph_base_metabolites = set(nodes.loc[nodes.node_type == "metabolite", "base_metabolite_id"].dropna())
        label_table = annotations.assign(
            direction=[mapping_by_term[term][0].direction for term in annotations.hpo_id],
            fluid=[mapping_by_term[term][0].fluid for term in annotations.hpo_id],
            route=[mapping_by_term[term][2] for term in annotations.hpo_id],
            metabolites_in_graph=[mapping_by_term[term][1] & graph_base_metabolites for term in annotations.hpo_id],
        )
        sections.append((label, coverage_rows(label_table, gene_symbols_of_evidence(evidence_path), graph_base_metabolites)))

    direction_counts = Counter(definition.direction for definition in definitions)
    route_counts = Counter(route for _, _, route in mapping_by_term.values())
    lines = [
        "# Laboratory-abnormality label coverage (generated by experiments/check_laboratory_label_coverage.py)",
        "",
        f"HPO terms whose logical definition names a ChEBI entity (hp-base.owl): {len(definitions)} "
        f"({', '.join(f'{DIRECTION_NAMES[direction]} {count}' for direction, count in sorted(direction_counts.items(), reverse=True))}). "
        f"Terms mapped to a Human-GEM metabolite: {route_counts['equivalent']} through an equivalent form (conjugate acid or base, tautomer), "
        f"{route_counts['specific_form']} through a specific form (an is_a child named after the compound, such as L-phenylalanine), "
        f"{route_counts['class_term']} left out as class terms (more than three metabolites), {route_counts['unmapped']} unmapped. "
        "Generic pseudo-metabolites ('[protein]') are excluded; an unannotated gene-metabolite pair is unlabelled, not normal.",
        "",
        "| | " + " | ".join(label for label, _ in sections) + " |",
        "|---|" + "---|" * len(sections),
    ]
    for key in sections[0][1]:
        lines.append(f"| {key} | " + " | ".join(str(values[key]) for _, values in sections) + " |")
    report = "\n".join(lines) + "\n"
    print(report)
    if arguments.markdown_output:
        arguments.markdown_output.write_text(report)


if __name__ == "__main__":
    main()
