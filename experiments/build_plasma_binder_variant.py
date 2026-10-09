"""Build the plasma-binder graph variant (mechanistic_pathway_learning/graph/plasma_binding.py) of a graph directory:
<graph-dir>_binders, with binds edges, sign 0, from each extracellular cargo metabolite to its binder's gene node.

Stage 1 of docs/plasma_protein_binding.md (the user's decision of 9 October 2026). It writes a new directory, so no
preregistered input changes; folding the edges into the confirmatory graph is a separate decision, because it would
need a dated amendment and the development pilots re-run.

Inputs:
  - data/raw/uniprot/uniprot_plasma_binders.tsv.gz, pinned by experiments/fetch_plasma_binder_annotations.py;
  - docs/curated_plasma_carriage.csv, one row per (binder, cargo) with its source, PMID, DOI and a quoted sentence;
  - data/raw/chebi/chebi_core.obo.gz and the Human-GEM metabolite table, for the ChEBI route that the Reactome import
    also uses (mechanistic_pathway_learning/evidence/laboratory_abnormality_labels.map_definition_to_metabolites).

No node is added. Every binder is already a gene node and every cargo that reaches plasma already has an extracellular
copy, so a row that finds neither is reported in the summary rather than invented. Node descriptors are unchanged
(the node set is the same), so the source graph's node_descriptors.parquet is copied when it exists.

Usage:
  PYTHONPATH=. python experiments/build_plasma_binder_variant.py --graph-dirs data/processed/graph_full_neuronal
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd

from experiments.fetch_plasma_binder_annotations import BINDER_GENES
from mechanistic_pathway_learning.evaluation.experiment_data import GENE_TO_PROTEIN_FILE
from mechanistic_pathway_learning.evidence.laboratory_abnormality_labels import (
    ChemicalDefinition, human_gem_metabolites_by_chebi, map_definition_to_metabolites, read_chebi_relations)
from mechanistic_pathway_learning.graph.plasma_binding import add_plasma_carriage, curated_carriage_rows, uniprot_carriage_rows

# Copied unchanged, which is sound only because carriage adds nothing these files read. Two model inputs follow the
# wiring rather than an annotation: the reaction_brain descriptor block and the reactions' cell-class weights both carry
# gene expression onto reactions through catalyzed_by (graph/rewired_expression_features.py). Carriage adds `binds`
# edges only, so neither moves. A later variant that touched catalysis edges would have to recompute both, as the
# rewired runs do, not copy them.
COPIED_FILES = ("node_descriptors.parquet", "node_descriptors_summary.json", "gene_to_protein.parquet")


def binder_nodes_on_split(graph_directory: Path, binder_genes) -> dict[str, list[str]] | None:
    """On a gene/protein split graph, the protein nodes of each binder, because binding is the protein's property and
    the gene node there holds only expression; None on a merged graph, where the gene node is the binder."""
    mapping_file = graph_directory / GENE_TO_PROTEIN_FILE
    if not mapping_file.exists():
        return None
    mapping = pd.read_parquet(mapping_file)
    nodes_of_gene: dict[str, list[str]] = {}
    for gene_node, protein_node in zip(mapping.gene_node_id, mapping.protein_node_id):
        symbol = str(gene_node).removeprefix("GENE:")
        if symbol in set(binder_genes):
            nodes_of_gene.setdefault(symbol, []).append(str(protein_node))
    return {symbol: sorted(set(ids)) for symbol, ids in nodes_of_gene.items()}


def carriage_report(summary: dict, graph_directory: Path, output_directory: Path, curated: pd.DataFrame,
                    other_outputs: list[Path] | None = None) -> str:
    lines = ["# Plasma substrate binders in the graph", "",
             f"Built by experiments/build_plasma_binder_variant.py from `{graph_directory}` into `{output_directory}`. "
             f"{summary['edges_added']} `binds` edges, sign 0, each from the extracellular copy of a cargo metabolite to its "
             f"binder, carried on {summary['carriage_on']}, which is the direction and sign the graph's own small-molecule "
             "edges use."]
    if other_outputs:
        lines.append("")
        lines.append("The same run built " + ", ".join(f"`{output}`" for output in other_outputs)
                     + " from the matching source graphs. On a gene/protein split graph the carriage targets the protein "
                       "node, because binding is the protein's property and the gene node there holds only expression, so "
                       "the counts above hold but the edge endpoints differ.")
    lines += ["",
              "## Cargo connected, by binder", "", "| binder | cargo | edges |", "|---|---|---|"]
    for gene, cargo in sorted(summary["cargo_by_binder"].items()):
        lines.append(f"| {gene} | {', '.join(cargo)} | {summary['edges_by_binder'].get(gene, 0)} |")
    lines += ["", f"By source: {summary['edges_by_source']}. Rows read: {summary['carriage_rows_read']} "
                  f"({len(curated)} curated, {summary['carriage_rows_read'] - len(curated)} from UniProt binding-site features).", ""]
    for key, heading in (("binder_genes_absent_from_graph", "Binder genes absent from the graph"),
                         ("cargo_not_mapped_to_human_gem", "Cargo with no Human-GEM metabolite"),
                         ("cargo_without_an_extracellular_node", "Cargo with no extracellular copy in this graph"),
                         ("edges_already_in_the_graph", "Edges the graph already held")):
        lines += [f"## {heading}", "", ("\n".join(f"- {item}" for item in summary[key]) if summary[key] else "None.") , ""]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dirs", type=Path, nargs="+", default=[Path("data/processed/graph_full_neuronal")])
    parser.add_argument("--binder-table", type=Path, default=Path("data/raw/uniprot/uniprot_plasma_binders.tsv.gz"))
    parser.add_argument("--curated-table", type=Path, default=Path("docs/curated_plasma_carriage.csv"))
    parser.add_argument("--chebi-obo", type=Path, default=Path("data/raw/chebi/chebi_core.obo.gz"))
    parser.add_argument("--human-gem-metabolites", type=Path, default=Path("data/raw/Human-GEM/model/metabolites.tsv"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/plasma_binder_graph.md"),
                        help="report of the first graph built; pass none to skip")
    arguments = parser.parse_args()

    binder_table = pd.read_csv(arguments.binder_table, sep="\t")
    curated = pd.read_csv(arguments.curated_table)
    carriage_rows = pd.concat([uniprot_carriage_rows(binder_table, BINDER_GENES), curated_carriage_rows(curated)], ignore_index=True)
    names, neighbours, children = read_chebi_relations(arguments.chebi_obo)
    metabolites_by_chebi = human_gem_metabolites_by_chebi(pd.read_csv(arguments.human_gem_metabolites, sep="\t"))

    def bases_of_chebi(chebi_id: str) -> set[str]:
        mapped, _ = map_definition_to_metabolites(ChemicalDefinition("", (chebi_id,), 0, ""), names, neighbours, children, metabolites_by_chebi)
        return mapped

    for position, graph_directory in enumerate(arguments.graph_dirs):
        output_directory = graph_directory.parent / f"{graph_directory.name}_binders"
        output_directory.mkdir(parents=True, exist_ok=True)
        nodes, edges = pd.read_parquet(graph_directory / "nodes.parquet"), pd.read_parquet(graph_directory / "edges.parquet")
        relation_types = json.loads((graph_directory / "relation_types.json").read_text())
        binder_nodes_of_gene = binder_nodes_on_split(graph_directory, BINDER_GENES)
        new_nodes, new_edges, new_relations, summary = add_plasma_carriage(nodes, edges, relation_types, carriage_rows, bases_of_chebi,
                                                                          binder_nodes_of_gene)
        new_nodes.to_parquet(output_directory / "nodes.parquet")
        new_edges.to_parquet(output_directory / "edges.parquet")
        (output_directory / "relation_types.json").write_text(json.dumps(new_relations, indent=1) + "\n")
        for file_name in COPIED_FILES:
            if (graph_directory / file_name).exists():
                shutil.copyfile(graph_directory / file_name, output_directory / file_name)
        summary.update(source_graph=str(graph_directory), binder_genes_asked_for=list(BINDER_GENES),
                       carriage_on=("protein nodes (gene/protein split graph)" if binder_nodes_of_gene else "gene nodes"),
                       curated_table=str(arguments.curated_table), binder_table=str(arguments.binder_table),
                       nodes=len(new_nodes), edges_before=len(edges), edges_after=len(new_edges),
                       binder_degrees={node: int(new_nodes.loc[new_nodes.node_id == node, "degree"].iloc[0])
                                       for node in summary["binder_nodes"]})
        (output_directory / "plasma_binder_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
        print(output_directory, json.dumps({key: value for key, value in summary.items()
                                            if key not in ("cargo_by_binder", "binder_genes_asked_for")}, indent=1))
        if position == 0 and str(arguments.markdown_output) != "none":
            other_outputs = [other.parent / f"{other.name}_binders" for other in arguments.graph_dirs[1:]]
            arguments.markdown_output.write_text(carriage_report(summary, graph_directory, output_directory, curated, other_outputs))
            print(f"wrote {arguments.markdown_output}")


if __name__ == "__main__":
    main()
