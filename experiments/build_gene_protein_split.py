"""Write a split copy of a graph: each gene node gets one protein node per reviewed UniProt entry, joined by encodes
edges (mechanistic_pathway_learning/graph/gene_protein_split.py; docs/gene_protein_split.md).

Reads a graph directory and writes a new one. It refuses to write into a directory that already holds a graph unless
--overwrite is given, so a graph a running job reads is never regenerated. It writes:
- nodes.parquet, edges.parquet and relation_types.json, with encodes appended last and degree recomputed;
- gene_to_protein.parquet, which tells experiment_data.load_experiment_data to seed a drug's target proteins (the
  drug_target rows) and a knockout's gene;
- node_descriptors.parquet, if the source graph has one, with the protein block on the protein nodes (one row per
  entry, from --protein-descriptors-per-entry) and the gene brain expression block left on the gene nodes;
- split_summary.json.
Each --descriptor-table and --cell-class-weights pair (source:destination) is placed the same way into a new file.
Each --merged-descriptor-table pair writes a merged-graph copy of the source in which genes the split resolved through
HGNC or an accession take their entries' protein descriptors, so a merged twin carries the same descriptors.

A gene with several entries gets a protein node per group of entries the data tell apart (mechanistic_pathway_learning/
graph/gene_protein_split.py, entry_groups_of_gene): the records are the graph's edge ends whose source names some of
the gene's entries and not all (the accession columns of the OmniPath and CollecTRI files the graph was built from, the
member accessions of the Reactome entities) and the ChEMBL targets that do the same. An edge end at a split gene goes to
the groups its source names, and a drug acting on it seeds the groups its ChEMBL target components name. Gene symbols
no reviewed entry lists are resolved through the HGNC table's previous symbols.

Usage:
  python experiments/build_gene_protein_split.py --graph-dir data/processed/graph --output-dir data/processed/graph_split \\
      --descriptor-table data/processed/node_descriptors/slice_descriptors_brain_expression.parquet:data/processed/node_descriptors/slice_split_descriptors_brain_expression.parquet \\
      --cell-class-weights data/processed/cell_class_weights/slice_cell_class_weights.parquet:data/processed/cell_class_weights/slice_split_cell_class_weights.parquet
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import GENE_TO_PROTEIN_FILE
from mechanistic_pathway_learning.graph.build_physiology_graph import interaction_sign
from mechanistic_pathway_learning.graph.gene_protein_split import (
    ENTRY_SEPARATOR,
    PROTEIN_DESCRIPTOR_FLAG,
    PROTEIN_NODE_TYPE,
    UNIPROT_ENTRY_COLUMN,
    UNIPROT_GENE_NAMES_COLUMN,
    cell_class_weights_for_split,
    entries_named_by_interaction_rows,
    entries_of_gene_symbol,
    entry_groups_of_gene,
    place_descriptors_on_split,
    protein_node_assignment,
    resolved_protein_descriptors_on_merged,
    resolve_entries_of_gene_symbols,
    split_graph,
)
from mechanistic_pathway_learning.graph.reactome_import import parse_reactome_sbml

REACTOME_MEMBER_RELATION = "member_of"
REACTOME_MEMBER_SOURCE = "Reactome entity member"
REACTOME_ENTITY_PREFIX = "RCTE_"
ONE_TO_ONE_TOLERANCE = 1e-9


def source_and_destination(pair: str) -> tuple[Path, Path]:
    source, separator, destination = pair.partition(":")
    if not separator or not destination:
        raise argparse.ArgumentTypeError(f"expected source:destination, got {pair!r}")
    return Path(source), Path(destination)


def refuse_existing(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise SystemExit(f"{path} exists; a split is written to a new place (pass --overwrite only when no job reads it)")


def entries_named_by_sources(edges: pd.DataFrame, entries_of_symbol: dict[str, list[str]], omnipath_path: Path, collectri_path: Path,
                             reactome_directory: Path) -> dict[tuple[str, str, str, str, str], set[str]]:
    """split_graph's named_entries from the raw files the graph was built from (only ends at genes with several entries)."""
    named: dict[tuple[str, str, str, str, str], set[str]] = {}
    for path, evidence_source, relation_of_row in ((omnipath_path, "OmniPath", lambda row: interaction_sign(row)[0]),
                                                   (collectri_path, "CollecTRI", lambda row: "regulates_transcription_of")):
        if path.exists():
            with open(path, encoding="utf-8") as handle:
                named.update(entries_named_by_interaction_rows(csv.DictReader(handle, delimiter="\t"), evidence_source, entries_of_symbol, relation_of_row))
    member_edges = edges[(edges.relation_type == REACTOME_MEMBER_RELATION) & edges.target_id.astype(str).str.startswith(REACTOME_ENTITY_PREFIX)]
    if len(member_edges) and reactome_directory.exists():
        accessions_of_entity: dict[str, set[str]] = defaultdict(set)
        for path in sorted(reactome_directory.glob("*.sbml")):
            species, _ = parse_reactome_sbml(path.read_text())
            for species_id, entity in species.items():
                accessions_of_entity[f"{REACTOME_ENTITY_PREFIX}{species_id.split('_')[-1]}"].update(entity.uniprot_ids)
        for source_id, target_id, evidence_source in zip(member_edges.source_id, member_edges.target_id, member_edges.evidence_source):
            symbol = str(source_id).removeprefix("GENE:")
            if len(entries_of_symbol.get(symbol, [])) > 1:
                named[(source_id, target_id, REACTOME_MEMBER_RELATION, evidence_source, "source")] = accessions_of_entity.get(target_id, set()) & set(entries_of_symbol[symbol])
    return named


def drug_target_entries_from_chembl(targets_path: Path, entries_of_symbol: dict[str, list[str]]) -> tuple[dict[str, set[str]], list[str]]:
    """gene symbol -> the accessions ChEMBL target components name for it, for genes with several entries; and the genes
    whose targets name different entries (each drug of such a gene would need its own reading)."""
    if not targets_path.exists():
        return {}, []
    symbols_of_entry: dict[str, set[str]] = defaultdict(set)
    for symbol, entries in entries_of_symbol.items():
        for entry in entries:
            symbols_of_entry[entry].add(symbol)
    named: dict[str, set[str]] = defaultdict(set)
    per_target: dict[str, set[frozenset[str]]] = defaultdict(set)
    for target in json.loads(targets_path.read_text()).values():
        by_symbol: dict[str, set[str]] = defaultdict(set)
        for accession in target.get("accessions", []):
            for symbol in symbols_of_entry.get(accession, ()):
                if len(entries_of_symbol[symbol]) > 1:
                    by_symbol[symbol].add(accession)
        for symbol, accessions in by_symbol.items():
            named[symbol] |= accessions
            per_target[symbol].add(frozenset(accessions))
    return dict(named), sorted(symbol for symbol, sets in per_target.items() if len(sets) > 1)


def chembl_target_records(targets_path: Path, entries_of_symbol: dict[str, list[str]]) -> dict[str, set]:
    """accession -> the ChEMBL targets whose components name it and not all of its gene's entries."""
    records: dict[str, set] = defaultdict(set)
    if not targets_path.exists():
        return records
    gene_entry_sets = [set(entries) for entries in entries_of_symbol.values() if len(entries) > 1]
    for target_id, target in json.loads(targets_path.read_text()).items():
        accessions = set(target.get("accessions", []))
        for entries in gene_entry_sets:
            subset = accessions & entries
            if subset and subset != entries:
                for accession in subset:
                    records[accession].add(("ChEMBL target", target_id))
    return records


def edge_end_records(edges: pd.DataFrame, named_entries: dict, entries_of_symbol: dict[str, list[str]], symbol_of_node: dict[str, str]) -> dict[str, set]:
    """accession -> the graph's edge ends whose source names it and not all of its gene's entries. Ends of rows the
    graph does not hold are not counted."""
    records: dict[str, set] = defaultdict(set)
    edge_keys = set(zip(edges.source_id, edges.target_id, edges.relation_type, edges.evidence_source.astype(str)))
    for key, accessions in named_entries.items():
        if key[:4] not in edge_keys:
            continue
        entries = set(entries_of_symbol.get(symbol_of_node.get(key[0] if key[4] == "source" else key[1], ""), []))
        subset = set(accessions) & entries
        if subset and subset != entries:
            for accession in subset:
                records[accession].add(("edge end",) + tuple(key))
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--uniprot-table", type=Path, default=Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz"))
    parser.add_argument("--protein-descriptors-per-entry", type=Path, default=Path("data/processed/node_descriptors/protein_descriptors_per_entry.parquet"),
                        help="per-entry protein descriptors (experiments/build_protein_descriptors.py --per-entry-only)")
    parser.add_argument("--omnipath-interactions", type=Path, default=Path("data/raw/omnipath/omnipath_interactions.tsv"))
    parser.add_argument("--collectri-interactions", type=Path, default=Path("data/raw/omnipath/collectri_interactions.tsv"))
    parser.add_argument("--reactome-dir", type=Path, default=Path("data/raw/reactome/v97"))
    parser.add_argument("--chembl-targets", type=Path, default=Path("data/raw/chembl/targets.json"))
    parser.add_argument("--hgnc-table", type=Path, default=Path("data/raw/hgnc/hgnc_complete_set.txt"),
                        help="resolves gene symbols no reviewed entry lists through their previous HGNC symbols")
    parser.add_argument("--descriptor-table", type=source_and_destination, action="append", default=[])
    parser.add_argument("--cell-class-weights", type=source_and_destination, action="append", default=[])
    parser.add_argument("--merged-descriptor-table", type=source_and_destination, action="append", default=[],
                        help="source:destination; a merged-graph copy with the protein descriptors of the resolved genes filled")
    parser.add_argument("--overwrite", action="store_true")
    arguments = parser.parse_args()

    refuse_existing(arguments.output_dir / "nodes.parquet", arguments.overwrite)
    for _, destination in arguments.descriptor_table + arguments.cell_class_weights + arguments.merged_descriptor_table:
        refuse_existing(destination, arguments.overwrite)

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet")
    relation_types = json.loads((arguments.graph_dir / "relation_types.json").read_text())
    uniprot_entries = pd.read_csv(arguments.uniprot_table, sep="\t", usecols=[UNIPROT_ENTRY_COLUMN, UNIPROT_GENE_NAMES_COLUMN])
    hgnc_table = pd.read_csv(arguments.hgnc_table, sep="\t", usecols=["symbol", "prev_symbol"], dtype=str) if arguments.hgnc_table.exists() else None
    genes = nodes[nodes.node_type == "gene"]
    symbol_of_node = dict(zip(genes.node_id, genes.gene_symbol.astype(str)))
    resolved, how_resolved = resolve_entries_of_gene_symbols(set(symbol_of_node.values()), uniprot_entries, hgnc_table)
    entries_of_symbol = {**entries_of_gene_symbol(uniprot_entries), **{symbol: entries for symbol, entries in resolved.items() if entries}}
    drug_target_entries, conflicting_drug_target_genes = drug_target_entries_from_chembl(arguments.chembl_targets, entries_of_symbol)
    named_entries = entries_named_by_sources(edges, entries_of_symbol, arguments.omnipath_interactions, arguments.collectri_interactions, arguments.reactome_dir)
    records_of_entry = edge_end_records(edges, named_entries, entries_of_symbol, symbol_of_node)
    for accession, records in chembl_target_records(arguments.chembl_targets, entries_of_symbol).items():
        records_of_entry.setdefault(accession, set()).update(records)
    assignment = protein_node_assignment(nodes, uniprot_entries, drug_target_entries, records_of_entry, hgnc_table)
    split_nodes, split_edges, split_relations, summary = split_graph(nodes, edges, relation_types, assignment, named_entries)

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    split_nodes.to_parquet(arguments.output_dir / "nodes.parquet", index=False)
    split_edges.to_parquet(arguments.output_dir / "edges.parquet", index=False)
    (arguments.output_dir / "relation_types.json").write_text(json.dumps(split_relations, indent=1) + "\n")
    assignment.to_parquet(arguments.output_dir / GENE_TO_PROTEIN_FILE, index=False)

    per_entry = pd.read_parquet(arguments.protein_descriptors_per_entry) if arguments.protein_descriptors_per_entry.exists() else None
    tables = list(arguments.descriptor_table)
    if (arguments.graph_dir / "node_descriptors.parquet").exists():
        tables.insert(0, (arguments.graph_dir / "node_descriptors.parquet", arguments.output_dir / "node_descriptors.parquet"))
    node_type = split_nodes.set_index("node_id").node_type
    single_protein = assignment[~assignment.several_proteins & ~assignment.shared].set_index("protein_node_id").gene_node_id
    summary["descriptor_tables"] = {}
    for source, destination in tables:
        table = pd.read_parquet(source)
        placed = place_descriptors_on_split(table, split_nodes, assignment, per_entry)
        destination.parent.mkdir(parents=True, exist_ok=True)
        placed.to_parquet(destination)
        is_protein = (node_type.reindex(placed.index) == PROTEIN_NODE_TYPE).to_numpy()
        is_gene = (node_type.reindex(placed.index) == "gene").to_numpy()
        protein_columns = [column for column in placed.columns if column.startswith("protein_")]
        other_columns = [column for column in placed.columns if not column.startswith("protein_")]
        record = {"source": str(source), "rows": len(placed), "columns": len(placed.columns),
                  "gene_rows_with_a_nonzero_protein_column": int((placed.loc[is_gene, protein_columns] != 0).any(axis=1).sum()),
                  "protein_rows_with_a_nonzero_other_column": int((placed.loc[is_protein, other_columns] != 0).any(axis=1).sum()),
                  "gene_rows_with_a_nonzero_other_column": int((placed.loc[is_gene, other_columns] != 0).any(axis=1).sum())}
        if PROTEIN_DESCRIPTOR_FLAG in placed.columns:
            record["protein_nodes_with_protein_descriptors"] = int(placed.loc[is_protein, PROTEIN_DESCRIPTOR_FLAG].sum())
            record["gene_nodes_with_protein_descriptors_before"] = int(table.loc[table.index.intersection(assignment.gene_node_id), PROTEIN_DESCRIPTOR_FLAG].sum())
            one_to_one = single_protein[single_protein.isin(table.index)]
            record["one_to_one_protein_rows_differing_from_their_gene_row_before"] = int(  # the per-entry refit agrees to about 1e-12
                (np.abs(placed.loc[one_to_one.index, protein_columns].to_numpy() - table.loc[one_to_one.to_numpy(), protein_columns].to_numpy())
                 > ONE_TO_ONE_TOLERANCE).any(axis=1).sum())
        unchanged = table.index.intersection(placed.index)
        record["other_columns_unchanged_on_existing_nodes"] = bool(np.array_equal(placed.loc[unchanged, other_columns].to_numpy(dtype=float),
                                                                                  table.loc[unchanged, other_columns].fillna(0.0).to_numpy(dtype=float)))
        summary["descriptor_tables"][str(destination)] = record
    summary["merged_descriptor_tables"] = {}
    for source, destination in arguments.merged_descriptor_table:
        table = pd.read_parquet(source)
        filled, genes_filled = resolved_protein_descriptors_on_merged(table, assignment, place_descriptors_on_split(table, split_nodes, assignment, per_entry))
        destination.parent.mkdir(parents=True, exist_ok=True)
        filled.to_parquet(destination)
        other_columns = [column for column in table.columns if not column.startswith("protein_")]
        summary["merged_descriptor_tables"][str(destination)] = {
            "source": str(source), "gene_nodes_filled": genes_filled,
            "rows_changed": int((filled.fillna(0.0) != table.fillna(0.0)).any(axis=1).sum()),
            "other_columns_unchanged": bool(filled[other_columns].fillna(0.0).equals(table[other_columns].fillna(0.0)))}
    summary["cell_class_weights"] = {}
    for source, destination in arguments.cell_class_weights:
        weights = cell_class_weights_for_split(pd.read_parquet(source), split_nodes, assignment)
        destination.parent.mkdir(parents=True, exist_ok=True)
        weights.to_parquet(destination)
        summary["cell_class_weights"][str(destination)] = {"source": str(source), "rows": len(weights), "rows_with_a_missing_value": int(weights.isna().any(axis=1).sum())}

    is_gene = split_nodes.node_type == "gene"
    is_protein = split_nodes.node_type == PROTEIN_NODE_TYPE
    several = assignment[assignment.several_proteins]
    summary.update({
        "graph_dir": str(arguments.graph_dir), "output_dir": str(arguments.output_dir),
        "nodes_before": len(nodes), "nodes_after": len(split_nodes),
        "gene_degree_quantiles": {str(q): float(np.quantile(split_nodes.degree[is_gene], q)) for q in (0.25, 0.5, 0.75, 0.95)},
        "protein_degree_quantiles": {str(q): float(np.quantile(split_nodes.degree[is_protein], q)) for q in (0.25, 0.5, 0.75, 0.95)},
        "gene_nodes_of_degree_one": int((split_nodes.degree[is_gene] == 1).sum()),
        "shared_protein_node_names": sorted(assignment[assignment.shared].protein_display_name.unique().tolist()),
        "protein_nodes_of_each_gene_with_several": {gene: group.protein_node_id.tolist() for gene, group in several.groupby("gene_node_id")},
        "entries_of_each_protein_node_of_several_entries": {protein: entries for protein, entries in assignment.drop_duplicates("protein_node_id")[["protein_node_id", "uniprot_entry"]].itertuples(index=False) if ENTRY_SEPARATOR in entries},
        "entries_left_out": {symbol: entry_groups_of_gene(resolved[symbol], records_of_entry)[1] for symbol in sorted(symbol_of_node.values()) if entry_groups_of_gene(resolved[symbol], records_of_entry)[1]},
        "records_telling_entries_apart": {accession: sorted(" | ".join(map(str, record)) for record in records) for accession, records in sorted(records_of_entry.items())},
        "gene_symbols_resolved_without_a_listing_entry": dict(sorted(how_resolved.items())),
        "gene_nodes_kept_without_a_protein_node": sorted(set(symbol_of_node.values()) - set(symbol_of_node[gene] for gene in assignment.gene_node_id)),
        "drug_target_protein_nodes_of_genes_with_several": several[several.drug_target & (several.groupby("gene_node_id").drug_target.transform("sum") < several.groupby("gene_node_id").protein_node_id.transform("size"))].protein_node_id.tolist(),
        "genes_whose_drug_targets_name_different_entries": conflicting_drug_target_genes,
    })
    (arguments.output_dir / "split_summary.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key not in ("shared_protein_node_names", "edges_per_relation", "protein_nodes_of_each_gene_with_several",
                                                                                  "records_telling_entries_apart", "gene_nodes_kept_without_a_protein_node",
                                                                                  "entries_of_each_protein_node_of_several_entries")},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
