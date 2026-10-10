"""Build the drug-entry-node variant of a graph directory (mechanistic_pathway_learning/graph/drug_entry_nodes.py;
docs/drug_entry_nodes.md): <graph-dir>_drugs, with one node per drug perturbation, its mechanism edges and, for a
drug whose plasma carrier is known, the carrier's sequestration edge.

It writes a new directory and changes no existing one, so no registered input moves. Whether a confirmatory
configuration reads the variant is the user's decision and needs a dated amendment (docs/drug_entry_nodes.md).

Inputs:
  - the source graph directory, with gene_to_protein.parquet when it splits genes from proteins, so a drug's mechanism
    edges land on the protein nodes its seeds would land on;
  - the evidence table, whose drug rows name each drug perturbation and its seeds (perturbation_nodes);
  - configs/drug_plasma_carriers.csv, the carrier of each drug the labels or the fraction-unbound database name one for
    (experiments/scope_label_plasma_binder_coverage.py); pass --carrier-table none to add no sequestration edge;
  - configs/drugs_acting_as_graph_compounds.csv, the drugs that are themselves a graph metabolite and get no node.

Outputs in the new directory: nodes.parquet, edges.parquet and relation_types.json with rows appended and none changed;
drug_entry_nodes.parquet, one row per mechanism edge with its sign and magnitude, which the loader reads
(experiment_data.load_experiment_data, drug_entry); drug_entry_summary.json; and a copy of every other file of the
source directory. Each --descriptor-table and --cell-class-weights pair (source:destination) writes the per-node table
with a row per drug node: zeros for descriptors, ones for cell-class weights.

Usage:
  PYTHONPATH=. python experiments/build_drug_entry_node_variant.py --graph-dir data/processed/graph_full_neuronal_split_binders \\
      --evidence-dir data/processed/evidence_full_v3_parkinsonism
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import drug_seeds_on_proteins, read_protein_of_gene
from mechanistic_pathway_learning.graph.drug_entry_nodes import (
    DRUG_ENTRY_FILE,
    DRUG_NODE_TYPE,
    add_drug_entry_nodes,
    extend_table_for_drug_nodes,
)
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH

SUMMARY_FILE = "drug_entry_summary.json"
REWRITTEN_FILES = ("nodes.parquet", "edges.parquet", "relation_types.json", DRUG_ENTRY_FILE, SUMMARY_FILE)
DESCRIPTOR_FILL, CELL_CLASS_WEIGHT_FILL = 0.0, 1.0


def source_and_destination(value: str) -> tuple[Path, Path]:
    source, separator, destination = value.partition(":")
    if not separator or not source or not destination:
        raise argparse.ArgumentTypeError("expected source:destination")
    return Path(source), Path(destination)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def drug_seeds_of_evidence(evidence: pd.DataFrame, node_ids: set[str], protein_of_gene: dict[str, list[str]]) -> dict[str, tuple[str, list[tuple[str, float, float]]]]:
    """perturbation id -> (drug name, the seeds the loader places for it on this graph).

    The same two steps as experiment_data.load_experiment_data: keep the triples whose node is in the graph, then on a
    split graph move each target gene's triple onto its protein nodes."""
    drug_rows = evidence[evidence.perturbation_type == "drug"].drop_duplicates("perturbation_id")
    seeds = {}
    for row in drug_rows.itertuples(index=False):
        triples = [triple for triple in json.loads(row.perturbation_nodes) if triple[0] in node_ids]
        if protein_of_gene:
            triples = drug_seeds_on_proteins(triples, protein_of_gene)
        seeds[str(row.perturbation_id)] = (str(row.perturbation_label), [(str(node_id), float(sign), float(magnitude)) for node_id, sign, magnitude in triples])
    return seeds


def graph_compound_base_ids(path: Path) -> list[frozenset[str]]:
    if not path.exists():
        return []
    table = pd.read_csv(path, dtype=str).fillna("")
    return [frozenset(base_id for base_id in cell.split(";") if base_id) for cell in table.base_metabolite_ids if cell]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=None, help="default: <graph-dir>_drugs")
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--carrier-table", default="configs/drug_plasma_carriers.csv", help="per-drug plasma carrier; none adds no sequestration edge")
    parser.add_argument("--drugs-acting-as-graph-compounds", type=Path, default=DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH)
    parser.add_argument("--descriptor-table", type=source_and_destination, action="append", default=[])
    parser.add_argument("--cell-class-weights", type=source_and_destination, action="append", default=[])
    parser.add_argument("--overwrite", action="store_true")
    arguments = parser.parse_args()

    output_directory = arguments.output_dir or arguments.graph_dir.parent / f"{arguments.graph_dir.name}_drugs"
    if (output_directory / "nodes.parquet").exists() and not arguments.overwrite:
        raise SystemExit(f"{output_directory} already holds a graph; pass --overwrite to replace it")
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet")
    relation_types = json.loads((arguments.graph_dir / "relation_types.json").read_text())
    if (nodes.node_type == DRUG_NODE_TYPE).any():
        raise SystemExit(f"{arguments.graph_dir} already has drug nodes")
    protein_of_gene = read_protein_of_gene(arguments.graph_dir)
    evidence = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet")
    drug_seeds = drug_seeds_of_evidence(evidence, set(nodes.node_id), protein_of_gene)
    carriers = None
    if str(arguments.carrier_table).lower() != "none":
        carriers = pd.read_csv(arguments.carrier_table, dtype=str).fillna("")
        carriers = carriers.assign(evidence=carriers.source)
    carrier_nodes_of_gene = None
    if protein_of_gene:  # binding is the protein's property on a split graph, as for the plasma carriage edges
        carrier_nodes_of_gene = {gene_node.removeprefix("GENE:"): proteins for gene_node, proteins in protein_of_gene.items()}

    new_nodes, new_edges, new_relations, entry_table, summary = add_drug_entry_nodes(
        nodes, edges, relation_types, drug_seeds, graph_compound_base_ids(arguments.drugs_acting_as_graph_compounds), carriers, carrier_nodes_of_gene)
    if not new_nodes.iloc[: len(nodes)].reset_index(drop=True).equals(nodes) or not new_edges.iloc[: len(edges)].reset_index(drop=True).equals(edges):
        raise SystemExit("an existing node or edge row changed; the variant must only append")

    output_directory.mkdir(parents=True, exist_ok=True)
    for source_file in sorted(arguments.graph_dir.iterdir()):  # gene_to_protein.parquet and the source's own summaries
        if source_file.is_file() and source_file.name not in REWRITTEN_FILES:
            shutil.copy2(source_file, output_directory / source_file.name)
    new_nodes.to_parquet(output_directory / "nodes.parquet")
    new_edges.to_parquet(output_directory / "edges.parquet")
    (output_directory / "relation_types.json").write_text(json.dumps(new_relations, indent=1) + "\n")
    entry_table.to_parquet(output_directory / DRUG_ENTRY_FILE, index=False)
    drug_node_ids = list(new_nodes.node_id[new_nodes.node_type == DRUG_NODE_TYPE])
    written_tables = {}
    for (source, destination), fill_value in [(pair, DESCRIPTOR_FILL) for pair in arguments.descriptor_table] + [(pair, CELL_CLASS_WEIGHT_FILL) for pair in arguments.cell_class_weights]:
        if destination.exists() and not arguments.overwrite:
            raise SystemExit(f"{destination} exists; pass --overwrite to replace it")
        table = pd.read_parquet(source)
        missing = set(nodes.node_id) - set(table.index)
        if missing:
            raise SystemExit(f"{source} has no row for {len(missing)} nodes of {arguments.graph_dir}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        extend_table_for_drug_nodes(table, drug_node_ids, fill_value).to_parquet(destination)
        written_tables[str(destination)] = {"source": str(source), "source_sha256": file_sha256(source), "sha256": file_sha256(destination), "fill_value": fill_value}
    summary.update(
        source_graph=str(arguments.graph_dir), evidence=str(arguments.evidence_dir),
        evidence_records_sha256=file_sha256(arguments.evidence_dir / "evidence_records.parquet"),
        carrier_table=None if carriers is None else str(arguments.carrier_table),
        carrier_table_sha256=None if carriers is None else file_sha256(Path(arguments.carrier_table)),
        nodes=int(len(new_nodes)), edges=int(len(new_edges)),
        source_nodes_sha256=file_sha256(arguments.graph_dir / "nodes.parquet"), source_edges_sha256=file_sha256(arguments.graph_dir / "edges.parquet"),
        nodes_sha256=file_sha256(output_directory / "nodes.parquet"), edges_sha256=file_sha256(output_directory / "edges.parquet"),
        per_node_tables=written_tables,
    )
    (output_directory / SUMMARY_FILE).write_text(json.dumps(summary, indent=1) + "\n")
    print(output_directory, json.dumps({key: value for key, value in summary.items() if not isinstance(value, list) or len(value) < 12}, indent=1))


if __name__ == "__main__":
    main()
