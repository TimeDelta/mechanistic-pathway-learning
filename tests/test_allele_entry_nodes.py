"""A neomorphic allele as a node present only where it is seeded (graph/drug_entry_nodes.add_allele_entry_nodes,
decision 9.54 (d)): appended without changing a row, its cited mechanism edges to protein or gene nodes, seeded on its own
node by the loader under drug_entry nodes and on its targets under targets, absent from every other perturbation's
response, and the per-node tables extended for it without rewriting the graph."""
from __future__ import annotations

import json
from argparse import Namespace

import numpy as np
import pandas as pd
import torch

from experiments.build_drug_entry_node_variant import extend_tables_for_entry_nodes
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evidence.allele_mechanisms import AlleleMechanismEdge, NeomorphicAllele
from mechanistic_pathway_learning.graph.drug_entry_nodes import ALLELE_NODE_TYPE, DRUG_ENTRY_FILE, add_allele_entry_nodes
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder
from tests.test_drug_entry_nodes import RELATIONS, evidence_row, toy_graph, write_graph

TOXIC = NeomorphicAllele("B|toxic_expansion", "B", "B (toxic expansion)", (
    AlleleMechanismEdge("C", "protein", -1.0, 1.0, "1", "the mutant product represses C"),
    AlleleMechanismEdge("A", "gene", 1.0, 0.5, "2", "the mutant product raises A's transcription"),
    AlleleMechanismEdge("NOPE", "protein", -1.0, 1.0, "3", "a target the graph lacks")))


def test_allele_nodes_are_appended_with_their_cited_edges_and_change_no_row():
    nodes, edges = toy_graph()
    new_nodes, new_edges, relations, entries, summary = add_allele_entry_nodes(nodes, edges, RELATIONS, {"B": TOXIC})
    pd.testing.assert_frame_equal(new_nodes.iloc[: len(nodes)].reset_index(drop=True), nodes)
    pd.testing.assert_frame_equal(new_edges.iloc[: len(edges)].reset_index(drop=True), edges)
    added = new_nodes.iloc[len(nodes):]
    assert added.node_id.tolist() == ["ALLELE:B_toxic_expansion"] and added.node_type.tolist() == [ALLELE_NODE_TYPE] and added.degree.tolist() == [2]
    mechanism = new_edges.iloc[len(edges):]
    assert sorted(zip(mechanism.target_id, mechanism.relation_type, mechanism.sign)) == [("GENE:A", "targets", 1.0), ("GENE:C", "targets", -1.0)]
    assert set(mechanism.evidence_source) == {"neomorphic allele mechanism (PMID 1)", "neomorphic allele mechanism (PMID 2)"}
    assert sorted(zip(entries.perturbation_id, entries.target_node_id, entries.magnitude)) == [("B|toxic_expansion", "GENE:A", 0.5), ("B|toxic_expansion", "GENE:C", 1.0)]
    assert summary["allele_targets_absent_from_this_graph"] == ["B|toxic_expansion -> NOPE"] and relations == RELATIONS
    split = add_allele_entry_nodes(pd.concat([nodes, pd.DataFrame({"node_id": ["PROTEIN:C"], "node_type": ["protein"]})], ignore_index=True), edges, RELATIONS,
                                   {"B": TOXIC}, protein_nodes_of_gene={"GENE:C": ["PROTEIN:C"]})[1]
    assert "PROTEIN:C" in set(split.target_id) and "GENE:C" not in set(split.target_id[split.relation_type == "targets"]), "a protein mechanism goes to the protein on a split graph"


def written_allele_graph(tmp_path):
    nodes, edges = toy_graph()
    new_nodes, new_edges, relations, entries, _ = add_allele_entry_nodes(nodes, edges, RELATIONS, {"B": TOXIC})
    write_graph(tmp_path / "graph", new_nodes, new_edges, relations, entries)
    (tmp_path / "evidence").mkdir()
    rows = [evidence_row("B", "gene", "cluster:B", [["GENE:B", -1.0, 1.0]]),
            {**evidence_row("B|toxic_expansion", "gene", "cluster:B", [["GENE:C", -1.0, 1.0], ["GENE:A", 1.0, 0.5]]), "group_id": "B"},
            evidence_row("C", "gene", "cluster:C", [["GENE:C", -1.0, 1.0]])]
    pd.DataFrame(rows).to_parquet(tmp_path / "evidence/evidence_records.parquet", index=False)
    return tmp_path / "graph", tmp_path / "evidence"


def test_the_loader_seeds_an_allele_on_its_node_or_on_its_targets(tmp_path):
    graph, evidence = written_allele_graph(tmp_path)
    as_nodes = load_experiment_data(graph, evidence, drug_entry="nodes")
    seeds = {perturbation: [as_nodes.node_ids[index] for index in indices] for perturbation, indices in zip(as_nodes.perturbation_ids, as_nodes.perturbation_seeds)}
    assert seeds["B|toxic_expansion"] == ["ALLELE:B_toxic_expansion"] and seeds["B"] == ["GENE:B"], "the loss-of-function allele keeps the gene's own seed"
    assert [as_nodes.node_ids[index] for index in np.flatnonzero(as_nodes.entry_node_mask)] == ["ALLELE:B_toxic_expansion"]
    mechanism = np.asarray(as_nodes.edge_relation) == as_nodes.relation_types.index("targets")
    assert sorted(as_nodes.edge_weight[mechanism].tolist()) == [0.5, 1.0], "an allele's mechanism edge carries its stated magnitude"
    assert as_nodes.group_ids[as_nodes.perturbation_ids.index("B|toxic_expansion")] == as_nodes.group_ids[as_nodes.perturbation_ids.index("B")] == "B"
    as_targets = load_experiment_data(graph, evidence, drug_entry="targets")
    position = as_targets.perturbation_ids.index("B|toxic_expansion")
    assert [as_targets.node_ids[index] for index in as_targets.perturbation_seeds[position]] == ["GENE:C", "GENE:A"]
    assert as_targets.perturbation_signs[position].tolist() == [-1.0, 1.0] and as_targets.entry_node_mask is None


def test_a_knockouts_response_does_not_depend_on_the_allele_node(tmp_path):
    graph, evidence = written_allele_graph(tmp_path)
    data = load_experiment_data(graph, evidence, drug_entry="nodes")
    relation_types = data.relation_types
    tensors = (torch.as_tensor(data.edge_source), torch.as_tensor(data.edge_target), torch.as_tensor(data.edge_relation), torch.as_tensor(data.edge_sign, dtype=torch.float32))
    torch.manual_seed(0)
    encoder = LinearResponseEncoder(len(data.node_ids), relation_types, *tensors, torch.eye(len(data.node_ids)), 4, num_propagation_steps=20, fixed_point_tolerance=None,
                                    entry_node_mask=torch.as_tensor(data.entry_node_mask), entry_relation_names=("targets",), edge_weight=torch.as_tensor(data.edge_weight))
    knockout = encoder.response(torch.tensor([[data.node_index["GENE:C"]]]), torch.tensor([[[-1.0, 1.0]]]))
    allele_row = data.node_index["ALLELE:B_toxic_expansion"]
    assert float(knockout[0, allele_row].abs().max()) == 0.0, "an allele node that is not seeded holds nothing"
    allele = encoder.response(torch.tensor([[allele_row]]), torch.tensor([[[1.0, 1.0]]]))
    assert float(allele[0, data.node_index["GENE:C"]].abs().max()) > 0, "seeded, the allele acts on its targets"


def test_the_tables_are_extended_for_entry_nodes_without_writing_the_graph(tmp_path):
    graph, _ = written_allele_graph(tmp_path)
    source = tmp_path / "source"
    nodes, edges = toy_graph()
    write_graph(source, nodes, edges, RELATIONS)
    descriptors = pd.DataFrame({"x": np.arange(len(nodes), dtype=float)}, index=pd.Index(nodes.node_id, name="node_id"))
    descriptors.to_parquet(tmp_path / "descriptors.parquet")
    before = (graph / "nodes.parquet").read_bytes()
    extend_tables_for_entry_nodes(Namespace(graph_dir=source, output_dir=graph, descriptor_table=[(tmp_path / "descriptors.parquet", tmp_path / "extended.parquet")],
                                            cell_class_weights=[], overwrite=False))
    extended = pd.read_parquet(tmp_path / "extended.parquet")
    assert extended.loc["ALLELE:B_toxic_expansion", "x"] == 0.0 and len(extended) == len(nodes) + 1
    assert (graph / "nodes.parquet").read_bytes() == before and json.loads((graph / "relation_types.json").read_text()) == RELATIONS
    assert (graph / DRUG_ENTRY_FILE).exists()


def test_tables_only_records_each_table_in_the_variants_summary(tmp_path):
    graph, _ = written_allele_graph(tmp_path)
    (graph / "drug_entry_summary.json").write_text(json.dumps({"nodes": 4, "per_node_tables": {"earlier.parquet": {"sha256": "kept"}}}))
    source = tmp_path / "source"
    nodes, edges = toy_graph()
    write_graph(source, nodes, edges, RELATIONS)
    weights = pd.DataFrame({"neuron": np.full(len(nodes), 0.5)}, index=pd.Index(nodes.node_id, name="node_id"))
    weights.to_parquet(tmp_path / "weights.parquet")
    extend_tables_for_entry_nodes(Namespace(graph_dir=source, output_dir=graph, descriptor_table=[], cell_class_weights=[(tmp_path / "weights.parquet", tmp_path / "weights_extended.parquet")],
                                            overwrite=False))
    summary = json.loads((graph / "drug_entry_summary.json").read_text())
    recorded = summary["per_node_tables"][str(tmp_path / "weights_extended.parquet")]
    assert recorded["fill_value"] == 1.0 and recorded["written_by"] == "--tables-only" and len(recorded["sha256"]) == 64
    assert summary["per_node_tables"]["earlier.parquet"] == {"sha256": "kept"} and summary["nodes"] == 4, "the rest of the summary is kept"
    assert pd.read_parquet(tmp_path / "weights_extended.parquet").loc["ALLELE:B_toxic_expansion", "neuron"] == 1.0
