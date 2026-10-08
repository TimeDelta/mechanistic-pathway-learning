"""Tests: the gene and protein split (graph/gene_protein_split.py, docs/gene_protein_split.md) reroutes every relation
onto the protein node except the target end of transcription, gives genes of one shared UniProt entry one protein node,
moves the protein and gene brain blocks onto the protein nodes, makes the loader seed a drug's proteins and a knockout's
gene with the leakage groups of the merged graph, and leaves encodes out of the rewiring."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from experiments.run_main_model import rewiring_fixed_relations
from mechanistic_pathway_learning.evaluation.experiment_data import GENE_TO_PROTEIN_FILE, load_experiment_data
from mechanistic_pathway_learning.evaluation.negative_controls import fast_degree_preserving_rewiring
from mechanistic_pathway_learning.graph.gene_protein_split import (
    ENCODES_RELATION, PROTEIN_NODE_TYPE, cell_class_weights_for_split, move_gene_blocks_to_proteins, protein_node_assignment, split_graph)

RELATIONS = ["catalyzed_by", "activates", "regulates_transcription_of", "member_of"]


def toy_graph() -> tuple[pd.DataFrame, pd.DataFrame]:
    """TF regulates B; TF activates C; B and C (one shared entry) both activate D; D catalyses R1; D is in complex X."""
    nodes = pd.DataFrame({
        "node_id": ["GENE:TF", "GENE:B", "GENE:C", "GENE:D", "R1", "COMPLEX:X"],
        "node_type": ["gene", "gene", "gene", "gene", "reaction", "protein_entity"],
        "gene_symbol": ["TF", "B", "C", "D", None, None],
        "ensembl_gene_id": [None, "ENSG_B", None, "ENSG_D", None, None],
        "is_currency": [False] * 6, "degree": [0] * 6,
        "brain_median_tpm_max": [5.0, 2.0, 7.0, 1.0, 3.0, None],
        "brain_expressed": [True, True, True, False, True, None]})
    edges = pd.DataFrame({
        "source_id": ["GENE:TF", "GENE:TF", "GENE:B", "GENE:C", "GENE:D", "GENE:D"],
        "target_id": ["GENE:B", "GENE:C", "GENE:D", "GENE:D", "R1", "COMPLEX:X"],
        "relation_type": ["regulates_transcription_of", "activates", "activates", "activates", "catalyzed_by", "member_of"],
        "sign": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]})
    return nodes, edges


SHARED_ENTRY = pd.Series(["B; C", "TF", "D", "OTHER"])


def test_split_reroutes_relations_and_adds_one_encodes_edge_per_gene() -> None:
    nodes, edges = toy_graph()
    assignment = protein_node_assignment(nodes, SHARED_ENTRY)
    split_nodes, split_edges, relations, summary = split_graph(nodes, edges, RELATIONS, assignment)
    assert relations == RELATIONS + [ENCODES_RELATION]
    assert dict(zip(assignment.gene_node_id, assignment.protein_node_id)) == {
        "GENE:B": "PROTEIN:B", "GENE:C": "PROTEIN:B", "GENE:D": "PROTEIN:D", "GENE:TF": "PROTEIN:TF"}
    triples = set(zip(split_edges.source_id, split_edges.target_id, split_edges.relation_type))
    assert ("PROTEIN:TF", "GENE:B", "regulates_transcription_of") in triples  # the regulator's protein acts on the target's gene
    assert ("PROTEIN:TF", "PROTEIN:B", "activates") in triples
    assert ("PROTEIN:B", "PROTEIN:D", "activates") in triples  # B -> D and C -> D become one edge of the shared node
    assert ("PROTEIN:D", "R1", "catalyzed_by") in triples and ("PROTEIN:D", "COMPLEX:X", "member_of") in triples
    encodes = split_edges[split_edges.relation_type == ENCODES_RELATION]
    assert sorted(zip(encodes.source_id, encodes.target_id)) == [("GENE:B", "PROTEIN:B"), ("GENE:C", "PROTEIN:B"), ("GENE:D", "PROTEIN:D"), ("GENE:TF", "PROTEIN:TF")]
    assert summary["edges_merged_onto_a_shared_protein"] == 1 and summary["shared_protein_nodes"] == 1
    assert set(split_edges[split_edges.source_id.str.startswith("GENE:")].relation_type) == {ENCODES_RELATION}  # a gene acts only through its protein
    assert set(split_edges[split_edges.target_id.str.startswith("GENE:")].relation_type) == {"regulates_transcription_of"}
    degree = split_nodes.set_index("node_id").degree
    assert degree["GENE:B"] == 2 and degree["GENE:C"] == 1 and degree["GENE:TF"] == 1  # transcription in-edge plus encodes; encodes only; encodes only
    assert degree["PROTEIN:B"] == 4  # two encodes, the activation from TF, the activation of D
    protein = split_nodes.set_index("node_id").loc["PROTEIN:B"]
    assert protein.node_type == PROTEIN_NODE_TYPE and protein.display_name == "B/C" and protein.ensembl_gene_id == "ENSG_B"
    assert protein.brain_median_tpm_max == 7.0 and bool(protein.brain_expressed)  # the larger of its genes
    assert pd.isna(split_nodes.set_index("node_id").loc["GENE:B"].brain_median_tpm_max)  # moved off the gene
    assert split_nodes.set_index("node_id").loc["R1"].brain_median_tpm_max == 3.0  # other node types keep theirs


def test_descriptor_blocks_move_to_protein_nodes_and_shared_nodes_take_their_entry_row() -> None:
    nodes, edges = toy_graph()
    assignment = protein_node_assignment(nodes, SHARED_ENTRY)
    split_nodes, _, _, _ = split_graph(nodes, edges, RELATIONS, assignment)
    table = pd.DataFrame({"protein_rrr_1": [0.5, 0.0, 0.0, -1.0, 0.0, 0.0], "protein_has_protein_descriptors": [1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                          "gene_brain_region_x": [0.3, 1.2, -0.4, 0.0, 0.0, 0.0], "reaction_brain_region_x": [0.0, 0.0, 0.0, 0.0, 0.8, 0.0]},
                         index=pd.Index(nodes.node_id, name="node_id"))
    entry_descriptors = pd.DataFrame({"rrr_1": [2.5]}, index=pd.Index(["B; C"], name="gene_symbol"))
    moved = move_gene_blocks_to_proteins(table, split_nodes, assignment, entry_descriptors, SHARED_ENTRY)
    assert list(moved.index) == list(split_nodes.node_id)
    assert (moved.loc[["GENE:TF", "GENE:B", "GENE:C", "GENE:D"], ["protein_rrr_1", "protein_has_protein_descriptors", "gene_brain_region_x"]] == 0).all().all()
    assert moved.loc["PROTEIN:TF", "protein_rrr_1"] == 0.5 and moved.loc["PROTEIN:D", "protein_rrr_1"] == -1.0
    assert moved.loc["PROTEIN:B", "protein_rrr_1"] == 2.5 and moved.loc["PROTEIN:B", "protein_has_protein_descriptors"] == 1.0  # the shared entry's own row
    assert moved.loc["PROTEIN:B", "gene_brain_region_x"] == 1.2  # the larger of B and C
    assert moved.loc["R1", "reaction_brain_region_x"] == 0.8 and moved.loc["PROTEIN:B", "reaction_brain_region_x"] == 0.0
    weights = pd.DataFrame({"neuron": [0.9, 0.2, 0.6, 0.1, 0.5, 1.0]}, index=pd.Index(nodes.node_id, name="node_id"))
    split_weights = cell_class_weights_for_split(weights, split_nodes, assignment)
    assert split_weights.loc["GENE:B", "neuron"] == 0.2 and split_weights.loc["PROTEIN:B", "neuron"] == 0.6 and split_weights.loc["PROTEIN:TF", "neuron"] == 0.9
    assert not split_weights.isna().any().any()


def write_graph(directory: Path, nodes: pd.DataFrame, edges: pd.DataFrame, relations: list[str], assignment: pd.DataFrame | None = None) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    nodes.to_parquet(directory / "nodes.parquet", index=False)
    edges.to_parquet(directory / "edges.parquet", index=False)
    (directory / "relation_types.json").write_text(json.dumps(relations))
    if assignment is not None:
        assignment.to_parquet(directory / GENE_TO_PROTEIN_FILE, index=False)


def evidence_row(perturbation: str, kind: str, cluster: str, seeds: list) -> dict:
    return {"perturbation_id": perturbation, "perturbation_type": kind, "perturbation_label": perturbation, "group_id": perturbation,
            "disease_cluster_id": cluster, "symptom": "anxiety", "relation": "induces", "evidence_class": "monogenic", "grade": "A",
            "weight": 1.0, "label_frequency": None, "in_metabolic_layer": True, "evidence_date": None, "perturbation_nodes": json.dumps(seeds)}


def test_loader_seeds_a_drug_on_proteins_and_a_knockout_on_its_gene_with_the_merged_graph_groups(tmp_path: Path) -> None:
    nodes, edges = toy_graph()
    assignment = protein_node_assignment(nodes, SHARED_ENTRY)
    split_nodes, split_edges, relations, _ = split_graph(nodes, edges, RELATIONS, assignment)
    write_graph(tmp_path / "merged", nodes, edges, RELATIONS)
    write_graph(tmp_path / "split", split_nodes, split_edges, relations, assignment)
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    pd.DataFrame([evidence_row("knockout_D", "gene", "cluster:D", [["GENE:D", -1.0, 1.0]]),
                  evidence_row("knockout_TF", "gene", "cluster:TF", [["GENE:TF", -1.0, 1.0]]),
                  evidence_row("drug_BC", "drug", "drug_BC", [["GENE:B", -1.0, 0.5], ["GENE:C", -1.0, 0.5], ["GENE:D", -1.0, 0.5]])]
                 ).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    merged = load_experiment_data(tmp_path / "merged", evidence_directory, group_by="disease_cluster_and_targets")
    split = load_experiment_data(tmp_path / "split", evidence_directory, group_by="disease_cluster_and_targets")
    seeds = {p: [split.node_ids[i] for i in s] for p, s in zip(split.perturbation_ids, split.perturbation_seeds)}
    assert seeds == {"drug_BC": ["PROTEIN:B", "PROTEIN:D"], "knockout_D": ["GENE:D"], "knockout_TF": ["GENE:TF"]}  # B and C share one protein, seeded once
    assert split.perturbation_magnitudes[split.perturbation_ids.index("drug_BC")].tolist() == [0.5, 0.5]
    assert split.group_ids == merged.group_ids  # the drug joins the knockout of its target D on both graphs
    assert split.group_ids[split.perturbation_ids.index("drug_BC")] == split.group_ids[split.perturbation_ids.index("knockout_D")]
    degree = dict(zip(split.perturbation_ids, split.perturbation_degrees))
    assert degree["knockout_TF"] == 1.0  # the split degree of the seed itself (the user's decision)


def test_rewiring_leaves_encodes_fixed_and_the_other_relations_as_without_it() -> None:
    generator = np.random.default_rng(3)
    activation_edges = generator.choice(40, size=(2, 120))
    activation_edges = activation_edges[:, activation_edges[0] != activation_edges[1]]
    encodes_edges = np.stack([np.arange(40, 60), np.arange(60, 80)])
    edges = np.concatenate([activation_edges, encodes_edges], axis=1)
    relations = np.array([0] * activation_edges.shape[1] + [1] * encodes_edges.shape[1])
    fixed = fast_degree_preserving_rewiring(edges, relations, num_swaps_per_edge=10, random_seed=7, fixed_relations=[1])
    np.testing.assert_array_equal(fixed[:, relations == 1], encodes_edges)
    free = fast_degree_preserving_rewiring(edges, relations, num_swaps_per_edge=10, random_seed=7)
    np.testing.assert_array_equal(fixed[:, relations == 0], free[:, relations == 0])  # encodes comes last, so the stream before it is the same
    assert not np.array_equal(free[:, relations == 1], encodes_edges)
    assert rewiring_fixed_relations(RELATIONS, False) == []
    assert rewiring_fixed_relations(RELATIONS + [ENCODES_RELATION], False) == [len(RELATIONS)]
    assert rewiring_fixed_relations(RELATIONS + [ENCODES_RELATION], True) == []
