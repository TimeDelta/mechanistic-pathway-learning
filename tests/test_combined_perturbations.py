"""Tests: the combined perturbation type (evidence/combined_perturbations.py), which rows 4 to 6 of section 8 of
docs/label_source_review.md named and the owner asked for on 12 October 2026 ("Build the combined perturbation type").

What is checked here is what a combination has to get right to be safe to add after a held-out set is drawn:

  - a row of one perturbation reads as one part, so nothing changes for a table without combinations;
  - the parts of a combination are routed per part, which is the whole reason the parts are kept: on a split graph the
    drug part moves to proteins and the gene part does not, and on a graph with drug entry nodes the drug part moves to
    the drug's own node and the gene part stays where it is;
  - a node two parts seed carries the sum of their signed magnitudes;
  - a combination is in one leakage group with every part of it (section 8: "a combination inherits the groups of its
    parts");
  - a drug with a combination does not take two drugs' share of the loss under --equal-drug-shares;
  - the input-size controls read a combination as a combination and not as a drug with one more target.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from experiments.run_main_model import drug_input_normalised, drug_loss_slots, permuted_drug_targets
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.input_size_strata import COMBINATION_STRATUM, input_size_strata
from mechanistic_pathway_learning.evidence.combined_perturbations import (
    COMBINATION_PERTURBATION_TYPE,
    PARTS_COLUMN,
    PerturbationPart,
    combined_perturbation_id,
    nodes_json_of_parts,
    parts_json,
    parts_with_disjoint_nodes,
    read_parts,
)
from mechanistic_pathway_learning.graph.drug_entry_nodes import add_drug_entry_nodes
from mechanistic_pathway_learning.models.soft_constraint_losses import equal_drug_loss_shares
from tests.test_drug_entry_nodes import DRUG_SEEDS, GRAPH_COMPOUNDS, RELATIONS, evidence_row, toy_graph, write_graph

CLEARING_ENZYMES = pd.DataFrame({"perturbation_id": ["drug_one"], "enzyme_gene": ["C"], "evidence": ["table"]})


def combination_of_a_drug_and_a_lowered_gene(drug: str = "drug_one", gene_node: str = "GENE:C") -> list[PerturbationPart]:
    _, triples = DRUG_SEEDS[drug]
    return [PerturbationPart(drug, "drug", tuple((node_id, sign, magnitude) for node_id, sign, magnitude in triples)),
            PerturbationPart("C", "gene", ((gene_node, -1.0, 1.0),))]


def combination_row(parts: list[PerturbationPart], symptom: str = "anxiety") -> dict:
    row = evidence_row(combined_perturbation_id([part.part_id for part in parts]), COMBINATION_PERTURBATION_TYPE, "cluster:C", [])
    row.update({"symptom": symptom, "perturbation_nodes": nodes_json_of_parts(parts), PARTS_COLUMN: parts_json(parts),
                "group_id": "drug_one", "disease_cluster_id": "drug_one"})
    return row


def test_a_single_perturbation_reads_as_one_part_of_its_own_kind() -> None:
    parts = read_parts("ABAT", "gene", json.dumps([["GENE:ABAT", -1.0, 1.0]]))
    assert parts == [PerturbationPart("ABAT", "gene", (("GENE:ABAT", -1.0, 1.0),))]
    assert read_parts("CID1", "drug", "[]") == [PerturbationPart("CID1", "drug", ())]


def test_a_combination_round_trips_and_is_refused_without_its_parts() -> None:
    parts = combination_of_a_drug_and_a_lowered_gene()
    perturbation_id = combined_perturbation_id([part.part_id for part in parts])
    assert perturbation_id == "drug_one+C"
    assert read_parts(perturbation_id, COMBINATION_PERTURBATION_TYPE, nodes_json_of_parts(parts), parts_json(parts)) == parts
    with pytest.raises(ValueError, match="cannot be read from the union"):
        read_parts(perturbation_id, COMBINATION_PERTURBATION_TYPE, nodes_json_of_parts(parts))
    with pytest.raises(ValueError, match="a combination's id is its part ids joined"):
        read_parts("another_name", COMBINATION_PERTURBATION_TYPE, nodes_json_of_parts(parts), parts_json(parts))
    with pytest.raises(ValueError, match="at least two parts"):
        combined_perturbation_id(["drug_one"])
    with pytest.raises(ValueError, match="is not one of"):
        PerturbationPart("drug_one", "compound", ())


def test_a_node_two_parts_seed_carries_the_sum_of_their_signed_magnitudes() -> None:
    both_inhibit = [PerturbationPart("first", "drug", (("GENE:B", -1.0, 0.5), ("GENE:A", -1.0, 1.0))),
                    PerturbationPart("second", "drug", (("GENE:B", -1.0, 0.25),))]
    disjoint, shared = parts_with_disjoint_nodes(both_inhibit)
    assert shared == ["GENE:B"]
    assert disjoint[0].nodes == (("GENE:B", -1.0, 0.75), ("GENE:A", -1.0, 1.0))  # held by the first part that named it
    assert disjoint[1].nodes == ()
    opposite = [PerturbationPart("up", "drug", (("GENE:B", 1.0, 1.0),)), PerturbationPart("down", "gene", (("GENE:B", -1.0, 0.25),))]
    assert parts_with_disjoint_nodes(opposite)[0][0].nodes == (("GENE:B", 1.0, 0.75),)
    cancelling = [PerturbationPart("up", "drug", (("GENE:B", 1.0, 1.0),)), PerturbationPart("down", "gene", (("GENE:B", -1.0, 1.0),))]
    assert parts_with_disjoint_nodes(cancelling)[0][0].nodes == (("GENE:B", 1.0, 0.0),)  # an input of nothing, not of unknown direction


def written_graphs_and_evidence(tmp_path: Path, parts: list[PerturbationPart]) -> tuple[Path, Path, Path, Path]:
    """A toy graph with and without drug entry nodes, an evidence table of single perturbations, and a combination table."""
    nodes, edges = toy_graph()
    new_nodes, new_edges, relations, entry_table, _ = add_drug_entry_nodes(nodes, edges, RELATIONS, DRUG_SEEDS, GRAPH_COMPOUNDS, None, None, CLEARING_ENZYMES)
    write_graph(tmp_path / "source", nodes, edges, RELATIONS)
    write_graph(tmp_path / "drugs", new_nodes, new_edges, relations, entry_table)
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    rows = [evidence_row("knockout_C", "gene", "cluster:C", [["GENE:C", -1.0, 1.0]])]
    rows += [evidence_row(perturbation, "drug", perturbation, [list(triple) for triple in triples]) for perturbation, (_, triples) in DRUG_SEEDS.items()]
    pd.DataFrame(rows).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    combined = tmp_path / "combined.parquet"
    pd.DataFrame([combination_row(parts)]).to_parquet(combined, index=False)
    return tmp_path / "source", tmp_path / "drugs", evidence_directory, combined


def test_the_drug_part_is_seeded_on_the_drug_node_and_the_gene_part_stays_where_it_is(tmp_path: Path) -> None:
    parts = combination_of_a_drug_and_a_lowered_gene()
    source, drugs, evidence_directory, combined = written_graphs_and_evidence(tmp_path, parts)
    without = load_experiment_data(drugs, evidence_directory, drug_entry="nodes", group_by="disease_cluster_and_targets")
    data = load_experiment_data(drugs, evidence_directory, drug_entry="nodes", group_by="disease_cluster_and_targets", combined_perturbations=combined)
    assert without.perturbation_parts is None and data.perturbation_parts is not None
    assert len(data.perturbation_ids) == len(without.perturbation_ids) + 1
    position = data.perturbation_ids.index("drug_one+C")
    seeded = [data.node_ids[node] for node in data.perturbation_seeds[position]]
    assert seeded == ["DRUG:drug_one", "GENE:C"]  # the drug's own node, then the lowered gene's
    np.testing.assert_array_equal(data.perturbation_signs[position], [1.0, -1.0])
    np.testing.assert_array_equal(data.perturbation_seed_parts[position], [0, 1])
    # without the entry nodes the drug part is its targets and the gene part is unchanged
    as_targets = load_experiment_data(drugs, evidence_directory, drug_entry="targets", group_by="disease_cluster_and_targets", combined_perturbations=combined)
    position = as_targets.perturbation_ids.index("drug_one+C")
    assert [as_targets.node_ids[node] for node in as_targets.perturbation_seeds[position]] == ["GENE:A", "GENE:B", "GENE:C"]
    np.testing.assert_array_equal(as_targets.perturbation_seed_parts[position], [0, 0, 1])
    assert load_experiment_data(source, evidence_directory, group_by="disease_cluster_and_targets", combined_perturbations=combined).perturbation_ids == as_targets.perturbation_ids


def test_a_combination_is_in_one_leakage_group_with_every_part_of_it(tmp_path: Path) -> None:
    parts = combination_of_a_drug_and_a_lowered_gene()
    _, drugs, evidence_directory, combined = written_graphs_and_evidence(tmp_path, parts)
    data = load_experiment_data(drugs, evidence_directory, drug_entry="nodes", group_by="disease_cluster_and_targets", combined_perturbations=combined)
    group_of = dict(zip(data.perturbation_ids, data.group_ids))
    # GENE:C is the gene knockout_C perturbs and GENE:B is a target of drug_one and of drug_two, so one group holds them all
    assert group_of["drug_one+C"] == group_of["drug_one"] == group_of["knockout_C"] == group_of["drug_two"]
    assert len(set(data.group_ids)) < len(data.perturbation_ids)


def test_a_combination_is_refused_when_its_id_is_already_a_perturbation(tmp_path: Path) -> None:
    parts = combination_of_a_drug_and_a_lowered_gene()
    _, drugs, evidence_directory, combined = written_graphs_and_evidence(tmp_path, parts)
    clashing = pd.read_parquet(combined).assign(perturbation_id="drug_one", perturbation_type=COMBINATION_PERTURBATION_TYPE)
    clashing.to_parquet(tmp_path / "clashing.parquet", index=False)
    with pytest.raises(ValueError, match="already in the evidence table"):
        load_experiment_data(drugs, evidence_directory, drug_entry="nodes", combined_perturbations=tmp_path / "clashing.parquet")
    not_a_combination = pd.read_parquet(combined).assign(perturbation_type="drug")
    not_a_combination.to_parquet(tmp_path / "plain.parquet", index=False)
    with pytest.raises(ValueError, match="every row of it is a"):
        load_experiment_data(drugs, evidence_directory, drug_entry="nodes", combined_perturbations=tmp_path / "plain.parquet")


def test_a_drug_with_a_combination_does_not_take_two_drugs_share_of_the_loss() -> None:
    # three drugs with 2, 2 and 2 positive pairs; the third also appears in a combination carrying its own two pairs
    pair_weights = np.array([[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 0.0]])
    is_positive = pair_weights > 0
    is_drug = np.array([True, True, True, True])
    fit_indices = np.arange(4)
    slots = np.array(["drug_a", "drug_b", "drug_c", "drug_c"], dtype=object)
    on_its_own = equal_drug_loss_shares(pair_weights, is_positive, is_drug, fit_indices)
    assert on_its_own[2, 0] == pytest.approx(on_its_own[3, 0])  # four separate drugs: the combination is a drug of its own
    shared = equal_drug_loss_shares(pair_weights, is_positive, is_drug, fit_indices, slots)
    total_of = {"drug_a": shared[0].sum(), "drug_c": shared[2].sum() + shared[3].sum()}
    assert total_of["drug_c"] == pytest.approx(total_of["drug_a"])  # one slot, one share, split between the two rows
    np.testing.assert_allclose(shared[0], shared[1])


def test_the_loss_slot_of_a_combination_is_its_drug(tmp_path: Path) -> None:
    parts = combination_of_a_drug_and_a_lowered_gene()
    _, drugs, evidence_directory, combined = written_graphs_and_evidence(tmp_path, parts)
    data = load_experiment_data(drugs, evidence_directory, drug_entry="nodes", combined_perturbations=combined)
    slots = dict(zip(data.perturbation_ids, drug_loss_slots(data)))
    assert slots["drug_one+C"] == "drug_one" == slots["drug_one"]
    assert slots["knockout_C"] == "knockout_C"
    plain = load_experiment_data(drugs, evidence_directory, drug_entry="nodes")
    assert drug_loss_slots(plain) is None  # no combination: every row is its own drug, the rule the shares were registered with


def test_the_input_size_controls_read_a_combination_as_one(tmp_path: Path) -> None:
    labels, order = input_size_strata(["gene", "drug", COMBINATION_PERTURBATION_TYPE], np.array([1, 2, 3]))
    assert labels[2] == COMBINATION_STRATUM and order[-1] == COMBINATION_STRATUM
    parts = combination_of_a_drug_and_a_lowered_gene()
    source, _, evidence_directory, combined = written_graphs_and_evidence(tmp_path, parts)
    data = load_experiment_data(source, evidence_directory, combined_perturbations=combined)
    position = data.perturbation_ids.index("drug_one+C")
    # the drug part is rescaled to sum to 1 among its own seeds and the lowered gene keeps its magnitude
    normalised = drug_input_normalised(data.perturbation_types, data.perturbation_magnitudes, data.perturbation_parts, data.perturbation_seed_parts)
    np.testing.assert_allclose(normalised[position], [0.5, 0.5, 1.0])
    # the permuted-target control leaves a combination alone
    seeds, _, _, moved = permuted_drug_targets(data.perturbation_types, data.perturbation_seeds, data.perturbation_signs, data.perturbation_magnitudes, random_seed=0)
    np.testing.assert_array_equal(seeds[position], data.perturbation_seeds[position])
    assert moved == sum(kind == "drug" for kind in data.perturbation_types) - 1  # the one-target drug cannot be exchanged with a two-target drug


def test_the_restriction_to_a_development_pool_keeps_the_parts(tmp_path: Path) -> None:
    from mechanistic_pathway_learning.evaluation.experiment_data import restrict_to_perturbations

    parts = combination_of_a_drug_and_a_lowered_gene()
    _, drugs, evidence_directory, combined = written_graphs_and_evidence(tmp_path, parts)
    data = load_experiment_data(drugs, evidence_directory, drug_entry="nodes", combined_perturbations=combined)
    keep = np.array([perturbation_id != "drug_two" for perturbation_id in data.perturbation_ids])
    smaller = restrict_to_perturbations(data, keep)
    position = smaller.perturbation_ids.index("drug_one+C")
    assert [part.part_id for part in smaller.perturbation_parts[position]] == ["drug_one", "C"]
    np.testing.assert_array_equal(smaller.perturbation_seed_parts[position], [0, 1])
    assert len(smaller.perturbation_parts) == len(smaller.perturbation_ids) == len(data.perturbation_ids) - 1


def test_the_builder_writes_a_row_per_kept_positive_symptom_and_skips_a_negative_control(tmp_path: Path) -> None:
    from experiments.build_combined_perturbations import parts_of_check
    from mechanistic_pathway_learning.evaluation.comparative_checks import ComparativeCheck

    drug_nodes = {"drug_one": json.dumps([list(triple) for triple in DRUG_SEEDS["drug_one"][1]])}
    node_ids = {"GENE:A", "GENE:B", "GENE:C"}
    check = ComparativeCheck("enzyme_drug_one_C", "enzyme_activity", "drug_one", lowered_genes=("C",))
    parts, reason = parts_of_check(check, drug_nodes, {}, node_ids)
    assert not reason and [(part.part_id, part.part_kind) for part in parts] == [("drug_one", "drug"), ("C", "gene")]
    assert parts[1].nodes == (("GENE:C", -1.0, 1.0),)
    missing, reason = parts_of_check(ComparativeCheck("x", "enzyme_activity", "drug_one", lowered_genes=("NOT_A_GENE",)), drug_nodes, {}, node_ids)
    assert not missing and reason == "no node for NOT_A_GENE"
    unknown, reason = parts_of_check(ComparativeCheck("y", "enzyme_activity", "drug_nine", lowered_genes=("C",)), drug_nodes, {}, node_ids)
    assert not unknown and "not a drug of the evidence table" in reason
    # a gene part is seeded on the protein nodes of a split graph, as a check seeds it
    split, _ = parts_of_check(check, drug_nodes, {"GENE:C": ["PROTEIN:C1", "PROTEIN:C2"]}, {*node_ids, "PROTEIN:C1", "PROTEIN:C2"})
    assert split[1].nodes == (("PROTEIN:C1", -1.0, 0.5), ("PROTEIN:C2", -1.0, 0.5))
