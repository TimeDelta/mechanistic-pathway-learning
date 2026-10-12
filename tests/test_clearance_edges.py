"""Tests: clearance edges (graph/drug_entry_nodes.py, relation clears; the user, 11 October 2026: "Add the enzyme labels
as the carrier twin instead then"). An enzyme's edge to the drug it clears is the carrier's edge with another relation
name. The builder appends it after every other row and leaves the node table alone, so the graph without it is the
graph it was; the loader's "targets" mode still gives the source graph back; and the message-passing encoder carries
the enzyme's change as it carries the carrier's, each relation through its own weights, silent while both are at rest.
The table of clearing enzymes takes a row only when the regulator's sentence gives the direction."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import torch

from experiments.build_drug_clearing_enzymes import edge_decision
from experiments.run_main_model import rewiring_fixed_relations, sequestration_inputs
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.graph.drug_entry_nodes import (
    CHANGE_CARRYING_RELATIONS,
    CLEARS_RELATION,
    DRUG_MECHANISM_RELATION,
    DRUG_NODE_TYPE,
    SEQUESTERS_RELATION,
    add_drug_entry_nodes,
)
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder
from tests.test_drug_entry_nodes import (
    CARRIERS_OF_DRUG_NODES,
    DRUG_SEEDS,
    GRAPH_COMPOUNDS,
    RELATIONS,
    evidence_row,
    toy_graph,
    typed_features,
    write_graph,
)

# GENE:C stands for the enzyme: it clears drug_one (which also has a carrier) and drug_two (which has none)
CLEARING_ENZYMES = pd.DataFrame({"perturbation_id": ["drug_one", "drug_two"], "enzyme_gene": ["C", "C"], "evidence": ["table", "table"]})


def variant_with_enzymes(clearing_enzymes: pd.DataFrame | None = CLEARING_ENZYMES):
    nodes, edges = toy_graph()
    return (nodes, edges, *add_drug_entry_nodes(nodes, edges, RELATIONS, DRUG_SEEDS, GRAPH_COMPOUNDS, CARRIERS_OF_DRUG_NODES, None, clearing_enzymes))


def test_clearance_edges_are_appended_last_and_leave_the_graph_without_them_as_it_was() -> None:
    _, _, plain_nodes, plain_edges, plain_relations, plain_entry, plain_summary = variant_with_enzymes(None)
    _, _, new_nodes, new_edges, relations, entry_table, summary = variant_with_enzymes()
    pd.testing.assert_frame_equal(new_nodes, plain_nodes)  # a drug node's degree does not count its clearance edges
    pd.testing.assert_frame_equal(new_edges.iloc[: len(plain_edges)].reset_index(drop=True), plain_edges)
    pd.testing.assert_frame_equal(entry_table, plain_entry)
    clearance = new_edges.iloc[len(plain_edges):]
    assert sorted(zip(clearance.source_id, clearance.target_id, clearance.relation_type, clearance.sign)) == [
        ("GENE:C", "DRUG:drug_one", CLEARS_RELATION, -1.0), ("GENE:C", "DRUG:drug_two", CLEARS_RELATION, -1.0)]
    assert clearance.evidence_source.tolist() == ["enzyme that clears the drug (C; table)"] * 2
    assert relations == plain_relations + [CLEARS_RELATION] and plain_relations[-1] == SEQUESTERS_RELATION
    assert summary["clearance_edges_added"] == 2 and summary["drugs_given_a_clearance_edge"] == 2 and summary["enzyme_genes_without_a_node"] == []
    assert not any("clearance" in key or "enzyme" in key for key in plain_summary)  # a build without the table has the summary it always had
    absent = add_drug_entry_nodes(*toy_graph(), RELATIONS, DRUG_SEEDS, GRAPH_COMPOUNDS, None, None,
                                  pd.DataFrame({"perturbation_id": ["drug_one"], "enzyme_gene": ["NOT_IN_THE_GRAPH"]}))[-1]
    assert absent["clearance_edges_added"] == 0 and absent["enzyme_genes_without_a_node"] == ["NOT_IN_THE_GRAPH"]


def test_the_loader_gives_the_source_graph_back_from_a_graph_with_clearance_edges(tmp_path: Path) -> None:
    nodes, edges, new_nodes, new_edges, relations, entry_table, _ = variant_with_enzymes()
    write_graph(tmp_path / "source", nodes, edges, RELATIONS)
    write_graph(tmp_path / "enzymes", new_nodes, new_edges, relations, entry_table)
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    rows = [evidence_row("knockout_C", "gene", "cluster:C", [["GENE:C", -1.0, 1.0]])]
    rows += [evidence_row(perturbation, "drug", perturbation, [list(triple) for triple in triples]) for perturbation, (_, triples) in DRUG_SEEDS.items()]
    pd.DataFrame(rows).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    source = load_experiment_data(tmp_path / "source", evidence_directory)
    as_targets = load_experiment_data(tmp_path / "enzymes", evidence_directory, drug_entry="targets")
    as_nodes = load_experiment_data(tmp_path / "enzymes", evidence_directory, drug_entry="nodes")
    assert as_targets.node_ids == source.node_ids and as_targets.relation_types == source.relation_types
    for name in ("edge_source", "edge_target", "edge_relation", "edge_sign"):
        np.testing.assert_array_equal(getattr(as_targets, name), getattr(source, name))
    assert as_nodes.relation_types[-2:] == list(CHANGE_CARRYING_RELATIONS)
    fixed = rewiring_fixed_relations(as_nodes.relation_types, False, has_entry_nodes=True)
    assert fixed == [as_nodes.relation_types.index(relation) for relation in (DRUG_MECHANISM_RELATION, SEQUESTERS_RELATION, CLEARS_RELATION)]

    # the encoder keywords: one relation index per edge with both relations, and the keywords they always were with sequestration alone
    arguments = SimpleNamespace(sequestration_carries=None)
    keywords = sequestration_inputs(as_nodes, arguments)
    assert set(keywords) == {"change_only_edge_index", "change_only_edge_relation"}
    assert sorted(keywords["change_only_edge_relation"].tolist()) == sorted([as_nodes.relation_types.index(SEQUESTERS_RELATION)] + [as_nodes.relation_types.index(CLEARS_RELATION)] * 2)
    assert sequestration_inputs(as_nodes, SimpleNamespace(sequestration_carries="presence")) == {}
    plain_nodes, plain_edges, plain_relations, plain_entry = variant_with_enzymes(None)[2:6]
    write_graph(tmp_path / "plain", plain_nodes, plain_edges, plain_relations, plain_entry)
    plain = load_experiment_data(tmp_path / "plain", evidence_directory, drug_entry="nodes")
    assert set(sequestration_inputs(plain, arguments)) == {"change_only_edge_index", "change_only_relation_index"}
    assert sequestration_inputs(source, arguments) == {}


def encoder_on(new_nodes: pd.DataFrame, new_edges: pd.DataFrame, relations: list[str], entry_table: pd.DataFrame, change_carrying: tuple[str, ...]):
    """A message-passing encoder on a drug-node graph, with the edges of the named relations handed over as change-only edges."""
    all_types = sorted(set(new_nodes.node_type))
    index_of = {node_id: index for index, node_id in enumerate(new_nodes.node_id)}
    edge_index = torch.tensor([[index_of[node] for node in new_edges.source_id], [index_of[node] for node in new_edges.target_id]])
    relation = torch.tensor([relations.index(name) for name in new_edges.relation_type])
    is_entry = relation == relations.index(DRUG_MECHANISM_RELATION)
    magnitude = {(source, target): weight for source, target, weight in zip(entry_table.drug_node_id, entry_table.target_node_id, entry_table.magnitude)}
    entry_weight = torch.tensor([magnitude[pair] for pair in zip(new_edges.source_id[is_entry.numpy()], new_edges.target_id[is_entry.numpy()])])
    carries_change = torch.isin(relation, torch.tensor([relations.index(name) for name in change_carrying]))
    if len(change_carrying) == 1:
        change_keywords = {"change_only_edge_index": edge_index[:, carries_change], "change_only_relation_index": relations.index(change_carrying[0])}
    else:
        change_keywords = {"change_only_edge_index": edge_index[:, carries_change], "change_only_edge_relation": relation[carries_change]}
    torch.manual_seed(0)
    encoder = RelationalMessagePassingEncoder(
        len(new_nodes), len(relations), 6, 2, node_features=typed_features(list(new_nodes.node_type), all_types),
        entry_node_mask=torch.tensor(list(new_nodes.node_type == DRUG_NODE_TYPE)), entry_edge_index=edge_index[:, is_entry],
        entry_edge_sign=torch.tensor(list(new_edges.sign[is_entry.numpy()]), dtype=torch.float32), entry_edge_weight=entry_weight,
        entry_relation_indices=(relations.index(DRUG_MECHANISM_RELATION),), **change_keywords)
    with torch.no_grad():
        encoder.layer_bias.copy_(torch.randn_like(encoder.layer_bias))
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, relation, len(new_nodes), len(relations))
    return encoder, adjacencies, index_of


def test_a_clearance_edge_carries_the_enzyme_s_change_and_leaves_the_carrier_s_edge_as_it_was() -> None:
    """GENE:C clears drug_one and drug_two; ALB sequesters drug_one. With both relations carrying change: a drug taken
    alone hears neither edge; the enzyme lowered with a drug moves that drug's node; and with the clearance weights at
    zero the encoder is the one that knows sequestration alone."""
    _, _, new_nodes, new_edges, relations, entry_table, _ = variant_with_enzymes()
    both, adjacencies, index_of = encoder_on(new_nodes, new_edges, relations, entry_table, CHANGE_CARRYING_RELATIONS)
    clears_index, sequesters_index = relations.index(CLEARS_RELATION), relations.index(SEQUESTERS_RELATION)
    assert both.change_only_relation_indices == (sequesters_index, clears_index)
    drug_one, drug_two, enzyme, carrier = (index_of[name] for name in ("DRUG:drug_one", "DRUG:drug_two", "GENE:C", "GENE:ALB"))
    # two change-only edges of one relation end on different drugs, one edge of the other relation on drug_one: each is a mean within its relation
    assert both.change_only_edge_scale.tolist() == [1.0, 1.0, 1.0]

    def field(encoder, seeds: list[int], signs: list[float]) -> torch.Tensor:
        sign_and_magnitude = torch.tensor([[[sign, 1.0] for sign in signs]])
        with torch.no_grad():
            return encoder.perturbation_difference_field(torch.tensor([seeds]), sign_and_magnitude, adjacencies)

    def with_weights_zeroed(*relation_indices: int):
        silenced, _, _ = encoder_on(new_nodes, new_edges, relations, entry_table, CHANGE_CARRYING_RELATIONS)
        with torch.no_grad():
            silenced.load_state_dict(both.state_dict())
            for relation_index in relation_indices:
                silenced.relation_weight[:, relation_index] = 0.0
        return silenced

    without_either = with_weights_zeroed(clears_index, sequesters_index)
    # taken alone, each drug's field is the field of the graph without either edge
    for drug in (drug_one, drug_two):
        torch.testing.assert_close(field(both, [drug], [1.0]), field(without_either, [drug], [1.0]), atol=1e-6, rtol=0)
    # the enzyme lowered with drug_two: the clearance edge sends, and the drug's node and its target move
    moved = (field(both, [drug_two, enzyme], [1.0, -1.0]) - field(without_either, [drug_two, enzyme], [1.0, -1.0])).abs().amax(dim=2)[0]
    assert float(moved[drug_two]) > 1e-4 and float(moved[index_of["GENE:B"]]) > 1e-4
    assert float(field(both, [enzyme], [-1.0])[0, [drug_one, drug_two]].abs().max()) == 0.0  # no drug is taken: no drug node holds anything
    # the carrier lowered with drug_one: the sequestration edge sends. In the toy graph ALB also activates C, the enzyme,
    # so after one layer the enzyme has moved too and the clearance edge carries that into drug_one's node in the second
    # layer, and nowhere else: an enzyme that a perturbation moves indirectly is heard as one that is seeded
    without_clearance = with_weights_zeroed(clears_index)
    through_clearance = (field(both, [drug_one, carrier], [1.0, -1.0]) - field(without_clearance, [drug_one, carrier], [1.0, -1.0])).abs().amax(dim=2)[0]
    assert float(through_clearance[drug_one]) > 1e-4
    assert float(through_clearance[[index for index in range(len(new_nodes)) if index != drug_one]].max()) < 1e-6
    assert float((field(without_clearance, [drug_one, carrier], [1.0, -1.0]) - field(without_either, [drug_one, carrier], [1.0, -1.0])).abs().max()) > 1e-4

    # sequestration alone, built the way it always was, gives the same field as both relations with the clearance weights at zero
    plain_nodes, plain_edges, plain_relations, plain_entry = variant_with_enzymes(None)[2:6]
    alone, plain_adjacencies, _ = encoder_on(plain_nodes, plain_edges, plain_relations, plain_entry, (SEQUESTERS_RELATION,))
    with torch.no_grad():
        state = both.state_dict()
        state["relation_weight"] = state["relation_weight"][:, : len(plain_relations)]
        alone.load_state_dict(state)
        expected = alone.perturbation_difference_field(torch.tensor([[drug_one, carrier]]), torch.tensor([[[1.0, 1.0], [-1.0, 1.0]]]), plain_adjacencies)
    torch.testing.assert_close(field(without_clearance, [drug_one, carrier], [1.0, -1.0]), expected, atol=1e-6, rtol=0)

    with pytest.raises(ValueError, match="give the relation of the change-only edges once"):
        RelationalMessagePassingEncoder(5, 3, 4, 2, change_only_edge_index=torch.tensor([[0], [1]]), change_only_relation_index=1, change_only_edge_relation=torch.tensor([1]))
    with pytest.raises(ValueError, match="one entry per change-only edge"):
        RelationalMessagePassingEncoder(5, 3, 4, 2, change_only_edge_index=torch.tensor([[0, 2], [1, 1]]), change_only_edge_relation=torch.tensor([1]))


def test_a_row_of_the_regulator_s_table_gives_an_edge_only_when_it_states_the_direction() -> None:
    higher = "Results in higher systemic concentrations and higher adverse reaction risk. Dosage adjustment is recommended."
    assert edge_decision("aripiprazole", "CYP2D6", "poor metabolizers", higher)[0]
    assert edge_decision("doxepin", "CYP2D6", "ultrarapid, intermediate, or poor metabolizers", "May alter systemic concentrations. Poor metabolizers have higher systemic concentrations and adverse reaction risk.")[0]
    assert edge_decision("clobazam", "CYP2C19", "intermediate or poor metabolizers", "Results in higher systemic active metabolite concentrations. Poor metabolism results in higher adverse reaction risk.")[0]
    assert edge_decision("cevimeline", "CYP2D6", "poor metabolizers", "May result in higher adverse reaction risk. Use with caution.")[0]
    assert edge_decision("amitriptyline", "CYP2D6", "ultrarapid, intermediate, or poor metabolizers", "May alter systemic concentrations.") == (False, "no direction stated")
    # the enzyme forms tramadol's active metabolite: less enzyme, less effect, so a clearance edge would have the wrong sign
    assert not edge_decision("tramadol", "CYP2D6", "poor metabolizers", "Results in lower systemic active metabolite concentrations and may result in reduced efficacy.")[0]
    assert not edge_decision("tramadol", "CYP2D6", "Ultrarapid metabolizers3", "Results in higher systemic and breast milk active metabolite concentrations, which may result in respiratory depression and death.")[0]
    assert not edge_decision("carbamazepine", "HLA-B", "*15:02 allele positive", "Results in higher adverse reaction risk (severe skin reactions).")[0]
    assert edge_decision("bupropion", "CYP2D6", "poor metabolizers", "Results in higher systemic concentrations of dextromethorphan and may result in higher adverse reaction risk.") == (
        False, "the row is about dextromethorphan, another ingredient")
    assert json.dumps(edge_decision("phenytoin", "CYP2C9", "intermediate or poor metabolizers", "May result in higher systemic concentrations and higher adverse reaction risk (central nervous system toxicity).")[0]) == "true"
