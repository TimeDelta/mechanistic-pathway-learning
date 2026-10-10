"""Tests: drug entry nodes (graph/drug_entry_nodes.py, docs/drug_entry_nodes.md). The builder appends a node per drug
with its mechanism edges and its carriers' sequestration edges and changes no existing row, a drug that is a graph
compound gets no node, the loader's "targets" mode gives the source graph back and its "nodes" mode seeds each drug on
its own node with the leakage groups and degrees of the targets, and both encoders treat a drug node as present only
where it is seeded: a gene knockout's field does not depend on the drug nodes, a drug's field does not depend on the
other drugs, and the linear response of a drug seeded on its node is the response of the drug seeded on its targets
scaled by the entry gain."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import torch

from experiments.run_main_model import perturbations_with_a_zero_signed_input, rewiring_fixed_relations, seed_mask_partner_pairs
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.graph.drug_entry_nodes import (
    DRUG_ENTRY_FILE,
    DRUG_MECHANISM_RELATION,
    DRUG_NODE_TYPE,
    ENTRY_RELATIONS,
    SEQUESTERS_RELATION,
    add_drug_entry_nodes,
    drug_node_id,
    extend_table_for_drug_nodes,
    is_graph_compound,
)
from mechanistic_pathway_learning.models.descriptor_treatments import seed_unit_mask
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder

RELATIONS = ["activates", "inhibits", "targets", "catalyzed_by"]
NODE_IDS = ["GENE:A", "GENE:B", "GENE:C", "GENE:ALB", "R1", "MAM1e", "MAM1c"]
DRUG_SEEDS = {
    "drug_one": ("one", [("GENE:A", -1.0, 0.5), ("GENE:B", -1.0, 0.5)]),
    "drug_two": ("two", [("GENE:B", 0.5, 1.0)]),
    "compound": ("compound", [("MAM1e", 1.0, 0.5), ("MAM1c", 1.0, 0.5)]),
}
CARRIERS = pd.DataFrame({"perturbation_id": ["drug_one", "compound"], "carrier": ["albumin", "albumin"], "evidence": ["label", "database"]})
GRAPH_COMPOUNDS = [frozenset({"MAM1"})]
# for the encoder tests: a carrier on a drug node only. The compound's carriage edge joins two nodes the source graph
# already has (the carrier and the compound), so it is a real edge for every perturbation and is tested on its own.
CARRIERS_OF_DRUG_NODES = CARRIERS[CARRIERS.perturbation_id == "drug_one"]


def toy_graph() -> tuple[pd.DataFrame, pd.DataFrame]:
    """A activates B, B inhibits C, C catalyses R1, ALB activates C; one metabolite in two compartments."""
    nodes = pd.DataFrame({
        "node_id": NODE_IDS,
        "node_type": ["gene", "gene", "gene", "gene", "reaction", "metabolite", "metabolite"],
        "display_name": NODE_IDS,
        "compartment": [None, None, None, None, "c", "e", "c"],
        "base_metabolite_id": [None, None, None, None, None, "MAM1", "MAM1"],
        "is_currency": [False] * 7,
        "degree": [1, 2, 3, 1, 1, 0, 0]})
    edges = pd.DataFrame({
        "source_id": ["GENE:A", "GENE:B", "GENE:C", "GENE:ALB"],
        "target_id": ["GENE:B", "GENE:C", "R1", "GENE:C"],
        "relation_type": ["activates", "inhibits", "catalyzed_by", "activates"],
        "sign": [1.0, -1.0, 1.0, 1.0],
        "evidence_source": ["toy"] * 4})
    return nodes, edges


def variant(drug_seeds=None, carriers=CARRIERS):
    nodes, edges = toy_graph()
    return (nodes, edges, *add_drug_entry_nodes(nodes, edges, RELATIONS, DRUG_SEEDS if drug_seeds is None else drug_seeds, GRAPH_COMPOUNDS, carriers))


def test_builder_appends_drug_nodes_mechanism_and_sequestration_edges_and_changes_no_existing_row() -> None:
    nodes, edges, new_nodes, new_edges, relations, entry_table, summary = variant()
    pd.testing.assert_frame_equal(new_nodes.iloc[: len(nodes)].reset_index(drop=True), nodes)  # degrees of existing nodes included
    pd.testing.assert_frame_equal(new_edges.iloc[: len(edges)].reset_index(drop=True), edges)
    assert list(new_nodes.node_id[len(nodes):]) == ["DRUG:drug_one", "DRUG:drug_two"]  # the graph compound gets no node
    assert set(new_nodes.node_type[len(nodes):]) == {DRUG_NODE_TYPE}
    added = new_edges.iloc[len(edges):]
    mechanism = added[added.relation_type == DRUG_MECHANISM_RELATION]
    assert sorted(zip(mechanism.source_id, mechanism.target_id, mechanism.sign)) == [
        ("DRUG:drug_one", "GENE:A", -1.0), ("DRUG:drug_one", "GENE:B", -1.0), ("DRUG:drug_two", "GENE:B", 0.5)]
    carriage = added[added.relation_type == SEQUESTERS_RELATION]
    # the carrier sequesters the drug node, and for the graph compound its extracellular node only
    assert sorted(zip(carriage.source_id, carriage.target_id, carriage.sign)) == [("GENE:ALB", "DRUG:drug_one", -1.0), ("GENE:ALB", "MAM1e", -1.0)]
    assert relations == RELATIONS + [SEQUESTERS_RELATION]  # targets was reserved by the first graph build, so it is not added again
    assert dict(zip(new_nodes.node_id, new_nodes.degree))["DRUG:drug_one"] == 3  # two mechanism edges and one carrier
    assert sorted(zip(entry_table.perturbation_id, entry_table.target_node_id, entry_table.magnitude)) == [
        ("drug_one", "GENE:A", 0.5), ("drug_one", "GENE:B", 0.5), ("drug_two", "GENE:B", 1.0)]
    assert summary["drug_nodes_added"] == 2 and summary["drugs_attached_to_a_graph_compound"] == ["compound"]
    assert summary["mechanism_edges_added"] == 3 and summary["sequestration_edges_added"] == 2


def test_a_drug_whose_target_is_a_metabolite_is_not_a_graph_compound() -> None:
    base_of_node = {"MAM1e": "MAM1", "MAM1c": "MAM1", "MAM2e": "MAM2"}
    assert is_graph_compound(["MAM1e", "MAM1c"], base_of_node, GRAPH_COMPOUNDS)
    assert not is_graph_compound(["MAM2e"], base_of_node, GRAPH_COMPOUNDS)  # a metabolite target the compound table does not list
    assert not is_graph_compound(["MAM1e", "GENE:A"], base_of_node, GRAPH_COMPOUNDS)
    assert not is_graph_compound([], base_of_node, GRAPH_COMPOUNDS)
    with pytest.raises(ValueError, match="carrier of drug_one"):
        variant(carriers=pd.DataFrame({"perturbation_id": ["drug_one"], "carrier": ["transferrin"]}))


def test_per_node_tables_get_one_row_per_drug_node() -> None:
    table = pd.DataFrame({"first": [0.25, 0.5], "second": [1.0, 0.0]}, index=pd.Index(["GENE:A", "GENE:B"], name="node_id"))
    extended = extend_table_for_drug_nodes(table, ["DRUG:x", "DRUG:y"], fill_value=1.0)
    assert list(extended.index) == ["GENE:A", "GENE:B", "DRUG:x", "DRUG:y"]
    pd.testing.assert_frame_equal(extended.iloc[:2], table)
    assert extended.loc["DRUG:y"].tolist() == [1.0, 1.0]
    with pytest.raises(ValueError, match="already have a row"):
        extend_table_for_drug_nodes(extended, ["DRUG:x"], fill_value=0.0)


def write_graph(directory: Path, nodes: pd.DataFrame, edges: pd.DataFrame, relations: list[str], entry_table: pd.DataFrame | None = None) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    nodes.to_parquet(directory / "nodes.parquet", index=False)
    edges.to_parquet(directory / "edges.parquet", index=False)
    (directory / "relation_types.json").write_text(json.dumps(relations))
    if entry_table is not None:
        entry_table.to_parquet(directory / DRUG_ENTRY_FILE, index=False)


def evidence_row(perturbation: str, kind: str, cluster: str, seeds: list) -> dict:
    return {"perturbation_id": perturbation, "perturbation_type": kind, "perturbation_label": perturbation, "group_id": perturbation,
            "disease_cluster_id": cluster, "symptom": "anxiety", "relation": "induces", "evidence_class": "monogenic", "grade": "A",
            "weight": 1.0, "label_frequency": None, "in_metabolic_layer": True, "evidence_date": None, "perturbation_nodes": json.dumps(seeds)}


def loaded(tmp_path: Path):
    nodes, edges, new_nodes, new_edges, relations, entry_table, _ = variant()
    write_graph(tmp_path / "source", nodes, edges, RELATIONS)
    write_graph(tmp_path / "drugs", new_nodes, new_edges, relations, entry_table)
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    rows = [evidence_row("knockout_B", "gene", "cluster:B", [["GENE:B", -1.0, 1.0]]), evidence_row("knockout_ALB", "gene", "cluster:ALB", [["GENE:ALB", -1.0, 1.0]])]
    rows += [evidence_row(perturbation, "drug", perturbation, [list(triple) for triple in triples]) for perturbation, (_, triples) in DRUG_SEEDS.items()]
    pd.DataFrame(rows).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    keywords = {"group_by": "disease_cluster_and_targets"}
    return (load_experiment_data(tmp_path / "source", evidence_directory, **keywords),
            load_experiment_data(tmp_path / "drugs", evidence_directory, drug_entry="targets", **keywords),
            load_experiment_data(tmp_path / "drugs", evidence_directory, drug_entry="nodes", **keywords))


def test_loader_targets_mode_is_the_source_graph_and_nodes_mode_seeds_each_drug_on_its_own_node(tmp_path: Path) -> None:
    source, as_targets, as_nodes = loaded(tmp_path)
    assert as_targets.node_ids == source.node_ids and as_targets.relation_types == source.relation_types
    for name in ("edge_source", "edge_target", "edge_relation", "edge_sign", "node_degree", "outcomes"):
        np.testing.assert_array_equal(getattr(as_targets, name), getattr(source, name))
    assert all(np.array_equal(first, second) for first, second in zip(as_targets.perturbation_seeds, source.perturbation_seeds))
    assert as_targets.entry_node_mask is None and as_targets.edge_weight is None and as_targets.drug_entry is None
    np.testing.assert_array_equal(as_targets.structural_node_features(), source.structural_node_features())

    seeds = {perturbation: [as_nodes.node_ids[index] for index in indices] for perturbation, indices in zip(as_nodes.perturbation_ids, as_nodes.perturbation_seeds)}
    assert seeds == {"compound": ["MAM1e", "MAM1c"], "drug_one": ["DRUG:drug_one"], "drug_two": ["DRUG:drug_two"], "knockout_ALB": ["GENE:ALB"], "knockout_B": ["GENE:B"]}
    position = as_nodes.perturbation_ids.index("drug_one")
    assert as_nodes.perturbation_signs[position].tolist() == [1.0] and as_nodes.perturbation_magnitudes[position].tolist() == [1.0]
    assert as_nodes.perturbation_magnitudes[as_nodes.perturbation_ids.index("compound")].tolist() == [0.5, 0.5]  # the compound keeps its own seeds
    assert as_nodes.group_ids == source.group_ids  # groups come from the targets: drug_one and drug_two share GENE:B with its knockout
    assert len({as_nodes.group_ids[as_nodes.perturbation_ids.index(name)] for name in ("drug_one", "drug_two", "knockout_B")}) == 1
    np.testing.assert_array_equal(as_nodes.perturbation_degrees, source.perturbation_degrees)  # degrees come from the targets too
    np.testing.assert_array_equal(as_nodes.perturbation_degrees_for_strata, source.perturbation_degrees_for_strata)
    np.testing.assert_array_equal(as_nodes.node_degree[: len(source.node_ids)], source.node_degree)
    mechanism = np.asarray(as_nodes.edge_relation) == as_nodes.relation_types.index(DRUG_MECHANISM_RELATION)
    weights = {(as_nodes.node_ids[source_index], as_nodes.node_ids[target_index]): weight
               for source_index, target_index, weight in zip(as_nodes.edge_source[mechanism], as_nodes.edge_target[mechanism], as_nodes.edge_weight[mechanism])}
    assert weights == {("DRUG:drug_one", "GENE:A"): 0.5, ("DRUG:drug_one", "GENE:B"): 0.5, ("DRUG:drug_two", "GENE:B"): 1.0}
    assert (as_nodes.edge_weight[~mechanism] == 1.0).all()
    assert [as_nodes.node_ids[index] for index in np.flatnonzero(as_nodes.entry_node_mask)] == ["DRUG:drug_one", "DRUG:drug_two"]
    with pytest.raises(ValueError, match="drug_entry must be one of"):
        load_experiment_data(tmp_path / "drugs", tmp_path / "evidence", drug_entry="both")


def test_rewiring_holds_a_drug_s_edges_fixed_and_a_binding_agent_has_no_signed_input(tmp_path: Path) -> None:
    _, _, as_nodes = loaded(tmp_path)
    relations = as_nodes.relation_types
    assert rewiring_fixed_relations(relations, False, has_entry_nodes=True) == [relations.index(DRUG_MECHANISM_RELATION), relations.index(SEQUESTERS_RELATION)]
    assert rewiring_fixed_relations(relations, False) == []  # without entry nodes nothing is added, as before
    assert perturbations_with_a_zero_signed_input(as_nodes) == []
    binding_agent = {"drug_one": ("one", [("GENE:A", 0.0, 1.0)])}  # sign 0, as ChEMBL's BINDING AGENT maps
    nodes, edges, new_nodes, new_edges, relations, entry_table, _ = variant(drug_seeds=binding_agent, carriers=None)
    write_graph(tmp_path / "binding", new_nodes, new_edges, relations, entry_table)
    evidence_directory = tmp_path / "binding_evidence"
    evidence_directory.mkdir()
    pd.DataFrame([evidence_row("drug_one", "drug", "drug_one", [["GENE:A", 0.0, 1.0]]), evidence_row("knockout_B", "gene", "cluster:B", [["GENE:B", -1.0, 1.0]])]
                 ).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    data = load_experiment_data(tmp_path / "binding", evidence_directory, drug_entry="nodes")
    assert perturbations_with_a_zero_signed_input(data) == ["drug_one"]  # its seed has sign 1, but its one mechanism edge has sign 0


def typed_features(node_types: list[str], all_types: list[str]) -> torch.Tensor:
    return torch.tensor([[1.0 if node_type == name else 0.0 for name in all_types] for node_type in node_types])


def message_passing_pair(num_layers: int = 2, carriers: pd.DataFrame = CARRIERS_OF_DRUG_NODES, sequestration_carries: str = "presence"):
    """An encoder on the toy graph and one on its drug-node variant sharing every parameter the two have in common.
    sequestration_carries "change" hands the variant's sequestration edges to the encoder as change-only edges."""
    nodes, edges, new_nodes, new_edges, relations, entry_table, _ = variant(carriers=carriers)
    all_types = sorted(set(new_nodes.node_type))
    index_of = {node_id: index for index, node_id in enumerate(new_nodes.node_id)}

    def arrays(edge_table: pd.DataFrame, relation_list: list[str]):
        return (torch.tensor([[index_of[node] for node in edge_table.source_id], [index_of[node] for node in edge_table.target_id]]),
                torch.tensor([relation_list.index(relation) for relation in edge_table.relation_type]))

    torch.manual_seed(0)
    source_edge_index, source_relation = arrays(edges, RELATIONS)
    source_encoder = RelationalMessagePassingEncoder(len(nodes), len(RELATIONS), 6, num_layers, node_features=typed_features(list(nodes.node_type), all_types))
    source_adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(source_edge_index, source_relation, len(nodes), len(RELATIONS))
    edge_index, relation = arrays(new_edges, relations)
    is_entry = relation == relations.index(DRUG_MECHANISM_RELATION)
    magnitude = {(source, target): weight for source, target, weight in zip(entry_table.drug_node_id, entry_table.target_node_id, entry_table.magnitude)}
    entry_weight = torch.tensor([magnitude[pair] for pair in zip(new_edges.source_id[is_entry.numpy()], new_edges.target_id[is_entry.numpy()])])
    encoder = RelationalMessagePassingEncoder(
        len(new_nodes), len(relations), 6, num_layers, node_features=typed_features(list(new_nodes.node_type), all_types),
        entry_node_mask=torch.tensor(list(new_nodes.node_type == DRUG_NODE_TYPE)), entry_edge_index=edge_index[:, is_entry],
        entry_edge_sign=torch.tensor(list(new_edges.sign[is_entry.numpy()]), dtype=torch.float32), entry_edge_weight=entry_weight,
        entry_relation_indices=(relations.index(DRUG_MECHANISM_RELATION),),
        **({"change_only_edge_index": edge_index[:, relation == relations.index(SEQUESTERS_RELATION)], "change_only_relation_index": relations.index(SEQUESTERS_RELATION)}
           if sequestration_carries == "change" else {}))
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, relation, len(new_nodes), len(relations))
    with torch.no_grad():  # the variant has one relation more (sequesters, appended last) and the same feature columns
        encoder.feature_projection.load_state_dict(source_encoder.feature_projection.state_dict())
        encoder.perturbation_injection.load_state_dict(source_encoder.perturbation_injection.state_dict())
        encoder.relation_weight[:, : len(RELATIONS)] = source_encoder.relation_weight
        encoder.self_weight.copy_(source_encoder.self_weight)
        encoder.layer_bias.copy_(torch.randn_like(encoder.layer_bias))  # a nonzero bias is what would give an absent node a state
        source_encoder.layer_bias.copy_(encoder.layer_bias)
    return source_encoder, source_adjacencies, encoder, adjacencies, index_of, len(nodes)


def test_message_passing_field_of_a_gene_knockout_does_not_depend_on_the_drug_nodes() -> None:
    source_encoder, source_adjacencies, encoder, adjacencies, index_of, num_source_nodes = message_passing_pair()
    for gene in ("GENE:B", "GENE:ALB"):  # ALB is the carrier, so its knockout is the one a sequestration edge could leak through
        seeds, sign_and_magnitude = torch.tensor([[index_of[gene]]]), torch.tensor([[[-1.0, 1.0]]])
        expected = source_encoder.perturbation_difference_field(seeds, sign_and_magnitude, source_adjacencies)
        field = encoder.perturbation_difference_field(seeds, sign_and_magnitude, adjacencies)
        torch.testing.assert_close(field[:, :num_source_nodes], expected, atol=1e-6, rtol=0)
        assert float(field[:, num_source_nodes:].abs().max()) == 0.0  # every drug node is absent: exactly zero


def test_a_graph_compound_s_carriage_edge_is_a_real_edge_for_every_perturbation() -> None:
    """The carrier of a drug that is a graph compound sequesters the compound's own node, which the body also makes, so
    a knockout of the carrier moves that node whether or not the drug is taken. This is the one way the variant
    changes a gene knockout's field, and it is meant."""
    source_encoder, source_adjacencies, encoder, adjacencies, index_of, num_source_nodes = message_passing_pair(carriers=CARRIERS)
    seeds, sign_and_magnitude = torch.tensor([[index_of["GENE:ALB"]]]), torch.tensor([[[-1.0, 1.0]]])
    expected = source_encoder.perturbation_difference_field(seeds, sign_and_magnitude, source_adjacencies)
    field = encoder.perturbation_difference_field(seeds, sign_and_magnitude, adjacencies)
    changed = (field[0, :num_source_nodes] - expected[0]).abs().amax(dim=1) > 1e-6
    assert [node for node, index in index_of.items() if index < num_source_nodes and bool(changed[index])] == ["MAM1e"]


def test_message_passing_drug_node_is_present_only_where_seeded_and_sends_its_mechanism_with_its_sign() -> None:
    _, _, encoder, adjacencies, index_of, num_source_nodes = message_passing_pair(num_layers=1)
    drug_one, drug_two = index_of["DRUG:drug_one"], index_of["DRUG:drug_two"]
    seeds, sign_and_magnitude = torch.tensor([[drug_one], [drug_two]]), torch.tensor([[[1.0, 1.0]], [[1.0, 1.0]]])
    field = encoder.perturbation_difference_field(seeds, sign_and_magnitude, adjacencies)
    assert float(field[0, drug_two].abs().max()) == 0.0 and float(field[1, drug_one].abs().max()) == 0.0  # the other drug is absent
    assert float(field[0, drug_one].abs().max()) > 0.0 and float(field[1, drug_two].abs().max()) > 0.0
    assert float(field[0, index_of["GENE:A"]].abs().max()) > 0.0  # one layer reaches the drug's targets
    assert float(field[1, index_of["GENE:A"]].abs().max()) == 0.0  # drug_two does not target A
    # the message is weight * (x W_unsigned) + sign * weight * (x W_signed), added without dividing by the two drugs that target B
    drug_state = encoder.base_state()[drug_two] + encoder.perturbation_injection(sign_and_magnitude[1, 0])
    target_before = encoder.base_state()[index_of["GENE:B"]]
    incoming = encoder.base_state()[index_of["GENE:A"]] @ encoder.relation_weight[0, RELATIONS.index("activates")]  # A activates B, in both passes
    message = 1.0 * (drug_state @ encoder.entry_unsigned_weight[0]) + 0.5 * 1.0 * (drug_state @ encoder.entry_signed_weight[0])
    with_drug = torch.relu(target_before @ encoder.self_weight[0] + incoming + message + encoder.layer_bias[0])
    without_drug = torch.relu(target_before @ encoder.self_weight[0] + incoming + encoder.layer_bias[0])
    torch.testing.assert_close(field[1, index_of["GENE:B"]], with_drug - without_drug, atol=1e-6, rtol=0)


def test_message_passing_field_of_a_drug_does_not_depend_on_the_other_drugs_in_the_table() -> None:
    _, _, encoder, adjacencies, index_of, num_source_nodes = message_passing_pair()
    seeds, sign_and_magnitude = torch.tensor([[index_of["DRUG:drug_two"]]]), torch.tensor([[[1.0, 1.0]]])
    field = encoder.perturbation_difference_field(seeds, sign_and_magnitude, adjacencies)
    # the same graph without drug_one's node: drug_two's field on the shared nodes must be the same
    nodes, edges, alone_nodes, alone_edges, relations, entry_table, _ = variant(drug_seeds={"drug_two": DRUG_SEEDS["drug_two"]}, carriers=None)
    alone_index = {node_id: index for index, node_id in enumerate(alone_nodes.node_id)}
    edge_index = torch.tensor([[alone_index[node] for node in alone_edges.source_id], [alone_index[node] for node in alone_edges.target_id]])
    relation = torch.tensor([relations.index(name) for name in alone_edges.relation_type])
    is_entry = relation == relations.index(DRUG_MECHANISM_RELATION)
    all_types = sorted({*alone_nodes.node_type})
    alone = RelationalMessagePassingEncoder(
        len(alone_nodes), len(relations), 6, 2, node_features=typed_features(list(alone_nodes.node_type), all_types),
        entry_node_mask=torch.tensor(list(alone_nodes.node_type == DRUG_NODE_TYPE)), entry_edge_index=edge_index[:, is_entry],
        entry_edge_sign=torch.tensor([0.5]), entry_edge_weight=torch.tensor([1.0]), entry_relation_indices=(relations.index(DRUG_MECHANISM_RELATION),))
    with torch.no_grad():
        alone.feature_projection.load_state_dict(encoder.feature_projection.state_dict())
        alone.perturbation_injection.load_state_dict(encoder.perturbation_injection.state_dict())
        alone.relation_weight.copy_(encoder.relation_weight[:, : len(relations)])
        for name in ("self_weight", "layer_bias", "entry_unsigned_weight", "entry_signed_weight"):
            getattr(alone, name).copy_(getattr(encoder, name))
    alone_field = alone.perturbation_difference_field(torch.tensor([[alone_index["DRUG:drug_two"]]]), sign_and_magnitude,
                                                      RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, relation, len(alone_nodes), len(relations)))
    torch.testing.assert_close(field[:, :num_source_nodes], alone_field[:, :num_source_nodes], atol=1e-6, rtol=0)


def test_encoders_without_entry_nodes_have_the_parameters_they_always_had() -> None:
    plain = RelationalMessagePassingEncoder(5, 2, 4, 2)
    assert sorted(plain.state_dict()) == ["base_node_state.weight", "layer_bias", "perturbation_injection.bias", "perturbation_injection.weight", "relation_weight", "self_weight"]
    assert plain.present_nodes(torch.tensor([[0]])) is None
    with pytest.raises(ValueError, match="entry edges need entry_node_mask"):
        RelationalMessagePassingEncoder(5, 2, 4, 2, entry_relation_indices=(1,))
    with pytest.raises(ValueError, match="every entry edge must leave an entry node"):
        RelationalMessagePassingEncoder(5, 2, 4, 2, entry_node_mask=torch.tensor([False, False, False, False, True]), entry_edge_index=torch.tensor([[0], [1]]))
    linear = LinearResponseEncoder(3, ["activates"], torch.tensor([0]), torch.tensor([1]), torch.tensor([0]), torch.tensor([1.0]), torch.ones(3, 2), 4)
    assert "entry_gain_logit" not in linear.state_dict() and linear.present_nodes(torch.tensor([[0]])) is None
    with pytest.raises(ValueError, match="entry relations need entry_node_mask"):
        LinearResponseEncoder(3, ["activates"], torch.tensor([0]), torch.tensor([1]), torch.tensor([0]), torch.tensor([1.0]), torch.ones(3, 2), 4, entry_relation_names=("activates",))


def linear_response_pair(num_propagation_steps: int = 60):
    nodes, edges, new_nodes, new_edges, relations, entry_table, _ = variant(carriers=CARRIERS_OF_DRUG_NODES)
    index_of = {node_id: index for index, node_id in enumerate(new_nodes.node_id)}
    all_types = sorted(set(new_nodes.node_type))

    def tensors(edge_table: pd.DataFrame, relation_list: list[str]):
        return (torch.tensor([index_of[node] for node in edge_table.source_id]), torch.tensor([index_of[node] for node in edge_table.target_id]),
                torch.tensor([relation_list.index(relation) for relation in edge_table.relation_type]), torch.tensor(list(edge_table.sign), dtype=torch.float32))

    torch.manual_seed(1)
    source_encoder = LinearResponseEncoder(len(nodes), RELATIONS, *tensors(edges, RELATIONS), typed_features(list(nodes.node_type), all_types), 6,
                                           num_propagation_steps=num_propagation_steps, propagation_channels=3)
    magnitude = {(source, target): weight for source, target, weight in zip(entry_table.drug_node_id, entry_table.target_node_id, entry_table.magnitude)}
    edge_weight = torch.tensor([magnitude.get(pair, 1.0) for pair in zip(new_edges.source_id, new_edges.target_id)])
    encoder = LinearResponseEncoder(len(new_nodes), relations, *tensors(new_edges, relations), typed_features(list(new_nodes.node_type), all_types), 6,
                                    num_propagation_steps=num_propagation_steps, propagation_channels=3,
                                    entry_node_mask=torch.tensor(list(new_nodes.node_type == DRUG_NODE_TYPE)), entry_relation_names=ENTRY_RELATIONS, edge_weight=edge_weight)
    with torch.no_grad():  # the gains of the relations both graphs have, and the input weights
        encoder.gain_logit[: len(RELATIONS)] = source_encoder.gain_logit
        encoder.input_weight.copy_(source_encoder.input_weight)
    return source_encoder, encoder, index_of, len(nodes)


def test_linear_response_of_a_gene_knockout_does_not_depend_on_the_drug_nodes() -> None:
    source_encoder, encoder, index_of, num_source_nodes = linear_response_pair()
    for gene in ("GENE:A", "GENE:ALB"):  # ALB sequesters drug_one: its knockout must not reach drug_one's targets through the absent node
        seeds, sign_and_magnitude = torch.tensor([[index_of[gene]]]), torch.tensor([[[-1.0, 1.0]]])
        response = encoder.response(seeds, sign_and_magnitude)
        torch.testing.assert_close(response[:, :num_source_nodes], source_encoder.response(seeds, sign_and_magnitude), atol=1e-7, rtol=0)
        assert float(response[:, num_source_nodes:].abs().max()) == 0.0


def test_linear_response_of_a_drug_on_its_node_is_its_response_on_its_targets_scaled_by_the_entry_gain() -> None:
    source_encoder, encoder, index_of, num_source_nodes = linear_response_pair(num_propagation_steps=200)
    on_node = encoder.response(torch.tensor([[index_of["DRUG:drug_one"]]]), torch.tensor([[[1.0, 1.0]]]))
    on_targets = source_encoder.response(torch.tensor([[index_of["GENE:A"], index_of["GENE:B"]]]), torch.tensor([[[-1.0, 0.5], [-1.0, 0.5]]]))
    entry_gain = torch.sigmoid(encoder.entry_gain_logit)[0]  # one positive gain per channel
    torch.testing.assert_close(on_node[:, :num_source_nodes], on_targets * entry_gain, atol=1e-6, rtol=0)
    torch.testing.assert_close(on_node[0, index_of["DRUG:drug_one"]], encoder.input_weight, atol=1e-6, rtol=0)  # the node holds its own input
    assert float(on_node[0, index_of["DRUG:drug_two"]].abs().max()) == 0.0  # the other drug is absent
    # drug_two alone targets B with sign 0.5 and weight 1: B is not divided by the two drugs that share it
    two_on_node = encoder.response(torch.tensor([[index_of["DRUG:drug_two"]]]), torch.tensor([[[1.0, 1.0]]]))
    two_on_targets = source_encoder.response(torch.tensor([[index_of["GENE:B"]]]), torch.tensor([[[0.5, 1.0]]]))
    torch.testing.assert_close(two_on_node[:, :num_source_nodes], two_on_targets * entry_gain, atol=1e-6, rtol=0)
    assert drug_node_id("drug_one") == "DRUG:drug_one"


def test_mechanism_delivered_before_the_first_layer_reaches_as_far_as_a_seed_on_the_targets() -> None:
    """drug_two targets B, B inhibits C, C catalyses R1. With one layer a seed on B moves C; a seed on drug_two's node
    moves B and stops there, one edge short, unless the mechanism is also delivered before the first layer."""
    nodes, edges, new_nodes, new_edges, relations, entry_table, _ = variant(drug_seeds={"drug_two": DRUG_SEEDS["drug_two"]}, carriers=None)
    index_of = {node_id: index for index, node_id in enumerate(new_nodes.node_id)}
    edge_index = torch.tensor([[index_of[node] for node in new_edges.source_id], [index_of[node] for node in new_edges.target_id]])
    relation = torch.tensor([relations.index(name) for name in new_edges.relation_type])
    is_entry = relation == relations.index(DRUG_MECHANISM_RELATION)
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, relation, len(new_nodes), len(relations))
    features = typed_features(list(new_nodes.node_type), sorted(set(new_nodes.node_type)))
    seeds, sign_and_magnitude = torch.tensor([[index_of["DRUG:drug_two"]]]), torch.tensor([[[1.0, 1.0]]])

    def moved_nodes(before_first_layer: bool) -> set[str]:
        torch.manual_seed(4)
        encoder = RelationalMessagePassingEncoder(
            len(new_nodes), len(relations), 6, 1, node_features=features, entry_node_mask=torch.tensor(list(new_nodes.node_type == DRUG_NODE_TYPE)),
            entry_edge_index=edge_index[:, is_entry], entry_edge_sign=torch.tensor([0.5]), entry_edge_weight=torch.tensor([1.0]),
            entry_relation_indices=(relations.index(DRUG_MECHANISM_RELATION),), entry_before_first_layer=before_first_layer)
        field = encoder.perturbation_difference_field(seeds, sign_and_magnitude, adjacencies)[0]
        return {node for node, index in index_of.items() if float(field[index].abs().max()) > 0}

    assert moved_nodes(before_first_layer=False) == {"DRUG:drug_two", "GENE:B"}
    assert moved_nodes(before_first_layer=True) == {"DRUG:drug_two", "GENE:B", "GENE:C"}
    with pytest.raises(ValueError, match="entry edges need entry_node_mask"):
        RelationalMessagePassingEncoder(5, 2, 4, 2, entry_before_first_layer=True)


def entry_partner_pairs(new_edges: pd.DataFrame, index_of: dict[str, int]) -> torch.Tensor:
    """(drug node, target) pairs of the mechanism edges, as the trainer passes them for seed_masked on a graph without
    encodes edges."""
    mechanism = new_edges[new_edges.relation_type == DRUG_MECHANISM_RELATION]
    return torch.tensor([[index_of[node] for node in mechanism.source_id], [index_of[node] for node in mechanism.target_id]])


def test_seed_mask_partners_of_a_drug_node_are_its_targets_and_their_encodes_partners_one_way(tmp_path: Path) -> None:
    _, _, as_nodes = loaded(tmp_path)
    names = as_nodes.node_ids
    pairs = seed_mask_partner_pairs(as_nodes)
    assert sorted((names[node], names[partner]) for node, partner in pairs.T.tolist()) == [
        ("DRUG:drug_one", "GENE:A"), ("DRUG:drug_one", "GENE:B"), ("DRUG:drug_two", "GENE:B")]
    # a split graph: gene 0 encodes protein 1, the drug on node 2 targets the protein, gene 3 encodes protein 4 and has no drug
    split = SimpleNamespace(relation_types=["encodes", "targets"], edge_source=np.array([0, 2, 3]), edge_target=np.array([1, 1, 4]),
                            edge_relation=np.array([0, 1, 0]), edge_sign=np.array([1.0, -1.0, 1.0]), edge_weight=None,
                            entry_node_mask=np.array([False, False, True, False, False]))
    split_pairs = torch.as_tensor(seed_mask_partner_pairs(split))
    assert split_pairs.T.tolist() == [[0, 1], [3, 4], [1, 0], [4, 3], [2, 0], [2, 1]]
    masks = seed_unit_mask(torch.tensor([[2], [0], [1], [3]]), 5, split_pairs)
    assert masks.tolist() == [[True, True, True, False, False],     # the drug: its node, its target protein and that protein's gene
                              [True, True, False, False, False],    # the gene's knockout: the gene and its protein, as before, and no drug
                              [True, True, False, False, False],
                              [False, False, False, True, True]]
    without_entry_nodes = SimpleNamespace(**{**vars(split), "entry_node_mask": None})
    assert seed_mask_partner_pairs(without_entry_nodes).T.tolist() == [[0, 1], [3, 4], [1, 0], [4, 3]]  # the pairs a split graph always had


def test_seed_masked_hides_the_descriptors_of_a_drug_s_targets_under_message_passing() -> None:
    """Under seed_masked a drug seeded on its node reads its targets on their structural columns only, as a drug seeded
    on its targets does: its field does not depend on the descriptor values at the targets, it does depend on those of a
    node further on, and the reference of the difference field holds no drug node."""
    nodes, edges, new_nodes, new_edges, relations, entry_table, _ = variant(carriers=CARRIERS_OF_DRUG_NODES)
    index_of = {node_id: index for index, node_id in enumerate(new_nodes.node_id)}
    edge_index = torch.tensor([[index_of[node] for node in new_edges.source_id], [index_of[node] for node in new_edges.target_id]])
    relation = torch.tensor([relations.index(name) for name in new_edges.relation_type])
    is_entry = relation == relations.index(DRUG_MECHANISM_RELATION)
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, relation, len(new_nodes), len(relations))
    structural = typed_features(list(new_nodes.node_type), sorted(set(new_nodes.node_type)))
    descriptors = torch.randn(len(new_nodes), 2, generator=torch.Generator().manual_seed(7))
    descriptors[len(nodes):] = 0.0  # a drug node has no descriptor

    def drug_one_fields(changed_nodes: tuple[str, ...]) -> tuple[torch.Tensor, torch.Tensor]:
        """(the perturbed field, the difference field) of drug_one, with 3 added to the descriptors of changed_nodes."""
        changed = descriptors.clone()
        for node in changed_nodes:
            changed[index_of[node]] += 3.0
        torch.manual_seed(5)
        encoder = RelationalMessagePassingEncoder(
            len(new_nodes), len(relations), 6, 2, node_features=torch.cat([structural, changed], dim=1), num_descriptor_columns=2,
            descriptor_treatment="seed_masked", seed_mask_partner_index=entry_partner_pairs(new_edges, index_of),
            entry_node_mask=torch.tensor(list(new_nodes.node_type == DRUG_NODE_TYPE)), entry_edge_index=edge_index[:, is_entry],
            entry_edge_sign=torch.tensor(list(new_edges.sign[is_entry.numpy()]), dtype=torch.float32),
            entry_relation_indices=(relations.index(DRUG_MECHANISM_RELATION),))
        with torch.no_grad():
            encoder.layer_bias.copy_(torch.randn_like(encoder.layer_bias))
        seeds, sign_and_magnitude = torch.tensor([[index_of["DRUG:drug_one"]]]), torch.tensor([[[1.0, 1.0]]])
        with torch.no_grad():
            return (encoder(seeds, sign_and_magnitude, relation_adjacencies=adjacencies)[0],
                    encoder.perturbation_difference_field(seeds, sign_and_magnitude, adjacencies)[0])

    # the perturbed field holds each node's own base state, so it shows which nodes were read with their descriptors; the
    # difference field cancels a base state wherever the perturbation flips no rectifier
    perturbed, difference = drug_one_fields(())
    perturbed_targets_changed, difference_targets_changed = drug_one_fields(("GENE:A", "GENE:B"))
    torch.testing.assert_close(perturbed_targets_changed, perturbed, atol=1e-6, rtol=0)  # drug_one's targets: masked
    torch.testing.assert_close(difference_targets_changed, difference, atol=1e-6, rtol=0)
    assert float((drug_one_fields(("GENE:C",))[0] - perturbed).abs().max()) > 1e-4  # one edge further on: read with its descriptors
    assert float(difference[index_of["DRUG:drug_one"]].abs().max()) > 0.0  # perturbed pass minus a reference without the node
    assert float(difference[index_of["DRUG:drug_two"]].abs().max()) == 0.0  # absent in both passes


def test_seed_masked_gates_a_drug_s_targets_on_their_structural_columns_under_linear_response() -> None:
    nodes, edges, new_nodes, new_edges, relations, entry_table, _ = variant(carriers=CARRIERS_OF_DRUG_NODES)
    index_of = {node_id: index for index, node_id in enumerate(new_nodes.node_id)}
    structural = typed_features(list(new_nodes.node_type), sorted(set(new_nodes.node_type)))
    descriptors = torch.randn(len(new_nodes), 2, generator=torch.Generator().manual_seed(8))
    descriptors[len(nodes):] = 0.0
    torch.manual_seed(6)
    encoder = LinearResponseEncoder(
        len(new_nodes), relations, torch.tensor([index_of[node] for node in new_edges.source_id]), torch.tensor([index_of[node] for node in new_edges.target_id]),
        torch.tensor([relations.index(name) for name in new_edges.relation_type]), torch.tensor(list(new_edges.sign), dtype=torch.float32),
        torch.cat([structural, descriptors], dim=1), 6, num_propagation_steps=20, propagation_channels=3, num_descriptor_columns=2,
        descriptor_treatment="seed_masked", seed_mask_partner_index=entry_partner_pairs(new_edges, index_of),
        entry_node_mask=torch.tensor(list(new_nodes.node_type == DRUG_NODE_TYPE)), entry_relation_names=ENTRY_RELATIONS)
    with torch.no_grad():
        torch.nn.init.normal_(encoder.output_gate.weight)  # so that a descriptor column moves the gate
    with_descriptors = torch.sigmoid(encoder.output_gate(encoder.node_features))
    without_descriptors = torch.sigmoid(encoder.output_gate(encoder.node_features_without_descriptors))
    gate = encoder.node_gate(torch.tensor([[index_of["DRUG:drug_one"]], [index_of["GENE:B"]]]))
    for row, masked_nodes in ((0, {"DRUG:drug_one", "GENE:A", "GENE:B"}), (1, {"GENE:B"})):
        for node, index in index_of.items():
            torch.testing.assert_close(gate[row, index], (without_descriptors if node in masked_nodes else with_descriptors)[index], atol=1e-7, rtol=0)
    assert float((with_descriptors - without_descriptors)[index_of["GENE:A"]].abs().max()) > 1e-4  # the mask changes something
    field = encoder(torch.tensor([[index_of["DRUG:drug_one"]]]), torch.tensor([[[1.0, 1.0]]]))
    assert field.shape == (1, len(new_nodes), 6) and bool(torch.isfinite(field).all())


def test_a_sequestration_edge_that_carries_change_is_silent_while_the_carrier_is_at_rest() -> None:
    """ALB sequesters drug_one. Carrying the carrier's state, the edge gives drug_one an input that drug_two, which has
    no known carrier, does not get. Carrying the carrier's change from its unperturbed state, it sends nothing when the
    drug is taken alone, sends when the carrier is lowered with the drug, and never reaches a drug that is not taken."""
    _, _, with_presence, adjacencies, index_of, num_source_nodes = message_passing_pair(sequestration_carries="presence")
    _, _, with_change, _, _, _ = message_passing_pair(sequestration_carries="change")
    with torch.no_grad():  # one set of parameters for the two readings of the edge
        with_change.load_state_dict(with_presence.state_dict())
    drug_one, carrier = index_of["DRUG:drug_one"], index_of["GENE:ALB"]
    drug_alone = (torch.tensor([[drug_one, -1]]), torch.tensor([[[1.0, 1.0], [0.0, 0.0]]]))
    drug_with_carrier_lowered = (torch.tensor([[drug_one, carrier]]), torch.tensor([[[1.0, 1.0], [-1.0, 1.0]]]))
    carrier_lowered = (torch.tensor([[carrier, -1]]), torch.tensor([[[-1.0, 1.0], [0.0, 0.0]]]))

    def without_the_edge(seeds: torch.Tensor, sign_and_magnitude: torch.Tensor) -> torch.Tensor:
        """The field with the sequestration relation's maps at zero, which is the graph without the edge."""
        _, _, silenced, _, _, _ = message_passing_pair(sequestration_carries="presence")
        with torch.no_grad():
            silenced.load_state_dict(with_presence.state_dict())
            silenced.relation_weight[:, -1] = 0.0  # sequesters is the relation appended last
        return silenced.perturbation_difference_field(seeds, sign_and_magnitude, adjacencies)

    with torch.no_grad():
        # taken alone: the change-carrying edge is the graph without the edge, exactly, and the state-carrying edge is not
        torch.testing.assert_close(with_change.perturbation_difference_field(*drug_alone, adjacencies), without_the_edge(*drug_alone), atol=1e-6, rtol=0)
        assert float((with_presence.perturbation_difference_field(*drug_alone, adjacencies) - without_the_edge(*drug_alone)).abs().max()) > 1e-3
        # the carrier lowered with the drug: the edge sends, and the drug's targets move with it
        moved = (with_change.perturbation_difference_field(*drug_with_carrier_lowered, adjacencies) - without_the_edge(*drug_with_carrier_lowered)).abs().amax(dim=2)[0]
        assert float(moved[drug_one]) > 1e-4 and float(moved[index_of["GENE:A"]]) > 1e-4
        # the carrier lowered alone: no drug is taken, so no drug node holds anything and the rest is the graph without the edge
        alone = with_change.perturbation_difference_field(*carrier_lowered, adjacencies)
        assert float(alone[0, num_source_nodes:].abs().max()) == 0.0
        torch.testing.assert_close(alone, without_the_edge(*carrier_lowered), atol=1e-6, rtol=0)
        # a batch is its rows: the three perturbations together give the three fields
        together = with_change.perturbation_difference_field(torch.cat([drug_alone[0], drug_with_carrier_lowered[0], carrier_lowered[0]]),
                                                             torch.cat([drug_alone[1], drug_with_carrier_lowered[1], carrier_lowered[1]]), adjacencies)
        torch.testing.assert_close(together[2:3], alone, atol=1e-6, rtol=0)
        torch.testing.assert_close(together[0:1], with_change.perturbation_difference_field(*drug_alone, adjacencies), atol=1e-6, rtol=0)
    with pytest.raises(ValueError, match="change-only edges need both"):
        RelationalMessagePassingEncoder(5, 2, 4, 2, change_only_relation_index=1)


def test_lowering_a_carrier_with_its_drug_raises_the_drug_s_linear_response() -> None:
    """The linear response is a deviation from rest, so the sequestration edge carries the carrier's change by
    construction: with ALB lowered beside drug_one, the part of the response that passes through the drug node is the
    drug's own response times a positive factor per channel, and nothing passes when the drug is not taken."""
    _, encoder, index_of, num_source_nodes = linear_response_pair(num_propagation_steps=400)
    drug_one, carrier = index_of["DRUG:drug_one"], index_of["GENE:ALB"]
    with torch.no_grad():
        drug = encoder.response(torch.tensor([[drug_one]]), torch.tensor([[[1.0, 1.0]]]))
        both = encoder.response(torch.tensor([[drug_one, carrier]]), torch.tensor([[[1.0, 1.0], [-1.0, 1.0]]]))
        carrier_alone = encoder.response(torch.tensor([[carrier]]), torch.tensor([[[-1.0, 1.0]]]))
    assert float(carrier_alone[0, num_source_nodes:].abs().max()) == 0.0  # no drug node without the drug
    through_the_drug = (both - drug - carrier_alone)[0]  # [nodes, channels]
    factor = through_the_drug[drug_one] / drug[0, drug_one]  # one per channel
    assert bool((factor > 1e-3).all())  # less carrier, more of the drug
    for target in ("GENE:A", "GENE:B", "GENE:C"):
        torch.testing.assert_close(through_the_drug[index_of[target]], factor * drug[0, index_of[target]], atol=1e-5, rtol=1e-4)
    with torch.no_grad():  # the carrier raised: the same factor with the other sign
        raised = encoder.response(torch.tensor([[drug_one, carrier]]), torch.tensor([[[1.0, 1.0], [1.0, 1.0]]]))
        raised_alone = encoder.response(torch.tensor([[carrier]]), torch.tensor([[[1.0, 1.0]]]))
    torch.testing.assert_close((raised - drug - raised_alone)[0, drug_one], -factor * drug[0, drug_one], atol=1e-6, rtol=1e-4)

