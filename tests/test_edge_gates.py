"""Edge gates pulled toward evidence priors (models/edge_gates.py, decision 9.57) and how the two encoders read them:
a gate starts at its prior and costs nothing there; the message-passing encoder with gates is the plain encoder run on
the gated adjacencies, and reads no gate on a conjunction relation; in the linear-response encoder a gate is a share at
its target, so equal gates on a node's inputs cancel; with drug entry nodes the gates still index the graph's edges."""
from __future__ import annotations

import numpy as np
import pytest
import torch

from mechanistic_pathway_learning.models.edge_gates import (
    PROBABILITY_FLOOR,
    EdgeGates,
    gated_relation_adjacencies,
    penalty_denominator,
    relation_adjacency_parts,
    sparse_layout,
    sparse_values_product,
)
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder
from mechanistic_pathway_learning.models.soft_constraint_losses import bernoulli_kl_divergence

# six nodes, relations activates (0), inhibits (1) and member_of (2, a conjunction relation); edge 1 repeats edge 0, so
# their two values share one entry of the adjacency; edge 4 has no tier and is ungated
NUM_NODES, NUM_RELATIONS = 6, 3
EDGE_SOURCE = [0, 0, 2, 1, 3, 4, 5, 4]
EDGE_TARGET = [1, 1, 1, 3, 4, 5, 2, 2]
EDGE_RELATION = [0, 0, 0, 1, 0, 1, 2, 2]
EDGE_SIGN = [1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0]
EDGE_WEIGHT = [1.0, 1.0, 0.5, 1.0, 2.0, 1.0, 1.0, 1.0]
GATED = np.array([True, True, True, True, False, True, True, True])
PRIOR = np.array([0.9, 0.6, 0.7, 0.55, np.nan, 0.8, 0.9, 0.9])
COUNT = np.array([1.0, 0.5, 0.5, 1.0, 0.0, 2.0, 1.0, 1.0])


def graph_tensors():
    return (torch.tensor([EDGE_SOURCE, EDGE_TARGET]), torch.tensor(EDGE_RELATION), torch.tensor(EDGE_SIGN), torch.tensor(EDGE_WEIGHT))


def test_a_gate_starts_at_its_prior_costs_nothing_there_and_pays_its_pseudo_count_to_move():
    gates = EdgeGates(GATED, PRIOR, COUNT)
    assert gates.num_edges == 8 and gates.gated_edges.tolist() == [0, 1, 2, 3, 5, 6, 7]
    torch.testing.assert_close(gates.gate(), torch.as_tensor(PRIOR[GATED], dtype=torch.float32), atol=1e-6, rtol=0)
    values = gates.edge_values()
    assert float(values[4]) == 1.0, "an ungated edge always carries its full weight"
    assert float(gates.kl_penalty()) < 1e-9
    with torch.no_grad():
        gates.logit[0] = torch.logit(torch.tensor(0.5))
    expected = float(COUNT[0] * bernoulli_kl_divergence(torch.tensor(0.9), torch.tensor(0.5)))
    assert float(gates.kl_penalty()) == pytest.approx(expected, rel=1e-5), "n KL(prior || gate), with the prior first"
    gates.kl_penalty().backward()
    assert float(gates.logit.grad[0]) == pytest.approx(COUNT[0] * (0.5 - 0.9), abs=1e-5), "d/d logit of n KL(p || sigmoid(logit)) = n (g - p)"
    assert torch.all(gates.logit.grad[1:].abs() < 1e-6)
    per_edge = gates.values_per_edge()
    assert np.isnan(per_edge[4]) and per_edge[0] == pytest.approx(0.5) and per_edge[5] == pytest.approx(0.8)
    clipped = EdgeGates(np.array([True, True]), np.array([1.0, 0.0]), np.array([1.0, 1.0]))
    assert clipped.prior_probability.tolist() == pytest.approx([1.0 - PROBABILITY_FLOOR, PROBABILITY_FLOOR]), "a gate never reaches 0 or 1"
    with pytest.raises(ValueError, match="pseudo-count"):
        EdgeGates(np.array([True]), np.array([0.5]), np.array([-1.0]))


def test_the_summary_counts_only_the_gates_the_encoder_reads():
    gates = EdgeGates(GATED, PRIOR, COUNT)
    gates.set_read_edges(torch.tensor([0, 1, 2, 3, 4, 5]))  # the member_of edges 6 and 7 are not read
    summary = gates.summary()
    assert summary["gated_edges"] == 7 and summary["gated_edges_the_encoder_reads"] == 5
    assert summary["mean_prior"] == pytest.approx(float(np.mean([0.9, 0.6, 0.7, 0.55, 0.8])), abs=1e-6) and summary["mean_absolute_move"] < 1e-6
    assert summary["prior_below_0_5"] == 0 and summary["closed_below_0_5"] == 0


def test_the_penalty_denominator_is_the_fit_rows_total_loss_weight():
    weights = np.array([[1.0, 0.2], [0.2, 0.2], [0.5, 1.0]])
    assert penalty_denominator(weights, np.array([0, 2])) == pytest.approx(2.7)
    assert penalty_denominator(np.zeros((2, 2)), np.array([0])) == 1.0, "never below one"


@pytest.mark.parametrize("signed", [True, False])
def test_the_parts_rebuild_the_message_passing_adjacencies_entry_for_entry(signed):
    edge_index, relation, sign, weight = graph_tensors()
    sign, weight = (sign, weight) if signed else (None, None)
    plain = RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, relation, NUM_NODES, NUM_RELATIONS, edge_sign=sign, edge_weight=weight)
    parts = relation_adjacency_parts(edge_index, relation, NUM_NODES, NUM_RELATIONS, sign, weight)
    rebuilt = gated_relation_adjacencies(parts, torch.ones(len(EDGE_SOURCE)), NUM_NODES)
    for first, second in zip(plain, rebuilt):
        torch.testing.assert_close(first.to_dense(), second.to_dense())
        assert torch.equal(first.indices(), second.indices()), "the same coalesced layout"
    assert parts[0]["slot"].tolist()[:2] == [parts[0]["slot"].tolist()[0]] * 2, "the repeated edge shares its entry"


def message_passing_pair(seed: int = 0):
    torch.manual_seed(seed)
    plain = RelationalMessagePassingEncoder(NUM_NODES, NUM_RELATIONS, 4, 2, conjunction_relation_indices=(2,))
    torch.manual_seed(seed)
    gated = RelationalMessagePassingEncoder(NUM_NODES, NUM_RELATIONS, 4, 2, conjunction_relation_indices=(2,))
    gated.load_state_dict(plain.state_dict())
    edge_index, relation, sign, weight = graph_tensors()
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, relation, NUM_NODES, NUM_RELATIONS, edge_sign=sign, edge_weight=weight)
    gated.attach_edge_gates(EdgeGates(GATED, PRIOR, COUNT), edge_index, relation, adjacencies, edge_sign=sign, edge_weight=weight)
    return plain, gated, adjacencies


PERTURBATION = (torch.tensor([[0], [3]]), torch.tensor([[[-1.0, 1.0]], [[1.0, 0.5]]]))


def test_a_gated_message_passing_encoder_is_the_plain_one_on_the_gated_adjacencies():
    plain, gated, adjacencies = message_passing_pair()
    with torch.no_grad():
        gated.edge_gates.logit.copy_(torch.tensor([2.0, -1.0, 0.3, -2.0, 1.5, 0.0, 0.0]))
    edge_index, relation, sign, weight = graph_tensors()
    parts = relation_adjacency_parts(edge_index, relation, NUM_NODES, NUM_RELATIONS, sign, weight)
    by_hand = gated_relation_adjacencies(parts, gated.edge_gates.edge_values().detach(), NUM_NODES)
    field = gated.perturbation_difference_field(*PERTURBATION, adjacencies)
    torch.testing.assert_close(field, plain.perturbation_difference_field(*PERTURBATION, by_hand), atol=1e-6, rtol=0)
    assert not torch.allclose(field, plain.perturbation_difference_field(*PERTURBATION, adjacencies)), "gates away from 1 change the field"
    field.square().sum().backward()
    gradient = gated.edge_gates.logit.grad  # gated edges 0, 1, 2, 3, 5, 6, 7 in that order
    assert torch.all(gradient[[0, 1, 3, 4]] != 0), "a gate on a stacked relation whose source the perturbations move gets a gradient"
    assert float(gradient[2]) == 0.0, "node 2 moves only at the last layer's output, so edge 2 sends the same message with and without the perturbation"
    assert torch.all(gradient[5:] == 0), "the conjunction relation's members take a minimum, not a weighted mean, and are not gated"
    assert gated.edge_gates.summary()["gated_edges_the_encoder_reads"] == 5


def test_gates_at_one_give_the_plain_message_passing_field_and_the_state_carries_the_gates():
    plain, gated, adjacencies = message_passing_pair()
    with torch.no_grad():
        gated.edge_gates.logit.fill_(30.0)  # sigmoid(30) is 1 in single precision
    torch.testing.assert_close(gated.perturbation_difference_field(*PERTURBATION, adjacencies), plain.perturbation_difference_field(*PERTURBATION, adjacencies),
                               atol=1e-6, rtol=0)
    assert {"edge_gates.logit", "edge_gates.prior_probability", "edge_gates.pseudo_count"} <= set(gated.state_dict())
    assert not any(key.startswith("edge_gates") for key in plain.state_dict()), "an encoder without gates has the parameters it always had"


def test_the_message_passing_gates_refuse_other_adjacencies_other_signs_and_a_second_attachment():
    plain, gated, adjacencies = message_passing_pair()
    edge_index, relation, sign, weight = graph_tensors()
    other = RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, relation, NUM_NODES, NUM_RELATIONS, edge_sign=sign, edge_weight=weight)
    with pytest.raises(ValueError, match="not the ones the edge gates were attached to"):
        gated.perturbation_difference_field(*PERTURBATION, other)
    with pytest.raises(ValueError, match="attached once"):
        gated.attach_edge_gates(EdgeGates(GATED, PRIOR, COUNT), edge_index, relation, adjacencies, edge_sign=sign, edge_weight=weight)
    with pytest.raises(ValueError, match="differs from the static one"):
        plain.attach_edge_gates(EdgeGates(GATED, PRIOR, COUNT), edge_index, relation, adjacencies)  # unsigned parts against signed adjacencies
    with pytest.raises(ValueError, match="edges and the graph"):
        RelationalMessagePassingEncoder(NUM_NODES, NUM_RELATIONS, 4, 2).attach_edge_gates(EdgeGates(GATED[:3], PRIOR[:3], COUNT[:3]), edge_index, relation, adjacencies)


# the linear-response graph: node 1 has two inputs (0 and 2), so its shares can move
LINEAR_RELATIONS = ["activates", "inhibits"]
LINEAR_EDGES = [(0, 1, 0, 1.0), (2, 1, 0, 1.0), (1, 2, 1, -1.0), (1, 3, 0, 1.0), (3, 4, 1, -1.0), (2, 4, 0, 1.0)]


def linear_response_pair(normalisation: str = "total_in_degree", edges=LINEAR_EDGES, relations=LINEAR_RELATIONS, num_nodes: int = 5, gates: EdgeGates | None = None, **keywords):
    tensors = (torch.tensor([edge[0] for edge in edges]), torch.tensor([edge[1] for edge in edges]), torch.tensor([edge[2] for edge in edges]),
               torch.tensor([edge[3] for edge in edges]))
    features = torch.nn.functional.one_hot(torch.arange(num_nodes) % 2, num_classes=2).float()
    common = dict(num_propagation_steps=12, fixed_point_tolerance=None, normalisation=normalisation, **keywords)
    torch.manual_seed(0)
    plain = LinearResponseEncoder(num_nodes, relations, *tensors, features, 4, **common)
    if gates is None:
        gates = EdgeGates(np.ones(len(edges), dtype=bool), np.full(len(edges), 0.7), np.ones(len(edges)))
    torch.manual_seed(0)
    gated = LinearResponseEncoder(num_nodes, relations, *tensors, features, 4, edge_gates=gates, **common)
    missing, unexpected = gated.load_state_dict(plain.state_dict(), strict=False)
    assert not unexpected and all(key.startswith("edge_gates.") for key in missing)
    return plain, gated


@pytest.mark.parametrize("normalisation", ["total_in_degree", "in_degree"])
def test_equal_gates_on_a_nodes_inputs_cancel_in_the_linear_response_share(normalisation):
    plain, gated = linear_response_pair(normalisation)
    seeds, injection = torch.tensor([[0]]), torch.tensor([[[1.0, 1.0]]])
    for logit in (30.0, 0.0, -2.0):
        with torch.no_grad():
            gated.edge_gates.logit.fill_(logit)
        torch.testing.assert_close(gated.response(seeds, injection), plain.response(seeds, injection), atol=1e-6, rtol=0)
    with torch.no_grad():
        gated.edge_gates.logit.copy_(torch.tensor([2.0, -2.0, 0.0, 0.0, 0.0, 0.0]))  # node 1 now listens mostly to node 0
    response = gated.response(seeds, injection)
    assert not torch.allclose(response, plain.response(seeds, injection))
    response.square().sum().backward()
    assert float(gated.edge_gates.logit.grad[0].abs()) > 0 and float(gated.edge_gates.logit.grad[1].abs()) > 0


def test_the_linear_response_gates_refuse_the_spectral_normalisation_and_another_edge_count():
    with pytest.raises(ValueError, match="spectral"):
        linear_response_pair("spectral")
    with pytest.raises(ValueError, match="edges and the graph"):
        linear_response_pair(gates=EdgeGates(np.ones(3, dtype=bool), np.full(3, 0.7), np.ones(3)))


def test_with_drug_entry_nodes_the_linear_response_gates_still_index_the_graphs_edges():
    """The entry edges leave the arrays the adjacency is laid out on; the gates must not shift with them. The drug's
    mechanism edge comes first and its gate is nearly closed, while node 1's two inputs have equal gates: read through
    the right edges, node 1's shares are those of the plain encoder."""
    relations = ["activates", "inhibits", "targets"]
    edges = [(4, 0, 2, -1.0), (0, 1, 0, 1.0), (2, 1, 0, 1.0), (1, 2, 1, -1.0), (1, 3, 0, 1.0)]
    gates = EdgeGates(np.ones(5, dtype=bool), np.array([0.01, 0.99, 0.99, 0.5, 0.5]), np.ones(5))
    entry = dict(entry_node_mask=torch.tensor([False, False, False, False, True]), entry_relation_names=("targets",), edge_weight=torch.tensor([0.5, 1.0, 1.0, 1.0, 1.0]))
    plain, gated = linear_response_pair(edges=edges, relations=relations, num_nodes=5, gates=gates, **entry)
    assert sorted(set(gated.gate_entry_origin.tolist())) == [1, 2, 3, 4], "graph edge indices; the mechanism edge 0 is delivered apart"
    assert gated.edge_gates.read_by_encoder.tolist() == [False, True, True, True, True]
    seeds, injection = torch.tensor([[4]]), torch.tensor([[[1.0, 1.0]]])
    torch.testing.assert_close(gated.response(seeds, injection)[:, :2], plain.response(seeds, injection)[:, :2], atol=1e-6, rtol=0)


def test_the_gated_product_has_the_dense_products_gradients_and_forms_no_dense_matrix():
    """PyTorch's own backward of a sparse product in its values forms grad @ dense.T, [rows, columns]; on the v4 graph
    that asked for 15.8 GB per linear-response round. The gated product gathers each value's gradient instead."""
    torch.manual_seed(0)
    num_rows, num_columns = 7, 5
    matrix = (torch.rand(num_rows, num_columns) < 0.4).double() * torch.randn(num_rows, num_columns, dtype=torch.float64)
    indices = matrix.to_sparse().coalesce().indices()
    values = matrix[indices[0], indices[1]].clone().requires_grad_(True)
    dense = torch.randn(num_columns, 3, dtype=torch.float64, requires_grad=True)
    transposed_indices, order = sparse_layout(indices, num_rows)
    torch.testing.assert_close(sparse_values_product(indices, values, dense, num_rows, transposed_indices, order), matrix @ dense)
    assert torch.autograd.gradcheck(lambda v, d: sparse_values_product(indices, v, d, num_rows, transposed_indices, order), (values, dense))

    size = 100_000  # a dense [size, size] gradient would need 40 GB
    indices = torch.tensor([[0, 5, 99_999], [3, 99_998, 7]])
    values = torch.ones(3, requires_grad=True)
    dense = torch.randn(size, 2, requires_grad=True)
    transposed_indices, order = sparse_layout(indices, size)
    sparse_values_product(indices, values, dense, size, transposed_indices, order).sum().backward()
    assert values.grad.tolist() == pytest.approx(dense.detach()[indices[1]].sum(dim=1).tolist())
    assert dense.grad[indices[1]].tolist() == [[1.0, 1.0]] * 3 and float(dense.grad.abs().sum()) == 6.0

