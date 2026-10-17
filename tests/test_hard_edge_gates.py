"""Hard edge gates (decision 9.62 (e), models/edge_gates.HardConcreteEdgeGates): started at the binary prior
sigmoid(n logit p); sampled as hard concrete in the stochastic phase and read out as exactly 0 or 1 otherwise; a read-out
is the same as deleting the closed edges, in both encoders and through the limiting-step bottleneck; the penalty is the
KL divergence to Bernoulli(prior) in the stochastic phase and -log Bernoulli(z; prior) in the binary one; and the flip
estimates of the binary pass are exact where the loss is linear in the adjacency and see a sole input's deletion,
which the local slope cannot."""
from __future__ import annotations

import numpy as np
import pytest
import torch

from experiments.measure_edge_gate_readings import hard_concrete_masses
from mechanistic_pathway_learning.models.edge_gates import HardConcreteEdgeGates, binary_edge_prior, flip_effect_estimates
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder

# node 1 has two inputs (0 and 2); node 4 has two (3 and 2); node 3's only input is 1 (edge 3), node 2's only input is 1 (edge 2)
LINEAR_RELATIONS = ["activates", "inhibits"]
LINEAR_EDGES = [(0, 1, 0, 1.0), (2, 1, 0, 1.0), (1, 2, 1, -1.0), (1, 3, 0, 1.0), (3, 4, 1, -1.0), (2, 4, 0, 1.0)]
NUM_NODES = 5


def hard_gates(num_edges: int, prior=0.8, count=1.0) -> HardConcreteEdgeGates:
    return HardConcreteEdgeGates(np.ones(num_edges, dtype=bool), np.full(num_edges, prior), np.full(num_edges, count))


def set_read_out(gates: HardConcreteEdgeGates, open_edges: list[bool]) -> None:
    with torch.no_grad():
        gates.log_alpha.copy_(torch.tensor([20.0 if is_open else -20.0 for is_open in open_edges]))
    gates.set_phase("binary")


def linear_response(edges, gates=None, **keywords) -> LinearResponseEncoder:
    tensors = (torch.tensor([edge[0] for edge in edges], dtype=torch.long), torch.tensor([edge[1] for edge in edges], dtype=torch.long),
               torch.tensor([edge[2] for edge in edges], dtype=torch.long), torch.tensor([edge[3] for edge in edges]))
    features = torch.nn.functional.one_hot(torch.arange(NUM_NODES) % 2, num_classes=2).float()
    torch.manual_seed(0)
    return LinearResponseEncoder(NUM_NODES, LINEAR_RELATIONS, *tensors, features, 4, edge_gates=gates, **keywords)


def test_the_gates_start_at_the_binary_prior_and_sample_hard_concrete():
    gates = HardConcreteEdgeGates(np.array([True, True, False]), np.array([0.9, 0.7, np.nan]), np.array([1.0, 0.5, 0.0]))
    assert np.allclose(gates.prior_open.numpy(), binary_edge_prior([0.9, 0.7], [1.0, 0.5]), atol=1e-6)
    assert abs(float(gates.prior_open[1]) - 1.0 / (1.0 + np.exp(-0.5 * np.log(0.7 / 0.3)))) < 1e-6, "n = 0.5 halves the log-odds"
    torch.testing.assert_close(gates.open_probability().detach(), gates.prior_open)
    assert gates.read_out().tolist() == [1.0, 1.0] and gates.edge_values().detach().tolist()[2] == 1.0, "an ungated edge is always open"
    many = hard_gates(20000, prior=0.9)
    many.train()
    sample = many.resample(torch.Generator().manual_seed(0))
    masses = hard_concrete_masses(float(many.log_alpha[0]))
    assert abs(float((sample == 0).float().mean()) - float(masses["closed"])) < 0.01
    assert abs(float((sample == 1).float().mean()) - float(masses["open"])) < 0.01
    many.eval()
    assert set(many.edge_values().tolist()) == {1.0}, "evaluation reads the binary read-out, not a sample"


def test_the_penalty_is_the_kl_to_the_prior_and_in_the_binary_phase_the_cost_of_each_departure():
    gates = HardConcreteEdgeGates(np.ones(3, dtype=bool), np.array([0.9, 0.6, 0.3]), np.ones(3))
    assert float(gates.kl_penalty()) < 1e-6, "zero at the prior"
    with torch.no_grad():
        gates.log_alpha.add_(torch.tensor([-1.0, 0.5, 0.0]))
    q, prior = gates.open_probability().detach(), gates.prior_open
    expected = (q * torch.log(q / prior) + (1 - q) * torch.log((1 - q) / (1 - prior))).sum()
    torch.testing.assert_close(gates.kl_penalty().detach(), expected)
    set_read_out(gates, [False, True, True])  # close the 0.9 edge and open the 0.3 one
    expected_binary = float(torch.log(prior[0] / (1 - prior[0])) + torch.log((1 - prior[2]) / prior[2]))
    assert abs(float(gates.kl_penalty()) - expected_binary) < 1e-5, "logit(0.9) to close a likely edge, logit(0.7) to open an unlikely one"
    change = gates.flip_prior_change()
    assert float(change[0]) < 0 and float(change[1]) > 0, "reopening the closed 0.9 edge pays back its cost; closing the 0.6 edge costs logit(0.6)"


def test_the_read_out_version_follows_the_open_edges_and_is_none_while_a_sample_is_read():
    gates = hard_gates(4)
    gates.train()
    assert gates.binary_state_version() is None
    set_read_out(gates, [True, True, False, True])
    first = gates.binary_state_version()
    set_read_out(gates, [True, False, False, True])
    assert first is not None and gates.binary_state_version() != first


@pytest.mark.parametrize("open_edges", [[True, False, True, True, False, True], [False, True, True, False, True, True], [True] * 6])
def test_a_linear_response_read_out_is_the_graph_without_the_closed_edges(open_edges):
    gates = hard_gates(len(LINEAR_EDGES))
    gated = linear_response(LINEAR_EDGES, gates, num_propagation_steps=12, fixed_point_tolerance=None)
    set_read_out(gates, open_edges)
    pruned = linear_response([edge for edge, is_open in zip(LINEAR_EDGES, open_edges) if is_open], num_propagation_steps=12, fixed_point_tolerance=None)
    missing, unexpected = pruned.load_state_dict({key: value for key, value in gated.state_dict().items() if not key.startswith("edge_gates.")}, strict=False)
    assert not unexpected
    seeds, injection = torch.tensor([[0], [3]]), torch.tensor([[[1.0, 1.0]], [[-1.0, 0.5]]])
    torch.testing.assert_close(gated.response(seeds, injection), pruned.response(seeds, injection), atol=1e-6, rtol=0)


def test_the_limiting_step_bottleneck_walks_only_the_open_edges():
    """Class B expresses node 2 at 0.1 and everything else at 1. Each read-out must give what the graph without its
    closed edges gives, bottleneck included; with edge 0 (0 -> 1) closed the seed at 0 reaches nothing in either class."""
    weights = torch.ones(NUM_NODES, 2)
    weights[2, 1] = 0.1
    gates = hard_gates(len(LINEAR_EDGES))
    common = dict(fixed_point_tolerance=1e-3, cell_class_weights=weights, channels_per_cell_class=1)
    gated = linear_response(LINEAR_EDGES, gates, **common)
    gated.eval()
    seeds, injection = torch.tensor([[0]]), torch.tensor([[[1.0, 1.0]]])
    for open_edges in ([True] * 6, [True, False, True, True, True, False], [False, True, True, True, True, True]):
        set_read_out(gates, open_edges)
        pruned = linear_response([edge for edge, is_open in zip(LINEAR_EDGES, open_edges) if is_open], **common)
        pruned.load_state_dict({key: value for key, value in gated.state_dict().items() if not key.startswith("edge_gates.")}, strict=False)
        pruned.eval()
        torch.testing.assert_close(gated.response(seeds, injection), pruned.response(seeds, injection), atol=1e-6, rtol=0)
    assert float(gated.response(seeds, injection)[0, 1:].abs().max()) == 0.0, "with edge 0 closed the seed reaches nothing"


def message_passing(edges, gates=None):
    edge_index = torch.tensor([[edge[0] for edge in edges], [edge[1] for edge in edges]], dtype=torch.long)
    relation = torch.tensor([edge[2] for edge in edges], dtype=torch.long)
    sign = torch.tensor([edge[3] for edge in edges])
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, relation, NUM_NODES, len(LINEAR_RELATIONS), edge_sign=sign)
    torch.manual_seed(0)
    encoder = RelationalMessagePassingEncoder(NUM_NODES, len(LINEAR_RELATIONS), 4, 2)
    if gates is not None:
        encoder.attach_edge_gates(gates, edge_index, relation, adjacencies, edge_sign=sign)
    return encoder, adjacencies


@pytest.mark.parametrize("open_edges", [[True, False, True, True, False, True], [False, True, False, False, True, True]])
def test_a_message_passing_read_out_is_the_graph_without_the_closed_edges(open_edges):
    gates = hard_gates(len(LINEAR_EDGES))
    gated, adjacencies = message_passing(LINEAR_EDGES, gates)
    set_read_out(gates, open_edges)
    pruned, pruned_adjacencies = message_passing([edge for edge, is_open in zip(LINEAR_EDGES, open_edges) if is_open])
    pruned.load_state_dict({key: value for key, value in gated.state_dict().items() if not key.startswith("edge_gates.")})
    seeds, injection = torch.tensor([[0], [3]]), torch.tensor([[[1.0, 1.0]], [[-1.0, 0.5]]])
    torch.testing.assert_close(gated.perturbation_difference_field(seeds, injection, adjacencies),
                               pruned.perturbation_difference_field(seeds, injection, pruned_adjacencies), atol=1e-6, rtol=0)


def test_the_flip_estimate_is_exact_when_the_loss_is_linear_in_the_adjacency():
    """One round of linear response with damping 1 is h = sum_r g_r S_r u + u, linear in the adjacency's values, and a
    linear loss keeps it so: each estimate must equal the exact change in the loss when that edge alone flips."""
    open_edges = [True, True, True, True, False, True]
    gates = hard_gates(len(LINEAR_EDGES))
    common = dict(num_propagation_steps=1, fixed_point_tolerance=None, damping=1.0)
    encoder = linear_response(LINEAR_EDGES, gates, **common)
    set_read_out(gates, open_edges)
    seeds, injection = torch.tensor([[0], [1], [2]]), torch.tensor([[[1.0, 1.0]], [[1.0, 1.0]], [[-1.0, 1.0]]])
    readout = torch.randn(NUM_NODES, encoder.propagation_channels, generator=torch.Generator().manual_seed(1))

    def loss() -> torch.Tensor:
        return (encoder.response(seeds, injection) * readout).sum()

    encoder.collect_flip_gradients = True
    encoder.reset_flip_collection()
    base = loss()
    base.backward()
    estimates = encoder.flip_effect_estimates().detach()
    encoder.collect_flip_gradients = False
    for edge in range(len(LINEAR_EDGES)):
        flipped = list(open_edges)
        flipped[edge] = not flipped[edge]
        set_read_out(gates, flipped)
        with torch.no_grad():
            exact = float(loss() - base)
        assert abs(float(estimates[edge]) - exact) < 1e-5, (edge, float(estimates[edge]), exact)
        set_read_out(gates, open_edges)


def test_a_sole_inputs_deletion_is_seen_by_the_flip_estimate_and_not_by_the_slope():
    """Edge 3 (1 -> 3) is node 3's only input: its share is 1 at every gate above 0, so dL/dz is 0 there, while closing it
    removes all of node 3's input. The estimate must be the exact change of the closure in a loss linear in v."""
    gates = hard_gates(len(LINEAR_EDGES))
    encoder = linear_response(LINEAR_EDGES, gates, num_propagation_steps=1, fixed_point_tolerance=None, damping=1.0)
    gates.set_phase("stochastic")
    gates.train()
    with torch.no_grad():
        gates.log_alpha.fill_(20.0)  # every sample at exactly 1
    seeds, injection = torch.tensor([[1]]), torch.tensor([[[1.0, 1.0]]])
    gates.resample(torch.Generator().manual_seed(0))
    encoder.response(seeds, injection)[0, 3].sum().backward()
    assert float(gates.log_alpha.grad[3].abs()) < 1e-6, "the local slope of a sole input is zero"
    set_read_out(gates, [True] * 6)
    encoder.collect_flip_gradients = True
    encoder.reset_flip_collection()
    encoder.response(seeds, injection)[0, 3].sum().backward()
    estimate = float(encoder.flip_effect_estimates()[3])
    with torch.no_grad():
        before = float(encoder.response(seeds, injection)[0, 3].sum())
        set_read_out(gates, [True, True, True, False, True, True])
        after = float(encoder.response(seeds, injection)[0, 3].sum())
    assert abs(estimate - (after - before)) < 1e-5 and abs(estimate) > 1e-3


def test_the_flip_estimates_of_the_function_on_a_hand_built_case():
    """Two groups: group 0 holds edges 0 and 1 (weights 1 and 3), group 1 holds edge 2 alone; edge 1 is closed.
    Values: v0 = a0 / 1, v1 = 0, v2 = a2 / 1. Closing edge 0 removes group 0's only open input; opening edge 1 makes
    v0 = a0 / 4 and v1 = a1 / 4; closing edge 2 empties group 1."""
    group, origin = torch.tensor([0, 0, 1]), torch.tensor([0, 1, 2])
    numerator, weight = torch.tensor([2.0, -1.0, 5.0]), torch.tensor([1.0, 3.0, 1.0])
    gradient = torch.tensor([0.5, 0.25, -2.0])
    state = torch.tensor([True, False, True])
    estimate = flip_effect_estimates(group, origin, numerator, weight, gradient, state, 3)
    expected_close_0 = -0.5 * 2.0  # v0 from 2 to 0
    expected_open_1 = 0.5 * (2.0 / 4.0 - 2.0) + 0.25 * (-1.0 / 4.0)
    expected_close_2 = -(-2.0 * 5.0)  # v2 from 5 to 0: -g v
    torch.testing.assert_close(estimate, torch.tensor([expected_close_0, expected_open_1, expected_close_2]))


def test_the_proxy_gradient_pushes_each_gate_toward_the_side_that_lowers_the_objective():
    gates = hard_gates(3)
    set_read_out(gates, [True, True, False])
    with torch.no_grad():
        gates.log_alpha.copy_(torch.tensor([0.3, 0.3, -0.3]) + gates.beta * float(np.log(-gates.gamma / gates.zeta)))  # near the threshold
    gates.log_alpha.grad = None
    # closing edge 0 would lower the objective by 1; closing edge 1 would raise it by 1; opening edge 2 would lower it
    gates.apply_flip_gradient(torch.tensor([-1.0, 1.0, -1.0]), surrogate_width=0.5)
    gradient = gates.log_alpha.grad
    assert float(gradient[0]) > 0 and float(gradient[1]) < 0 and float(gradient[2]) < 0, "descent closes 0, keeps 1 open and opens 2"


def test_the_estimated_closure_evidence_is_the_exact_one_when_the_loss_is_linear_in_the_adjacency():
    """Decision 9.63 (c) lists every closure; the first-order estimate ranks them in one pass. Where the estimate is exact
    (one round of linear response, a linear loss), it must give each closed edge the exact closure evidence and P(false)."""
    from mechanistic_pathway_learning.models.edge_gate_passes import closure_evidence, estimated_closure_evidence

    open_edges = [True, True, False, True, False, True]
    gates = hard_gates(len(LINEAR_EDGES), prior=0.7)
    encoder = linear_response(LINEAR_EDGES, gates, num_propagation_steps=1, fixed_point_tolerance=None, damping=1.0)
    set_read_out(gates, open_edges)
    seeds, injection = torch.tensor([[0], [1], [2]]), torch.tensor([[[1.0, 1.0]], [[1.0, 1.0]], [[-1.0, 1.0]]])
    readout = torch.randn(NUM_NODES, encoder.propagation_channels, generator=torch.Generator().manual_seed(1))

    def data_loss(batch) -> torch.Tensor:
        return (encoder.response(seeds, injection) * readout).sum()

    estimated = estimated_closure_evidence(encoder, gates, data_loss, [None], penalty_denominator=10.0)
    exact = {record["edge"]: record for record in closure_evidence(encoder, gates, data_loss, [None], penalty_denominator=10.0)["records"]}
    assert set(estimated) == set(exact) == {2, 4}, "every closed edge the encoder reads, keyed by graph edge"
    for edge, record in exact.items():
        assert estimated[edge]["estimated_evidence_nats"] == pytest.approx(record["evidence_nats"], abs=1e-4)
        assert estimated[edge]["estimated_posterior_false"] == pytest.approx(record["posterior_false"], abs=1e-5)
    set_read_out(gates, [True] * len(LINEAR_EDGES))
    assert estimated_closure_evidence(encoder, gates, data_loss, [None], penalty_denominator=10.0) == {}, "no closure, nothing to list"


def test_a_drugs_mechanism_edge_stays_a_route_of_the_bottleneck_under_hard_gates():
    """The bottleneck's open entries are the gated adjacency's open entries followed by the drug's mechanism edges,
    which carry no gate; closing the gate of the edge past the target cuts the route there and nowhere else."""
    gates = HardConcreteEdgeGates(np.array([False, True]), np.array([np.nan, 0.8]), np.array([0.0, 1.0]))
    torch.manual_seed(0)
    encoder = LinearResponseEncoder(
        3, ["activates", "targets"], edge_source=torch.tensor([0, 1]), edge_target=torch.tensor([1, 2]),
        edge_relation=torch.tensor([1, 0]), edge_sign=torch.tensor([1.0, 1.0]), node_features=torch.eye(3), node_state_dim=4,
        cell_class_weights=torch.tensor([[1.0], [0.6], [0.3]]), class_weight_rule="limiting_step",
        entry_node_mask=torch.tensor([True, False, False]), entry_relation_names=("targets",), edge_weight=torch.tensor([1.0, 1.0]), edge_gates=gates)
    set_read_out(gates, [True])
    with torch.no_grad():
        assert encoder.bottleneck_open_entries().tolist() == [True, True]
        assert encoder.class_bottleneck(torch.tensor([[0]]))[:, 0, 0].tolist() == pytest.approx([1.0, 0.6, 0.3])
        set_read_out(gates, [False])
        assert encoder.bottleneck_open_entries().tolist() == [False, True], "the mechanism edge has no gate to close"
        assert encoder.class_bottleneck(torch.tensor([[0]]))[:, 0, 0].tolist() == pytest.approx([1.0, 0.6, 0.0])
