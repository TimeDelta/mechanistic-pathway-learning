"""Amendment item 8: one set of relation transforms applied at every round, diagonal plus a rank-2 factorisation for
the relations enough training examples reach, with an explicit norm bound (the owner, 12 October 2026: "set the
reduced rank bases with the example share flow through each edge type and 1 repeated layer regardless of depth")."""
import torch

import pytest

from mechanistic_pathway_learning.models.relational_message_passing_encoder import (
    DEFAULT_RELATION_NORM_BOUND,
    PERMITTED_RELATION_FACTOR_RANKS,
    RELATION_PARAMETERISATIONS,
    RelationalMessagePassingEncoder,
)

REGISTERED_RELATIONS = ["encodes", "catalyzed_by", "substrate_of", "activates", "product_of", "binds",
                        "regulates_transcription_of", "inhibits", "member_of"] + [f"other_{index}" for index in range(11)]
FACTORISED = tuple(range(9))  # the nine the registration names


def build(parameterisation="shared_low_rank", layers=3, rank=2, node_state_dim=32, seed=0):
    torch.manual_seed(seed)
    return RelationalMessagePassingEncoder(
        num_graph_nodes=40, num_relation_types=len(REGISTERED_RELATIONS), node_state_dim=node_state_dim,
        num_message_passing_layers=layers, relation_parameterisation=parameterisation,
        factorised_relation_indices=FACTORISED if parameterisation == "shared_low_rank" else None,
        relation_factor_rank=rank)


def relation_parameter_count(encoder):
    held = ("relation_weight", "relation_diagonal", "relation_factor_left", "relation_factor_right")
    return sum(parameter.numel() for name, parameter in encoder.named_parameters() if name in held)


def test_the_parameter_counts_are_the_registered_ones():
    """61,440 against 1,792, which is what item 8 states, and the arithmetic it states them from."""
    assert relation_parameter_count(build("per_layer_dense")) == 61_440
    assert relation_parameter_count(build()) == 1_792 == 20 * 32 + 9 * 2 * 32 * 2
    assert relation_parameter_count(build(rank=4)) == 20 * 32 + 9 * 2 * 32 * 4  # the only permitted step up
    assert PERMITTED_RELATION_FACTOR_RANKS == (2, 4)
    with pytest.raises(ValueError):
        build(rank=3)
    with pytest.raises(ValueError):
        build(parameterisation="something_else")
    assert set(RELATION_PARAMETERISATIONS) == {"per_layer_dense", "shared_low_rank"}


def test_depth_stops_being_a_parameter_count():
    """"1 repeated layer regardless of depth": the point of the change, and what lets the rounds run to a fixed point."""
    assert relation_parameter_count(build(layers=3)) == relation_parameter_count(build(layers=30))
    assert relation_parameter_count(build("per_layer_dense", layers=30)) == 10 * relation_parameter_count(build("per_layer_dense", layers=3))
    shared, per_layer = build(layers=4), build("per_layer_dense", layers=4)
    assert shared.self_weight.shape[0] == 1 and per_layer.self_weight.shape[0] == 4
    assert shared.layer_bias.shape[0] == 1
    assert all(shared.parameter_layer(round_index) == 0 for round_index in range(4))
    assert [per_layer.parameter_layer(round_index) for round_index in range(4)] == [0, 1, 2, 3]


def test_a_relation_without_a_factorisation_is_its_diagonal_alone():
    encoder = build()
    messages = torch.randn(5, 32)
    plain = encoder.relation_transform(15, 0)(messages)  # one of the eleven
    diagonal = DEFAULT_RELATION_NORM_BOUND * torch.tanh(encoder.relation_diagonal[15])
    assert torch.allclose(plain, messages * diagonal)
    factorised = encoder.relation_transform(0, 0)(messages)  # one of the nine
    assert not torch.allclose(factorised, messages * (DEFAULT_RELATION_NORM_BOUND * torch.tanh(encoder.relation_diagonal[0])))


def test_every_transform_stays_inside_the_norm_bound_even_at_large_parameters():
    """The bound is what makes the repeated map a contraction, so the iteration of item 7 has a fixed point to reach.
    It must hold for parameters far outside their initial scale, since training moves them."""
    encoder = build()
    with torch.no_grad():
        encoder.relation_diagonal.mul_(50.0)
        encoder.relation_factor_left.mul_(50.0)
        encoder.relation_factor_right.mul_(50.0)
    identity = torch.eye(32)
    for relation_index in (0, 4, 8, 15, 19):
        matrix = encoder.relation_transform(relation_index, 0)(identity)  # rows are the transform's action
        induced = float(torch.linalg.matrix_norm(matrix, ord=2))
        assert induced <= DEFAULT_RELATION_NORM_BOUND + 1e-5, (relation_index, induced)
    assert encoder.relation_transform_norm_bound() == DEFAULT_RELATION_NORM_BOUND


def test_the_shared_transform_runs_and_carries_gradients_to_both_factors():
    encoder = build(layers=3)
    source, target = torch.randint(0, 40, (120,)), torch.randint(0, 40, (120,))
    relation = torch.randint(0, len(REGISTERED_RELATIONS), (120,))
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(
        torch.stack([source, target]), relation, 40, len(REGISTERED_RELATIONS))
    field = encoder(torch.tensor([[0], [3]]), torch.tensor([[[-1.0, 1.0]], [[1.0, 1.0]]]),
                    relation_adjacencies=adjacencies)
    assert field.shape == (2, 40, 32) and bool(torch.isfinite(field).all())
    field.pow(2).mean().backward()
    assert encoder.relation_diagonal.grad is not None and float(encoder.relation_diagonal.grad.abs().sum()) > 0
    assert encoder.relation_factor_left.grad is not None and float(encoder.relation_factor_left.grad.abs().sum()) > 0
    assert encoder.relation_factor_right.grad is not None and float(encoder.relation_factor_right.grad.abs().sum()) > 0


def test_the_dense_parameterisation_is_untouched():
    """Every run before 13 October 2026 built per_layer_dense, and it must still build the same parameters."""
    encoder = build("per_layer_dense", layers=2)
    assert encoder.relation_weight.shape == (2, 20, 32, 32)
    assert encoder.relation_diagonal is None and encoder.relation_factor_left is None
    assert encoder.parameter_layers == 2


# amendment item 7 for the message-passing encoder, which item 8's shared weights are what made possible
def toy_graph(num_nodes=40, num_edges=160, seed=3):
    torch.manual_seed(seed)
    source, target = torch.randint(0, num_nodes, (num_edges,)), torch.randint(0, num_nodes, (num_edges,))
    relation = torch.randint(0, len(REGISTERED_RELATIONS), (num_edges,))
    return RelationalMessagePassingEncoder.build_relation_adjacencies(
        torch.stack([source, target]), relation, num_nodes, len(REGISTERED_RELATIONS))


def build_settling(tolerance=1e-3, cap=60, layers=3, parameterisation="shared_low_rank"):
    torch.manual_seed(0)
    return RelationalMessagePassingEncoder(
        num_graph_nodes=40, num_relation_types=len(REGISTERED_RELATIONS), node_state_dim=32,
        num_message_passing_layers=layers, relation_parameterisation=parameterisation,
        factorised_relation_indices=FACTORISED if parameterisation == "shared_low_rank" else None,
        fixed_point_tolerance=tolerance, max_message_passing_rounds=cap)


def test_the_rounds_run_to_a_fixed_point_and_report_what_they_took():
    adjacencies = toy_graph()
    seeds, sign_and_magnitude = torch.tensor([[0], [7]]), torch.tensor([[[-1.0, 1.0]], [[1.0, 1.0]]])
    settling = build_settling()
    with torch.no_grad():
        settled = settling(seeds, sign_and_magnitude, relation_adjacencies=adjacencies)
    assert 1 < settling.rounds_taken < 60  # it settles, and before the cap
    assert bool(torch.isfinite(settled).all())
    # another round moves the state by less than the tolerance, which is what settled means
    once_more = build_settling(tolerance=None, layers=settling.rounds_taken + 1)
    once_more.load_state_dict(settling.state_dict())
    with torch.no_grad():
        further = once_more(seeds, sign_and_magnitude, relation_adjacencies=adjacencies)
    change = torch.linalg.vector_norm(further - settled) / torch.linalg.vector_norm(settled).clamp_min(1e-12)
    assert float(change) < 1e-3


def test_a_fixed_round_count_is_unchanged_and_still_reports_its_rounds():
    adjacencies = toy_graph()
    fixed = build_settling(tolerance=None, layers=3)
    with torch.no_grad():
        fixed(torch.tensor([[0]]), torch.tensor([[[-1.0, 1.0]]]), relation_adjacencies=adjacencies)
    assert fixed.rounds_taken == 3


def test_convergence_is_refused_without_shared_weights():
    """With one set of weights per layer the depth is the parameter count, so a run to convergence would be a
    different model at every round. The two amendments depend on each other and the constructor says so."""
    with pytest.raises(ValueError, match="shared_low_rank"):
        build_settling(parameterisation="per_layer_dense")
    with pytest.raises(ValueError):
        build_settling(tolerance=-1.0)
    with pytest.raises(ValueError):
        build_settling(cap=0)


def test_the_cap_ends_a_run_that_does_not_settle():
    encoder = build_settling(tolerance=1e-12, cap=5)  # a tolerance it cannot reach
    with torch.no_grad():
        encoder(torch.tensor([[0]]), torch.tensor([[[-1.0, 1.0]]]), relation_adjacencies=toy_graph())
    assert encoder.rounds_taken == 5
