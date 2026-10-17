"""Tests for the linear-response encoder: stoichiometric signs, linearity, reach per step and convergence."""
import pytest
import torch

from mechanistic_pathway_learning.models.linear_response_encoder import (
    MIXTURE_STATISTICS,
    MIXTURE_WEIGHTINGS,
    NORMALISATIONS,
    CrossRelationAggregator,
    LinearResponseEncoder,
    project_onto_simplex,
)

# gene 0 catalyses reaction 1, which turns substrate 2 into product 3; product 3 is the substrate of reaction 4,
# which makes product 5
RELATION_TYPES = ["substrate_of", "product_of", "catalyzed_by"]
EDGES = [  # (source, target, relation, sign)
    (0, 1, "catalyzed_by", 1.0),
    (2, 1, "substrate_of", 1.0),
    (1, 3, "product_of", 1.0),
    (3, 4, "substrate_of", 1.0),
    (4, 5, "product_of", 1.0),
]
NODE_TYPES = [0, 1, 2, 2, 1, 2]  # gene, reaction, metabolite


def build_encoder(edges=EDGES, relation_types=RELATION_TYPES, num_nodes=6, node_types=NODE_TYPES, steps=8) -> LinearResponseEncoder:
    torch.manual_seed(0)
    node_features = torch.nn.functional.one_hot(torch.tensor(node_types), num_classes=max(node_types) + 1).float()
    return LinearResponseEncoder(
        num_nodes, relation_types,
        edge_source=torch.tensor([edge[0] for edge in edges]), edge_target=torch.tensor([edge[1] for edge in edges]),
        edge_relation=torch.tensor([relation_types.index(edge[2]) for edge in edges]), edge_sign=torch.tensor([edge[3] for edge in edges]),
        node_features=node_features, node_state_dim=4, num_propagation_steps=steps,
    )


def perturb(node: int, sign: float) -> tuple[torch.Tensor, torch.Tensor]:
    return torch.tensor([[node]]), torch.tensor([[[sign, 1.0]]])


def test_loss_of_function_lowers_the_product_and_raises_the_substrate() -> None:
    encoder = build_encoder()
    with torch.no_grad():
        response = encoder.response(*perturb(0, -1.0))[0]  # [nodes, channels]
    # every channel carries the input with its own scale, so compare signs within a channel
    assert torch.all(response[3] * response[1] > 0)  # product moves with its reaction's flux
    assert torch.all(response[2] * response[1] < 0)  # substrate moves against it (depletes_substrate, sign -1)
    assert torch.all(response[0] * response[1] > 0)  # the sustained input at the gene drives its reaction


def test_response_is_linear_and_zero_without_input() -> None:
    encoder = build_encoder()
    with torch.no_grad():
        loss_of_function = encoder.response(*perturb(0, -1.0))
        gain_of_function = encoder.response(*perturb(0, 1.0))
        padding_only = encoder.response(torch.tensor([[-1]]), torch.tensor([[[1.0, 1.0]]]))
    assert torch.allclose(loss_of_function, -gain_of_function)
    assert torch.count_nonzero(padding_only) == 0


def test_reach_grows_by_one_edge_per_step() -> None:
    encoder = build_encoder()
    with torch.no_grad():
        two_steps = encoder.response(*perturb(0, -1.0), num_steps=2)[0]
        three_steps = encoder.response(*perturb(0, -1.0), num_steps=3)[0]
    assert torch.count_nonzero(two_steps[4]) == 0  # reaction 4 is three edges from the gene
    assert torch.count_nonzero(three_steps[4]) > 0


def test_iteration_converges_to_a_fixed_point() -> None:
    encoder = build_encoder()
    with torch.no_grad():
        many_steps = encoder.response(*perturb(0, -1.0), num_steps=80)
        one_more = encoder.response(*perturb(0, -1.0), num_steps=81)
    assert torch.max(torch.abs(many_steps - one_more)) < 1e-5 * torch.max(torch.abs(many_steps))


def test_repressing_transcription_edge_inverts_the_target_response() -> None:
    relation_types = ["regulates_transcription_of"]
    activating = build_encoder(edges=[(0, 1, "regulates_transcription_of", 1.0)], relation_types=relation_types, num_nodes=2, node_types=[0, 0])
    repressing = build_encoder(edges=[(0, 1, "regulates_transcription_of", -1.0)], relation_types=relation_types, num_nodes=2, node_types=[0, 0])
    with torch.no_grad():
        activated = activating.response(*perturb(0, 1.0))[0]
        repressed = repressing.response(*perturb(0, 1.0))[0]
    assert torch.all(activated[1] * activated[0] > 0)
    assert torch.all(repressed[1] * repressed[0] < 0)


def test_gains_bound_the_transition_and_gradients_reach_every_parameter() -> None:
    encoder = build_encoder()
    assert float(encoder.relation_gain().detach().abs().max()) < encoder.maximum_gain + 1e-9
    field = encoder.perturbation_difference_field(*perturb(0, -1.0))
    field.pow(2).sum().backward()
    for name, parameter in encoder.named_parameters():
        assert parameter.grad is not None and torch.any(parameter.grad != 0), name


def test_signed_log_scale_keeps_sign_and_odd_symmetry_and_lifts_small_changes() -> None:
    torch.manual_seed(0)
    node_features = torch.nn.functional.one_hot(torch.tensor(NODE_TYPES), num_classes=3).float()
    arguments = dict(edge_source=torch.tensor([edge[0] for edge in EDGES]), edge_target=torch.tensor([edge[1] for edge in EDGES]),
                     edge_relation=torch.tensor([RELATION_TYPES.index(edge[2]) for edge in EDGES]), edge_sign=torch.tensor([edge[3] for edge in EDGES]),
                     node_features=node_features, node_state_dim=4)
    encoder = LinearResponseEncoder(6, RELATION_TYPES, response_scale="signed_log", **arguments)
    with torch.no_grad():
        loss_of_function = encoder(*perturb(0, -1.0))
        gain_of_function = encoder(*perturb(0, 1.0))
        raw = encoder.response(*perturb(0, -1.0))[0]
        scaled = torch.sign(raw) * torch.log1p(raw.abs() / torch.exp(encoder.log_response_scale))
    assert torch.allclose(loss_of_function, -gain_of_function, atol=1e-6)
    assert torch.equal(torch.sign(scaled), torch.sign(raw))
    nonzero = raw.abs() > 0
    smallest, largest = raw.abs()[nonzero].min(), raw.abs()[nonzero].max()
    assert scaled.abs()[nonzero].min() / scaled.abs()[nonzero].max() > smallest / largest  # the scale compresses the range


def _chain_encoder(normalisation, steps=64):
    """A chain 0 -> 1 -> 2 -> 3 -> 4 whose every link destination also has three distractor sources."""
    chain = [(0, 1), (1, 2), (2, 3), (3, 4)]
    edge_list, distractor = list(chain), 5
    for _, target in chain:
        for _ in range(3):
            edge_list.append((distractor, target))
            distractor += 1
    source = torch.tensor([edge[0] for edge in edge_list])
    target = torch.tensor([edge[1] for edge in edge_list])
    torch.manual_seed(0)
    return LinearResponseEncoder(
        distractor, ["activates"], source, target,
        torch.zeros(len(edge_list), dtype=torch.long), torch.ones(len(edge_list)),
        node_features=torch.eye(distractor), node_state_dim=2, num_propagation_steps=steps,
        propagation_channels=1, normalisation=normalisation,
    )


@pytest.mark.parametrize("normalisation", NORMALISATIONS)
def test_every_normalisation_reaches_a_fixed_point(normalisation):
    """What the encoder needs is a fixed point, which needs the spectral radius below one. A row sum of sum_r |S_r|
    below one is sufficient for that but not necessary, and the spectral arm bounds the radius directly while leaving
    rows above one, so the property to assert for every arm is that the field stops moving as the steps grow."""
    with torch.no_grad():
        at_32 = _chain_encoder(normalisation, steps=32)(torch.tensor([[0]]), torch.tensor([[[1.0, 1.0]]]))
        at_128 = _chain_encoder(normalisation, steps=128)(torch.tensor([[0]]), torch.tensor([[[1.0, 1.0]]]))
    assert torch.isfinite(at_32).all() and torch.isfinite(at_128).all()
    assert torch.allclose(at_32, at_128, atol=1e-6), normalisation


@pytest.mark.parametrize("normalisation", ["in_degree", "total_in_degree"])
def test_the_per_node_arms_keep_every_row_sum_below_one(normalisation):
    """Sufficient for the contraction, though not necessary: the spectral arm legitimately exceeds it and is
    asserted on its fixed point instead."""
    encoder = _chain_encoder(normalisation)
    unsigned_aggregate = encoder.stacked_adjacency.to_dense().abs().view(len(encoder.relation_names), encoder.num_graph_nodes, -1).sum(dim=0)
    largest_row_sum = float(unsigned_aggregate.sum(dim=1).max())
    assert largest_row_sum <= 1.0 + 1e-6, largest_row_sum


def _metabolite_row_sum(normalisation: str, num_producing: int, num_consuming: int, edge_signs: str = "graph") -> float:
    """Signed row sum, over every relation, at metabolite 0, made by num_producing reactions and used by num_consuming."""
    producing = list(range(1, num_producing + 1))
    consuming = list(range(num_producing + 1, num_producing + num_consuming + 1))
    edges = [(reaction, 0, "product_of") for reaction in producing] + [(0, reaction, "substrate_of") for reaction in consuming]
    num_nodes = 1 + num_producing + num_consuming
    stacked, relation_names = LinearResponseEncoder.signed_stacked_adjacency(
        num_nodes, RELATION_TYPES, torch.tensor([edge[0] for edge in edges]), torch.tensor([edge[1] for edge in edges]),
        torch.tensor([RELATION_TYPES.index(edge[2]) for edge in edges]), torch.ones(len(edges)), normalisation=normalisation,
        edge_signs=edge_signs)
    return float(stacked.to_dense().view(len(relation_names), num_nodes, num_nodes).sum(dim=0)[0].sum())


def test_in_degree_erases_the_stoichiometric_counts_and_total_in_degree_keeps_them():
    """Five producing reactions against one consuming reaction: averaging inside each relation leaves
    mean(+1) / 2 + mean(-1) / 2 = 0 whatever the counts, while one divisor over all relations leaves (5 - 1) / 6."""
    assert _metabolite_row_sum("in_degree", 5, 1) == pytest.approx(0.0, abs=1e-6)
    assert _metabolite_row_sum("in_degree", 1, 5) == pytest.approx(0.0, abs=1e-6)
    assert _metabolite_row_sum("total_in_degree", 5, 1) == pytest.approx(4 / 6)
    assert _metabolite_row_sum("total_in_degree", 1, 5) == pytest.approx(-4 / 6)
    assert _metabolite_row_sum("total_in_degree", 2, 2) == pytest.approx(0.0, abs=1e-6)  # equal counts still cancel, which is correct


def test_the_sign_ablation_turns_the_depletion_edge_positive_so_nothing_cancels():
    """One producing and one consuming reaction: the depletion edge's -1 cancels the product's +1 under the graph's
    signs, and with every sign +1 the two relation means add to one."""
    assert _metabolite_row_sum("in_degree", 1, 1) == pytest.approx(0.0, abs=1e-6)
    assert _metabolite_row_sum("in_degree", 1, 1, edge_signs="all_positive") == pytest.approx(1.0)
    encoder = LinearResponseEncoder(
        6, RELATION_TYPES, torch.tensor([edge[0] for edge in EDGES]), torch.tensor([edge[1] for edge in EDGES]),
        torch.tensor([RELATION_TYPES.index(edge[2]) for edge in EDGES]), torch.tensor([edge[3] for edge in EDGES]),
        torch.zeros(6, 3), node_state_dim=4, edge_signs="all_positive")
    assert bool((encoder.stacked_adjacency.values() > 0).all())
    assert not bool(encoder.relation_is_unsigned.any())


def test_an_unknown_edge_sign_rule_is_refused():
    with pytest.raises(ValueError, match="edge_signs must be one of"):
        _metabolite_row_sum("in_degree", 1, 1, edge_signs="learned")


def test_an_unknown_normalisation_is_refused():
    with pytest.raises(ValueError, match="normalisation must be one of"):
        _chain_encoder("row_stochastic")


def build_mixed_sign_encoder(cross_relation_aggregator, mixture_statistics=None, dominant_statistic=None):
    """One metabolite fed by a producing reaction and a consuming one, so it receives both signs.

    Node 0 is the metabolite, nodes 1 and 2 are reactions, node 3 is the perturbed gene. The
    substrate_of edge 0 -> 2 makes the encoder derive depletes_substrate 2 -> 0 with sign -1, while
    product_of 1 -> 0 carries +1, so node 0 is exactly the case the cross-relation mean annihilates.
    """
    num_nodes = 4
    relation_types = ["substrate_of", "product_of", "catalyzed_by"]
    edge_source = torch.tensor([0, 1, 3, 3])
    edge_target = torch.tensor([2, 0, 1, 2])
    edge_relation = torch.tensor([0, 1, 2, 2])
    edge_sign = torch.tensor([1.0, 1.0, 1.0, 1.0])
    node_type_index = torch.tensor([0, 1, 1, 2])
    torch.manual_seed(0)
    encoder = LinearResponseEncoder(num_nodes, relation_types, edge_source, edge_target, edge_relation, edge_sign,
                                    torch.eye(num_nodes), node_state_dim=3, propagation_channels=1,
                                    num_propagation_steps=6, cross_relation_aggregator=cross_relation_aggregator,
                                    node_type_index=node_type_index)
    if dominant_statistic is not None:
        with torch.no_grad():
            statistics = encoder.cross_relation_mixture.statistics
            encoder.cross_relation_mixture.mixture_logit.fill_(-8.0)
            encoder.cross_relation_mixture.mixture_logit[:, statistics.index(dominant_statistic), :] = 8.0
    return encoder


def test_zero_input_gives_exactly_zero_field_under_both_aggregators():
    """The difference-field property: the mixture is not linear, but zero still maps to zero, so the
    unperturbed field is zero and the field is the difference."""
    for aggregator in ("mean", "softmax_mixture"):
        encoder = build_mixed_sign_encoder(aggregator)
        response = encoder.response(torch.tensor([[3]]), torch.tensor([[[0.0, 1.0]]]))
        assert float(response.abs().max()) == 0.0


def test_mean_cancels_opposite_messages_where_the_signed_maximum_magnitude_keeps_one():
    """The aggregator's own behaviour on the case the real graph presents at 88% of metabolites.

    This is a statement about the operator, not about the field: a node whose incoming signs cancel has
    an aggregate row sum of zero, so it is blind to any uniform change upstream and passes on only the
    difference between its production and consumption inputs. The field there is not zero in the
    dynamics, because those inputs carry different upstream states; what is zero is the response to the
    part of the input they share.
    """
    aggregator = CrossRelationAggregator(num_node_types=1, propagation_channels=1)
    opposite_messages = torch.tensor([[[[1.0]]], [[[-1.0]]]])  # [R=2, N=1, B=1, C=1]
    both_relations_feed = torch.ones(2, 1, dtype=torch.bool)
    node_type_index = torch.zeros(1, dtype=torch.long)
    statistics = aggregator.statistics

    with torch.no_grad():
        aggregator.mixture_logit.fill_(-8.0)
        aggregator.mixture_logit[:, statistics.index("mean"), :] = 8.0
    assert abs(float(aggregator(opposite_messages, both_relations_feed, node_type_index))) < 1e-3

    with torch.no_grad():
        aggregator.mixture_logit.fill_(-8.0)
        aggregator.mixture_logit[:, statistics.index("signed_maximum_magnitude"), :] = 8.0
    assert float(aggregator(opposite_messages, both_relations_feed, node_type_index)) > 0.9


def test_a_relation_that_does_not_feed_a_node_is_left_out_of_the_statistics():
    """Otherwise a relation with no edge into the node contributes a zero that the median reads as a
    real message, pulling the median of one positive message towards zero."""
    aggregator = CrossRelationAggregator(num_node_types=1, propagation_channels=1)
    with torch.no_grad():
        aggregator.mixture_logit.fill_(-8.0)
        aggregator.mixture_logit[:, aggregator.statistics.index("median"), :] = 8.0
    one_message_and_one_absent_relation = torch.tensor([[[[0.6]]], [[[0.0]]]])
    only_the_first_feeds = torch.tensor([[True], [False]])
    node_type_index = torch.zeros(1, dtype=torch.long)
    median_over_feeding_only = float(aggregator(one_message_and_one_absent_relation, only_the_first_feeds, node_type_index))
    assert abs(median_over_feeding_only - 0.6) < 1e-3, f"expected 0.6, got {median_over_feeding_only}"


def test_the_mixture_carries_more_field_than_the_mean_at_a_mixed_sign_node():
    """In the dynamics the cancellation is partial rather than exact, so the claim is an inequality."""
    mean_encoder = build_mixed_sign_encoder("mean")
    maximum_encoder = build_mixed_sign_encoder("softmax_mixture", dominant_statistic="signed_maximum_magnitude")
    for encoder in (mean_encoder, maximum_encoder):
        with torch.no_grad():
            encoder.gain_logit.fill_(0.0)  # equal gains, so the two signs arrive with equal weight
    perturbation_index, perturbation_value = torch.tensor([[3]]), torch.tensor([[[1.0, 1.0]]])
    mean_field = abs(float(mean_encoder.response(perturbation_index, perturbation_value)[0, 0, 0]))
    maximum_field = abs(float(maximum_encoder.response(perturbation_index, perturbation_value)[0, 0, 0]))
    assert maximum_field > mean_field, f"maximum {maximum_field} did not exceed mean {mean_field}"


def test_mixture_stays_within_the_contraction_bound():
    """Each statistic is bounded by the largest per-relation message and softmax weights are convex,
    so the response cannot grow without bound however many steps are taken."""
    for dominant_statistic in MIXTURE_STATISTICS:
        encoder = build_mixed_sign_encoder("softmax_mixture", dominant_statistic=dominant_statistic)
        perturbation_index, perturbation_value = torch.tensor([[3]]), torch.tensor([[[1.0, 1.0]]])
        short_run = float(encoder.response(perturbation_index, perturbation_value, num_steps=4).abs().max())
        long_run = float(encoder.response(perturbation_index, perturbation_value, num_steps=200).abs().max())
        assert long_run < 1.0 / (1.0 - encoder.maximum_gain), f"{dominant_statistic} grew to {long_run}"
        assert abs(long_run - float(encoder.response(perturbation_index, perturbation_value, num_steps=240).abs().max())) < 1e-5, (
            f"{dominant_statistic} had not converged by 200 steps")
        assert short_run > 0.0


def test_mixture_rejects_an_unknown_statistic_and_a_mismatched_node_type_vector():
    with pytest.raises(ValueError, match="unknown mixture statistics"):
        CrossRelationAggregator(num_node_types=2, propagation_channels=1, statistics=("mean", "skewness"))
    with pytest.raises(ValueError, match="at least one mixture statistic"):
        CrossRelationAggregator(num_node_types=2, propagation_channels=1, statistics=())
    with pytest.raises(ValueError, match="node_type_index must have one entry per graph node"):
        LinearResponseEncoder(3, ["substrate_of"], torch.tensor([0]), torch.tensor([1]), torch.tensor([0]),
                              torch.tensor([1.0]), torch.eye(3), node_state_dim=2,
                              cross_relation_aggregator="softmax_mixture", node_type_index=torch.tensor([0, 1]))


def test_unknown_cross_relation_aggregator_is_rejected():
    with pytest.raises(ValueError, match="cross_relation_aggregator must be one of"):
        LinearResponseEncoder(3, ["substrate_of"], torch.tensor([0]), torch.tensor([1]), torch.tensor([0]),
                              torch.tensor([1.0]), torch.eye(3), node_state_dim=2,
                              cross_relation_aggregator="median")


def test_every_statistic_is_odd_so_the_response_mirrors_under_a_sign_flip():
    """The property that ruled out the plain maximum, minimum and standard deviation: a gain of
    function must produce the negative of the field a loss of function produces."""
    for dominant_statistic in MIXTURE_STATISTICS:
        encoder = build_mixed_sign_encoder("softmax_mixture", dominant_statistic=dominant_statistic)
        perturbation_index = torch.tensor([[3]])
        rising = encoder.response(perturbation_index, torch.tensor([[[1.0, 1.0]]]))
        falling = encoder.response(perturbation_index, torch.tensor([[[-1.0, 1.0]]]))
        largest_departure = float((rising + falling).abs().max())
        assert largest_departure < 1e-6, f"{dominant_statistic} is not odd, mirror error {largest_departure}"
        assert float(rising.abs().max()) > 1e-9, f"{dominant_statistic} produced no field at all"


def test_sparsemax_gives_exact_zeros_where_softmax_leaves_a_residue():
    """Why sparsemax is the default: the statistics a node type does not use can be read off the
    weights rather than inferred from a small share."""
    logit = torch.tensor([[[2.0], [0.5], [0.4], [0.3], [0.2]]])  # [1, statistics, channels]
    sparse_weight = project_onto_simplex(logit, dim=1)
    dense_weight = torch.softmax(logit, dim=1)
    assert abs(float(sparse_weight.sum()) - 1.0) < 1e-6
    assert float(sparse_weight.min()) >= 0.0
    assert int((sparse_weight == 0).sum()) == 4
    assert int((dense_weight == 0).sum()) == 0


def test_sparsemax_leaves_a_uniform_slice_uniform():
    """The projection of an already-uniform point is itself, so an untrained mixture starts spread."""
    uniform_weight = project_onto_simplex(torch.zeros(2, 4, 3), dim=1)
    assert torch.allclose(uniform_weight, torch.full((2, 4, 3), 0.25), atol=1e-6)


def test_power_mean_spans_from_the_mean_to_the_dominant_message():
    """One learned exponent per node type covers the span a discrete mixture only samples."""
    aggregator = CrossRelationAggregator(num_node_types=1, propagation_channels=1, statistics=("power_mean",))
    messages = torch.tensor([[[[0.9]]], [[[-0.3]]]])  # [relations, nodes, batch, channels]
    both_feed = torch.ones(2, 1, dtype=torch.bool)
    node_type_index = torch.zeros(1, dtype=torch.long)

    with torch.no_grad():
        aggregator.log_power_mean_exponent.fill_(0.0)  # exponent 1, the plain mean
    at_the_mean = float(aggregator(messages, both_feed, node_type_index))
    assert abs(at_the_mean - 0.3) < 1e-4, at_the_mean

    with torch.no_grad():
        aggregator.log_power_mean_exponent.fill_(4.0)  # a large exponent, so the dominant message
    at_the_dominant = float(aggregator(messages, both_feed, node_type_index))
    assert abs(at_the_dominant - 0.9) < 1e-3, at_the_dominant


def test_the_median_equals_the_mean_at_two_messages():
    """Recorded as a test because it is why the statistic list is short: no destination on
    data/processed/graph is fed by three or more relations, and 73.1% are fed by exactly two, so on
    that graph the median carries nothing the mean does not."""
    median_only = CrossRelationAggregator(num_node_types=1, propagation_channels=1, statistics=("median",))
    mean_only = CrossRelationAggregator(num_node_types=1, propagation_channels=1, statistics=("mean",))
    messages = torch.tensor([[[[0.7]]], [[[-0.2]]]])
    both_feed = torch.ones(2, 1, dtype=torch.bool)
    node_type_index = torch.zeros(1, dtype=torch.long)
    assert abs(float(median_only(messages, both_feed, node_type_index))
               - float(mean_only(messages, both_feed, node_type_index))) < 1e-6


def test_signed_standard_deviation_keeps_the_magnitude_when_one_side_dominates():
    """At two messages of opposite sign it is sign(m1 + m2) * (|m1| + |m2|) / 2, so it keeps the whole
    magnitude where the mean keeps only the difference, and takes its direction from the dominant side."""
    aggregator = CrossRelationAggregator(num_node_types=1, propagation_channels=1,
                                         statistics=("signed_standard_deviation",))
    both_feed = torch.ones(2, 1, dtype=torch.bool)
    node_type_index = torch.zeros(1, dtype=torch.long)

    production_wins = float(aggregator(torch.tensor([[[[0.6]]], [[[-0.4]]]]), both_feed, node_type_index))
    consumption_wins = float(aggregator(torch.tensor([[[[0.4]]], [[[-0.6]]]]), both_feed, node_type_index))
    assert abs(production_wins - 0.5) < 1e-4, production_wins     # (0.6 + 0.4) / 2, positive
    assert abs(consumption_wins + 0.5) < 1e-4, consumption_wins   # the same magnitude, negative


PERMUTATION_INVARIANT_STATISTICS = tuple(statistic for statistic in MIXTURE_STATISTICS
                                         if statistic != "signed_maximum_magnitude")


@pytest.mark.parametrize("statistic", PERMUTATION_INVARIANT_STATISTICS)
def test_no_symmetric_odd_statistic_can_rescue_an_exact_cancellation(statistic):
    """An impossibility, recorded because a test written on the opposite assumption failed.

    Every statistic here is odd and permutation invariant. Negating a message multiset that is already
    symmetric under negation, such as {+0.5, -0.5}, returns the same multiset, so f = -f and therefore
    f = 0. No choice of symmetric statistic escapes this, and the mixture is a convex combination of
    them, so it cannot either. Only something that distinguishes the relations can, which the
    per-relation gains already do; what they lack is freedom per node, and that is the normalisation's
    business rather than the aggregator's (task 25, the spectral arm).

    The mixture is still worth having: an exactly symmetric multiset is the measure-zero case, and at
    the 73.1% of destinations fed by two relations the generic case is a mean that is small rather than
    zero, where the dispersion and the dominant message both carry more than the mean does.

    signed_maximum_magnitude is excluded because it is the one escape, and the way it escapes is worth
    knowing rather than relying on: it is odd but not permutation invariant, since argmax resolves the
    tie between two equal magnitudes by relation index. That index is the relation's position in the
    stack, which is an arbitrary convention, so the value it returns at an exactly cancelling node is
    an artifact of ordering rather than a mechanism. test_signed_maximum_magnitude_breaks_a_tie_by
    _relation_order records that behaviour so a change to the ordering is caught.
    """
    aggregator = CrossRelationAggregator(num_node_types=1, propagation_channels=1, statistics=(statistic,))
    symmetric_messages = torch.tensor([[[[0.5]]], [[[-0.5]]]])
    both_feed = torch.ones(2, 1, dtype=torch.bool)
    node_type_index = torch.zeros(1, dtype=torch.long)
    assert abs(float(aggregator(symmetric_messages, both_feed, node_type_index))) < 1e-6


def test_an_unknown_mixture_weighting_is_refused():
    with pytest.raises(ValueError, match="weighting must be one of"):
        CrossRelationAggregator(num_node_types=1, propagation_channels=1, weighting="argmax")


def test_signed_maximum_magnitude_breaks_a_tie_by_relation_order():
    """The one escape from the impossibility above, and the reason not to lean on it: at two equal
    magnitudes of opposite sign the statistic returns whichever relation comes first in the stack, so
    its value at an exactly cancelling node follows the ordering convention rather than the biology."""
    aggregator = CrossRelationAggregator(num_node_types=1, propagation_channels=1,
                                         statistics=("signed_maximum_magnitude",))
    both_feed = torch.ones(2, 1, dtype=torch.bool)
    node_type_index = torch.zeros(1, dtype=torch.long)
    positive_relation_first = float(aggregator(torch.tensor([[[[0.5]]], [[[-0.5]]]]), both_feed, node_type_index))
    negative_relation_first = float(aggregator(torch.tensor([[[[-0.5]]], [[[0.5]]]]), both_feed, node_type_index))
    assert positive_relation_first == 0.5
    assert negative_relation_first == -0.5


def test_median_is_taken_per_batch_row_and_channel_not_broadcast_from_the_first():
    """Regression: gather does not broadcast its index, so a [N, 1, 1] median index silently returned batch
    row 0 of channel 0 for every row. Every aggregator test used one batch row and one channel, which is the
    only shape where that is correct."""
    num_relations, num_nodes, batch_size, num_channels = 3, 4, 5, 2
    aggregator = CrossRelationAggregator(num_node_types=1, propagation_channels=num_channels, statistics=("median",))
    torch.manual_seed(0)
    relation_messages = torch.randn(num_relations, num_nodes, batch_size, num_channels)
    relation_feeds_node = torch.ones(num_relations, num_nodes, dtype=torch.bool)
    with torch.no_grad():
        aggregated = aggregator(relation_messages, relation_feeds_node, torch.zeros(num_nodes, dtype=torch.long))
    assert aggregated.shape == (num_nodes, batch_size, num_channels)
    expected_median = relation_messages.median(dim=0).values
    assert torch.allclose(aggregated, expected_median, atol=1e-6)
    # and the rows genuinely differ, so the comparison above is not vacuous
    assert not torch.allclose(expected_median[:, 0, :], expected_median[:, 1, :], atol=1e-3)


@pytest.mark.parametrize("weighting", MIXTURE_WEIGHTINGS)
def test_every_statistic_keeps_a_gradient_at_the_warm_start(weighting):
    """Regression: sparsemax has zero gradient outside its support, so a warm start that gave the other
    statistics exactly zero weight froze the whole mixture at the mean for the duration of training."""
    num_relations, num_nodes, batch_size, num_channels = 4, 6, 3, 2
    aggregator = CrossRelationAggregator(num_node_types=1, propagation_channels=num_channels, weighting=weighting)
    torch.manual_seed(0)
    relation_messages = torch.randn(num_relations, num_nodes, batch_size, num_channels)
    relation_feeds_node = torch.ones(num_relations, num_nodes, dtype=torch.bool)
    aggregator(relation_messages, relation_feeds_node, torch.zeros(num_nodes, dtype=torch.long)).pow(2).sum().backward()
    gradient_per_statistic = aggregator.mixture_logit.grad[0, :, 0]
    assert gradient_per_statistic.abs().min() > 0.0, f"a statistic cannot be learned under {weighting}"


@pytest.mark.parametrize("weighting", MIXTURE_WEIGHTINGS)
def test_the_warm_start_puts_the_requested_weight_on_the_mean_under_either_weighting(weighting):
    """Both weightings must start from the same mixture, or a comparison between them is confounded."""
    requested_mean_weight = 0.9
    aggregator = CrossRelationAggregator(num_node_types=2, propagation_channels=3, weighting=weighting,
                                         initial_mean_weight=requested_mean_weight)
    weight = (torch.softmax(aggregator.mixture_logit, dim=1) if weighting == "softmax"
              else project_onto_simplex(aggregator.mixture_logit, dim=1))
    mean_position = MIXTURE_STATISTICS.index("mean")
    assert torch.allclose(weight[:, mean_position, :], torch.full_like(weight[:, mean_position, :], requested_mean_weight), atol=1e-6)
    remaining = torch.cat([weight[:, :mean_position, :], weight[:, mean_position + 1:, :]], dim=1)
    expected_remaining = (1.0 - requested_mean_weight) / (len(MIXTURE_STATISTICS) - 1)
    assert torch.allclose(remaining, torch.full_like(remaining, expected_remaining), atol=1e-6)


def test_the_warm_start_rejects_a_weight_that_would_leave_a_statistic_without_a_gradient():
    with pytest.raises(ValueError, match="initial_mean_weight"):
        CrossRelationAggregator(num_node_types=1, propagation_channels=1, initial_mean_weight=1.0)
    with pytest.raises(ValueError, match="initial_mean_weight"):
        CrossRelationAggregator(num_node_types=1, propagation_channels=1, initial_mean_weight=0.1)


def test_a_shared_gain_gives_every_relation_the_same_gain():
    encoder = LinearResponseEncoder(
        6, RELATION_TYPES, torch.tensor([edge[0] for edge in EDGES]), torch.tensor([edge[1] for edge in EDGES]),
        torch.tensor([RELATION_TYPES.index(edge[2]) for edge in EDGES]), torch.tensor([edge[3] for edge in EDGES]),
        torch.zeros(6, 3), node_state_dim=4, relation_gains="shared")
    gain = encoder.relation_gain()
    assert encoder.gain_logit.shape[0] == 1
    assert gain.shape[0] == len(encoder.relation_names)
    assert torch.allclose(gain, gain[:1].expand_as(gain))


def test_an_unknown_relation_gain_rule_is_refused():
    with pytest.raises(ValueError, match="relation_gains must be one of"):
        LinearResponseEncoder(
            6, RELATION_TYPES, torch.tensor([edge[0] for edge in EDGES]), torch.tensor([edge[1] for edge in EDGES]),
            torch.tensor([RELATION_TYPES.index(edge[2]) for edge in EDGES]), torch.tensor([edge[3] for edge in EDGES]),
            torch.zeros(6, 3), node_state_dim=4, relation_gains="typed")


def test_the_sign_permutation_keeps_the_number_of_negative_entries_and_is_reproducible():
    def stacked_values(edge_signs):
        stacked, _ = LinearResponseEncoder.signed_stacked_adjacency(
            6, RELATION_TYPES, torch.tensor([edge[0] for edge in EDGES]), torch.tensor([edge[1] for edge in EDGES]),
            torch.tensor([RELATION_TYPES.index(edge[2]) for edge in EDGES]), torch.tensor([edge[3] for edge in EDGES]),
            edge_signs=edge_signs)
        return stacked.values()

    graph_signs, permuted_signs = stacked_values("graph"), stacked_values("permuted")
    assert int((permuted_signs < 0).sum()) == int((graph_signs < 0).sum())
    assert torch.equal(permuted_signs, stacked_values("permuted"))


def build_cell_class_encoder(cell_class_weights, channels_per_cell_class=1, extracellular_pool_nodes=None, shared_pool_classes=None, steps=8):
    torch.manual_seed(0)
    node_features = torch.nn.functional.one_hot(torch.tensor(NODE_TYPES), num_classes=3).float()
    return LinearResponseEncoder(
        6, RELATION_TYPES,
        edge_source=torch.tensor([edge[0] for edge in EDGES]), edge_target=torch.tensor([edge[1] for edge in EDGES]),
        edge_relation=torch.tensor([RELATION_TYPES.index(edge[2]) for edge in EDGES]), edge_sign=torch.tensor([edge[3] for edge in EDGES]),
        node_features=node_features, node_state_dim=4, num_propagation_steps=steps, cell_class_weights=cell_class_weights,
        channels_per_cell_class=channels_per_cell_class, extracellular_pool_nodes=extracellular_pool_nodes, shared_pool_classes=shared_pool_classes,
    )


def test_one_all_ones_cell_class_reproduces_the_encoder_without_classes():
    plain = build_encoder()  # four propagation channels by default
    with_one_class = build_cell_class_encoder(torch.ones(6, 1), channels_per_cell_class=4)
    with torch.no_grad():
        assert torch.allclose(plain.response(*perturb(0, -1.0)), with_one_class.response(*perturb(0, -1.0)))
        assert torch.allclose(plain(*perturb(0, -1.0)), with_one_class(*perturb(0, -1.0)))


def test_a_class_that_does_not_express_the_reaction_does_not_carry_the_response_past_it():
    weights = torch.ones(6, 2)
    weights[1, 1] = 0.0  # class 1 does not express reaction 1
    encoder = build_cell_class_encoder(weights)
    with torch.no_grad():
        response = encoder.response(*perturb(0, -1.0))[0]  # [nodes, classes * 1 channel]
    assert torch.all(response[[1, 2, 3, 4, 5], 0] != 0)  # the all-expressing class carries it to the end of the chain
    assert torch.count_nonzero(response[[1, 2, 3, 4, 5], 1]) == 0  # the other stops at the gene
    assert response[0, 1] != 0


def test_partial_expression_scales_the_response_down_and_keeps_its_sign():
    full, half = torch.ones(6, 1), torch.ones(6, 1)
    half[1, 0] = 0.5
    with torch.no_grad():
        full_response = build_cell_class_encoder(full).response(*perturb(0, -1.0))[0, :, 0]
        half_response = build_cell_class_encoder(half).response(*perturb(0, -1.0))[0, :, 0]
    downstream = [1, 2, 3, 4, 5]
    assert torch.all(half_response[downstream].abs() < full_response[downstream].abs())
    assert torch.all(torch.sign(half_response[downstream]) == torch.sign(full_response[downstream]))


def test_the_extracellular_pool_carries_a_product_into_a_class_that_cannot_make_it():
    weights = torch.ones(6, 2)
    weights[1, 1] = 0.0  # class 1 cannot run reaction 1, but uses product 3 through reaction 4
    pool = torch.tensor([False, False, False, True, False, False])
    with torch.no_grad():
        isolated = build_cell_class_encoder(weights).response(*perturb(0, -1.0))[0]
        pooled = build_cell_class_encoder(weights, extracellular_pool_nodes=pool).response(*perturb(0, -1.0))[0]
    assert isolated[3, 1] == 0 and isolated[4, 1] == 0
    assert pooled[3, 1] != 0 and pooled[4, 1] != 0
    assert torch.allclose(pooled[3, 0], pooled[3, 1])  # one shared pool


def test_cell_class_weights_outside_the_unit_interval_and_a_pool_without_classes_are_refused():
    with pytest.raises(ValueError):
        build_cell_class_encoder(torch.full((6, 2), 1.5))
    with pytest.raises(ValueError):
        build_cell_class_encoder(torch.ones(6, 2), shared_pool_classes=torch.tensor([True, False]), extracellular_pool_nodes=torch.zeros(6, dtype=torch.bool))


def test_edge_wise_propagation_equals_the_stacked_product_with_gradients():
    """The mean-across-relations path sums over edges (edge_*_index buffers); it must equal the stacked sparse product,
    with and without cell-class channels and the extracellular pool, in value and in the gradients of gains and inputs."""
    torch.manual_seed(3)
    num_nodes, relations = 40, ["activates", "inhibits", "substrate_of", "product_of", "binds"]
    source, target = torch.randint(0, num_nodes, (300,)), torch.randint(0, num_nodes, (300,))
    relation = torch.randint(0, len(relations), (300,))
    sign = torch.where(relation == 1, -1.0, 1.0)
    for cell_classes in (None, torch.rand(num_nodes, 3)):
        pool = None if cell_classes is None else torch.rand(num_nodes) < 0.2
        # the reference below is the product rule at a fixed step count, so the encoder is built in that regime
        encoder = LinearResponseEncoder(num_nodes, relations, source, target, relation, sign, torch.randn(num_nodes, 4), 8,
                                        cell_class_weights=cell_classes, extracellular_pool_nodes=pool, channels_per_cell_class=2,
                                        fixed_point_tolerance=None, class_weight_rule="product")
        node_index, sign_and_magnitude = torch.randint(0, num_nodes, (4, 2)), torch.randn(4, 2, 2)

        def stacked_product_response():
            sustained = encoder.sustained_input(node_index, sign_and_magnitude)
            gain = encoder.relation_gain()
            node_weight = None if encoder.channel_node_weight is None else encoder.channel_node_weight[:, None, :]
            node_major_input = sustained.permute(1, 0, 2).contiguous() * (1.0 if node_weight is None else node_weight)
            state = node_major_input
            for _ in range(encoder.num_propagation_steps):
                aggregated = torch.sparse.mm(encoder.stacked_adjacency, state.reshape(num_nodes, -1))
                messages = (aggregated.view(gain.shape[0], num_nodes, 4, encoder.propagation_channels) * gain[:, None, None, :]).sum(dim=0)
                messages = messages if node_weight is None else messages * node_weight
                state = (1.0 - encoder.damping) * state + encoder.damping * (messages + node_major_input)
                if encoder.extracellular_pool_nodes is not None:
                    state = encoder.share_extracellular_pool(state)
            return state.permute(1, 0, 2)

        parameters = [parameter for parameter in encoder.parameters() if parameter.requires_grad]
        edge_wise = encoder.response(node_index, sign_and_magnitude)
        edge_wise_gradients = torch.autograd.grad((edge_wise ** 2).sum(), parameters, allow_unused=True)
        reference = stacked_product_response()
        reference_gradients = torch.autograd.grad((reference ** 2).sum(), parameters, allow_unused=True)
        assert torch.allclose(edge_wise, reference, atol=1e-6, rtol=1e-5)
        for edge_wise_gradient, reference_gradient in zip(edge_wise_gradients, reference_gradients):
            if reference_gradient is not None:
                assert torch.allclose(edge_wise_gradient, reference_gradient, atol=1e-5, rtol=1e-4)


def test_the_propagation_runs_to_a_fixed_point_and_reports_the_rounds_it_took():
    """Amendment item 7 of 12 October 2026: both encoders iterate until the relative change falls below 1e-3, in place
    of the fixed eight steps (the owner: "8 steps is arbitrary, not physiological")."""
    from mechanistic_pathway_learning.models.linear_response_encoder import (
        DEFAULT_FIXED_POINT_TOLERANCE,
        LinearResponseEncoder,
    )

    def encoder_with(**overrides):
        torch.manual_seed(0)
        return LinearResponseEncoder(
            num_graph_nodes=6, relation_types=["activates"],
            edge_source=torch.tensor([0, 1, 2, 3, 4]), edge_target=torch.tensor([1, 2, 3, 4, 5]),
            edge_relation=torch.zeros(5, dtype=torch.long), edge_sign=torch.ones(5),
            node_features=torch.eye(6), node_state_dim=4, **overrides)

    node_index = torch.tensor([[0]])
    sign_and_magnitude = torch.tensor([[[1.0, 1.0]]])
    fixed = encoder_with(num_propagation_steps=8, fixed_point_tolerance=None)
    with torch.no_grad():
        fixed.response(node_index, sign_and_magnitude)
    assert fixed.propagation_steps_taken == 8  # a fixed count when the tolerance is turned off explicitly
    assert encoder_with().fixed_point_tolerance == DEFAULT_FIXED_POINT_TOLERANCE  # the fixed point is the default

    settling = encoder_with(fixed_point_tolerance=DEFAULT_FIXED_POINT_TOLERANCE, max_propagation_steps=60)
    with torch.no_grad():
        settled = settling.response(node_index, sign_and_magnitude)
    assert 1 < settling.propagation_steps_taken < 60  # it settles, and does so before the cap
    # another round changes the state by less than the tolerance, which is what settled means
    with torch.no_grad():
        one_more = settling.response(node_index, sign_and_magnitude, num_steps=settling.propagation_steps_taken + 1)
    change = torch.linalg.vector_norm(one_more - settled) / torch.linalg.vector_norm(settled)
    assert float(change) < DEFAULT_FIXED_POINT_TOLERANCE

    # a cap of one round stops after one, so a graph that does not settle still terminates
    capped = encoder_with(fixed_point_tolerance=DEFAULT_FIXED_POINT_TOLERANCE, max_propagation_steps=1)
    with torch.no_grad():
        capped.response(node_index, sign_and_magnitude)
    assert capped.propagation_steps_taken == 1
    # num_steps overrides both, for the measurements that read a chosen round
    with torch.no_grad():
        settling.response(node_index, sign_and_magnitude, num_steps=3)
    assert settling.propagation_steps_taken == 3
    with pytest.raises(ValueError, match="fixed_point_tolerance must be positive"):
        encoder_with(fixed_point_tolerance=0.0)
    with pytest.raises(ValueError, match="max_propagation_steps must be at least 1"):
        encoder_with(fixed_point_tolerance=1e-3, max_propagation_steps=0)


# the limiting-step cell-class rule (the owner, 13 October 2026: "min along the path makes WAY more sense to me")
DIAMOND_RELATIONS = ["activates"]
DIAMOND_SOURCE, DIAMOND_TARGET = [0, 1, 0, 2, 3], [1, 3, 2, 3, 4]  # 0 -> 1 -> 3 and 0 -> 2 -> 3, then 3 -> 4


def build_diamond_encoder(cell_class_weights, class_weight_rule="limiting_step", channels_per_cell_class=1, steps=8):
    torch.manual_seed(0)
    num_nodes = 5
    return LinearResponseEncoder(
        num_nodes, DIAMOND_RELATIONS,
        edge_source=torch.tensor(DIAMOND_SOURCE), edge_target=torch.tensor(DIAMOND_TARGET),
        edge_relation=torch.zeros(len(DIAMOND_SOURCE), dtype=torch.long), edge_sign=torch.ones(len(DIAMOND_SOURCE)),
        node_features=torch.eye(num_nodes), node_state_dim=4, num_propagation_steps=steps,
        cell_class_weights=cell_class_weights, channels_per_cell_class=channels_per_cell_class,
        class_weight_rule=class_weight_rule,
    )


def test_the_bottleneck_takes_the_widest_route_and_not_the_product_along_it():
    """b(t) = max over routes of min over nodes on the route. The thin route through node 1 at 0.1 must not set the
    ceiling at node 3 when the route through node 2 at 0.4 exists, and node 4 must inherit 0.4 rather than 0.4 * 0.5."""
    weights = torch.tensor([[1.0], [0.1], [0.4], [1.0], [0.5]])
    encoder = build_diamond_encoder(weights)
    with torch.no_grad():
        bottleneck = encoder.class_bottleneck(torch.tensor([[0]]), rounds=8)[:, 0, 0]
    assert torch.allclose(bottleneck, torch.tensor([1.0, 0.1, 0.4, 0.4, 0.4]))


def test_the_bottleneck_counts_no_route_longer_than_the_rounds_it_is_given():
    """Node 4 is three edges from the seed, so two rounds must leave it at zero: the attenuation may not be set by a
    route the propagation never travelled."""
    encoder = build_diamond_encoder(torch.tensor([[1.0], [0.1], [0.4], [1.0], [0.5]]))
    with torch.no_grad():
        after_two = encoder.class_bottleneck(torch.tensor([[0]]), rounds=2)[:, 0, 0]
        after_three = encoder.class_bottleneck(torch.tensor([[0]]), rounds=3)[:, 0, 0]
    assert after_two[4] == 0.0 and after_two[3] == pytest.approx(0.4)
    assert after_three[4] == pytest.approx(0.4)


def test_a_stoichiometric_demand_divides_the_bottleneck_past_the_step_that_consumes_it():
    """The owner, 13 October 2026: "ACh needing to bind two molecules to the receptor [...] would cut that expression
    in half". A demand of two on the wide route leaves it at 0.2, which is still above the thin route's 0.1."""
    encoder = build_diamond_encoder(torch.tensor([[1.0], [0.1], [0.4], [1.0], [0.5]]))
    wide_route_edge = [index for index, (source, target) in enumerate(zip(encoder.bottleneck_edge_source.tolist(),
                                                                         encoder.bottleneck_edge_target.tolist()))
                       if (source, target) == (2, 3)]
    demand = torch.ones(encoder.bottleneck_edge_source.shape[0])
    demand[wide_route_edge] = 2.0
    encoder.bottleneck_edge_demand = demand
    with torch.no_grad():
        bottleneck = encoder.class_bottleneck(torch.tensor([[0]]), rounds=8)[:, 0, 0]
    assert torch.allclose(bottleneck, torch.tensor([1.0, 0.1, 0.4, 0.2, 0.2]))


def test_the_bottleneck_of_two_seeds_is_the_better_of_their_own():
    encoder = build_diamond_encoder(torch.tensor([[1.0], [0.1], [0.4], [1.0], [0.5]]))
    with torch.no_grad():
        from_the_thin_node = encoder.class_bottleneck(torch.tensor([[1]]), rounds=8)[:, 0, 0]
        from_the_wide_node = encoder.class_bottleneck(torch.tensor([[2]]), rounds=8)[:, 0, 0]
        from_both = encoder.class_bottleneck(torch.tensor([[1, 2]]), rounds=8)[:, 0, 0]
    assert torch.allclose(from_both, torch.maximum(from_the_thin_node, from_the_wide_node))


def test_the_bottleneck_does_not_depend_on_how_the_gather_is_chunked():
    weights = torch.rand(5, 3)
    encoder = build_diamond_encoder(weights)
    with torch.no_grad():
        one_column_at_a_time = encoder.class_bottleneck(torch.tensor([[0], [2]]), rounds=8, chunk=1)
        all_at_once = encoder.class_bottleneck(torch.tensor([[0], [2]]), rounds=8, chunk=1024)
    assert torch.equal(one_column_at_a_time, all_at_once)


def test_the_bottleneck_is_per_class_and_does_not_change_with_the_channels_per_class():
    weights = torch.rand(5, 2)
    single = build_diamond_encoder(weights, channels_per_cell_class=1)
    tripled = build_diamond_encoder(weights, channels_per_cell_class=3)
    with torch.no_grad():
        assert torch.equal(single.class_bottleneck(torch.tensor([[0]]), rounds=8),
                           tripled.class_bottleneck(torch.tensor([[0]]), rounds=8))
    assert single.cell_class_node_weight.shape == (5, 2) and tripled.channel_node_weight.shape == (5, 6)


def test_the_limiting_step_rule_attenuates_by_the_bottleneck_and_not_by_the_product():
    """The two rules must differ where the product compounds: along the wide route the product charges 0.4 * 1.0 * 0.5
    at node 4 while the limiting step charges 0.4 once, so the limiting step's response there is the larger."""
    weights = torch.tensor([[1.0], [0.1], [0.4], [1.0], [0.5]])
    with torch.no_grad():
        unweighted = build_diamond_encoder(torch.ones(5, 1), class_weight_rule="limiting_step").response(*perturb(0, 1.0))[0, :, 0]
        limiting = build_diamond_encoder(weights, class_weight_rule="limiting_step").response(*perturb(0, 1.0))[0, :, 0]
        product = build_diamond_encoder(weights, class_weight_rule="product").response(*perturb(0, 1.0))[0, :, 0]
    expected_bottleneck = torch.tensor([1.0, 0.1, 0.4, 0.4, 0.4])
    assert torch.allclose(limiting, unweighted * expected_bottleneck, atol=1e-6)
    assert abs(limiting[4]) > abs(product[4])  # the compounding is what the rule removes


def test_the_limiting_step_rule_keeps_a_blocked_node_blocked_when_no_route_avoids_it():
    """Node 4 has one route in, through node 3. A class that does not express node 3 must hold nothing at node 4,
    even though the attenuation is applied after the propagation rather than inside it."""
    weights = torch.ones(5, 2)
    weights[3, 1] = 0.0
    with torch.no_grad():
        response = build_diamond_encoder(weights, class_weight_rule="limiting_step").response(*perturb(0, 1.0))[0]
    assert response[4, 0] != 0 and response[3, 0] != 0
    assert response[4, 1] == 0 and response[3, 1] == 0
    assert response[1, 1] != 0 and response[2, 1] != 0  # upstream of the block the class still answers


def test_the_limiting_step_rule_stays_linear_and_zero_without_input():
    encoder = build_diamond_encoder(torch.rand(5, 2))
    with torch.no_grad():
        raised = encoder.response(*perturb(0, 1.0))
        lowered = encoder.response(*perturb(0, -1.0))
        padding_only = encoder.response(torch.tensor([[-1]]), torch.tensor([[[1.0, 1.0]]]))
        doubled = encoder.response(torch.tensor([[0]]), torch.tensor([[[1.0, 2.0]]]))
    assert torch.allclose(raised, -lowered)
    assert torch.count_nonzero(padding_only) == 0
    assert torch.allclose(doubled, 2.0 * raised, atol=1e-6)


def test_the_limiting_step_rule_leaves_a_run_without_cell_classes_alone():
    """Without cell-class weights there is nothing to take a minimum over, so the rule must change no arithmetic."""
    torch.manual_seed(0)
    plain = build_encoder()
    torch.manual_seed(0)
    named = build_encoder()
    named.class_weight_rule = "limiting_step"
    with torch.no_grad():
        assert torch.equal(plain.response(*perturb(0, -1.0)), named.response(*perturb(0, -1.0)))


def test_an_unknown_class_weight_rule_is_refused_and_the_bottleneck_needs_the_weights():
    with pytest.raises(ValueError):
        build_diamond_encoder(torch.ones(5, 1), class_weight_rule="minimum")
    with pytest.raises(ValueError):
        build_encoder().class_bottleneck(torch.tensor([[0]]), rounds=4)


def test_the_bottleneck_shares_the_extracellular_pool_across_the_pooled_classes_only():
    """Class 1 cannot run reaction 1, so its only route to product 3 is through the shared pool. Without sharing, the
    rule would set its ceiling at the pool to zero and undo --extracellular-coupling; the all-cells class (the last,
    not pooled) keeps its own value."""
    weights = torch.ones(6, 3)
    weights[1, 1] = 0.0  # class 1 cannot run reaction 1
    weights[1, 2] = 0.0  # neither can class 2, which is not pooled
    pool = torch.tensor([False, False, False, True, False, False])
    encoder = build_cell_class_encoder(weights, extracellular_pool_nodes=pool, shared_pool_classes=torch.tensor([True, True, False]))
    with torch.no_grad():
        bottleneck = encoder.class_bottleneck(torch.tensor([[0]]))[:, 0, :]
    assert bottleneck[3, 0] == 1.0 and bottleneck[3, 1] == 1.0  # shared
    assert bottleneck[4, 1] == 1.0  # and it continues out of the pool in class 1
    assert bottleneck[3, 2] == 0.0  # the unpooled class has no route


def test_the_limiting_step_rule_keeps_the_pool_carrying_a_product_into_a_class_that_cannot_make_it():
    weights = torch.ones(6, 2)
    weights[1, 1] = 0.0
    pool = torch.tensor([False, False, False, True, False, False])
    encoder = build_cell_class_encoder(weights, extracellular_pool_nodes=pool)
    assert encoder.class_weight_rule == "limiting_step"  # the default since 13 October 2026
    with torch.no_grad():
        pooled = encoder.response(*perturb(0, -1.0))[0]
    assert pooled[3, 1] != 0 and pooled[4, 1] != 0
    assert torch.allclose(pooled[3, 0], pooled[3, 1])


def test_an_entry_node_the_perturbation_does_not_seed_carries_no_route_in_the_bottleneck():
    """Node 2 is an entry node standing between the seed and node 3: absent unless seeded, so the only route to node 3
    is the thin one through node 1."""
    weights = torch.tensor([[1.0], [0.1], [1.0], [1.0], [1.0]])
    torch.manual_seed(0)
    encoder = LinearResponseEncoder(
        5, DIAMOND_RELATIONS, edge_source=torch.tensor(DIAMOND_SOURCE), edge_target=torch.tensor(DIAMOND_TARGET),
        edge_relation=torch.zeros(len(DIAMOND_SOURCE), dtype=torch.long), edge_sign=torch.ones(len(DIAMOND_SOURCE)),
        node_features=torch.eye(5), node_state_dim=4, cell_class_weights=weights,
        entry_node_mask=torch.tensor([False, False, True, False, False]))
    with torch.no_grad():
        without_the_entry = encoder.class_bottleneck(torch.tensor([[0]]))[:, 0, 0]
        with_the_entry = encoder.class_bottleneck(torch.tensor([[0, 2]]))[:, 0, 0]
    assert without_the_entry[2] == 0.0 and without_the_entry[3] == pytest.approx(0.1)
    assert with_the_entry[3] == pytest.approx(1.0)


def test_the_steady_state_bottleneck_counts_every_route_and_does_not_depend_on_the_batch():
    """At the fixed point every route contributes, so a route longer than the round the propagation stopped at still
    sets the ceiling, and a perturbation's bottleneck is the same alone or batched with others."""
    weights = torch.rand(5, 2)
    encoder = build_diamond_encoder(weights)
    with torch.no_grad():
        settled = encoder.class_bottleneck(torch.tensor([[0]]))
        assert torch.equal(settled, encoder.class_bottleneck(torch.tensor([[0]]), rounds=50))
        alone = encoder.cached_class_bottleneck(torch.tensor([[0]]))
        encoder.clear_class_bottleneck_cache()
        batched = encoder.cached_class_bottleneck(torch.tensor([[2, -1], [0, -1], [1, 3]]))
    assert torch.equal(alone[:, 0], batched[:, 1])
    assert torch.equal(alone, settled)


def test_the_cache_is_reused_and_emptied_when_the_demands_change():
    encoder = build_diamond_encoder(torch.tensor([[1.0], [0.1], [0.4], [1.0], [0.5]]))
    with torch.no_grad():
        first = encoder.response(*perturb(0, 1.0))
        assert len(encoder.class_bottleneck_cache) == 1
        again = encoder.response(*perturb(0, 1.0))
        assert torch.equal(first, again) and len(encoder.class_bottleneck_cache) == 1
        demand = torch.ones(encoder.bottleneck_edge_source.shape[0])
        demand[[index for index, pair in enumerate(zip(encoder.bottleneck_edge_source.tolist(), encoder.bottleneck_edge_target.tolist()))
                if pair == (2, 3)]] = 2.0
        encoder.set_bottleneck_edge_demand(demand)
        assert len(encoder.class_bottleneck_cache) == 0
        halved = encoder.response(*perturb(0, 1.0))
    assert abs(halved[0, 4, 0]) < abs(first[0, 4, 0])
    with pytest.raises(ValueError):
        encoder.set_bottleneck_edge_demand(torch.full_like(demand, 0.5))  # below 1 would raise availability past a step


def test_the_input_weights_are_positive_and_start_equal_so_the_pool_does_not_cancel():
    """The owner, 13 October 2026: "each class should be consistently signed within itself". A positive perturbation on
    an all-activating chain must reach the shared pool positive in every class, at every seed; the former signed draw
    is kept only to reproduce earlier runs."""
    pool = torch.tensor([False, False, False, True, False, False])
    for seed in range(6):
        torch.manual_seed(seed)
        encoder = LinearResponseEncoder(
            6, ["activates"], edge_source=torch.tensor([0, 1, 2, 3, 4]), edge_target=torch.tensor([1, 2, 3, 4, 5]),
            edge_relation=torch.zeros(5, dtype=torch.long), edge_sign=torch.ones(5), node_features=torch.eye(6), node_state_dim=4,
            cell_class_weights=torch.ones(6, 6), extracellular_pool_nodes=pool)
        weights = encoder.injection_weights()
        assert bool((weights > 0).all()) and torch.allclose(weights, weights[0].expand_as(weights))
        with torch.no_grad():
            assert bool((encoder.response(*perturb(0, 1.0))[0, 3] > 0).all())
    torch.manual_seed(0)
    signed = build_encoder()
    assert "input_weight" not in dict(signed.named_parameters())  # the default is positive
    torch.manual_seed(0)
    legacy = LinearResponseEncoder(6, RELATION_TYPES, edge_source=torch.tensor([edge[0] for edge in EDGES]),
                                   edge_target=torch.tensor([edge[1] for edge in EDGES]),
                                   edge_relation=torch.tensor([RELATION_TYPES.index(edge[2]) for edge in EDGES]),
                                   edge_sign=torch.tensor([edge[3] for edge in EDGES]), node_features=torch.eye(6), node_state_dim=4,
                                   input_weights="signed")
    assert "input_weight" in dict(legacy.named_parameters())
    with pytest.raises(ValueError):
        build_diamond_encoder(torch.ones(5, 1)).__class__(5, DIAMOND_RELATIONS, torch.tensor(DIAMOND_SOURCE), torch.tensor(DIAMOND_TARGET),
                                                          torch.zeros(5, dtype=torch.long), torch.ones(5), torch.eye(5), 4, input_weights="free")


# compartment pH as a throughput ceiling (the owner, 13 October 2026: "perfect solution. implement it")
def build_ph_encoder(num_nodes=6, class_weights=None, sensitivity=None):
    """A chain 0 -> 1 -> 2 -> 3 -> 4 with node 5 the proton of compartment "c"; nodes 1 and 2 are in "c", the rest in
    no compartment, so only they take a pH weight."""
    torch.manual_seed(0)
    encoder = LinearResponseEncoder(
        num_nodes, ["activates"], edge_source=torch.tensor([0, 1, 2, 3, 1]), edge_target=torch.tensor([1, 2, 3, 4, 5]),
        edge_relation=torch.zeros(5, dtype=torch.long), edge_sign=torch.ones(5),
        node_features=torch.eye(num_nodes), node_state_dim=4,
        cell_class_weights=torch.ones(num_nodes, 1) if class_weights is None else class_weights,
        ph_proton_node_index=torch.tensor([5]), ph_compartment_of_node=torch.tensor([-1, 0, 0, -1, -1, -1]),
        ph_compartment_names=["c"])
    if sensitivity is not None:
        with torch.no_grad():
            encoder.log_ph_sensitivity.fill_(float(torch.log(torch.tensor(sensitivity))))
    return encoder


def test_the_ph_weight_is_one_where_a_node_has_no_compartment_and_falls_with_the_proton_deviation():
    encoder = build_ph_encoder()
    node_index, sign_and_magnitude = torch.tensor([[0]]), torch.tensor([[[1.0, 1.0]]])
    with torch.no_grad():
        settled = encoder.response(node_index, sign_and_magnitude).permute(1, 0, 2)  # [N, B, C]
        weight = encoder.compartment_ph_weight(settled)
    assert torch.allclose(weight[[0, 3, 4, 5]], torch.ones_like(weight[[0, 3, 4, 5]]))  # no compartment, no weight
    assert torch.equal(weight[1], weight[2])  # one compartment, one weight
    assert bool((weight[1] > 0).all()) and bool((weight[1] <= 1).all())
    # a larger sensitivity pushes the weight down, because the argument grows
    small, large = build_ph_encoder(sensitivity=1.0), build_ph_encoder(sensitivity=1e6)
    with torch.no_grad():
        small_weight = small.compartment_ph_weight(small.response(node_index, sign_and_magnitude).permute(1, 0, 2))
        large_weight = large.compartment_ph_weight(large.response(node_index, sign_and_magnitude).permute(1, 0, 2))
    assert float(large_weight[1].max()) < float(small_weight[1].max())


def test_the_log_scale_costs_the_exact_negative_and_that_is_the_physiology():
    """The owner, 13 October 2026: "why not use true pH on the log scale so it fits nicely into equations?" A pH shift
    is not even in the deviation - acidifying by a fraction and alkalinising by it move pH by different amounts - so a
    gain-of-function stops being the exact negative of a loss-of-function where the compartment moves, which the
    previous linear form preserved. Outside the compartment the two are still exact negatives."""
    encoder = build_ph_encoder(sensitivity=1.0)
    with torch.no_grad():
        raised = encoder.response(torch.tensor([[0]]), torch.tensor([[[1.0, 1.0]]]))
        lowered = encoder.response(torch.tensor([[0]]), torch.tensor([[[-1.0, 1.0]]]))
    assert torch.allclose(raised[0, 0], -lowered[0, 0], atol=1e-6)  # node 0 is in no compartment
    assert not torch.allclose(raised[0, 1], -lowered[0, 1], atol=1e-6)  # node 1 is, and its weight differs
    # acidifying costs more throughput than alkalinising by the same fraction, since -log10(1 + x) is steeper below 0
    acidified = encoder.compartment_ph_weight(raised.permute(1, 0, 2))[1]
    alkalinised = encoder.compartment_ph_weight(lowered.permute(1, 0, 2))[1]
    assert float(acidified.min()) != float(alkalinised.min())


def test_the_resting_ph_values_cover_every_compartment_the_graph_names():
    """The constants the log scale needs, with their provenance separated: a value read from the review and a value
    reasoned where the review has none."""
    from mechanistic_pathway_learning.graph.compartment_resting_ph import RESTING_PH_BY_COMPARTMENT, resting_ph

    assert resting_ph("c") == 7.2 and resting_ph("l") == 4.7 and resting_ph("m") == 8.0
    assert resting_ph("c;m") is None  # a transport reaction is in no single compartment
    assert resting_ph("not a compartment") is None
    assert all(4.0 <= value <= 8.5 for value, _ in RESTING_PH_BY_COMPARTMENT.values())
    reasoned = {letter for letter, (_, why) in RESTING_PH_BY_COMPARTMENT.items() if why.startswith("reasoned")}
    assert reasoned == {"e", "i", "v"}  # the three the review does not cover, each with its reasoning recorded
    assert all(why.startswith(("review", "reasoned")) for _, why in RESTING_PH_BY_COMPARTMENT.values())


def test_the_ph_weight_breaks_positive_homogeneity_which_is_the_registered_cost():
    """A sensitivity where the weight is neither 1 nor saturated: on this chain the proton deviation is of order 1, so
    sensitivity 1 puts the weight near 0.74 and leaves it responsive in both directions."""
    encoder = build_ph_encoder(sensitivity=1.0)
    with torch.no_grad():
        single = encoder.response(torch.tensor([[0]]), torch.tensor([[[1.0, 1.0]]]))
        doubled = encoder.response(torch.tensor([[0]]), torch.tensor([[[1.0, 2.0]]]))
    assert not torch.allclose(doubled, 2.0 * single, atol=1e-6)  # the weight reads the deviation's size
    # compared at node 1, which is in the compartment; a node outside it scales exactly, since its weight stays 1
    assert float(doubled[0, 1].abs().max()) < float(2.0 * single[0, 1].abs().max())
    assert torch.allclose(doubled[0, 0], 2.0 * single[0, 0], atol=1e-6)


def test_calibration_sets_the_argument_to_its_target_and_reports_what_it_set():
    encoder = build_ph_encoder()
    node_index, sign_and_magnitude = torch.tensor([[0], [0]]), torch.tensor([[[1.0, 1.0]], [[-1.0, 2.0]]])
    with torch.no_grad():
        settled = encoder.response(node_index, sign_and_magnitude).permute(1, 0, 2)
    assert not encoder.ph_sensitivity_calibrated
    report = encoder.calibrate_ph_sensitivity(settled)
    assert encoder.ph_sensitivity_calibrated and set(report) == {"c"}
    deviation = settled[encoder.ph_proton_node_index]
    argument = torch.exp(encoder.log_ph_sensitivity)[:, None, None] * deviation
    assert float(argument.pow(2).mean().sqrt()) == pytest.approx(0.5, rel=1e-5)
    assert float(encoder.compartment_ph_weight(settled)[1].max()) < 1.0  # the weight bites at the calibrated scale


def test_an_uncalibrated_sensitivity_leaves_the_weight_at_one_on_a_realistic_deviation():
    """The case calibration exists for: on the confirmatory graph the settled field is of order 1e-11, so at
    sensitivity exp(0) = 1 the argument squares to 1e-22 and the weight is exactly 1 with no gradient."""
    encoder = build_ph_encoder()
    tiny = torch.zeros(6, 1, encoder.propagation_channels)
    tiny[5] = 2.4e-11
    assert float(encoder.compartment_ph_weight(tiny)[1].min()) == 1.0
    assert encoder.calibrate_ph_sensitivity(tiny) == {"c": pytest.approx(0.5 / 2.4e-11, rel=1e-5)}
    assert float(encoder.compartment_ph_weight(tiny)[1].max()) < 1.0


def test_a_compartment_the_batch_does_not_reach_keeps_its_sensitivity():
    encoder = build_ph_encoder()
    before = float(encoder.log_ph_sensitivity[0])
    report = encoder.calibrate_ph_sensitivity(torch.zeros(6, 1, encoder.propagation_channels))
    assert report == {} and float(encoder.log_ph_sensitivity[0]) == before


def test_the_ph_inputs_are_refused_unless_they_agree():
    with pytest.raises(ValueError):
        LinearResponseEncoder(6, ["activates"], torch.tensor([0]), torch.tensor([1]), torch.zeros(1, dtype=torch.long),
                              torch.ones(1), torch.eye(6), 4, ph_proton_node_index=torch.tensor([5]))
    with pytest.raises(ValueError):
        LinearResponseEncoder(6, ["activates"], torch.tensor([0]), torch.tensor([1]), torch.zeros(1, dtype=torch.long),
                              torch.ones(1), torch.eye(6), 4, ph_proton_node_index=torch.tensor([5]),
                              ph_compartment_of_node=torch.tensor([0, 0]))  # not one entry per node
    with pytest.raises(ValueError):
        LinearResponseEncoder(6, ["activates"], torch.tensor([0]), torch.tensor([1]), torch.zeros(1, dtype=torch.long),
                              torch.ones(1), torch.eye(6), 4, ph_proton_node_index=torch.tensor([5]),
                              ph_compartment_of_node=torch.tensor([0, 1, 0, 0, 0, 0]))  # compartment 1 has no proton node


def test_the_governing_proton_node_is_the_pool_the_reactions_use():
    """Two pools describe one species: Human-GEM's node carries 50 to 100 times the reactions of Reactome's in every
    compartment both describe, so it governs; a compartment only Reactome names is governed by its node."""
    from mechanistic_pathway_learning.graph.compartment_ph import compartment_index_per_node, proton_node_of_compartment

    node_ids = ["MAM02039c", "CHEBI15378c", "CHEBI15378v", "MAM02039e", "GENE:A", "REACTION:1"]
    display = ["H+", "H+", "H+", "H+", "A", "r1"]
    compartments = ["c", "c", "v", "e", "", "c;e"]
    governing = proton_node_of_compartment(node_ids, display, compartments, [3329, 44, 14, 719, 2, 2])
    assert governing == {"c": 0, "v": 2, "e": 3}  # Human-GEM wins in c; the vacuole is Reactome's alone
    order = sorted(governing)
    assert list(compartment_index_per_node(compartments, order)) == [0, 0, 2, 1, -1, -1]  # the transport reaction has none


def test_the_trainers_compartment_ph_inputs_build_the_encoder_on_a_graph_with_proton_nodes():
    """The registered model turns compartment pH on; until 16 October the trainer passed numpy index arrays and the
    constructor failed on every graph with a proton node in a compartment, which no test graph had."""
    from types import SimpleNamespace

    import numpy as np

    from experiments.run_main_model import compartment_ph_encoder_inputs

    data = SimpleNamespace(node_ids=["MAM02039c", "R1", "MAM00001c", "MAM00002m", "MAM02039m"], node_display_name=["H+", "R1", "A", "B", "H+"],
                           node_compartment=np.array(["c", "c", "c", "m", "m"]), node_degree=np.array([1.0, 3.0, 1.0, 1.0, 0.0]))
    inputs = compartment_ph_encoder_inputs(data)
    assert inputs["ph_compartment_names"] == ["c", "m"] and inputs["ph_proton_node_index"].tolist() == [0, 4]
    assert inputs["ph_compartment_of_node"].tolist() == [0, 0, 0, 1, 1] and inputs["ph_proton_node_index"].dtype == torch.long
    relations = ["substrate_of", "product_of"]
    encoder = LinearResponseEncoder(5, relations, torch.tensor([2, 1, 1]), torch.tensor([1, 0, 3]), torch.tensor([0, 1, 1]), torch.ones(3), torch.eye(5), 4,
                                    num_propagation_steps=6, fixed_point_tolerance=None, **inputs)
    assert encoder.ph_proton_node_index.tolist() == [0, 4] and torch.isfinite(encoder.response(*perturb(2, -1.0))).all()



def test_a_drugs_mechanism_edges_are_routes_of_the_limiting_step_bottleneck():
    """A drug's mechanism edges (entry relations) leave the normalised adjacency and are added on their own, so the
    bottleneck pass did not walk them: every node past the entry node had a bottleneck of 0, and the limiting-step rule
    multiplied the drug's whole response by it (found on the v5 graph, 17 October 2026). Node 0 is the drug, 1 its
    target, 2 downstream of the target."""
    torch.manual_seed(0)
    encoder = LinearResponseEncoder(
        3, ["activates", "targets"], edge_source=torch.tensor([0, 1]), edge_target=torch.tensor([1, 2]),
        edge_relation=torch.tensor([1, 0]), edge_sign=torch.tensor([1.0, 1.0]), node_features=torch.eye(3), node_state_dim=4,
        cell_class_weights=torch.tensor([[1.0], [0.6], [0.3]]), class_weight_rule="limiting_step",
        entry_node_mask=torch.tensor([True, False, False]), entry_relation_names=("targets",), edge_weight=torch.tensor([1.0, 1.0]))
    with torch.no_grad():
        bottleneck = encoder.class_bottleneck(torch.tensor([[0]]))[:, 0, 0]
        response = encoder.response(*perturb(0, 1.0))[0]
    assert bottleneck.tolist() == pytest.approx([1.0, 0.6, 0.3]), "the widest route runs through the mechanism edge"
    assert bool((response[1].abs() > 0).any()) and bool((response[2].abs() > 0).any()), "the drug's response reaches past its target"
    demand = encoder.set_bottleneck_edge_demand_from_edges(torch.tensor([0, 1]), torch.tensor([1, 2]), torch.tensor([1, 0]), torch.tensor([1.0, 1.0]))
    assert demand["bottleneck_entries"] == 2 and demand["largest_demand"] == 1.0, "the mechanism edge takes demand 1"
