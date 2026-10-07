"""Tests for the linear-response encoder: stoichiometric signs, linearity, reach per step and convergence."""
import pytest
import torch

from mechanistic_pathway_learning.models.linear_response_encoder import (CrossRelationAggregator, LinearResponseEncoder,
                                                                           project_onto_simplex,
                                                                           MIXTURE_STATISTICS, MIXTURE_WEIGHTINGS,
                                                                           NORMALISATIONS)

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


def _metabolite_row_sum(normalisation: str, num_producing: int, num_consuming: int) -> float:
    """Signed row sum, over every relation, at metabolite 0, made by num_producing reactions and used by num_consuming."""
    producing = list(range(1, num_producing + 1))
    consuming = list(range(num_producing + 1, num_producing + num_consuming + 1))
    edges = [(reaction, 0, "product_of") for reaction in producing] + [(0, reaction, "substrate_of") for reaction in consuming]
    num_nodes = 1 + num_producing + num_consuming
    stacked, relation_names = LinearResponseEncoder.signed_stacked_adjacency(
        num_nodes, RELATION_TYPES, torch.tensor([edge[0] for edge in edges]), torch.tensor([edge[1] for edge in edges]),
        torch.tensor([RELATION_TYPES.index(edge[2]) for edge in edges]), torch.ones(len(edges)), normalisation=normalisation)
    return float(stacked.to_dense().view(len(relation_names), num_nodes, num_nodes).sum(dim=0)[0].sum())


def test_in_degree_erases_the_stoichiometric_counts_and_total_in_degree_keeps_them():
    """Five producing reactions against one consuming reaction: averaging inside each relation leaves
    mean(+1) / 2 + mean(-1) / 2 = 0 whatever the counts, while one divisor over all relations leaves (5 - 1) / 6."""
    assert _metabolite_row_sum("in_degree", 5, 1) == pytest.approx(0.0, abs=1e-6)
    assert _metabolite_row_sum("in_degree", 1, 5) == pytest.approx(0.0, abs=1e-6)
    assert _metabolite_row_sum("total_in_degree", 5, 1) == pytest.approx(4 / 6)
    assert _metabolite_row_sum("total_in_degree", 1, 5) == pytest.approx(-4 / 6)
    assert _metabolite_row_sum("total_in_degree", 2, 2) == pytest.approx(0.0, abs=1e-6)  # equal counts still cancel, which is correct


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
