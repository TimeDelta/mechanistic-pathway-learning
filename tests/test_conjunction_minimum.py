"""The obligate-complex relation aggregates by a per-channel minimum, not a mean.

The user, 9 October 2026: "i thought it would be multiple mins (one for each preserved quantity)" and "even if it had
costed compute, it would have to be a lot of compute to outweigh the true physiology". These tests pin the three
properties that make the minimum worth a second code path: it is the smallest member rather than their average, it is
taken channel by channel so one member can set one quantity and another member the next, and a node with no member is
untouched.
"""
import torch

from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder

MEMBER_RELATION = 0
NUM_NODES = 4
STATE_DIM = 3
# nodes 1 and 2 are the members of the complex at node 0; node 3 is in no complex
MEMBER_EDGE_INDEX = torch.tensor([[1, 2], [0, 0]])
MEMBER_RELATION_TYPE = torch.tensor([MEMBER_RELATION, MEMBER_RELATION])


def encoder_with(aggregation: str, temperature: float = 0.1) -> RelationalMessagePassingEncoder:
    """One layer, the member relation reduced as asked, and an identity relation weight so the message is readable."""
    encoder = RelationalMessagePassingEncoder(NUM_NODES, num_relation_types=1, node_state_dim=STATE_DIM,
                                             num_message_passing_layers=1,
                                             conjunction_relation_indices=(MEMBER_RELATION,),
                                             conjunction_aggregation=aggregation, soft_minimum_temperature=temperature)
    with torch.no_grad():
        encoder.relation_weight.copy_(torch.eye(STATE_DIM)[None, None].expand_as(encoder.relation_weight))
        encoder.self_weight.zero_()
        encoder.layer_bias.zero_()
    return encoder


def member_states() -> torch.Tensor:
    """[1, NUM_NODES, STATE_DIM]: member 1 is scarcest in channels 0 and 2, member 2 in channel 1."""
    field = torch.zeros(1, NUM_NODES, STATE_DIM)
    field[0, 1] = torch.tensor([0.2, 0.9, 0.3])
    field[0, 2] = torch.tensor([0.7, 0.4, 0.8])
    field[0, 3] = torch.tensor([0.5, 0.5, 0.5])
    return field


def adjacencies() -> list[torch.Tensor | None]:
    return RelationalMessagePassingEncoder.build_relation_adjacencies(MEMBER_EDGE_INDEX, MEMBER_RELATION_TYPE, NUM_NODES, 1)


def test_the_complex_takes_the_smallest_member_in_every_channel():
    propagated = encoder_with("minimum").propagate(member_states(), adjacencies())
    # relu of the identity-mapped minimum: channel 0 and 2 from member 1, channel 1 from member 2
    assert torch.allclose(propagated[0, 0], torch.tensor([0.2, 0.4, 0.3]), atol=1e-6)


def test_the_minimum_is_not_the_mean_of_the_members():
    minimum = encoder_with("minimum").propagate(member_states(), adjacencies())
    mean = encoder_with("mean").propagate(member_states(), adjacencies())
    assert torch.allclose(mean[0, 0], torch.tensor([0.45, 0.65, 0.55]), atol=1e-6)
    assert not torch.allclose(minimum[0, 0], mean[0, 0], atol=1e-3)


def test_one_scarce_channel_cannot_be_covered_by_a_plentiful_member():
    """Raising member 2 without limit leaves the channels member 1 is scarcest in exactly where they were."""
    field = member_states()
    propagated = encoder_with("minimum").propagate(field, adjacencies())
    field[0, 2] = torch.tensor([50.0, 50.0, 50.0])
    raised = encoder_with("minimum").propagate(field, adjacencies())
    assert torch.allclose(raised[0, 0, 0], propagated[0, 0, 0], atol=1e-6)
    assert torch.allclose(raised[0, 0, 2], propagated[0, 0, 2], atol=1e-6)
    assert raised[0, 0, 1] > propagated[0, 0, 1]  # channel 1's scarcest member was the one raised


def test_a_node_in_no_complex_gets_no_conjunction_message():
    propagated = encoder_with("minimum").propagate(member_states(), adjacencies())
    assert torch.allclose(propagated[0, 3], torch.zeros(STATE_DIM), atol=1e-6)
    assert torch.allclose(propagated[0, 1], torch.zeros(STATE_DIM), atol=1e-6)


def test_the_soft_minimum_approaches_the_hard_one_as_the_temperature_falls():
    hard = encoder_with("minimum").propagate(member_states(), adjacencies())
    warm = encoder_with("soft_minimum", temperature=0.5).propagate(member_states(), adjacencies())
    cold = encoder_with("soft_minimum", temperature=0.005).propagate(member_states(), adjacencies())
    assert torch.allclose(cold[0, 0], hard[0, 0], atol=1e-4)
    assert (warm[0, 0] - hard[0, 0]).abs().max() > (cold[0, 0] - hard[0, 0]).abs().max()


def test_the_soft_minimum_never_exceeds_the_hard_one():
    """-T log sum exp(-x/T) <= min(x), so the surrogate is a lower bound and cannot invent presence."""
    hard = encoder_with("minimum").propagate(member_states(), adjacencies())
    warm = encoder_with("soft_minimum", temperature=0.5).propagate(member_states(), adjacencies())
    assert bool((warm[0, 0] <= hard[0, 0] + 1e-6).all())


def test_the_conjunction_relation_is_not_also_mean_aggregated():
    """Its rows must leave the stacked product, or the member edges would be counted twice."""
    encoder = encoder_with("minimum")
    _, present_index = encoder.stacked_relation_adjacency(adjacencies())
    assert MEMBER_RELATION not in present_index.tolist()
    mean_encoder = encoder_with("mean")
    _, mean_present = mean_encoder.stacked_relation_adjacency(adjacencies())
    assert MEMBER_RELATION in mean_present.tolist()


def test_the_minimum_carries_a_gradient_to_the_member_that_set_it():
    encoder = encoder_with("minimum")
    field = member_states().requires_grad_(True)
    encoder.propagate(field, adjacencies())[0, 0, 0].backward()
    assert field.grad[0, 1, 0] != 0.0   # member 1 is scarcest in channel 0
    assert field.grad[0, 2, 0] == 0.0   # member 2 is not, so it has no say in that channel


def test_a_relation_absent_from_the_graph_is_skipped_rather_than_failing():
    encoder = RelationalMessagePassingEncoder(NUM_NODES, num_relation_types=2, node_state_dim=STATE_DIM,
                                              num_message_passing_layers=1, conjunction_relation_indices=(1,),
                                              conjunction_aggregation="minimum")
    propagated = encoder.propagate(member_states(), RelationalMessagePassingEncoder.build_relation_adjacencies(
        MEMBER_EDGE_INDEX, MEMBER_RELATION_TYPE, NUM_NODES, 2))
    assert propagated.shape == (1, NUM_NODES, STATE_DIM)


def test_mean_aggregation_keeps_the_encoder_it_was_before_the_option():
    """Asking for the mean must empty the conjunction list, whatever relations were named."""
    encoder = encoder_with("mean")
    assert encoder.conjunction_relation_indices == ()
