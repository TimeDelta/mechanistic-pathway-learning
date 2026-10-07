"""Tests for the linear-response encoder: stoichiometric signs, linearity, reach per step and convergence."""
import pytest
import torch

from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder

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


def _chain_encoder(normalisation, exponent=1.0, steps=64):
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
        propagation_channels=1, normalisation=normalisation, normalisation_exponent=exponent,
    )


def test_the_in_degree_power_exponent_of_one_reproduces_in_degree_normalisation():
    plain = _chain_encoder("in_degree").stacked_adjacency.to_dense()
    power = _chain_encoder("in_degree_power", exponent=1.0).stacked_adjacency.to_dense()
    assert torch.allclose(plain, power)


@pytest.mark.parametrize("normalisation,exponent", [("in_degree", 1.0), ("in_degree_power", 0.5),
                                                    ("in_degree_power", 0.25), ("spectral", 1.0)])
def test_every_normalisation_reaches_a_fixed_point(normalisation, exponent):
    """What the encoder needs is a fixed point, which needs the spectral radius below one. A row sum of sum_r |S_r|
    below one is sufficient for that but not necessary, and the spectral arm bounds the radius directly while leaving
    rows above one, so the property to assert for every arm is that the field stops moving as the steps grow."""
    with torch.no_grad():
        at_32 = _chain_encoder(normalisation, exponent=exponent, steps=32)(torch.tensor([[0]]), torch.tensor([[[1.0, 1.0]]]))
        at_128 = _chain_encoder(normalisation, exponent=exponent, steps=128)(torch.tensor([[0]]), torch.tensor([[[1.0, 1.0]]]))
    assert torch.isfinite(at_32).all() and torch.isfinite(at_128).all()
    assert torch.allclose(at_32, at_128, atol=1e-6), (normalisation, exponent)


@pytest.mark.parametrize("normalisation,exponent", [("in_degree", 1.0), ("in_degree_power", 0.5), ("in_degree_power", 0.25)])
def test_the_in_degree_arms_keep_every_row_sum_below_one(normalisation, exponent):
    """Softening the divisor raises the row sums above one, which is why that arm rescales once globally."""
    largest_row_sum = float(_chain_encoder(normalisation, exponent=exponent).stacked_adjacency.to_dense().abs().sum(dim=1).max())
    assert largest_row_sum <= 1.0 + 1e-6, (normalisation, exponent, largest_row_sum)


def test_the_exponent_redistributes_division_rather_than_reducing_it():
    """On a graph whose destinations all share one in-degree the global rescale exactly undoes the softening, so the
    exponent changes nothing; it only acts where in-degrees differ, moving division off the hub and onto the sparse
    nodes. Stating it as a test because the first version of this test assumed the exponent reduced division overall
    and failed on the uniform case."""
    uniform = {exponent: _chain_encoder("in_degree_power", exponent=exponent).stacked_adjacency.to_dense()
               for exponent in (1.0, 0.25)}
    assert torch.allclose(uniform[1.0], uniform[0.25]), "a uniform in-degree leaves nothing for the exponent to move"

    # 0 -> 1 -> 2, where node 1 is a hub with a hundred other sources and node 2 has only one other
    edge_list = [(0, 1), (1, 2)] + [(source, 1) for source in range(3, 103)] + [(103, 2)]
    source = torch.tensor([edge[0] for edge in edge_list])
    target = torch.tensor([edge[1] for edge in edge_list])

    def share_of_the_hub_edge(exponent):
        torch.manual_seed(0)
        encoder = LinearResponseEncoder(104, ["activates"], source, target, torch.zeros(len(edge_list), dtype=torch.long),
                                        torch.ones(len(edge_list)), node_features=torch.eye(104), node_state_dim=2,
                                        num_propagation_steps=32, propagation_channels=1,
                                        normalisation="in_degree_power", normalisation_exponent=exponent)
        dense = encoder.stacked_adjacency.to_dense()
        return float(dense[1, 0].abs() / dense[2, 1].abs())  # the hub's incoming weight against the sparse node's

    assert share_of_the_hub_edge(0.25) > share_of_the_hub_edge(1.0)


def test_an_unknown_normalisation_is_refused():
    with pytest.raises(ValueError, match="normalisation must be one of"):
        _chain_encoder("row_stochastic")
