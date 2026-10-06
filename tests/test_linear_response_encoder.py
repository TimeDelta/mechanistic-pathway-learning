"""Tests for the linear-response encoder: stoichiometric signs, linearity, reach per step and convergence."""
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
