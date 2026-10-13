"""Tests for the descriptor treatments (models/descriptor_treatments.py): under seed_masked a perturbation's own
descriptors do not reach its field, under zero_init_slow no descriptor reaches it at initialisation, and plain keeps the
parameters it had before the treatments existed."""
from argparse import Namespace

import pytest
import torch

from experiments.run_main_model import optimizer_parameter_groups
from mechanistic_pathway_learning.graph.node_descriptors import descriptor_blocks
from mechanistic_pathway_learning.models.descriptor_treatments import seed_node_mask
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder

# gene 0 catalyses reaction 1, which turns substrate 2 into product 3; product 3 is the substrate of reaction 4,
# which makes product 5 (the graph of tests/test_linear_response_encoder.py)
RELATION_TYPES = ["substrate_of", "product_of", "catalyzed_by"]
EDGES = [(0, 1, "catalyzed_by", 1.0), (2, 1, "substrate_of", 1.0), (1, 3, "product_of", 1.0), (3, 4, "substrate_of", 1.0), (4, 5, "product_of", 1.0)]
NODE_TYPES = [0, 1, 2, 2, 1, 2]
NUM_NODES = 6
NUM_DESCRIPTOR_COLUMNS = 2
SEED_NODE, NEIGHBOUR_NODE = 0, 1


def node_features_with_descriptors() -> torch.Tensor:
    structural = torch.nn.functional.one_hot(torch.tensor(NODE_TYPES), num_classes=3).float()
    descriptors = torch.randn(NUM_NODES, NUM_DESCRIPTOR_COLUMNS, generator=torch.Generator().manual_seed(1))
    return torch.cat([structural, descriptors], dim=1)


def message_passing_encoder(treatment: str) -> tuple[RelationalMessagePassingEncoder, list]:
    torch.manual_seed(0)
    encoder = RelationalMessagePassingEncoder(NUM_NODES, len(RELATION_TYPES), node_state_dim=8, num_message_passing_layers=2,
                                              node_features=node_features_with_descriptors(),
                                              num_descriptor_columns=0 if treatment == "plain" else NUM_DESCRIPTOR_COLUMNS, descriptor_treatment=treatment)
    edge_index = torch.tensor([[edge[0] for edge in EDGES], [edge[1] for edge in EDGES]])
    edge_relation = torch.tensor([RELATION_TYPES.index(edge[2]) for edge in EDGES])
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, edge_relation, NUM_NODES, len(RELATION_TYPES))
    return encoder, adjacencies


def linear_response_encoder(treatment: str) -> LinearResponseEncoder:
    torch.manual_seed(0)
    return LinearResponseEncoder(
        NUM_NODES, RELATION_TYPES,
        edge_source=torch.tensor([edge[0] for edge in EDGES]), edge_target=torch.tensor([edge[1] for edge in EDGES]),
        edge_relation=torch.tensor([RELATION_TYPES.index(edge[2]) for edge in EDGES]), edge_sign=torch.tensor([edge[3] for edge in EDGES]),
        node_features=node_features_with_descriptors(), node_state_dim=4, num_propagation_steps=4,
        num_descriptor_columns=0 if treatment == "plain" else NUM_DESCRIPTOR_COLUMNS, descriptor_treatment=treatment,
    )


def perturbation_at(node: int) -> tuple[torch.Tensor, torch.Tensor]:
    return torch.tensor([[node]]), torch.tensor([[[-1.0, 1.0]]])


def change_descriptors(encoder, node: int) -> None:
    with torch.no_grad():
        encoder.node_features[node, -NUM_DESCRIPTOR_COLUMNS:] += 3.0


def message_passing_field(encoder, adjacencies) -> torch.Tensor:
    with torch.no_grad():
        return encoder.perturbation_difference_field(*perturbation_at(SEED_NODE), adjacencies)


def linear_response_field(encoder) -> torch.Tensor:
    with torch.no_grad():
        return encoder.perturbation_difference_field(*perturbation_at(SEED_NODE))


def test_seed_mask_marks_the_perturbed_nodes_and_ignores_padding() -> None:
    mask = seed_node_mask(torch.tensor([[0, -1], [-1, 3]]), 5)
    assert mask.tolist() == [[True, False, False, False, False], [False, False, False, True, False]]


def test_plain_keeps_the_parameters_it_had_before_the_treatments() -> None:
    encoder, _ = message_passing_encoder("plain")
    assert set(encoder.state_dict()) == {"node_features", "feature_projection.weight", "feature_projection.bias", "perturbation_injection.weight",
                                         "perturbation_injection.bias", "relation_weight", "self_weight", "layer_bias"}
    # the input weight is exp of log_input_weight under the default positive form (13 October 2026); plain adds nothing to either
    assert set(linear_response_encoder("plain").state_dict()) == {"gain_logit", "log_input_weight", "channel_expansion", "output_gate.weight", "output_gate.bias"}


@pytest.mark.parametrize("build, field", [(message_passing_encoder, message_passing_field), (linear_response_encoder, linear_response_field)])
def test_plain_lets_the_seed_descriptors_reach_the_field(build, field) -> None:
    built = build("plain")
    encoder, extra = (built if isinstance(built, tuple) else (built, None))
    before = field(encoder, extra) if extra is not None else field(encoder)
    change_descriptors(encoder, SEED_NODE)
    after = field(encoder, extra) if extra is not None else field(encoder)
    assert not torch.allclose(before, after)


@pytest.mark.parametrize("build, field", [(message_passing_encoder, message_passing_field), (linear_response_encoder, linear_response_field)])
def test_seed_masked_hides_the_seed_descriptors_and_keeps_the_others(build, field) -> None:
    built = build("seed_masked")
    encoder, extra = (built if isinstance(built, tuple) else (built, None))
    read = (lambda: field(encoder, extra)) if extra is not None else (lambda: field(encoder))
    before = read()
    change_descriptors(encoder, SEED_NODE)
    assert torch.allclose(read(), before)
    change_descriptors(encoder, NEIGHBOUR_NODE)
    assert not torch.allclose(read(), before)


@pytest.mark.parametrize("build, field", [(message_passing_encoder, message_passing_field), (linear_response_encoder, linear_response_field)])
def test_zero_init_slow_starts_without_descriptors(build, field) -> None:
    built = build("zero_init_slow")
    encoder, extra = (built if isinstance(built, tuple) else (built, None))
    read = (lambda: field(encoder, extra)) if extra is not None else (lambda: field(encoder))
    before = read()
    for node in range(NUM_NODES):
        change_descriptors(encoder, node)
    assert torch.allclose(read(), before)
    (descriptor_weight,) = encoder.descriptor_parameters()
    assert descriptor_weight.shape[1] == NUM_DESCRIPTOR_COLUMNS and torch.count_nonzero(descriptor_weight) == 0


def test_seed_masked_padding_gives_a_zero_field() -> None:
    encoder, adjacencies = message_passing_encoder("seed_masked")
    with torch.no_grad():
        padded = encoder.perturbation_difference_field(torch.tensor([[-1]]), torch.zeros(1, 1, 2), adjacencies)
    assert padded.shape == (1, NUM_NODES, 8) and torch.count_nonzero(padded) == 0


def test_descriptor_map_gets_its_own_learning_rate() -> None:
    encoder = linear_response_encoder("zero_init_slow")
    head = torch.nn.Linear(4, 2)
    arguments = Namespace(link_learning_rate=0.0, leak_learning_rate=0.0, module_bias_learning_rate=0.0, gate_learning_rate=0.0, descriptor_learning_rate=2e-4)
    groups = optimizer_parameter_groups(encoder, head, arguments)
    (descriptor_weight,) = encoder.descriptor_parameters()
    assert len(groups) == 2 and groups[1]["lr"] == 2e-4 and groups[1]["params"] == [descriptor_weight]
    assert all(parameter is not descriptor_weight for parameter in groups[0]["params"])
    plain_groups = optimizer_parameter_groups(linear_response_encoder("plain"), head, arguments)
    assert len(plain_groups) == 1



def test_descriptor_map_and_noisy_or_time_scales_get_one_group_each() -> None:
    """The confirmatory noisy-OR configurations set link, leak and gate rates; under zero_init_slow every parameter must
    sit in exactly one optimizer group, the descriptor map at its own rate."""
    from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import NoisyOrPathwayModuleHead
    encoder = linear_response_encoder("zero_init_slow")
    head = NoisyOrPathwayModuleHead(NUM_NODES, 4, 3, 2)
    arguments = Namespace(link_learning_rate=0.02, leak_learning_rate=2e-4, module_bias_learning_rate=0.0, gate_learning_rate=0.05, descriptor_learning_rate=2e-4)
    groups = optimizer_parameter_groups(encoder, head, arguments)
    grouped_ids = [id(parameter) for group in groups for parameter in group["params"]]
    every_parameter = list(encoder.parameters()) + list(head.parameters())
    assert sorted(grouped_ids) == sorted(id(parameter) for parameter in every_parameter)
    (descriptor_weight,) = encoder.descriptor_parameters()
    (descriptor_group,) = [group for group in groups if any(parameter is descriptor_weight for parameter in group["params"])]
    assert descriptor_group["lr"] == 2e-4 and len(descriptor_group["params"]) == 1
    assert {group.get("lr") for group in groups} >= {0.02, 2e-4, 0.05}


def test_seed_masked_message_passing_feeds_a_noisy_or_head_and_trains() -> None:
    """The slice arms ran the sigmoid head only; the confirmatory noisy-OR configurations would put the seed-masked
    field of a batch of perturbations into the noisy-OR head."""
    from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import NoisyOrPathwayModuleHead
    encoder, adjacencies = message_passing_encoder("seed_masked")
    head = NoisyOrPathwayModuleHead(NUM_NODES, 8, 3, 2, initial_readout_bias=-3.0)
    node_index = torch.tensor([[SEED_NODE, -1], [3, -1]])
    sign_and_magnitude = torch.tensor([[[-1.0, 1.0], [0.0, 0.0]], [[-1.0, 1.0], [0.0, 0.0]]])
    field = encoder.perturbation_difference_field(node_index, sign_and_magnitude, adjacencies)
    output = head(field)
    assert output.symptom_probability.shape == (2, 2) and torch.isfinite(output.symptom_probability).all()
    output.symptom_probability.sum().backward()
    assert encoder.feature_projection.weight.grad is not None and torch.isfinite(encoder.feature_projection.weight.grad).all()

def test_descriptor_blocks_cover_every_column_once() -> None:
    columns = ["metabolite_logp", "reaction_ec_class_1", "reaction_has_ec", "reaction_brain_region_pons", "protein_rrr_1", "gene_brain_has_expression"]
    assert descriptor_blocks(columns) == {"metabolite": ["metabolite_logp"], "reaction_ec": ["reaction_ec_class_1", "reaction_has_ec"],
                                          "reaction_brain": ["reaction_brain_region_pons"], "protein": ["protein_rrr_1"], "gene_brain": ["gene_brain_has_expression"]}
    with pytest.raises(ValueError):
        descriptor_blocks(["unknown_column"])


def test_a_perturbed_node_still_reads_its_regional_expression_but_not_its_identity():
    """The owner, 13 October 2026: "that brain gene expression masking not acceptable for perturbations". Whether the
    perturbed gene is expressed in a region decides whether its loss matters there, so those columns stay readable at
    the seed; the per-cell-class columns do not, being the same classes as the cell-class weight table."""
    from mechanistic_pathway_learning.graph.node_descriptors import columns_a_seed_still_reads
    from mechanistic_pathway_learning.models.descriptor_treatments import features_without_descriptors

    names = ["metabolite_logp", "protein_rrr_1", "gene_brain_region_amygdala", "gene_brain_gtex_brain_max",
             "gene_brain_has_expression", "gene_brain_class_astrocyte"]
    kept = columns_a_seed_still_reads(names)
    assert kept == [False, False, True, True, True, False]

    features = torch.arange(2 * (3 + len(names)), dtype=torch.float32).reshape(2, 3 + len(names))
    masked = features_without_descriptors(features, len(names), torch.tensor(kept))
    assert torch.equal(masked[:, :3], features[:, :3])  # structural columns untouched
    for position, name in enumerate(names):
        column = 3 + position
        if kept[position]:
            assert torch.equal(masked[:, column], features[:, column]), name
        else:
            assert torch.count_nonzero(masked[:, column]) == 0, name
    # with no mask given, every descriptor column goes, which is what every run before 13 October 2026 did
    assert torch.count_nonzero(features_without_descriptors(features, len(names))[:, 3:]) == 0
    with pytest.raises(ValueError):
        features_without_descriptors(features, len(names), torch.tensor(kept[:-1]))


def test_both_encoders_read_the_kept_columns_at_a_seed_and_the_rest_nowhere():
    """The gate (linear response) and the base state (message passing) must differ from the fully masked version by the
    kept columns alone, at the perturbed node and at no other."""
    from mechanistic_pathway_learning.models.descriptor_treatments import features_without_descriptors

    kept = torch.zeros(NUM_DESCRIPTOR_COLUMNS, dtype=torch.bool)
    kept[0] = True
    features = node_features_with_descriptors()
    fully_masked = features_without_descriptors(features, NUM_DESCRIPTOR_COLUMNS)
    partly_masked = features_without_descriptors(features, NUM_DESCRIPTOR_COLUMNS, kept)
    first = features.shape[1] - NUM_DESCRIPTOR_COLUMNS
    differing = (fully_masked != partly_masked).any(dim=0)
    assert bool(differing[first]) and not bool(differing[first + 1:].any()) and not bool(differing[:first].any())
