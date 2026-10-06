"""Tests for the noisy-OR pathway module head."""
import torch

from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import (
    HardConcreteNodeGate,
    NoisyOrPathwayModuleHead,
    noisy_or_combination,
)


def test_noisy_or_single_sufficient_module_explains_symptom() -> None:
    module_activation = torch.tensor([[1.0, 0.0, 0.0]])
    link_probability = torch.tensor([[0.999, 0.0], [0.0, 0.0], [0.0, 0.0]])
    leak_probability = torch.tensor([0.0, 0.0])
    symptom_probability = noisy_or_combination(module_activation, link_probability, leak_probability)
    assert symptom_probability[0, 0].item() > 0.99
    assert symptom_probability[0, 1].item() < 1e-5


def test_noisy_or_no_active_module_returns_leak() -> None:
    module_activation = torch.zeros(1, 3)
    link_probability = torch.full((3, 2), 0.9)
    leak_probability = torch.tensor([0.05, 0.2])
    symptom_probability = noisy_or_combination(module_activation, link_probability, leak_probability)
    assert torch.allclose(symptom_probability, leak_probability[None, :], atol=1e-5)


def test_noisy_or_modules_combine_without_interaction_term() -> None:
    # two half-sufficient modules: 1 - (1 - 0.5)(1 - 0.5) = 0.75
    module_activation = torch.tensor([[1.0, 1.0]])
    link_probability = torch.tensor([[0.5], [0.5]])
    leak_probability = torch.tensor([0.0])
    symptom_probability = noisy_or_combination(module_activation, link_probability, leak_probability)
    assert abs(symptom_probability.item() - 0.75) < 1e-5


def test_hard_concrete_gate_shapes_and_range() -> None:
    gate = HardConcreteNodeGate(num_pathway_modules=4, num_graph_nodes=10)
    gate.train()
    sampled = gate()
    assert sampled.shape == (4, 10)
    assert torch.all(sampled >= 0.0) and torch.all(sampled <= 1.0)
    gate.eval()
    evaluation_gate = gate()
    assert torch.all(evaluation_gate >= 0.0) and torch.all(evaluation_gate <= 1.0)
    assert gate.expected_active_node_count().shape == (4,)


def test_head_forward_shapes_and_gradients() -> None:
    batch_size, num_graph_nodes, node_state_dim, num_modules, num_symptoms = 3, 12, 8, 4, 5
    head = NoisyOrPathwayModuleHead(num_graph_nodes, node_state_dim, num_modules, num_symptoms)
    node_state_field = torch.randn(batch_size, num_graph_nodes, node_state_dim, requires_grad=True)
    output = head(node_state_field, relation_index=0)
    assert output.symptom_probability.shape == (batch_size, num_symptoms)
    assert output.module_activation.shape == (batch_size, num_modules)
    assert output.link_probability.shape == (num_modules, num_symptoms)
    assert output.module_support.shape == (num_modules, num_graph_nodes)
    assert torch.all(output.symptom_probability > 0.0) and torch.all(output.symptom_probability < 1.0)
    loss = output.symptom_probability.sum() + head.description_length_penalty()
    loss.backward()
    assert node_state_field.grad is not None
    assert head.support_gate.log_alpha.grad is not None
    assert output.dominant_module_per_perturbation(symptom_index=0).shape == (batch_size,)


def test_relation_index_selects_separate_link_matrix() -> None:
    head = NoisyOrPathwayModuleHead(num_graph_nodes=6, node_state_dim=4, num_pathway_modules=2, num_symptoms=3)
    with torch.no_grad():
        head.module_symptom_link_logit[0].fill_(5.0)
        head.module_symptom_link_logit[1].fill_(-5.0)
    assert head.link_probability(0).mean() > 0.9
    assert head.link_probability(1).mean() < 0.1


def test_leak_initialization_reproduces_base_rates_and_off_by_default_bias_silences_modules() -> None:
    torch.manual_seed(0)  # the symmetry-breaking check below compares a sample standard deviation of 18 gates with 0.3
    head = NoisyOrPathwayModuleHead(num_graph_nodes=6, node_state_dim=4, num_pathway_modules=3, num_symptoms=2, initial_readout_bias=-3.0, gate_initial_log_alpha_noise=0.5)
    base_rates = torch.tensor([0.25, 0.4])
    head.initialize_leak_from_base_rates(base_rates)
    assert torch.allclose(head.leak_probability(0), base_rates, atol=1e-6)
    head.eval()
    with torch.no_grad():
        output = head(torch.zeros(2, 6, 4), relation_index=0)  # no perturbation signal in the field
    assert torch.all(output.module_activation < 0.06)  # sigmoid(-3) with nothing to read
    assert torch.allclose(output.symptom_probability, 1.0 - (1.0 - base_rates) * torch.prod(1.0 - output.module_activation[:, :, None] * output.link_probability[None], dim=1), atol=1e-6)
    assert head.support_gate.log_alpha.std() > 0.3  # symmetry between modules is broken at initialization


def test_expected_gate_matches_the_mean_of_training_samples_and_exceeds_the_louizos_gate_when_fractional() -> None:
    torch.manual_seed(0)
    gate = HardConcreteNodeGate(num_pathway_modules=1, num_graph_nodes=3, initial_log_alpha=0.0, initial_log_alpha_noise=0.0)
    with torch.no_grad():
        gate.log_alpha.copy_(torch.tensor([[-2.3, 0.0, 2.0]]))
    gate.train()
    monte_carlo_mean = torch.stack([gate() for _ in range(20_000)]).mean(dim=0)
    expected = gate.expected_gate()
    assert torch.allclose(expected, monte_carlo_mean, atol=0.01)
    # at log_alpha -2.3 the training gate averages about 0.11 while the Louizos test-time estimator gives about 0.009
    assert expected[0, 0].item() > 5 * gate.louizos_test_time_gate()[0, 0].item()


def test_eval_mode_gate_is_the_expected_training_gate_and_time_scale_groups_are_links_leaks_and_biases() -> None:
    gate = HardConcreteNodeGate(num_pathway_modules=2, num_graph_nodes=5)
    gate.eval()
    assert torch.equal(gate(), gate.expected_gate())
    head = NoisyOrPathwayModuleHead(num_graph_nodes=5, node_state_dim=3, num_pathway_modules=2, num_symptoms=2)
    groups = head.time_scale_parameter_groups()
    assert [id(parameter) for parameter in groups["links"]] == [id(head.module_symptom_link_logit)]
    assert [id(parameter) for parameter in groups["leaks"]] == [id(head.symptom_leak_logit)]
    assert [id(parameter) for parameter in groups["module_biases"]] == [id(head.module_readout_bias)]


def test_degree_offset_makes_the_leak_depend_on_the_covariate_and_starts_neutral() -> None:
    torch.manual_seed(0)
    head = NoisyOrPathwayModuleHead(num_graph_nodes=6, node_state_dim=4, num_pathway_modules=2, num_symptoms=3, degree_offset=True)
    head.initialize_leak_from_base_rates(torch.tensor([0.2, 0.3, 0.4]))
    head.eval()  # deterministic gates, so two passes differ only through the covariate
    field = torch.randn(2, 6, 4)
    covariate = torch.tensor([-1.0, 2.0])
    with torch.no_grad():
        without_covariate = head(field).symptom_probability
        with_covariate = head(field, perturbation_covariate=covariate).symptom_probability
        assert torch.allclose(without_covariate, with_covariate)  # zero slope at initialisation
        head.leak_covariate_slope.fill_(1.0)
        leaks = head(field, perturbation_covariate=covariate).leak_probability
    assert leaks.shape == (2, 3)
    assert torch.all(leaks[1] > leaks[0])  # higher degree, higher leak once the slope is positive
    assert head.leak_covariate_slope is not None
    assert NoisyOrPathwayModuleHead(num_graph_nodes=6, node_state_dim=4, num_pathway_modules=2, num_symptoms=3).leak_covariate_slope is None
