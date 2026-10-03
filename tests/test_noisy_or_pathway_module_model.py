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
    deterministic = gate()
    assert torch.all(deterministic >= 0.0) and torch.all(deterministic <= 1.0)
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
