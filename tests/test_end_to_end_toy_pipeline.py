"""Smoke test: encoder + noisy-OR head train on a toy typed graph with two disjoint mechanisms."""
import torch

from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import NoisyOrPathwayModuleHead
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder
from mechanistic_pathway_learning.models.soft_constraint_losses import evidence_weighted_binary_cross_entropy


def build_toy_graph():
    # two linear chains of reactions (nodes 0-4 and 5-9) that never touch, one symptom reached by either
    edges = []
    relation_types = []
    for chain_start in (0, 5):
        for node in range(chain_start, chain_start + 4):
            edges.append((node, node + 1))
            relation_types.append(0)  # product-of
            edges.append((node + 1, node))
            relation_types.append(1)  # substrate-of (reverse view)
    edge_index = torch.tensor(edges, dtype=torch.long).t()
    edge_relation_type = torch.tensor(relation_types, dtype=torch.long)
    return edge_index, edge_relation_type


def test_toy_pipeline_learns_two_independent_mechanisms_for_one_symptom() -> None:
    torch.manual_seed(0)
    num_graph_nodes, node_state_dim, num_modules, num_symptoms = 10, 16, 2, 2
    edge_index, edge_relation_type = build_toy_graph()
    encoder = RelationalMessagePassingEncoder(num_graph_nodes, num_relation_types=2, node_state_dim=node_state_dim, num_message_passing_layers=2)
    head = NoisyOrPathwayModuleHead(num_graph_nodes, node_state_dim, num_modules, num_symptoms, gate_initial_log_alpha=0.0)
    # perturbations: knocking out any node in chain A (0-4) or chain B (5-9) induces symptom 0; symptom 1 is never induced
    perturbed_nodes = torch.arange(num_graph_nodes)[:, None]
    sign_and_magnitude = torch.tensor([[[-1.0, 1.0]]] * num_graph_nodes)
    observed_outcome = torch.zeros(num_graph_nodes, num_symptoms)
    observed_outcome[:, 0] = 1.0
    evidence_weight = torch.ones(num_graph_nodes, num_symptoms)
    optimizer = torch.optim.Adam(list(encoder.parameters()) + list(head.parameters()), lr=0.05)
    for _ in range(150):
        optimizer.zero_grad()
        node_state_field = encoder(perturbed_nodes, sign_and_magnitude, edge_index, edge_relation_type)
        output = head(node_state_field, relation_index=0)
        loss = evidence_weighted_binary_cross_entropy(output.symptom_probability, observed_outcome, evidence_weight)
        loss = loss + 1e-3 * head.description_length_penalty()
        loss.backward()
        optimizer.step()
    encoder.eval()
    head.eval()
    with torch.no_grad():
        output = head(encoder(perturbed_nodes, sign_and_magnitude, edge_index, edge_relation_type), relation_index=0)
    assert output.symptom_probability[:, 0].mean() > 0.9
    assert output.symptom_probability[:, 1].mean() < 0.1
