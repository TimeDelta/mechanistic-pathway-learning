"""The node-major, stacked-adjacency propagation must equal the batch-major per-relation loop it replaced."""
from __future__ import annotations

import torch

from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder


def reference_propagate(encoder: RelationalMessagePassingEncoder, node_state_field: torch.Tensor, relation_adjacencies: list) -> torch.Tensor:
    """The previous implementation, kept here as the numerical reference."""
    batch_size = node_state_field.shape[0]
    for layer_index in range(encoder.num_message_passing_layers):
        aggregated_messages = torch.zeros_like(node_state_field)
        for relation_index, adjacency in enumerate(relation_adjacencies):
            if adjacency is None:
                continue
            messages = node_state_field @ encoder.relation_weight[layer_index, relation_index]
            flattened = messages.permute(1, 0, 2).reshape(encoder.num_graph_nodes, batch_size * encoder.node_state_dim)
            aggregated = torch.sparse.mm(adjacency, flattened)
            aggregated_messages = aggregated_messages + aggregated.reshape(encoder.num_graph_nodes, batch_size, encoder.node_state_dim).permute(1, 0, 2)
        self_messages = node_state_field @ encoder.self_weight[layer_index]
        node_state_field = torch.relu(self_messages + aggregated_messages + encoder.layer_bias[layer_index])
    return node_state_field


def random_graph(num_nodes: int, num_edges: int, num_relation_types: int, generator: torch.Generator) -> tuple[torch.Tensor, torch.Tensor]:
    edge_index = torch.randint(0, num_nodes, (2, num_edges), generator=generator)
    edge_relation_type = torch.randint(0, num_relation_types, (num_edges,), generator=generator)
    edge_relation_type[edge_relation_type == 2] = 0  # relation 2 absent: exercises the None entry
    return edge_index, edge_relation_type


def test_stacked_propagation_matches_reference_and_gradients() -> None:
    generator = torch.Generator().manual_seed(7)
    num_nodes, num_relation_types, state_dim, batch_size = 60, 4, 8, 5
    edge_index, edge_relation_type = random_graph(num_nodes, 300, num_relation_types, generator)
    torch.manual_seed(7)
    encoder = RelationalMessagePassingEncoder(num_nodes, num_relation_types, state_dim, num_message_passing_layers=3)
    adjacencies = encoder.build_relation_adjacencies(edge_index, edge_relation_type, num_nodes, num_relation_types)
    assert adjacencies[2] is None
    perturbed_nodes = torch.randint(0, num_nodes, (batch_size, 3), generator=generator)
    perturbed_nodes[0, 1:] = -1
    sign_and_magnitude = torch.randn(batch_size, 3, 2, generator=generator)
    initial_field = encoder.initial_node_state_field(perturbed_nodes, sign_and_magnitude)

    reference = reference_propagate(encoder, initial_field, adjacencies)
    stacked = encoder.propagate(initial_field, adjacencies)
    assert stacked.shape == (batch_size, num_nodes, state_dim)
    assert torch.allclose(stacked, reference, atol=1e-5, rtol=1e-5)

    reference_loss = reference_propagate(encoder, encoder.initial_node_state_field(perturbed_nodes, sign_and_magnitude), adjacencies).pow(2).sum()
    reference_gradients = torch.autograd.grad(reference_loss, [encoder.relation_weight, encoder.self_weight, encoder.layer_bias])
    stacked_loss = encoder.forward(perturbed_nodes, sign_and_magnitude, relation_adjacencies=adjacencies).pow(2).sum()
    stacked_gradients = torch.autograd.grad(stacked_loss, [encoder.relation_weight, encoder.self_weight, encoder.layer_bias])
    for reference_gradient, stacked_gradient in zip(reference_gradients, stacked_gradients):
        assert torch.allclose(stacked_gradient, reference_gradient, atol=1e-4, rtol=1e-4)


def test_difference_field_is_zero_without_perturbation_and_cache_follows_the_adjacency_list() -> None:
    generator = torch.Generator().manual_seed(3)
    encoder = RelationalMessagePassingEncoder(30, 3, 4, num_message_passing_layers=2)
    edge_index, edge_relation_type = random_graph(30, 90, 3, generator)
    adjacencies = encoder.build_relation_adjacencies(edge_index, edge_relation_type, 30, 3)
    no_perturbation = encoder.perturbation_difference_field(torch.full((1, 1), -1), torch.zeros(1, 1, 2), adjacencies)
    assert torch.allclose(no_perturbation, torch.zeros_like(no_perturbation), atol=1e-6)
    first_stack, _ = encoder.stacked_relation_adjacency(adjacencies)
    other_edges, other_relations = random_graph(30, 40, 3, generator)
    other_adjacencies = encoder.build_relation_adjacencies(other_edges, other_relations, 30, 3)
    second_stack, _ = encoder.stacked_relation_adjacency(other_adjacencies)
    assert first_stack._nnz() != second_stack._nnz() or not torch.equal(first_stack.indices(), second_stack.indices())
