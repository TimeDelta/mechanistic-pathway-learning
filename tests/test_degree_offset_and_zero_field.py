"""The degree-only control: a zero field leaves the sigmoid head a logistic regression on the covariate."""
import torch

from mechanistic_pathway_learning.models.baselines.relational_gnn_sigmoid_baseline import RelationalGnnSigmoidHead
from mechanistic_pathway_learning.models.baselines.zero_field_encoder import ZeroFieldEncoder


def test_zero_field_head_depends_only_on_the_covariate() -> None:
    torch.manual_seed(0)
    encoder = ZeroFieldEncoder(num_graph_nodes=5, node_state_dim=3)
    head = RelationalGnnSigmoidHead(node_state_dim=3, num_symptoms=2, degree_offset=True)
    field = encoder.perturbation_difference_field(torch.tensor([[0], [3]]), torch.tensor([[[-1.0, 1.0]], [[1.0, 1.0]]]))
    assert torch.count_nonzero(field) == 0
    with torch.no_grad():
        head.covariate_slope.fill_(1.0)
        probabilities = head(field, perturbation_covariate=torch.tensor([0.0, 1.0])).symptom_probability
    assert torch.all(probabilities[1] > probabilities[0])  # one shared slope: higher degree raises every symptom
    same_covariate = head(field, perturbation_covariate=torch.tensor([0.5, 0.5])).symptom_probability
    assert torch.allclose(same_covariate[0], same_covariate[1])  # different perturbations, same degree: same prediction


def test_local_descriptor_encoder_writes_only_at_the_perturbed_node_and_depends_on_its_features() -> None:
    from mechanistic_pathway_learning.models.baselines.local_descriptor_encoder import LocalDescriptorEncoder

    torch.manual_seed(0)
    features = torch.randn(6, 5)
    encoder = LocalDescriptorEncoder(6, features, node_state_dim=4)
    node_index = torch.tensor([[2, -1], [4, -1]])
    sign_and_magnitude = torch.tensor([[[-1.0, 1.0], [0.0, 0.0]], [[-1.0, 1.0], [0.0, 0.0]]])
    field = encoder.perturbation_difference_field(node_index, sign_and_magnitude)
    assert field.shape == (2, 6, 4)
    other_nodes = [node for node in range(6) if node != 2]
    assert torch.count_nonzero(field[0, other_nodes]) == 0 and torch.count_nonzero(field[1, [0, 1, 2, 3, 5]]) == 0
    # the same perturbation at two nodes with different descriptors gives different fields
    assert not torch.allclose(field[0, 2], field[1, 4])
    field.sum().backward()
    assert encoder.feature_projection.weight.grad is not None
