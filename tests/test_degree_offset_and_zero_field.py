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
