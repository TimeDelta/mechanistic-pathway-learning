"""Tests for the auxiliary laboratory-label readout: label indexing, pooling over compartment copies and the loss."""
import pandas as pd
import torch

from mechanistic_pathway_learning.models.laboratory_readout import LaboratoryLabelIndex, LaboratoryReadout, laboratory_sign_loss

NODE_BASE_METABOLITES = [None, "MAM_PHE", "MAM_PHE", "MAM_TYR", "MAM_X"]  # node 0 is a gene; phenylalanine in two compartments
LABELS = pd.DataFrame({
    "gene_symbol": ["PAH", "PAH", "PAH", "PAH", "OTHER"],
    "base_metabolite_id": ["MAM_PHE", "MAM_TYR", "MAM_X", "MAM_X", "MAM_TYR"],
    "direction": [1, -1, 1, -1, 0],  # MAM_X conflicts for PAH; OTHER's label has no direction
})


def test_index_pools_compartment_copies_and_drops_conflicts_and_undirected_labels() -> None:
    index = LaboratoryLabelIndex(LABELS, ["PAH", "OTHER", "NONE"], NODE_BASE_METABOLITES)
    assert index.count([0]) == 2 and index.count([1, 2]) == 0
    positions, nodes, slots, directions = index.batch([2, 0])
    assert positions.tolist() == [1, 1, 1] and nodes.tolist() == [1, 2, 3] and slots.tolist() == [0, 0, 1]
    assert directions.tolist() == [1.0, -1.0]


def test_readout_sums_copies_predicts_nothing_where_the_field_is_zero_and_the_loss_prefers_the_measured_sign() -> None:
    torch.manual_seed(0)
    readout = LaboratoryReadout(node_state_dim=3)
    index = LaboratoryLabelIndex(LABELS, ["PAH"], NODE_BASE_METABOLITES)
    positions, nodes, slots, directions = index.batch([0])
    field = torch.zeros(1, 5, 3)
    assert torch.equal(readout.label_scores(field, positions, nodes, slots, 2), torch.zeros(2))
    field[0, 1] = field[0, 2] = torch.tensor([1.0, 0.0, 0.0])
    scores = readout.label_scores(field, positions, nodes, slots, 2)
    assert torch.isclose(scores[0], 2 * readout.linear.weight[0, 0])
    assert laboratory_sign_loss(torch.tensor([2.0, -2.0]), directions) < laboratory_sign_loss(torch.tensor([-2.0, 2.0]), directions)
