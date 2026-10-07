"""Auxiliary supervision of the perturbation field with measured metabolite changes (design section 4, auxiliary labels;
mechanistic_pathway_learning/evidence/laboratory_abnormality_labels.py).

A training gene's laboratory-abnormality labels say which metabolites rise or fall in its disorder (PAH: phenylalanine
up). One linear readout without a bias, shared by every node, maps the field at a node to a predicted change, and the
predicted change of a metabolite is the sum over its compartment copies (the measured pool); a logistic loss on the
measured direction adds to the symptom loss with weight --laboratory-label-weight. With no bias, a node the field does not
reach predicts no change. Only labels of the training perturbations enter the loss; those of held-out genes are used
only to score the held-out field (sign agreement and AUROC). Labels without a direction are not used, and a metabolite
whose direction conflicts between a gene's labels (raised in urine, lowered in blood) is left out for that gene.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
import torch
from torch import Tensor, nn


class LaboratoryLabelIndex:
    def __init__(self, labels: pd.DataFrame, perturbation_ids: list[str], node_base_metabolite_ids) -> None:
        copies_of_base: dict[str, list[int]] = defaultdict(list)
        for node_index, base_id in enumerate(node_base_metabolite_ids):
            if isinstance(base_id, str) and base_id:
                copies_of_base[base_id].append(node_index)
        signed = labels[labels.direction != 0]
        directions_by_gene: dict[str, dict[str, set[int]]] = defaultdict(lambda: defaultdict(set))
        for gene, base_id, direction in zip(signed.gene_symbol, signed.base_metabolite_id, signed.direction):
            if base_id in copies_of_base:
                directions_by_gene[gene][base_id].add(int(direction))
        self.label_nodes: list[list[np.ndarray]] = []
        self.label_directions: list[list[int]] = []
        for perturbation_id in perturbation_ids:
            consistent = {base: next(iter(directions)) for base, directions in directions_by_gene.get(perturbation_id, {}).items() if len(directions) == 1}
            self.label_nodes.append([np.asarray(copies_of_base[base], dtype=np.int64) for base in sorted(consistent)])
            self.label_directions.append([consistent[base] for base in sorted(consistent)])

    def count(self, indices) -> int:
        return sum(len(self.label_directions[index]) for index in indices)

    def batch(self, batch_indices, device=None) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        """(batch position, node index, label slot) per labelled node copy, and the direction per label slot."""
        positions, nodes, slots, directions = [], [], [], []
        for position, index in enumerate(batch_indices):
            for node_array, direction in zip(self.label_nodes[index], self.label_directions[index]):
                slot = len(directions)
                directions.append(direction)
                positions += [position] * len(node_array)
                nodes += node_array.tolist()
                slots += [slot] * len(node_array)
        as_long = lambda values: torch.as_tensor(values, dtype=torch.long, device=device)  # noqa: E731
        return as_long(positions), as_long(nodes), as_long(slots), torch.as_tensor(directions, dtype=torch.float32, device=device)


class LaboratoryReadout(nn.Module):
    def __init__(self, node_state_dim: int) -> None:
        super().__init__()
        self.linear = nn.Linear(node_state_dim, 1, bias=False)

    def label_scores(self, field: Tensor, positions: Tensor, nodes: Tensor, slots: Tensor, num_labels: int) -> Tensor:
        """Predicted change per label: the readout of the field summed over the metabolite's compartment copies."""
        values = self.linear(field[positions, nodes]).squeeze(-1)
        return torch.zeros(num_labels, device=field.device, dtype=values.dtype).index_add(0, slots, values)


def laboratory_sign_loss(scores: Tensor, directions: Tensor) -> Tensor:
    """Logistic loss of the measured direction given the predicted change."""
    return nn.functional.softplus(-directions * scores).mean()
