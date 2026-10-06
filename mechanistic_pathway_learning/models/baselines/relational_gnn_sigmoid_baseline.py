"""Baseline B3: relational GNN with a sigmoid head and no modules (design section 5.6).

The same relational message passing encoder as the proposed model produces the node-state
field; this head pools it over all graph nodes (sum of the field, which for the difference
field is the total propagated change), passes the pooled vector through one hidden layer
and emits an independent sigmoid per symptom. There are no pathway modules, no sparse
supports and no noisy-OR, so the comparison with B6 isolates the module structure. The head
exposes the attributes the training harness reads from the noisy-OR head (symptom_probability
on the output, a zero description-length penalty) so one loop trains both.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import Tensor, nn


@dataclass
class SigmoidHeadOutput:
    symptom_probability: Tensor  # [batch_size, num_symptoms]


class RelationalGnnSigmoidHead(nn.Module):
    def __init__(self, node_state_dim: int, num_symptoms: int, hidden_dim: int = 64, pooling: str = "sum", degree_offset: bool = False) -> None:
        super().__init__()
        # degree_offset: a per-symptom slope on a per-perturbation covariate (standardised log degree) added to each logit,
        # the sigmoid-head counterpart of the noisy-OR head's degree-dependent leak; starts at zero
        self.covariate_slope = nn.Parameter(torch.zeros(num_symptoms)) if degree_offset else None
        if pooling not in ("sum", "mean"):
            raise ValueError("pooling must be 'sum' or 'mean'")
        self.pooling = pooling
        self.readout = nn.Sequential(nn.Linear(node_state_dim, hidden_dim), nn.ReLU(), nn.Linear(hidden_dim, num_symptoms))

    def forward(self, node_state_field: Tensor, relation_index: int = 0, perturbation_covariate: Tensor | None = None) -> SigmoidHeadOutput:  # noqa: ARG002 - one link matrix only
        pooled = node_state_field.sum(dim=1) if self.pooling == "sum" else node_state_field.mean(dim=1)
        logits = self.readout(pooled)
        if self.covariate_slope is not None and perturbation_covariate is not None:
            logits = logits + perturbation_covariate[:, None] * self.covariate_slope[None, :]
        return SigmoidHeadOutput(symptom_probability=torch.sigmoid(logits))

    def description_length_penalty(self, node_cost: float = 1.0, link_cost: float = 1.0) -> Tensor:  # noqa: ARG002
        return torch.zeros((), device=self.readout[0].weight.device)


class RelationalGnnSigmoidBaseline:
    """Thin fit/predict wrapper kept for the baseline interface; training lives in experiments/run_main_model.py (--head sigmoid)."""

    def fit(self, training_observations, graph):
        raise NotImplementedError("train with experiments/run_main_model.py --head sigmoid")

    def predict(self, perturbations):
        raise NotImplementedError("train with experiments/run_main_model.py --head sigmoid")
