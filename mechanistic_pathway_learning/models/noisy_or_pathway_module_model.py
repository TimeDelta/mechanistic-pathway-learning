"""Noisy-OR pathway module head (experiment design, section 5.3).

A perturbation is encoded upstream as a node-state field of shape
[batch_size, num_graph_nodes, node_state_dim]. Each of the K pathway
modules owns a sparse support over graph nodes (one hard-concrete gate
per node), reads the gated field through a gated pool and produces an
activation in [0, 1]. Two poolings are available: "sum" (default) adds the gated node
states, so a perturbation that reaches a few support nodes is read at full strength and
a support that nothing reaches reads as zero; "mean" divides by the support mass, which
with the absolute field makes the activation a property of the support rather than of
the perturbation (the base states of thousands of gated nodes swamp a change at a few).
The sum pooling is meant for the perturbation difference field (encoder docstring). Module activations combine with learned
module-to-symptom link probabilities through a noisy-OR with a
per-symptom leak term:

    P(symptom s | perturbation p)
        = 1 - (1 - leak_s) * prod_k (1 - link_{k,s} * activation_k(p))

Any single module whose link times activation is near 1 explains the
symptom on its own, so the functional form encodes equifinality without a
cross-module interaction term. The relation index selects a separate link
matrix for "induces" (0) and "relieves" (1).

The hard-concrete gate follows Louizos, Welling and Kingma (2018), which
gives a differentiable relaxation of an L0 penalty so that the expected
number of active support nodes is available for the description-length
criterion in design section 5.5.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import torch
from torch import Tensor, nn

HARD_CONCRETE_STRETCH_LOW = -0.1
HARD_CONCRETE_STRETCH_HIGH = 1.1
PROBABILITY_EPSILON = 1e-6


class HardConcreteNodeGate(nn.Module):
    """One learnable [0, 1] gate per (pathway module, graph node)."""

    def __init__(
        self,
        num_pathway_modules: int,
        num_graph_nodes: int,
        temperature: float = 2.0 / 3.0,
        initial_log_alpha: float = -1.0,
        initial_log_alpha_noise: float = 0.01,
    ) -> None:
        super().__init__()
        self.temperature = temperature
        initial_values = torch.full((num_pathway_modules, num_graph_nodes), initial_log_alpha)
        initial_values = initial_values + initial_log_alpha_noise * torch.randn_like(initial_values)
        self.log_alpha = nn.Parameter(initial_values)

    def forward(self) -> Tensor:
        """Return gate values of shape [num_pathway_modules, num_graph_nodes]."""
        if self.training:
            uniform_noise = torch.rand_like(self.log_alpha).clamp(PROBABILITY_EPSILON, 1 - PROBABILITY_EPSILON)
            logistic_noise = torch.log(uniform_noise) - torch.log1p(-uniform_noise)
            relaxed_gate = torch.sigmoid((logistic_noise + self.log_alpha) / self.temperature)
        else:
            relaxed_gate = torch.sigmoid(self.log_alpha)
        stretched_gate = relaxed_gate * (HARD_CONCRETE_STRETCH_HIGH - HARD_CONCRETE_STRETCH_LOW) + HARD_CONCRETE_STRETCH_LOW
        return stretched_gate.clamp(0.0, 1.0)

    def expected_active_node_count(self) -> Tensor:
        """Expected number of nonzero gates per module (the L0 surrogate), shape [num_pathway_modules]."""
        shift = self.temperature * math.log(-HARD_CONCRETE_STRETCH_LOW / HARD_CONCRETE_STRETCH_HIGH)
        probability_nonzero = torch.sigmoid(self.log_alpha - shift)
        return probability_nonzero.sum(dim=1)


@dataclass
class NoisyOrOutput:
    symptom_probability: Tensor  # [batch_size, num_symptoms]
    module_activation: Tensor  # [batch_size, num_pathway_modules]
    link_probability: Tensor  # [num_pathway_modules, num_symptoms] for the requested relation
    leak_probability: Tensor  # [num_symptoms] for the requested relation
    module_support: Tensor  # [num_pathway_modules, num_graph_nodes]

    def module_contribution(self) -> Tensor:
        """Per-sample contribution link_{k,s} * activation_k, shape [batch_size, K, S].

        This is the natural attribution of a noisy-OR model and is what the
        equifinality sufficiency test (design section 6.5) uses to assign a
        held-out perturbation to its dominant module.
        """
        return self.module_activation[:, :, None] * self.link_probability[None, :, :]

    def dominant_module_per_perturbation(self, symptom_index: int) -> Tensor:
        """Index of the module with the largest contribution to one symptom, shape [batch_size]."""
        return self.module_contribution()[:, :, symptom_index].argmax(dim=1)


class NoisyOrPathwayModuleHead(nn.Module):
    """K sparse pathway modules combined by a noisy-OR into symptom probabilities."""

    def __init__(
        self,
        num_graph_nodes: int,
        node_state_dim: int,
        num_pathway_modules: int,
        num_symptoms: int,
        num_relation_types: int = 2,
        initial_link_logit: float = -3.0,
        initial_leak_logit: float = -4.0,
        gate_temperature: float = 2.0 / 3.0,
        gate_initial_log_alpha: float = -1.0,
        pooling: str = "sum",
        gate_initial_log_alpha_noise: float = 0.01,
        initial_readout_bias: float = 0.0,
    ) -> None:
        """gate_initial_log_alpha_noise breaks the symmetry between modules (all gates start at the same
        log-alpha otherwise, and identical modules stay identical); initial_readout_bias below zero makes a
        module silent unless the pooled perturbation signal drives it, the off-by-default reading of a pathway."""
        super().__init__()
        if pooling not in ("sum", "mean"):
            raise ValueError("pooling must be 'sum' or 'mean'")
        self.pooling = pooling
        self.num_graph_nodes = num_graph_nodes
        self.num_pathway_modules = num_pathway_modules
        self.num_symptoms = num_symptoms
        self.num_relation_types = num_relation_types
        self.support_gate = HardConcreteNodeGate(
            num_pathway_modules,
            num_graph_nodes,
            temperature=gate_temperature,
            initial_log_alpha=gate_initial_log_alpha,
            initial_log_alpha_noise=gate_initial_log_alpha_noise,
        )
        self.module_readout_weight = nn.Parameter(torch.randn(num_pathway_modules, node_state_dim) / math.sqrt(node_state_dim))
        self.module_readout_bias = nn.Parameter(torch.full((num_pathway_modules,), float(initial_readout_bias)))
        self.module_symptom_link_logit = nn.Parameter(
            torch.full((num_relation_types, num_pathway_modules, num_symptoms), initial_link_logit)
        )
        self.symptom_leak_logit = nn.Parameter(torch.full((num_relation_types, num_symptoms), initial_leak_logit))

    def module_support(self) -> Tensor:
        return self.support_gate()

    def module_activation(self, node_state_field: Tensor, module_support: Tensor | None = None) -> Tensor:
        """Activation of each module for each perturbation in the batch, shape [batch_size, K]."""
        if module_support is None:
            module_support = self.support_gate()
        gated_pool = torch.einsum("kn,bnd->bkd", module_support, node_state_field)
        if self.pooling == "mean":
            gated_pool = gated_pool / (module_support.sum(dim=1)[None, :, None] + PROBABILITY_EPSILON)
        activation_logit = (gated_pool * self.module_readout_weight[None, :, :]).sum(dim=-1) + self.module_readout_bias
        return torch.sigmoid(activation_logit)

    @torch.no_grad()
    def initialize_leak_from_base_rates(self, symptom_base_rates: Tensor, relation_index: int = 0) -> None:
        """Set the leak of each symptom to its training base rate, so the untrained model predicts the base rate.

        The leak absorbs unmodelled causes; starting it at the base rate instead of near zero removes the
        early epochs in which every positive is penalized against a probability of a few percent.
        """
        clamped = symptom_base_rates.clamp(PROBABILITY_EPSILON, 1.0 - PROBABILITY_EPSILON)
        self.symptom_leak_logit[relation_index] = torch.log(clamped) - torch.log1p(-clamped)

    def link_probability(self, relation_index: int = 0) -> Tensor:
        return torch.sigmoid(self.module_symptom_link_logit[relation_index])

    def leak_probability(self, relation_index: int = 0) -> Tensor:
        return torch.sigmoid(self.symptom_leak_logit[relation_index])

    def forward(self, node_state_field: Tensor, relation_index: int = 0) -> NoisyOrOutput:
        module_support = self.support_gate()
        module_activation = self.module_activation(node_state_field, module_support)
        link_probability = self.link_probability(relation_index)
        leak_probability = self.leak_probability(relation_index)
        symptom_probability = noisy_or_combination(module_activation, link_probability, leak_probability)
        return NoisyOrOutput(
            symptom_probability=symptom_probability,
            module_activation=module_activation,
            link_probability=link_probability,
            leak_probability=leak_probability,
            module_support=module_support,
        )

    def description_length_penalty(self, node_cost: float = 1.0, link_cost: float = 1.0) -> Tensor:
        """Model-cost half of the minimum description length criterion (design section 5.5).

        Returns node_cost * E[active support nodes] + link_cost * E[nonzero links], where
        the link count is relaxed to the sum of link probabilities. The data-cost half is
        the weighted negative log-likelihood computed by the training loop.
        """
        expected_active_nodes = self.support_gate.expected_active_node_count().sum()
        expected_nonzero_links = torch.sigmoid(self.module_symptom_link_logit).sum()
        return node_cost * expected_active_nodes + link_cost * expected_nonzero_links


def noisy_or_combination(module_activation: Tensor, link_probability: Tensor, leak_probability: Tensor) -> Tensor:
    """Combine module activations [B, K], links [K, S] and leaks [S] into symptom probabilities [B, S].

    Computed in log space: log(1 - P) = log(1 - leak) + sum_k log(1 - link * activation).
    """
    per_module_failure = 1.0 - (module_activation[:, :, None] * link_probability[None, :, :]).clamp(max=1.0 - PROBABILITY_EPSILON)
    log_probability_no_symptom = torch.log1p(-leak_probability.clamp(max=1.0 - PROBABILITY_EPSILON))[None, :] + torch.log(per_module_failure).sum(dim=1)
    return 1.0 - torch.exp(log_probability_no_symptom)
