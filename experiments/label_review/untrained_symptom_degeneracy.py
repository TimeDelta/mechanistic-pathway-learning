"""Why a per-symptom cross-validation fold is not available in the noisy-OR head: a symptom held out of training
keeps its links at their common initial value, and every such symptom then ranks the perturbations identically.

The head is

    P(symptom s | perturbation p) = 1 - (1 - leak_s) * prod_k (1 - link_{k,s} * activation_k(p)),

so with every link of symptom s left at the initial logit the link probability is one constant c for all k, and

    P(s | p) = 1 - (1 - leak_s) * Q(p),    Q(p) = prod_k (1 - c * activation_k(p)).

Q does not depend on the symptom. Two untrained symptoms therefore differ by a positive affine map of one another,
which fixes their Pearson correlation at 1 and their ranking of the perturbations identically, whatever their leaks:
the leak is the only per-symptom quantity left, and it scales rather than reorders. Q falls as any activation rises,
so the shared ranking is by total module activation, which is the degree-like quantity that the input-size and
degree controls of the registered design already account for.

This script states that identity and measures it, so the claim in preregistration items 25b2 and 26 is reproducible
without a trained run. It trains nothing and reads no data: the node-state field is random, because the identity
holds for any field.

Usage:
  python experiments/label_review/untrained_symptom_degeneracy.py
  python experiments/label_review/untrained_symptom_degeneracy.py --json-output docs/untrained_symptom_degeneracy.json
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import torch

from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import NoisyOrPathwayModuleHead

REGISTERED_PATHWAY_MODULES = 8
REGISTERED_SYMPTOMS = 24
RELATION_INDEX_INDUCES = 0


def smallest_correlation_between_columns(symptom_probabilities: np.ndarray) -> float:
    correlation = np.corrcoef(symptom_probabilities.T)
    off_diagonal = ~np.eye(correlation.shape[0], dtype=bool)
    return float(correlation[off_diagonal].min())

def largest_ranking_disagreement(symptom_probabilities: np.ndarray) -> int:
    """Largest number of perturbations placed differently by two symptoms' rankings."""
    orders = np.argsort(np.argsort(-symptom_probabilities, axis=0), axis=0)
    return int(max((orders[:, first] != orders[:, second]).sum()
                   for first in range(orders.shape[1]) for second in range(first + 1, orders.shape[1])))


def untrained_head(num_perturbations: int, num_graph_nodes: int, node_state_dim: int, seed: int,
                   base_rates: torch.Tensor) -> tuple[NoisyOrPathwayModuleHead, torch.Tensor]:
    torch.manual_seed(seed)
    head = NoisyOrPathwayModuleHead(num_graph_nodes=num_graph_nodes, node_state_dim=node_state_dim,
                                    num_pathway_modules=REGISTERED_PATHWAY_MODULES, num_symptoms=REGISTERED_SYMPTOMS)
    head.initialize_leak_from_base_rates(base_rates, relation_index=RELATION_INDEX_INDUCES)
    head.eval()
    return head, torch.randn(num_perturbations, num_graph_nodes, node_state_dim)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--perturbations", type=int, default=200)
    parser.add_argument("--graph-nodes", type=int, default=400)
    parser.add_argument("--node-state-dim", type=int, default=4)
    parser.add_argument("--seeds", type=int, nargs="*", default=[0, 1, 2])
    parser.add_argument("--json-output", type=str)
    arguments = parser.parse_args()

    # Leaks that differ per symptom, including the zero-positive case the clamp sends to 1e-6: if the ranking still
    # agrees across every pair, the leak is what the algebra says it is, a scale rather than an order.
    base_rates = torch.cat([torch.zeros(1), torch.linspace(0.01, 0.30, REGISTERED_SYMPTOMS - 1)])

    record = {"registered_pathway_modules": REGISTERED_PATHWAY_MODULES, "registered_symptoms": REGISTERED_SYMPTOMS,
              "perturbations": arguments.perturbations, "per_seed": []}
    for seed in arguments.seeds:
        head, node_state_field = untrained_head(arguments.perturbations, arguments.graph_nodes,
                                                arguments.node_state_dim, seed, base_rates)
        with torch.no_grad():
            untrained = head(node_state_field, relation_index=RELATION_INDEX_INDUCES).symptom_probability.numpy()
            # One symptom trained away from the initial value, to show the degeneracy is the untrained state and
            # not a property of the head at every setting.
            head.module_symptom_link_logit.data[RELATION_INDEX_INDUCES, :, 0] = torch.linspace(-4.0, 2.0, REGISTERED_PATHWAY_MODULES)
            one_trained = head(node_state_field, relation_index=RELATION_INDEX_INDUCES).symptom_probability.numpy()
        record["per_seed"].append({
            "seed": seed,
            "smallest_correlation_between_untrained_symptoms": smallest_correlation_between_columns(untrained),
            "largest_ranking_disagreement_between_untrained_symptoms": largest_ranking_disagreement(untrained),
            "smallest_correlation_once_one_symptom_is_trained": smallest_correlation_between_columns(one_trained)})

    for row in record["per_seed"]:
        print(f"seed {row['seed']}: smallest correlation between untrained symptoms "
              f"{row['smallest_correlation_between_untrained_symptoms']:.6f}, perturbations ranked differently by any "
              f"two of them {row['largest_ranking_disagreement_between_untrained_symptoms']}, smallest correlation once "
              f"one symptom is trained {row['smallest_correlation_once_one_symptom_is_trained']:.6f}")
    if arguments.json_output:
        with open(arguments.json_output, "w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=1)
            handle.write("\n")


if __name__ == "__main__":
    main()
