"""The numbers behind decision 9.62 (how a gate is read, and whether it is hard): hard concrete gates (Louizos, Welling
and Kingma 2018, arXiv 1712.01312) read at their initialisation, and the edge priors of decision 9.57 under three
readings of a prior for one yes-or-no fact.

Hard concrete. A gate is z = min(1, max(0, s (zeta - gamma) + gamma)) with s = sigmoid((logit(u) + log alpha) / beta),
u uniform; the paper sets gamma = -0.1, zeta = 1.1 and beta = 2/3 (section 4). P(z != 0) = sigmoid(log alpha - beta
log(-gamma / zeta)) is its eq. 12, and its test-time gate is eq. 13, min(1, max(0, sigmoid(log alpha) (zeta - gamma) +
gamma)). A gate is started here with P(z != 0) equal to its edge's prior.

The three binary priors of a (p, n) pair, n the pseudo-count of 9.57: p itself; the mean of the Beta(np + 1, n(1 - p) + 1)
density the soft gates' penalty is the negative log of; and sigmoid(n logit p), which gives p at n = 1.

The per-edge L0 cost. With q the probability the gate is open and pi its prior, KL(Bernoulli(q) || Bernoulli(pi)) =
q log((1 - pi) / pi) - H(q) - log(1 - pi): an L0 cost of log((1 - pi) / pi) nats per open edge, less the gate's entropy.

Usage:
  PYTHONPATH=. python experiments/measure_edge_gate_readings.py --graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v4
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.edge_evidence import edge_priors, read_edge_tiers

HARD_CONCRETE_GAMMA, HARD_CONCRETE_ZETA, HARD_CONCRETE_BETA = -0.1, 1.1, 2.0 / 3.0


def sigmoid(value):
    return 1.0 / (1.0 + np.exp(-np.asarray(value, dtype=float)))


def logit(probability):
    probability = np.asarray(probability, dtype=float)
    return np.log(probability / (1.0 - probability))


def log_alpha_for_open_probability(open_probability, beta: float = HARD_CONCRETE_BETA, gamma: float = HARD_CONCRETE_GAMMA, zeta: float = HARD_CONCRETE_ZETA):
    """The log alpha at which P(z != 0) (Louizos et al., eq. 12) equals open_probability."""
    return logit(open_probability) + beta * np.log(-gamma / zeta)


def hard_concrete_masses(log_alpha, beta: float = HARD_CONCRETE_BETA, gamma: float = HARD_CONCRETE_GAMMA, zeta: float = HARD_CONCRETE_ZETA) -> dict[str, np.ndarray]:
    """P(z = 0), P(0 < z < 1), P(z = 1) and the test-time gate (eq. 13) of a hard concrete gate."""
    closed_below = -gamma / (zeta - gamma)  # s at or below it gives z = 0
    open_above = (1.0 - gamma) / (zeta - gamma)  # s at or above it gives z = 1
    probability_closed = sigmoid(beta * logit(closed_below) - log_alpha)
    probability_open = sigmoid(log_alpha - beta * logit(open_above))
    test_time_gate = np.clip(sigmoid(log_alpha) * (zeta - gamma) + gamma, 0.0, 1.0)
    return {"closed": probability_closed, "fractional": 1.0 - probability_closed - probability_open, "open": probability_open, "test_time_gate": test_time_gate}


def sample_hard_concrete(log_alpha, num_samples: int, generator: np.random.Generator, beta: float = HARD_CONCRETE_BETA, gamma: float = HARD_CONCRETE_GAMMA,
                         zeta: float = HARD_CONCRETE_ZETA) -> np.ndarray:
    uniform = generator.uniform(1e-12, 1.0 - 1e-12, size=(num_samples,) + np.shape(log_alpha))
    stretched = sigmoid((logit(uniform) + log_alpha) / beta) * (zeta - gamma) + gamma
    return np.clip(stretched, 0.0, 1.0)


def binary_priors(probability, pseudo_count) -> dict[str, np.ndarray]:
    """The prior that an edge is real under the three readings of (p, n) put to the owner in decision 9.62."""
    probability, pseudo_count = np.asarray(probability, dtype=float), np.asarray(pseudo_count, dtype=float)
    return {"p": probability, "beta_mean": (pseudo_count * probability + 1.0) / (pseudo_count + 2.0), "log_odds_times_n": sigmoid(pseudo_count * logit(probability))}


def bernoulli_kl(open_probability, prior):
    open_probability, prior = np.asarray(open_probability, dtype=float), np.asarray(prior, dtype=float)
    return open_probability * np.log(open_probability / prior) + (1.0 - open_probability) * np.log((1.0 - open_probability) / (1.0 - prior))


def measure(graph_directory: Path, tiers_path: Path) -> dict:
    evidence = pd.read_parquet(graph_directory / "edge_evidence.parquet")
    gated, probability, pseudo_count = edge_priors(evidence, read_edge_tiers(tiers_path))
    probability, pseudo_count = probability[gated], pseudo_count[gated]
    readings = binary_priors(probability, pseudo_count)
    quantiles = [0.05, 0.25, 0.5, 0.75, 0.95]
    settings = pd.DataFrame({"p": np.round(probability, 3), "n": np.round(pseudo_count, 3)}).value_counts().rename("edges").reset_index()
    settings["closing_cost_nats_log_odds_times_n"] = settings.n * logit(settings.p)
    initial = {f"{prior:.2f}": {key: float(value) for key, value in hard_concrete_masses(log_alpha_for_open_probability(prior)).items()} for prior in (0.55, 0.6, 0.7, 0.85, 0.9)}
    return {
        "graph": str(graph_directory), "gated_edges": int(gated.sum()),
        "edges_with_p_below_one_half": int((probability < 0.5).sum()),
        "prior_settings": settings.to_dict(orient="records"),
        "binary_prior_quantiles": {name: dict(zip(map(str, quantiles), np.round(np.quantile(values, quantiles), 3).tolist())) for name, values in readings.items()},
        "beta_mean_log_odds_share_of_logit_p_median": float(np.median(logit(readings["beta_mean"]) / logit(probability))),
        "hard_concrete_at_initialisation": initial,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders_drugs_v4"))
    parser.add_argument("--edge-prior-tiers", type=Path, default=Path("configs/edge_prior_tiers.csv"))
    arguments = parser.parse_args()
    print(json.dumps(measure(arguments.graph_dir, arguments.edge_prior_tiers), indent=1))


if __name__ == "__main__":
    main()
