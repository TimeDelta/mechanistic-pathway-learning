"""Can the hard edge gates of decision 9.62 find the false edges of a graph? A simulation with planted truth.

Why a simulation. The registration lets pilots "run on the development set only", and no held-out set has been drawn
on the v4 or v5 tables, so no gate may be fitted to their labels yet. The method's settings (the hard concrete
temperature, the surrogate width, the step counts, the stability threshold) are therefore chosen here, on graphs whose
truth is known, and the real run takes them as fixed.

The world. A random signed graph of `--nodes` nodes, each node fed by a few sources over two relations (activates +1,
inhibits -1). Every edge carries an evidence tier with the v5 graph's proportions and values (configs/edge_prior_tiers
.csv: curated p 0.9 n 1, unreferenced p 0.7 n 0.5, text-mined only p 0.55 n 0.5), so its binary prior is
pi = sigmoid(n logit p), and the truth draws each edge real with probability pi: the priors are calibrated by
construction, which is the most favourable case for them and says nothing about whether the real tiers are. The labels
come from the true graph only: a linear-response encoder on the real edges, perturbations that seed one node each, and
symptoms that read a few nodes each through a logistic link, drawn as Bernoulli.

The fit. The same encoder on the whole database graph, its parameters set to the truth's (the oracle case, which
isolates what the gates can learn) or learned with them (--joint), with hard gates started at the priors; the stochastic
pass; the binary pass; the exact flip checks; then, with --stability-refits R, R refits on random halves of the
perturbations (Meinshausen and Buhlmann 2010) and each edge's closing frequency.

What is reported, against the truth: of the edges closed, the share that are false (precision); of the false edges,
the share closed (recall); the same for the stable closures at each threshold, with the Meinshausen and Buhlmann bound
on the expected number of false closures; and the held-out loss of the fitted graph against the database graph with
every edge open and against the true graph.

Usage:
  PYTHONPATH=. python experiments/simulate_edge_gate_recovery.py --seed 0
  PYTHONPATH=. python experiments/simulate_edge_gate_recovery.py --quick
"""
from __future__ import annotations

import argparse
import copy
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch

from mechanistic_pathway_learning.models.edge_gate_passes import (
    DEFAULT_SURROGATE_WIDTH,
    binary_pass_step,
    exact_flip_checks,
    objective_change_estimates,
    set_gate_state,
)
from mechanistic_pathway_learning.models.edge_gates import HardConcreteEdgeGates, binary_edge_prior
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder

RELATIONS = ["activates", "inhibits"]
# (name, p, n, share of the gated edges): the v5 graph's tiers grouped by their values (docs/figures.md, figure 12)
TIERS = (("curated", 0.9, 1.0, 0.68), ("unreferenced", 0.7, 0.5, 0.20), ("text_mined_only", 0.55, 0.5, 0.12))
TIER_P = {name: p for name, p, _, _ in TIERS}


@dataclass
class World:
    edges: np.ndarray  # [E, 4] source, target, relation, sign
    prior_p: np.ndarray
    prior_n: np.ndarray
    prior: np.ndarray  # pi
    real: np.ndarray  # bool
    seeds: torch.Tensor  # [P, 1]
    injection: torch.Tensor  # [P, 1, 2]
    readout: torch.Tensor  # [S, N]
    bias: torch.Tensor  # [S]
    labels: torch.Tensor  # [P, S]
    true_logits: torch.Tensor


def make_encoder(num_nodes: int, edges: np.ndarray, gates=None, seed: int = 0) -> LinearResponseEncoder:
    features = torch.nn.functional.one_hot(torch.arange(num_nodes) % 2, num_classes=2).float()
    torch.manual_seed(seed)
    return LinearResponseEncoder(num_nodes, RELATIONS, torch.as_tensor(edges[:, 0], dtype=torch.long), torch.as_tensor(edges[:, 1], dtype=torch.long),
                                 torch.as_tensor(edges[:, 2], dtype=torch.long), torch.as_tensor(edges[:, 3], dtype=torch.float32), features, 4,
                                 edge_gates=gates)


def make_world(num_nodes: int, mean_in_degree: float, num_perturbations: int, num_symptoms: int, nodes_per_symptom: int, seed: int,
               label_bias: float = -1.0, truth_real_probability: dict[str, float] | None = None, label_sensitivity: float = 1.0) -> World:
    """The planted world. By default each edge is real with its prior probability pi (calibrated priors) and every
    label is the true Bernoulli draw. truth_real_probability {tier name: probability} replaces pi as the truth's rate in
    the named tiers, so the priors are miscalibrated there while the fit still reads pi; label_sensitivity c observes
    each true positive as positive with probability c and as 0 otherwise (under-ascertained positives, as in a label
    table that records what was reported). Both draw from their own generators, so a world with the defaults is the
    same world the earlier batches used."""
    generator = np.random.default_rng(seed)
    in_degree = np.clip(generator.lognormal(np.log(mean_in_degree), 0.7, size=num_nodes).round().astype(int), 1, num_nodes - 1)
    rows = []
    for target in range(num_nodes):
        sources = generator.choice(np.delete(np.arange(num_nodes), target), size=in_degree[target], replace=False)
        for source in sources:
            relation = int(generator.random() < 0.3)
            rows.append((int(source), target, relation, -1.0 if relation else 1.0))
    edges = np.array(rows, dtype=float)
    tier = generator.choice(len(TIERS), size=len(edges), p=[share for *_, share in TIERS])
    prior_p = np.array([TIERS[index][1] for index in tier])
    prior_n = np.array([TIERS[index][2] for index in tier])
    prior = binary_edge_prior(prior_p, prior_n)
    truth_probability = prior.copy()
    for name, value in (truth_real_probability or {}).items():
        names = [tier_name for tier_name, *_ in TIERS]
        if name not in names or not 0.0 <= value <= 1.0:
            raise ValueError(f"truth_real_probability: {name}={value} needs a tier among {names} and a probability")
        truth_probability[tier == names.index(name)] = value
    real = generator.random(len(edges)) < truth_probability
    truth = make_encoder(num_nodes, edges[real], seed=seed)
    seeds = torch.as_tensor(generator.integers(0, num_nodes, size=(num_perturbations, 1)), dtype=torch.long)
    injection = torch.stack([torch.as_tensor(generator.choice([-1.0, 1.0], size=num_perturbations), dtype=torch.float32),
                             torch.ones(num_perturbations)], dim=1)[:, None, :]
    readout = torch.zeros(num_symptoms, num_nodes)
    for symptom in range(num_symptoms):
        readout[symptom, generator.choice(num_nodes, size=nodes_per_symptom, replace=False)] = torch.as_tensor(generator.normal(0.0, 1.0, nodes_per_symptom), dtype=torch.float32)
    with torch.no_grad():
        field = truth.response(seeds, injection)[:, :, 0]  # [P, N]
        raw = field @ readout.T
        scale = 2.0 / raw.std().clamp_min(1e-8)
        bias = torch.full((num_symptoms,), float(label_bias))
        readout = readout * scale
        logits = field @ readout.T + bias
        labels = torch.bernoulli(torch.sigmoid(logits), generator=torch.Generator().manual_seed(seed))
        if label_sensitivity < 1.0:
            observed = torch.bernoulli(torch.full_like(labels, float(label_sensitivity)), generator=torch.Generator().manual_seed(seed + 7919))
            labels = labels * observed
    return World(edges, prior_p, prior_n, prior, real, seeds, injection, readout, bias, labels, logits)


def batches_of(indices: np.ndarray, batch_size: int, generator: np.random.Generator | None = None) -> list[np.ndarray]:
    order = indices if generator is None else generator.permutation(indices)
    return [order[start:start + batch_size] for start in range(0, len(order), batch_size)]


def loss_of(encoder: LinearResponseEncoder, world: World, rows: np.ndarray, readout: torch.Tensor, bias: torch.Tensor) -> torch.Tensor:
    index = torch.as_tensor(rows, dtype=torch.long)
    field = encoder.response(world.seeds[index], world.injection[index])[:, :, 0]
    logits = field @ readout.T + bias
    return torch.nn.functional.binary_cross_entropy_with_logits(logits, world.labels[index])


def fit(world: World, train_rows: np.ndarray, arguments, seed: int) -> dict:
    """The fit of one data set. Oracle (the default): the encoder's parameters, the readout and the bias are the
    truth's and stay fixed, so only the gates learn. Joint (--joint): they learn with the gates, the way the trainer
    fits them, with AdamW at the trainer's weight decay and early stopping on an inner validation share of the training
    rows (patience as the trainer's), from the truth's values (--joint-start oracle) or from a fresh start (random);
    the readout is learned on every node (--joint-readout dense) or only on the nodes each symptom reads in the truth
    (support), which bounds how far a learned readout can make up for a false edge."""
    num_nodes = int(world.readout.shape[1])
    gated = np.ones(len(world.edges), dtype=bool)
    for name in getattr(arguments, "ungated_tiers", None) or []:  # a tier taken as known: its edges are always open, as an edge label would be
        gated &= world.prior_p != TIER_P[name]
    gates = HardConcreteEdgeGates(gated, world.prior_p, world.prior_n, beta=arguments.beta)
    encoder = make_encoder(num_nodes, world.edges, gates=gates, seed=arguments.world_seed)
    joint = bool(arguments.joint)
    random_start = joint and getattr(arguments, "joint_start", "oracle") == "random"
    if not random_start:
        truth_state = make_encoder(num_nodes, world.edges[world.real], seed=arguments.world_seed).state_dict()
        encoder.load_state_dict({key: value for key, value in truth_state.items()}, strict=False)  # the oracle's parameters
    support = (world.readout != 0).float()
    readout_on_support_only = joint and getattr(arguments, "joint_readout", "dense") == "support"
    if random_start:
        initial = torch.Generator().manual_seed(seed + 31)
        readout = 0.01 * torch.randn(world.readout.shape, generator=initial)
        if readout_on_support_only:
            readout = readout * support
        rate = float(world.labels[torch.as_tensor(train_rows, dtype=torch.long)].mean().clamp(1e-3, 1 - 1e-3))
        bias = torch.full_like(world.bias, float(np.log(rate / (1 - rate))))
    else:
        readout, bias = world.readout.clone(), world.bias.clone()
    readout.requires_grad_(joint)
    bias.requires_grad_(joint)
    if readout_on_support_only:
        readout.register_hook(lambda gradient: gradient * support)
    other_parameters = [parameter for name, parameter in encoder.named_parameters() if not name.startswith("edge_gates.")]
    for parameter in other_parameters:
        parameter.requires_grad_(joint)
    generator = np.random.default_rng(seed)
    torch.manual_seed(seed)
    stopping_rows = np.zeros(0, dtype=train_rows.dtype)
    if joint and arguments.joint_validation_share > 0:  # the inner validation share early stopping reads, as the trainer's
        shuffled = generator.permutation(train_rows)
        stopping_count = int(round(arguments.joint_validation_share * len(train_rows)))
        stopping_rows, train_rows = shuffled[:stopping_count], shuffled[stopping_count:]
    # W, the number of training labels: the batch loss is a mean per label, so with the prior's nats divided by W one
    # nat weighs as one label's log-loss, the scale of the evidence lower bound and of the trainer's W (models/edge_gates.py)
    penalty_denominator = float(len(train_rows) * world.labels.shape[1])
    gate_optimizer = torch.optim.Adam([gates.log_alpha], lr=arguments.gate_learning_rate)
    joint_optimizer = (torch.optim.AdamW(other_parameters + [readout, bias], lr=arguments.learning_rate, weight_decay=arguments.weight_decay)
                       if joint else None)
    history = []
    gates.set_phase("stochastic")
    best = None
    epochs_without_improvement = 0
    for epoch in range(arguments.stochastic_epochs):
        encoder.train()
        for rows in batches_of(train_rows, arguments.batch_size, generator):
            gates.resample()
            loss = loss_of(encoder, world, rows, readout, bias) + gates.kl_penalty() / penalty_denominator
            gate_optimizer.zero_grad(set_to_none=True)
            if joint_optimizer is not None:
                joint_optimizer.zero_grad(set_to_none=True)
            loss.backward()
            gate_optimizer.step()
            if joint_optimizer is not None:
                joint_optimizer.step()
        entry = {"phase": "stochastic", "epoch": epoch, "closed_at_read_out": int((gates.read_out() < 0.5).sum())}
        if len(stopping_rows):  # the read-out graph on the inner validation rows, as the trainer evaluates
            encoder.eval()
            with torch.no_grad():
                stopping_loss = float(loss_of(encoder, world, stopping_rows, readout, bias))
            entry["inner_validation_loss"] = stopping_loss
            if best is None or stopping_loss < best["loss"] - 1e-6:
                best = {"loss": stopping_loss, "epoch": epoch, "encoder": copy.deepcopy(encoder.state_dict()),
                        "readout": readout.detach().clone(), "bias": bias.detach().clone()}
                epochs_without_improvement = 0
            else:
                epochs_without_improvement += 1
        history.append(entry)
        if len(stopping_rows) and epochs_without_improvement >= arguments.patience:
            break
    if best is not None:  # the early-stopped state, gates included (they are in the encoder's state)
        encoder.load_state_dict(best["encoder"])
        with torch.no_grad():
            readout.copy_(best["readout"])
            bias.copy_(best["bias"])
        history.append({"phase": "early_stopping", "best_epoch": best["epoch"], "inner_validation_loss": best["loss"],
                        "closed_at_read_out": int((gates.read_out() < 0.5).sum())})
    for parameter in other_parameters + [readout, bias]:
        parameter.requires_grad_(False)
    binary_optimizer = torch.optim.Adam([gates.log_alpha], lr=arguments.binary_learning_rate)
    data_loss = lambda rows: loss_of(encoder, world, rows, readout, bias)  # noqa: E731
    history.append({"phase": "after_stochastic", "closed_at_read_out": int((gates.read_out() < 0.5).sum())})
    for epoch in range(arguments.binary_epochs):
        flipped = 0
        if arguments.binary_full_batch:  # one step per epoch on the estimate over every training batch: no minibatch noise
            flipped += binary_pass_step(encoder, gates, data_loss, None, penalty_denominator, binary_optimizer, arguments.surrogate_width,
                                        batches=batches_of(train_rows, arguments.batch_size))["flipped"]
        else:
            for rows in batches_of(train_rows, arguments.batch_size, generator):
                flipped += binary_pass_step(encoder, gates, data_loss, rows, penalty_denominator, binary_optimizer, arguments.surrogate_width)["flipped"]
        history.append({"phase": "binary", "epoch": epoch, "flips": flipped, "closed_at_read_out": int((gates.read_out() < 0.5).sum())})
    diagnostic = None
    if arguments.estimator_diagnostic:  # the first-order estimate against the exact change, on random edges
        diagnostic = estimator_quality(encoder, gates, data_loss, batches_of(train_rows, arguments.batch_size), penalty_denominator,
                                       arguments.estimator_diagnostic, arguments.surrogate_width, seed)
    checks = None
    if arguments.flip_checks:
        checks = exact_flip_checks(encoder, gates, data_loss, batches_of(train_rows, arguments.batch_size), penalty_denominator, arguments.flip_checks,
                                   arguments.surrogate_width)
        checks.pop("records")
    gates.set_phase("binary")
    encoder.eval()
    evidence = None
    if arguments.closure_evidence:  # per closed edge, the exact log-likelihood ratio for keeping it closed, in nats
        evidence = closure_evidence(encoder, gates, data_loss, batches_of(train_rows, arguments.batch_size), penalty_denominator, arguments.surrogate_width)
    return {"closed": (~gates.open_edges()).cpu().numpy(), "open_probability": np.nan_to_num(gates.values_per_edge(), nan=1.0), "closure_evidence": evidence,
            "history": history, "checks": checks, "diagnostic": diagnostic, "encoder": encoder, "readout": readout.detach(), "bias": bias.detach()}


def closure_evidence(encoder, gates, data_loss, batches, penalty_denominator, surrogate_width) -> dict:
    """{edge index: nats}: for each closed gate, how much the data part of the objective, times W (so in nats of
    log-likelihood over the training labels), rises when that edge alone is reopened; positive supports the closure.
    Keyed by the edge's index in the graph (gates.gated_edges), so it can be compared with the truth directly."""
    def data_objective() -> float:
        with torch.no_grad():
            return sum(float(data_loss(batch)) for batch in batches) / len(batches)

    base = data_objective()
    evidence = {}
    for gate_index in torch.nonzero(gates.read_out() < 0.5).flatten().tolist():
        previous = float(gates.log_alpha[gate_index])
        set_gate_state(gates, gate_index, True, 4.0 * surrogate_width)
        evidence[int(gates.gated_edges[gate_index])] = (data_objective() - base) * penalty_denominator
        with torch.no_grad():
            gates.log_alpha[gate_index] = previous
    return evidence


def estimator_quality(encoder, gates, data_loss, batches, penalty_denominator, sample_size, surrogate_width, seed) -> dict:
    """For sample_size random gated edges: the first-order estimate of the objective change if each flipped alone, and
    the exact change, with their correlation and how often their signs agree."""
    estimate = objective_change_estimates(encoder, gates, data_loss, batches, penalty_denominator)[gates.gated_edges]

    def objective() -> float:
        with torch.no_grad():
            return sum(float(data_loss(batch)) for batch in batches) / len(batches) + float(gates.kl_penalty()) / penalty_denominator

    base = objective()
    chosen = np.random.default_rng(seed).choice(len(gates.gated_edges), size=min(sample_size, len(gates.gated_edges)), replace=False)
    estimates, exact = [], []
    for gate_index in chosen:
        gate_index = int(gate_index)
        was_open = bool(gates.read_out()[gate_index] > 0.5)
        previous = float(gates.log_alpha[gate_index])
        set_gate_state(gates, gate_index, not was_open, 4.0 * surrogate_width)
        exact.append(objective() - base)
        estimates.append(float(estimate[gate_index]))
        with torch.no_grad():
            gates.log_alpha[gate_index] = previous
    estimates, exact = np.array(estimates), np.array(exact)
    data_part_estimates = estimates - np.array([float(gates.flip_prior_change()[gates.gated_edges[int(index)]]) / penalty_denominator for index in chosen])
    data_part_exact = exact - (estimates - data_part_estimates)
    return {"edges": int(len(chosen)), "correlation": float(np.corrcoef(estimates, exact)[0, 1]),
            "data_part_correlation": float(np.corrcoef(data_part_estimates, data_part_exact)[0, 1]),
            "sign_agreement": float(np.mean(np.sign(estimates) == np.sign(exact))),
            "median_absolute_error_over_median_absolute_exact": float(np.median(np.abs(estimates - exact)) / max(np.median(np.abs(exact)), 1e-12)),
            # per edge, in nats over the training labels (times W): the data part of the estimate and of the exact change
            "per_edge_data_part_nats": [[int(gates.gated_edges[int(index)]), float(estimate * penalty_denominator), float(exact_change * penalty_denominator)]
                                        for index, estimate, exact_change in zip(chosen, data_part_estimates, data_part_exact)]}


def closure_scores(closed: np.ndarray, real: np.ndarray) -> dict:
    false = ~real
    return {"closed": int(closed.sum()), "false_edges": int(false.sum()),
            "precision": float((closed & false).sum() / closed.sum()) if closed.sum() else None,
            "recall": float((closed & false).sum() / false.sum()) if false.sum() else None}


def held_out_loss(world: World, rows: np.ndarray, closed: np.ndarray, encoder_parameters: dict, readout, bias, seed: int) -> float:
    num_nodes = int(world.readout.shape[1])
    encoder = make_encoder(num_nodes, world.edges[~closed], seed=seed)
    encoder.load_state_dict({key: value for key, value in encoder_parameters.items() if not key.startswith("edge_gates.")}, strict=False)
    encoder.eval()
    with torch.no_grad():
        return float(loss_of(encoder, world, rows, readout, bias))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--world-seed", type=int, default=0)
    parser.add_argument("--nodes", type=int, default=300)
    parser.add_argument("--mean-in-degree", type=float, default=6.0)
    parser.add_argument("--perturbations", type=int, default=800)
    parser.add_argument("--symptoms", type=int, default=12)
    parser.add_argument("--nodes-per-symptom", type=int, default=6)
    parser.add_argument("--held-out-share", type=float, default=0.25)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--beta", type=float, default=2.0 / 3.0)
    parser.add_argument("--gate-learning-rate", type=float, default=0.05)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--binary-learning-rate", type=float, default=0.1)
    parser.add_argument("--surrogate-width", type=float, default=DEFAULT_SURROGATE_WIDTH)
    parser.add_argument("--stochastic-epochs", type=int, default=20)
    parser.add_argument("--binary-epochs", type=int, default=8)
    parser.add_argument("--flip-checks", type=int, default=40)
    parser.add_argument("--binary-full-batch", action="store_true", help="one binary-pass step per epoch on the estimate over all training batches")
    parser.add_argument("--estimator-diagnostic", type=int, default=0, help="compare the flip estimate with the exact change on this many random edges")
    parser.add_argument("--closure-evidence", action="store_true", help="for each closed edge, the exact log-likelihood ratio for keeping it closed, and precision by it")
    parser.add_argument("--per-edge-detail", action="store_true", help="write each closed edge (tier, prior, truth, evidence) and each stability frequency above 0")
    parser.add_argument("--joint", action="store_true", help="learn the encoder's parameters, the readout and the bias with the gates (see fit)")
    parser.add_argument("--joint-start", choices=("oracle", "random"), default="oracle")
    parser.add_argument("--joint-readout", choices=("dense", "support"), default="dense")
    parser.add_argument("--joint-validation-share", type=float, default=0.2, help="share of the training rows early stopping reads in a joint fit")
    parser.add_argument("--weight-decay", type=float, default=1e-4, help="AdamW weight decay of a joint fit (the trainer's default)")
    parser.add_argument("--patience", type=int, default=8, help="epochs without an improvement before a joint fit stops (the trainer's default)")
    parser.add_argument("--label-bias", type=float, default=-1.0, help="the symptoms' logit bias, which sets the positive rate")
    parser.add_argument("--truth-real-probability", nargs="*", default=[], metavar="TIER=P",
                        help="the truth's rate of real edges in a tier, in place of its prior pi (miscalibrated priors)")
    parser.add_argument("--label-sensitivity", type=float, default=1.0, help="each true positive is observed as positive with this probability")
    parser.add_argument("--ungated-tiers", nargs="*", default=[], choices=[name for name, *_ in TIERS],
                        help="tiers whose edges carry no gate: taken as known and always open, as edge labels would be (the owner's question of 16 October)")
    parser.add_argument("--stability-scheme", choices=("halves", "complementary"), default="halves",
                        help="independent random halves (Meinshausen and Buhlmann 2010) or complementary pairs (Shah and Samworth 2013)")
    parser.add_argument("--stability-refits", type=int, default=0)
    parser.add_argument("--stability-thresholds", type=float, nargs="*", default=[0.6, 0.7, 0.8, 0.9])
    parser.add_argument("--quick", action="store_true", help="a small world and few epochs, for a smoke run")
    parser.add_argument("--output", type=Path, default=None)
    arguments = parser.parse_args()
    if arguments.quick:
        arguments.nodes, arguments.perturbations, arguments.symptoms = 60, 160, 6
        arguments.stochastic_epochs, arguments.binary_epochs, arguments.flip_checks = 3, 2, 5
    started = time.time()
    truth_real_probability = {}
    for entry in arguments.truth_real_probability:
        name, _, value = entry.partition("=")
        truth_real_probability[name.strip()] = float(value)
    world = make_world(arguments.nodes, arguments.mean_in_degree, arguments.perturbations, arguments.symptoms, arguments.nodes_per_symptom, arguments.world_seed,
                       label_bias=arguments.label_bias, truth_real_probability=truth_real_probability, label_sensitivity=arguments.label_sensitivity)
    rows = np.random.default_rng(arguments.world_seed + 1).permutation(arguments.perturbations)
    held_out_count = int(round(arguments.held_out_share * arguments.perturbations))
    held_out_rows, train_rows = rows[:held_out_count], rows[held_out_count:]
    result = fit(world, train_rows, arguments, seed=arguments.world_seed)
    summary = {
        "world": {"nodes": arguments.nodes, "edges": int(len(world.edges)), "false_edges": int((~world.real).sum()),
                  "false_share_by_tier": {name: float((~world.real)[world.prior_p == p].mean()) for name, p, _, _ in TIERS if (world.prior_p == p).any()},
                  "perturbations": arguments.perturbations, "symptoms": arguments.symptoms, "positive_rate": float(world.labels.mean()),
                  "expected_true_positive_rate": float(torch.sigmoid(world.true_logits).mean()),
                  "training_labels_per_edge": float(len(train_rows) * arguments.symptoms / len(world.edges)),
                  "training_positives_per_edge": float(world.labels[torch.as_tensor(train_rows, dtype=torch.long)].sum() / len(world.edges))},
        "settings": {key: value for key, value in vars(arguments).items() if key not in ("output",)},
        "fit": closure_scores(result["closed"], world.real),
        "fit_by_tier": {name: closure_scores(result["closed"][world.prior_p == p], world.real[world.prior_p == p]) for name, p, _, _ in TIERS},
        "history": result["history"], "flip_checks": result["checks"], "estimator_diagnostic": result["diagnostic"],
        "closures_by_evidence": None if result["closure_evidence"] is None else {
            f"at least {threshold} nats": closure_scores(np.isin(np.arange(len(world.edges)), [index for index, nats in result["closure_evidence"].items() if nats >= threshold]),
                                                       world.real)
            for threshold in (0.0, 0.5, 1.0, 2.0, 3.0)},
    }
    encoder_parameters = result["encoder"].state_dict()
    summary["held_out_loss"] = {
        "fitted_graph": held_out_loss(world, held_out_rows, result["closed"], encoder_parameters, result["readout"], result["bias"], arguments.world_seed),
        "every_edge_open": held_out_loss(world, held_out_rows, np.zeros(len(world.edges), dtype=bool), encoder_parameters, result["readout"], result["bias"], arguments.world_seed),
        "true_graph": held_out_loss(world, held_out_rows, ~world.real, encoder_parameters, result["readout"], result["bias"], arguments.world_seed),
    }
    if arguments.stability_refits:
        closures = []
        for refit in range(arguments.stability_refits):
            if arguments.stability_scheme == "complementary":  # refit 2k and 2k+1 are the two halves of one split
                order = np.random.default_rng(1000 + refit // 2).permutation(train_rows)
                half = order[:len(order) // 2] if refit % 2 == 0 else order[len(order) // 2:2 * (len(order) // 2)]
            else:
                half = np.random.default_rng(1000 + refit).choice(train_rows, size=len(train_rows) // 2, replace=False)
            closures.append(fit(world, half, arguments, seed=1000 + refit)["closed"])
        frequency = np.mean(closures, axis=0)
        average_selected = float(np.mean([closed.sum() for closed in closures]))
        summary["stability"] = {"refits": arguments.stability_refits, "average_closed_per_refit": average_selected, "thresholds": {}}
        for threshold in arguments.stability_thresholds:
            stable = frequency >= threshold
            bound = average_selected ** 2 / ((2 * threshold - 1) * len(world.edges)) if threshold > 0.5 else None
            summary["stability"]["thresholds"][str(threshold)] = {**closure_scores(stable, world.real), "expected_false_closures_bound": bound}
    if arguments.per_edge_detail:  # what the figures and the record read edge by edge
        tier_name = {p: name for name, p, _, _ in TIERS}
        evidence = result["closure_evidence"] or {}
        summary["closed_edges"] = [{"edge": int(edge), "tier": tier_name[float(world.prior_p[edge])], "prior": float(world.prior[edge]), "real": bool(world.real[edge]),
                                    "evidence_nats": evidence.get(int(edge)), "open_probability": float(result["open_probability"][edge])}
                                   for edge in np.flatnonzero(result["closed"])]
        if arguments.stability_refits:
            summary["stability_frequencies"] = [{"edge": int(edge), "tier": tier_name[float(world.prior_p[edge])], "real": bool(world.real[edge]), "frequency": float(frequency[edge])}
                                                for edge in np.flatnonzero(frequency > 0)]
    summary["seconds"] = round(time.time() - started, 1)
    text = json.dumps(summary, indent=1, default=float)
    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
