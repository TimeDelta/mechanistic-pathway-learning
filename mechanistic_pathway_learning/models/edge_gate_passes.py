"""The binary pass over hard edge gates and its exact checks (decision 9.62 (e); models/edge_gates.py has the gates).

After the stochastic pass, the owner's second pass (16 October 2026: "Why not include that open probability threshold
as a second training pass over just that parameter or something? Can be given a proxy gradient that resembles the
cliff discontinuity") trains the gates alone. Every forward reads the binary read-out, so the graph that is fitted is
the graph that will be read. The proxy gradient is built from the cliff: per gated edge, D_e, the estimated change in
the objective if that edge alone flipped,

    D_e = (first-order change in the data loss, flip_effect_estimates) + (exact change in the prior term) / W,

and apply_flip_gradient turns D_e into a gradient on the gate's logit through a steep sigmoid at the threshold, so a
gate near the threshold moves toward the side that lowers the objective and a gate far from it barely moves.

Because the data part is first order, the pass ends with exact checks (exact_flip_checks): the edges whose estimate
promises the largest decrease are flipped one at a time, the objective is evaluated exactly with each, and a flip is
kept only when the exact objective falls. That leaves the read-out a local optimum under single flips among the
candidates checked, and the record of estimate against exact change says how far the first-order estimate can be
trusted on the graph at hand.

Both functions take the objective as a closure, so the trainer and the simulation study (experiments/
simulate_edge_gate_recovery.py) drive them the same way. data_loss(batch) must return the mean loss the trainer
minimises for that batch, without the gates' penalty; penalty_denominator is the W the gates' penalty is divided by.
"""
from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np
import torch
from torch import Tensor

from mechanistic_pathway_learning.models.edge_gates import HardConcreteEdgeGates

DEFAULT_SURROGATE_WIDTH = 0.5  # in open log-odds; chosen in experiments/simulate_edge_gate_recovery.py
DEFAULT_FLIP_MARGIN = 4.0  # a checked flip moves the open log-odds this far past the threshold, in surrogate widths


def objective_change_estimates(encoder, gates: HardConcreteEdgeGates, data_loss: Callable[[object], Tensor], batches: Iterable,
                               penalty_denominator: float) -> Tensor:
    """[num_edges] D_e over the given batches: the data part averaged over them (each batch's loss is a mean, so the
    average estimates the full-data change), plus the exact prior part divided by the penalty denominator."""
    gates.set_phase("binary")
    encoder.collect_flip_gradients = True
    total, count = None, 0
    try:
        for batch in batches:
            encoder.reset_flip_collection()
            loss = data_loss(batch)
            loss.backward()
            estimate = encoder.flip_effect_estimates().detach()
            total = estimate if total is None else total + estimate
            count += 1
    finally:
        encoder.collect_flip_gradients = False
        encoder.reset_flip_collection()
    if count == 0:
        raise ValueError("no batch to estimate the flips on")
    return total / count + gates.flip_prior_change().detach() / penalty_denominator


def binary_pass_step(encoder, gates: HardConcreteEdgeGates, data_loss: Callable[[object], Tensor], batch, penalty_denominator: float,
                     optimizer: torch.optim.Optimizer, surrogate_width: float = DEFAULT_SURROGATE_WIDTH, batches: list | None = None) -> dict:
    """One step of the binary pass: estimate D_e on the batch (or averaged over `batches`, which removes the minibatch
    noise from the estimate), set the gates' proxy gradient, step. Only the gates' parameters should be in the
    optimizer; the caller freezes the rest."""
    optimizer.zero_grad(set_to_none=True)
    before = gates.read_out().clone()
    change = objective_change_estimates(encoder, gates, data_loss, batches if batches is not None else [batch], penalty_denominator)
    gates.log_alpha.grad = None
    gates.apply_flip_gradient(change, surrogate_width)
    optimizer.step()
    after = gates.read_out()
    return {"flipped": int((before != after).sum()), "promised_decrease": float((-change[gates.gated_edges]).clamp_min(0).sum())}


def set_gate_state(gates: HardConcreteEdgeGates, gate_index: int, is_open: bool, margin: float) -> None:
    """Move one gate's open log-odds to +margin or -margin, so its read-out is open or closed."""
    with torch.no_grad():
        offset = float(gates.log_alpha[gate_index] - gates.open_log_odds()[gate_index])
        gates.log_alpha[gate_index] = (margin if is_open else -margin) + offset


def exact_flip_checks(encoder, gates: HardConcreteEdgeGates, data_loss: Callable[[object], Tensor], batches: list, penalty_denominator: float,
                      budget: int, surrogate_width: float = DEFAULT_SURROGATE_WIDTH, margin_in_widths: float = DEFAULT_FLIP_MARGIN) -> dict:
    """Flip, one at a time, the `budget` gated edges whose D_e promises the largest decrease; keep a flip when the exact
    objective (the mean of data_loss over all the batches plus the binary penalty over W) falls. Returns the accepted
    flips and, for every checked edge, its estimate and its exact change."""
    gates.set_phase("binary")

    def objective() -> float:
        with torch.no_grad():
            data = sum(float(data_loss(batch)) for batch in batches) / len(batches)
            return data + float(gates.kl_penalty()) / penalty_denominator

    change = objective_change_estimates(encoder, gates, data_loss, batches, penalty_denominator)
    gated_change = change[gates.gated_edges]
    read = gates.read_by_encoder
    candidates = [int(index) for index in torch.argsort(gated_change) if gated_change[index] < 0 and bool(read[index])][:budget]
    current = objective()
    records, accepted = [], []
    margin = margin_in_widths * surrogate_width
    for gate_index in candidates:
        was_open = bool(gates.read_out()[gate_index] > 0.5)
        previous_log_alpha = float(gates.log_alpha[gate_index])
        set_gate_state(gates, gate_index, not was_open, margin)
        trial = objective()
        exact = trial - current
        records.append({"gate": gate_index, "edge": int(gates.gated_edges[gate_index]), "was_open": was_open,
                        "estimated_change": float(gated_change[gate_index]), "exact_change": exact})
        if exact < 0:
            current = trial
            accepted.append(gate_index)
        else:
            with torch.no_grad():
                gates.log_alpha[gate_index] = previous_log_alpha
    estimated = torch.tensor([record["estimated_change"] for record in records]) if records else torch.zeros(0)
    exact_changes = torch.tensor([record["exact_change"] for record in records]) if records else torch.zeros(0)
    agreement = float(torch.corrcoef(torch.stack([estimated, exact_changes]))[0, 1]) if len(records) > 2 else None
    return {"checked": len(records), "accepted": len(accepted), "accepted_gates": accepted, "objective": current,
            "estimate_exact_correlation": agreement, "records": records}


def closure_evidence(encoder, gates: HardConcreteEdgeGates, data_loss: Callable[[object], Tensor], batches: list, penalty_denominator: float,
                     budget: int | None = None, surrogate_width: float = DEFAULT_SURROGATE_WIDTH, margin_in_widths: float = DEFAULT_FLIP_MARGIN) -> dict:
    """Each closed edge's exact evidence for its closure, and the posterior that it is false.

    For a gated edge the read-out closes (and the encoder reads), E is how much the data part of the objective (the mean
    of data_loss over the batches), times W, rises when that edge alone is reopened: nats of the fit's weighted
    log-likelihood, positive where the data support the closure. With pi the prior that the edge is real, the posterior
    log-odds that it is false, given every other edge as read out and the fitted parameters, is E - logit(pi), so
    P(false) = sigmoid(E - logit(pi)). At a local optimum under single flips every closed edge has E > logit(pi), so
    P(false) > 1/2; a closure below one half is one the passes left where reopening it alone would lower the objective.

    Why it is reported (docs/edge_gate_simulation.md): in planted worlds the closures with E of 1 nat or more were false
    far more often than the read-out's closures as a whole, so a closure's E says how far to believe it. Each closed edge
    costs one pass over the batches; with a budget, the closures against the strongest priors come first, since those
    are the ones a reader would most want checked. The gates are left as they were."""
    gates.set_phase("binary")

    def data_objective() -> float:
        with torch.no_grad():
            return sum(float(data_loss(batch)) for batch in batches) / len(batches)

    closed = [int(index) for index in torch.nonzero((gates.read_out() < 0.5) & gates.read_by_encoder).flatten()]
    closed.sort(key=lambda index: -float(gates.prior_open[index]))
    chosen = closed if budget is None else closed[:max(int(budget), 0)]
    base = data_objective() if chosen else None
    margin = margin_in_widths * surrogate_width
    records = []
    for gate_index in chosen:
        previous = float(gates.log_alpha[gate_index])
        set_gate_state(gates, gate_index, True, margin)
        evidence = (data_objective() - base) * penalty_denominator
        with torch.no_grad():
            gates.log_alpha[gate_index] = previous
        prior = float(gates.prior_open[gate_index])
        log_odds_false = evidence - float(np.log(prior / (1.0 - prior)))
        records.append({"gate": gate_index, "edge": int(gates.gated_edges[gate_index]), "prior": prior, "evidence_nats": evidence,
                        "posterior_false": float(1.0 / (1.0 + np.exp(-log_odds_false)))})
    return {"closed": len(closed), "evaluated": len(records), "records": records}
