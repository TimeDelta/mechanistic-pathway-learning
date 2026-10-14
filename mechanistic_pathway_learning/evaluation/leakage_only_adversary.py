"""The leakage-only adversary: a predictor given the leaked labels and nothing else (open decision 9.35).

Why it exists. The owner, 14 October 2026: "instead of doing the final approval on leakage after the draw, let's just
set an upper bound on the acceptable leakage. how can i do that in a non-arbitrary way?" A ceiling stated in counts of
overlapping targets is not checkable against anything, because nothing says what a count is worth. A ceiling stated in
the units of the test is: it asks how much of the registered reading a predictor could earn **by lookup alone**, and
compares that with the 0.05 difference H1 already requires and with the measured precision of the reading.

What the adversary is allowed to know. Exactly two things: which training perturbations share a seed node with the
held-out perturbation, and the training labels. No propagation, no edge, no descriptor, no distance beyond the shared
node itself. It is handed the leaked labels directly, so **its advantage is an upper bound on what any model could
extract from this channel by lookup**, which is the property that makes it usable as a ceiling.

Which channel this is, and which is deliberately excluded. Two perturbations sharing a seed node perturb the **same
node**, so a label transfers between them by identity and a lookup substitutes for mechanism. That is leakage. Two
perturbations whose seeds sit in one pathway are a different matter: the claim under test is that a mechanistic route
carries a symptom, so transfer along a route is the hypothesis rather than a leak, and scoring it here would bound the
study's own subject. Only the shared-node channel is measured. The excluded channel is named in the output so the
registered channel list says what it does not cover.

What item 25 closes and what it leaves open. `target_held_out_split` holds out every drug sharing a held-out target and
the gene perturbation of a held-out target, so neither reaches training. It does **not** hold out the gene perturbations
of a held-out drug's **other** targets, nor the drugs that share those: a drug with three targets of which one is held
out keeps two mechanisms in training, which the registration states as the price of keeping the label exact. So this
channel is open by design and the adversary is what measures its worth. `channel_breakdown` reports the two sides
separately rather than asserting either is empty.
"""
from __future__ import annotations

import numpy as np

# How a neighbour's labels are pooled. The strongest of the three is the bound, since the adversary is meant to be the
# best lookup available rather than a particular one. The default shrinks the neighbour mean towards the training base
# rate by the weight of evidence behind it, which matters more than it sounds: a raw neighbour mean over one or two
# neighbours is near-binary, and a near-binary score ranks badly whatever information it holds, so an unshrunk
# adversary understates the leak through its score scale rather than through the channel. A shrunk score is also
# continuous with the fallback, since a perturbation with no neighbour lands on the base rate exactly.
POOLING_RULES = ("base_rate_shrunk_mean", "shared_node_weighted_mean", "nearest_neighbour_maximum")
# weight given to the base rate in the shrunk mean, in units of one shared node: one neighbour sharing one node moves
# the score halfway from the base rate to that neighbour's label
BASE_RATE_PRIOR_STRENGTH = 1.0


def seed_nodes_of_perturbations(perturbation_ids, perturbation_seeds, node_ids) -> dict[str, frozenset[str]]:
    """perturbation id -> the node ids it seeds, for every perturbation regardless of type."""
    names = list(node_ids)
    return {str(perturbation_id): frozenset(names[int(index)] for index in seeds if 0 <= int(index) < len(names))
            for perturbation_id, seeds in zip(perturbation_ids, perturbation_seeds)}


def leaked_neighbours(held_out_ids, training_ids, seed_nodes_of_perturbation) -> dict[str, dict[str, int]]:
    """held-out perturbation id -> {training perturbation id: how many seed nodes they share}, shared nodes only.

    A held-out perturbation with no entry has no lookup available and the adversary falls back to the training base
    rate for it, which is what keeps the adversary from scoring *below* popularity on those rows and so keeps its
    advantage an upper bound rather than a mixture of a leak and a handicap.
    """
    training_by_node: dict[str, list[str]] = {}
    for training_id in training_ids:
        for node in seed_nodes_of_perturbation.get(str(training_id), ()):
            training_by_node.setdefault(node, []).append(str(training_id))
    neighbours: dict[str, dict[str, int]] = {}
    for held_out_id in held_out_ids:
        shared_node_count: dict[str, int] = {}
        for node in seed_nodes_of_perturbation.get(str(held_out_id), ()):
            for training_id in training_by_node.get(node, ()):
                shared_node_count[training_id] = shared_node_count.get(training_id, 0) + 1
        if shared_node_count:
            neighbours[str(held_out_id)] = shared_node_count
    return neighbours


def channel_breakdown(neighbours, perturbation_types_of_id: dict[str, str]) -> dict[str, int]:
    """Which side of the open channel each leaked pair comes through, so the registered channel list is counted rather
    than asserted: a training **drug** sharing a target, or the training **gene** perturbation of one of those targets."""
    counts = {"held_out_with_a_leaked_neighbour": len(neighbours), "pairs_through_a_training_drug": 0,
              "pairs_through_a_training_gene": 0, "pairs_through_another_type": 0}
    for shared_node_count in neighbours.values():
        for training_id in shared_node_count:
            kind = perturbation_types_of_id.get(training_id, "")
            key = {"drug": "pairs_through_a_training_drug", "gene": "pairs_through_a_training_gene"}.get(kind, "pairs_through_another_type")
            counts[key] += 1
    return counts


class LeakageOnlyAdversary:
    """Predicts a held-out perturbation's symptoms from the labels of the training perturbations sharing a seed node.

    Exposes fit/predict like the baselines of `models/baselines`, so the same reading scores it: the advantage it earns
    over the same baseline, in the same macro and micro AUPRC, is the leakage expressed in the units of the test.
    """

    def __init__(self, pooling_rule: str = "base_rate_shrunk_mean", base_rate_prior_strength: float = BASE_RATE_PRIOR_STRENGTH) -> None:
        if pooling_rule not in POOLING_RULES:
            raise ValueError(f"pooling_rule must be one of {POOLING_RULES}, not {pooling_rule!r}")
        self.pooling_rule = pooling_rule
        self.base_rate_prior_strength = float(base_rate_prior_strength)
        self.symptom_base_rate: np.ndarray | None = None
        self.training_kept_positive: np.ndarray | None = None
        self.training_row_of_id: dict[str, int] = {}

    def fit(self, training_ids, training_outcomes: np.ndarray, training_label_mask: np.ndarray | None = None) -> "LeakageOnlyAdversary":
        """training_outcomes [num_training, num_symptoms] in {0, 1}; with a label mask the base rate and the looked-up
        labels both run over the labelled pairs only, as in the trainer's loss and in PopularityBaseline."""
        labelled = np.ones(training_outcomes.shape, dtype=bool) if training_label_mask is None else np.asarray(training_label_mask, dtype=bool)
        self.training_kept_positive = np.asarray(training_outcomes, dtype=float) * labelled
        self.symptom_base_rate = self.training_kept_positive.sum(axis=0) / np.maximum(labelled.sum(axis=0), 1)
        self.training_row_of_id = {str(training_id): row for row, training_id in enumerate(training_ids)}
        return self

    def predict(self, held_out_ids, neighbours: dict[str, dict[str, int]]) -> np.ndarray:
        """[num_held_out, num_symptoms]. A held-out perturbation with no leaked neighbour takes the training base rate."""
        if self.symptom_base_rate is None or self.training_kept_positive is None:
            raise RuntimeError("fit before predict")
        scores = np.tile(self.symptom_base_rate, (len(held_out_ids), 1))
        for row, held_out_id in enumerate(held_out_ids):
            shared_node_count = neighbours.get(str(held_out_id))
            if not shared_node_count:
                continue
            rows_of_neighbours, weights = [], []
            for training_id, shared in shared_node_count.items():
                training_row = self.training_row_of_id.get(training_id)
                if training_row is not None:
                    rows_of_neighbours.append(training_row)
                    weights.append(float(shared))
            if not rows_of_neighbours:
                continue
            neighbour_labels = self.training_kept_positive[rows_of_neighbours]
            weight_vector = np.asarray(weights, dtype=float)
            if self.pooling_rule == "nearest_neighbour_maximum":
                scores[row] = neighbour_labels[int(np.argmax(weights))]
            elif self.pooling_rule == "shared_node_weighted_mean":
                scores[row] = (neighbour_labels * weight_vector[:, None]).sum(axis=0) / weight_vector.sum()
            else:  # base_rate_shrunk_mean
                looked_up = (neighbour_labels * weight_vector[:, None]).sum(axis=0)
                total_weight = weight_vector.sum() + self.base_rate_prior_strength
                scores[row] = (looked_up + self.base_rate_prior_strength * self.symptom_base_rate) / total_weight
        return scores


def seen_target_share(targets_of_drug: dict[str, set[str]], held_out_ids) -> float:
    """The label-free statistic the real draw is checked against: over the held-out drugs, the share of their mechanism
    targets that a training drug also has.

    It reads no label, so it can be computed on the real draw without opening the lockbox, which is what lets a ceiling
    calibrated on development pseudo-draws be applied to the draw that counts.
    """
    held = {str(held_out_id) for held_out_id in held_out_ids}
    trained_targets = {target for drug, own in targets_of_drug.items() if drug not in held for target in own}
    seen, total = 0, 0
    for drug, own in targets_of_drug.items():
        if drug in held:
            seen += len(own & trained_targets)
            total += len(own)
    return float(seen) / total if total else float("nan")
