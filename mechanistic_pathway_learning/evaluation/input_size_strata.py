"""Strata of how many mechanism targets a drug has, for the input-size reading (docs/label_source_review.md, 13.1).

A drug enters with magnitude 1 at each mechanism target, so its total input is the number of targets curated for it.
The user's decision of 11 October 2026 keeps that input rule and counters the confound in the evaluation instead. This
module is the second of the three controls: report the metric inside strata of the number of targets, so that an
advantage which lives only among the many-target drugs is visible rather than pooled away.

The first control is the input-size baseline (models/baselines/input_size_baseline.py), the third the permuted-target
control (experiments/run_main_model.py --permute-drug-targets).

A gene perturbation has one seed of magnitude 1 and no choice of target count, so it is its own stratum rather than
being pooled with the one-target drugs; otherwise the gene rows, which outnumber the drugs six to one, would decide
every stratum they fell in. A combination (evidence/combined_perturbations.py) is its own stratum for the same reason
read the other way: its seed count is a drug's targets plus a lowered gene, which is not a target count.
"""
from __future__ import annotations

import numpy as np

from mechanistic_pathway_learning.evidence.combined_perturbations import COMBINATION_PERTURBATION_TYPE

GENE_PERTURBATION_STRATUM = "gene_perturbation"
DRUG_STRATUM_PREFIX = "drug_targets_"
COMBINATION_STRATUM = "combined_perturbation"  # a drug carried with a lowered gene or a second drug: its own stratum


def target_count_of(perturbation_seeds, perturbation_types) -> np.ndarray:
    """The number of seeds of each perturbation, which for a drug is its number of mechanism target nodes.

    On a graph with drug entry nodes a drug is seeded on its own node, so the seed count is 1 for every drug; there the
    count has to come from the drug's targets, which the loader keeps as perturbation_target_seeds."""
    del perturbation_types
    return np.array([len(seeds) for seeds in perturbation_seeds], dtype=int)


def input_size_strata(perturbation_types: list[str], target_counts: np.ndarray, largest_named_count: int = 3) -> tuple[list[str], list[str]]:
    """One stratum label per perturbation, and the labels in reading order.

    Counts up to largest_named_count get a stratum each; above it they share one, because the drugs with more targets
    are few (8 of 230 without the limit of three) and a stratum of one or two drugs scores nothing."""
    labels: list[str] = []
    for kind, count in zip(perturbation_types, np.asarray(target_counts, dtype=int)):
        if kind == COMBINATION_PERTURBATION_TYPE:
            # Its input is a drug's targets together with a lowered gene, so it belongs in neither the gene stratum nor
            # a drug's target-count stratum: a count of 2 there would read as a two-target drug.
            labels.append(COMBINATION_STRATUM)
        elif kind != "drug":
            labels.append(GENE_PERTURBATION_STRATUM)
        elif count <= largest_named_count:
            labels.append(f"{DRUG_STRATUM_PREFIX}{count}")
        else:
            labels.append(f"{DRUG_STRATUM_PREFIX}{largest_named_count + 1}_or_more")
    present = set(labels)
    order = [GENE_PERTURBATION_STRATUM] + [f"{DRUG_STRATUM_PREFIX}{count}" for count in range(1, largest_named_count + 1)] + [f"{DRUG_STRATUM_PREFIX}{largest_named_count + 1}_or_more", COMBINATION_STRATUM]
    return labels, [label for label in order if label in present]


def strata_of_scored_rows(data, scored_rows: np.ndarray, largest_named_count: int = 3) -> tuple[list[str], list[str]]:
    """The strata of the rows a score was read on. Uses the drug's targets, so a graph with drug entry nodes, where a
    drug is seeded on its own node, gives the same strata as the graph without them."""
    seeds = data.perturbation_target_seeds if getattr(data, "perturbation_target_seeds", None) is not None else data.perturbation_seeds
    counts = target_count_of(seeds, data.perturbation_types)
    rows = np.asarray(scored_rows)
    rows = np.flatnonzero(rows) if rows.dtype == bool else rows
    return input_size_strata([data.perturbation_types[row] for row in rows], counts[rows], largest_named_count)
