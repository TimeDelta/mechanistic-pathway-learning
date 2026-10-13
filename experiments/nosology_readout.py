"""The nosology readout of preregistration item 26: cluster the symptoms by the mechanisms the model gives them, and
report the three readings that say whether the grouping means anything.

Descriptive. Nothing here is a held-out score, nothing enters H1, H2 or P1, and the output says so on its face. The
object is each symptom's column of module-to-symptom link probabilities from a fit on all the evidence, which is that
symptom's profile over mechanisms in the model's own basis.

What it computes, each fixed in item 26 before any confirmatory run:

- **the clustering** (26d): one minus the Pearson correlation between symptoms over the concatenated induces and
  relieves profile, average-linkage hierarchical clustering, and the number of clusters taken by the stability
  criterion of Ben-Hur, Elisseeff and Guyon (2002) rather than chosen after the solution is seen. That criterion also
  returns the absence of structure, so the readout can report that there is no grouping to report.
- **stability** (26 d2): over complementary-pair subsamples of the perturbations, compared at the level of the
  partition and not of the profile vectors - the items clustered are the same symptoms in every subsample, so two
  clusterings compare by one-to-one Jaccard over symptom membership and no coordinate has to mean the same thing in
  both. Reads the fitted runs it is given, one per subsample; with one run it reports the clustering and no stability.
- **the label co-occurrence baseline** (26e), which decides whether the readout says anything: the same procedure on
  the symptom columns of the raw outcome matrix. A cluster the baseline also returns is a grouping the labels already
  carry, so it is not reported as mechanism-derived. No threshold is drawn; the agreement is reported cluster by
  cluster.
- **link concentration** (26 d1): max_s link_ks / sum_s link_ks per module, beside the same statistic from the
  permuted-label runs, which is what concentration looks like when no symptom leans on anything. Reported as a
  continuous quantity with no line drawn, because this analysis makes no decision.

Usage:
  python experiments/nosology_readout.py --run-dirs runs/confirmatory/<model>/<split> ... \\
      --permuted-run-dirs runs/confirmatory/<model>/<split>_permuted ... \\
      --evidence-dir data/processed/evidence_full_v4 \\
      --label-selection data/processed/label_selection/better_v2_full_v4.parquet \\
      --markdown-output docs/nosology_readout.md
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

MINIMUM_POSITIVES_TO_CLUSTER = 5  # the registered condition; thin symptoms are placed afterwards and carry no claim
CANDIDATE_CLUSTER_COUNTS = range(2, 9)
STABILITY_SUBSAMPLE_FRACTION = 0.8  # of the symptoms, which are the objects clustered
STABILITY_DRAWS = 100
CHANCE_DRAWS = 50


def correlation_distances(symptom_profiles: np.ndarray) -> np.ndarray:
    """One minus the Pearson correlation between symptoms, as a square distance matrix.

    A symptom whose profile is constant has no correlation with anything, and the two ways of filling that in are not
    equivalent: treating it as correlation 1 puts every such symptom at distance 0 from every other, so the flat
    symptoms form a cluster of their own that reads as a nosological group and is an artefact of their having no
    profile (item 26f). They are given correlation 0, so distance 1, which leaves them far from everything including
    each other."""
    centred = symptom_profiles - symptom_profiles.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(centred, axis=0)
    usable = norms > 0
    correlation = np.zeros((symptom_profiles.shape[1], symptom_profiles.shape[1]))
    if usable.any():
        scaled = centred[:, usable] / norms[usable]
        correlation[np.ix_(usable, usable)] = scaled.T @ scaled
    np.fill_diagonal(correlation, 1.0)
    return np.clip(1.0 - correlation, 0.0, 2.0)


def cluster_at(distances: np.ndarray, num_clusters: int) -> np.ndarray:
    return fcluster(linkage(squareform(distances, checks=False), method="average"), num_clusters, criterion="maxclust")


def partition_agreement(first_labels: np.ndarray, second_labels: np.ndarray) -> float:
    """Mean over the first partition's clusters of the best Jaccard against any cluster of the second, matched one to
    one. The items are the same symptoms in both, so no coordinate has to be aligned first."""
    first_clusters = [set(np.flatnonzero(first_labels == label)) for label in np.unique(first_labels)]
    second_clusters = [set(np.flatnonzero(second_labels == label)) for label in np.unique(second_labels)]
    jaccard = np.array([[len(a & b) / len(a | b) for b in second_clusters] for a in first_clusters])
    from scipy.optimize import linear_sum_assignment
    rows, columns = linear_sum_assignment(-jaccard)
    return float(jaccard[rows, columns].sum() / max(len(first_clusters), len(second_clusters)))


def chance_agreement(partitions: list[np.ndarray], draws: int, generator: np.random.Generator) -> np.ndarray:
    """The null distribution of the mean pairwise agreement: these partitions with their items shuffled, redrawn.

    The null shuffles each observed partition's labels across the items, which keeps its cluster sizes exactly and
    destroys only which symptom is in which cluster. Both halves of that matter, and the second was learned the
    hard way. Dividing by the cluster count makes the statistic's chance level fall as k rises, so an uncorrected
    argmax reads that fall as evidence for k=2: measured, 20 noise matrices of the real artifact's shape chose k=2
    twenty times out of twenty. Correcting against *uniform* random labels then failed the other way, because
    average linkage returns unbalanced partitions - often one large cluster and some singletons - and two
    unbalanced partitions agree far more than two uniform ones, so noise still read as structure at k=7 or 8 in
    every seed. Only a null that keeps the observed sizes prices that in."""
    null_means = []
    for _ in range(draws):
        shuffled = [partition[generator.permutation(len(partition))] for partition in partitions]
        null_means.append(float(np.mean([partition_agreement(first, second) for first, second in combinations(shuffled, 2)])))
    return np.array(null_means)


def cluster_count_by_stability(profiles_per_fit: list[np.ndarray], seed: int = 0) -> tuple[int | None, dict[int, float]]:
    """The number of clusters most stable across **refits**, above chance, or None for no structure (Ben-Hur,
    Elisseeff and Guyon 2002: stability "provides a means of rationally defining an optimum number of clusters,
    and can also detect the lack of structure in data").

    The resampling has to be the one that refits the model, which is the perturbation subsampling of item 26d2.
    Neither resampling available inside a single fit works, and that is worth recording because both are the
    obvious thing to reach for. Subsampling the mechanism coordinates perturbs the distances but estimates a 24 by
    24 correlation matrix from four observations. Subsampling the symptoms perturbs nothing at all: the distance
    between two symptoms is computed over the mechanism coordinates, so dropping a third symptom leaves it
    unchanged, and two subsamples of pure noise then agree well enough to be read as structure - measured, that
    version picked seven or eight clusters on noise and never reported the absence of structure. Only a refit
    moves the profiles.

    So k is chosen by the same resampling that validates the clustering, and with fewer than two fits it cannot be
    chosen at all: the function returns None and the caller says so rather than printing a number."""
    generator = np.random.default_rng(seed)
    if len(profiles_per_fit) < 2:
        return None, {}
    num_symptoms = profiles_per_fit[0].shape[1]
    agreement_above_chance, beats_every_shuffle = {}, {}
    for num_clusters in CANDIDATE_CLUSTER_COUNTS:
        if num_clusters >= num_symptoms:
            continue
        partitions = [cluster_at(correlation_distances(profiles), num_clusters) for profiles in profiles_per_fit]
        observed = float(np.mean([partition_agreement(first, second)
                                  for first, second in combinations(partitions, 2)]))
        null_means = chance_agreement(partitions, CHANCE_DRAWS, generator)
        agreement_above_chance[num_clusters] = observed - float(null_means.mean())
        beats_every_shuffle[num_clusters] = bool(observed > null_means.max())
    # The bar is "above every shuffle", which is a permutation test at the resolution of the draws and not a
    # chosen level. It is needed because the selection takes a maximum over seven candidate counts: with the bar
    # at the null MEAN, a maximum of seven differences each centred on zero is positive almost always, and noise
    # was selected in nine seeds of ten.
    clearing = {count: value for count, value in agreement_above_chance.items() if beats_every_shuffle[count]}
    if not clearing:
        return None, agreement_above_chance
    return max(clearing, key=clearing.get), agreement_above_chance


def link_concentration(link_probabilities: np.ndarray) -> np.ndarray:
    """max_s link_ks / sum_s link_ks per module: 1 when a module's links are on one symptom, 1/num_symptoms when they
    are spread evenly. A module with no link mass is reported as NaN rather than as concentrated."""
    totals = link_probabilities.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(totals > 0, link_probabilities.max(axis=1) / totals, np.nan)


def symptom_profiles_of_run(run_directory: Path) -> tuple[np.ndarray, np.ndarray, int]:
    """(profiles [mechanism coordinate, symptom], link probabilities [module, symptom], relations recorded).

    Item 26a fixes the profile as the symptom's link column over 8 modules and 2 relations, stacked, so a symptom
    is described by what induces it and what relieves it together. A run that recorded only one relation - which
    every run does today, because run_main_model writes `head.link_probability(0)` - gives a profile of half that,
    and the count is returned so the report can say so instead of describing two relations it did not read."""
    results = json.loads((run_directory / "results.json").read_text())
    links = np.asarray(results["module_symptom_links"], dtype=float)
    if links.ndim == 2:
        links = links[None, :, :]
    return np.concatenate(list(links), axis=0), links.sum(axis=0), int(links.shape[0])


def kept_positive_counts(evidence_directory: Path, label_selection: Path | None, symptoms: list[str],
                         perturbations: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """(outcome matrix over every perturbation of the table, kept positives per symptom).

    The matrix keeps **every** perturbation, including those with no kept positive at all. An earlier version
    built its rows from the perturbations that survived the filter, which dropped 784 of 1,885 rows on
    `evidence_full_v4`; Pearson correlation between symptom columns reads those all-zero rows, so dropping them
    changes the baseline that 26e uses to decide whether the readout says anything.
    """
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet", columns=["perturbation_id", "symptom", "relation"])
    records = records[records.relation == "induces"]
    if label_selection is not None:
        selection = pd.read_parquet(label_selection)
        kept_pairs = set(zip(selection[selection.keep].perturbation_id, selection[selection.keep].symptom))
        records = records[[pair in kept_pairs for pair in zip(records.perturbation_id, records.symptom)]]
    row_of = {perturbation: row for row, perturbation in enumerate(perturbations)}
    column_of = {symptom: column for column, symptom in enumerate(symptoms)}
    outcomes = np.zeros((len(perturbations), len(symptoms)))
    for record in records.itertuples():
        if record.symptom in column_of and record.perturbation_id in row_of:
            outcomes[row_of[record.perturbation_id], column_of[record.symptom]] = 1.0
    return outcomes, outcomes.sum(axis=0)


def symptoms_thin_enough_to_exclude(positives_per_symptom: np.ndarray) -> np.ndarray:
    """Item 26f: the clustering is over the symptoms with at least the registered five kept positives. A symptom
    below that is placed afterwards at its nearest cluster and reported with its count and no claim, because a
    profile fitted on one or two positives is close to every other such profile and they would otherwise group
    together as an artefact of their counts."""
    return positives_per_symptom < MINIMUM_POSITIVES_TO_CLUSTER


def nearest_cluster_of(profiles: np.ndarray, partition: np.ndarray, clustered: np.ndarray, thin: np.ndarray) -> dict[int, int]:
    """The cluster each thin symptom is nearest to, by mean distance to that cluster's members. Reported beside
    the clustering, never inside it."""
    distances = correlation_distances(profiles)
    placements = {}
    for symptom in np.flatnonzero(thin):
        by_cluster = {label: float(np.mean(distances[symptom, clustered[partition == label]]))
                      for label in np.unique(partition)}
        placements[int(symptom)] = int(min(by_cluster, key=by_cluster.get))
    return placements


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dirs", type=Path, nargs="+", required=True,
                        help="the full-data fit first, then the perturbation-subsample refits of 26 d2; the cluster count and the "
                             "stability are both read across the refits, so with one run neither can be computed")
    parser.add_argument("--permuted-run-dirs", type=Path, nargs="*", default=[], help="permuted-label runs: the null for the link concentration")
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--label-selection", type=Path, required=True,
                        help="required: without it a pair the table records as an absence (no positive report) would count as a positive in the baseline")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    arguments = parser.parse_args()

    full_fit = json.loads((arguments.run_dirs[0] / "results.json").read_text())
    symptoms = list(full_fit["symptoms"])
    profiles, link_probabilities, relations_in_the_fit = symptom_profiles_of_run(arguments.run_dirs[0])
    profiles_per_fit = []
    for directory in arguments.run_dirs:
        its_symptoms = list(json.loads((directory / "results.json").read_text())["symptoms"])
        if its_symptoms != symptoms:
            raise SystemExit(f"{directory} scored symptoms {its_symptoms[:3]}... which are not those of {arguments.run_dirs[0]}; "
                             f"the clusterings would be compared over mismatched columns")
        profiles_per_fit.append(symptom_profiles_of_run(directory)[0])

    perturbations = sorted(set(pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet", columns=["perturbation_id"]).perturbation_id))
    baseline_profiles, positives_per_symptom = kept_positive_counts(arguments.evidence_dir, arguments.label_selection, symptoms, perturbations)
    thin = symptoms_thin_enough_to_exclude(positives_per_symptom)
    clustered = np.flatnonzero(~thin)
    if len(clustered) < 2:
        raise SystemExit(f"only {len(clustered)} symptoms reach {MINIMUM_POSITIVES_TO_CLUSTER} kept positives; there is nothing to cluster")

    num_clusters, stability_by_count = cluster_count_by_stability([fit[:, clustered] for fit in profiles_per_fit], arguments.seed)
    record = {"reading": "descriptive; no held-out score, in no hypothesis of the decision rule (preregistration item 26)",
              "relations_in_the_fit": relations_in_the_fit,
              "symptoms_clustered": [symptoms[i] for i in clustered],
              "symptoms_below_the_positive_condition": {symptoms[i]: int(positives_per_symptom[i]) for i in np.flatnonzero(thin)},
              "stability_above_chance_by_cluster_count": stability_by_count,
              "num_clusters": None if num_clusters is None else int(num_clusters)}
    lines = ["# Nosology readout (generated by experiments/nosology_readout.py)", "",
             "**Descriptive.** Fitted on all the evidence, so no reading here is a held-out score and none enters H1, H2 or "
             "P1 (preregistration item 26b). It carries no claim that the grouping predicts anything, none about a symptom "
             "outside the set and none that it would recur in another cohort.", ""]
    if relations_in_the_fit == 1:
        lines += ["**The profile covers one relation, not two.** Item 26a fixes it as 8 modules by 2 relations, and the run "
                  "recorded only `induces`, so the relieves direction is absent from everything below.", ""]
    if thin.any():
        lines += [f"Clustered over the {len(clustered)} symptoms with at least {MINIMUM_POSITIVES_TO_CLUSTER} kept positives "
                  f"(item 26f). Below that and placed afterwards with no claim: "
                  + ", ".join(f"{symptoms[i]} ({int(positives_per_symptom[i])})" for i in np.flatnonzero(thin)) + ".", ""]
    if num_clusters is None:
        lines += ["**No structure.** No cluster count was more stable across the refits than its own shuffled null, which is "
                  "the second half of the criterion of (benhur2002stability) and a finding rather than a failure. Nothing "
                  "below is reported as a grouping.", ""]
        if len(arguments.run_dirs) < 2:
            lines += ["One run was given, so the count could not be chosen at all: both it and the stability are read across "
                      "the subsample refits of item 26 d2.", ""]
    else:
        partition = cluster_at(correlation_distances(profiles[:, clustered]), num_clusters)
        baseline_partition = cluster_at(correlation_distances(baseline_profiles[:, clustered]), num_clusters)
        placements = nearest_cluster_of(profiles, partition, clustered, thin)
        record["clusters"] = {str(label): [symptoms[clustered[i]] for i in np.flatnonzero(partition == label)] for label in np.unique(partition)}
        record["label_co_occurrence_clusters"] = {str(label): [symptoms[clustered[i]] for i in np.flatnonzero(baseline_partition == label)] for label in np.unique(baseline_partition)}
        record["agreement_with_the_label_co_occurrence_baseline"] = partition_agreement(partition, baseline_partition)
        record["nearest_cluster_of_a_thin_symptom"] = {symptoms[symptom]: int(label) for symptom, label in placements.items()}
        lines += [f"Clusters: **{num_clusters}**, the count most stable across {len(profiles_per_fit)} fits above its shuffled "
                  "null; agreement above that null by count: "
                  + ", ".join(f"{count} ({value:+.3f})" for count, value in stability_by_count.items()) + ".", ""]
        for label, members in record["clusters"].items():
            lines.append(f"- **cluster {label}**: {', '.join(members)}")
        lines += ["", f"Agreement with the label co-occurrence baseline at the same count: "
                  f"**{record['agreement_with_the_label_co_occurrence_baseline']:.3f}**. A cluster the baseline also returns is "
                  "a grouping the labels already carry and is not reported as mechanism-derived; what the graph adds is the "
                  "difference between the two partitions, and that difference is this readout's only positive claim.", ""]
        if placements:
            lines.append("Thin symptoms at their nearest cluster, with no claim attached: "
                         + ", ".join(f"{symptoms[s]} -> {label}" for s, label in placements.items()) + ".")

    concentration = link_concentration(link_probabilities)
    permuted = np.concatenate([link_concentration(symptom_profiles_of_run(d)[1]) for d in arguments.permuted_run_dirs]) if arguments.permuted_run_dirs else np.array([])
    record["link_concentration"] = [None if np.isnan(value) else float(value) for value in concentration]
    record["link_concentration_permuted_null"] = {"runs": len(arguments.permuted_run_dirs),
                                                  "mean": float(np.nanmean(permuted)) if permuted.size else None,
                                                  "max": float(np.nanmax(permuted)) if permuted.size else None}
    lines += ["", "Link concentration per module (max share of one symptom; 1 means the module's links sit on one symptom, "
              f"{1 / max(len(symptoms), 1):.3f} means they are spread evenly): "
              + ", ".join("n/a" if value is None else f"{value:.3f}" for value in record["link_concentration"]) + ".",
              "", "No threshold is drawn on it: this analysis makes no decision, so the concentration is reported as the "
              "continuous quantity it is."]
    if permuted.size:
        lines.append(f"Permuted-label null over {len(arguments.permuted_run_dirs)} runs: mean "
                     f"{record['link_concentration_permuted_null']['mean']:.3f}, largest "
                     f"{record['link_concentration_permuted_null']['max']:.3f}.")
    table = "\n".join(lines)
    print(table)
    if arguments.json_output:
        arguments.json_output.write_text(json.dumps(record, indent=1) + "\n")
    if arguments.markdown_output:
        arguments.markdown_output.write_text(table + "\n")


if __name__ == "__main__":
    main()
