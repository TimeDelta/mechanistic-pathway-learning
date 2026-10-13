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
STABILITY_SUBSAMPLE_FRACTION = 0.5
STABILITY_DRAWS = 100


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


def cluster_count_by_stability(symptom_profiles: np.ndarray, seed: int = 0) -> tuple[int, dict[int, float]]:
    """The number of clusters whose solution is most stable under subsampling the profile's rows (Ben-Hur, Elisseeff
    and Guyon 2002: stability "provides a means of rationally defining an optimum number of clusters, and can also
    detect the lack of structure in data"). The rows are the mechanism coordinates, so a subsample asks whether the
    grouping rests on a few of them."""
    generator = np.random.default_rng(seed)
    num_coordinates = symptom_profiles.shape[0]
    take = max(2, int(round(num_coordinates * STABILITY_SUBSAMPLE_FRACTION)))
    agreement_at = {}
    for num_clusters in CANDIDATE_CLUSTER_COUNTS:
        if num_clusters >= symptom_profiles.shape[1]:
            continue
        agreements = []
        for _ in range(STABILITY_DRAWS):
            first = generator.choice(num_coordinates, take, replace=False)
            second = generator.choice(num_coordinates, take, replace=False)
            agreements.append(partition_agreement(cluster_at(correlation_distances(symptom_profiles[first]), num_clusters),
                                                  cluster_at(correlation_distances(symptom_profiles[second]), num_clusters)))
        agreement_at[num_clusters] = float(np.mean(agreements))
    best = max(agreement_at, key=agreement_at.get) if agreement_at else 1
    return best, agreement_at


def link_concentration(link_probabilities: np.ndarray) -> np.ndarray:
    """max_s link_ks / sum_s link_ks per module: 1 when a module's links are on one symptom, 1/num_symptoms when they
    are spread evenly. A module with no link mass is reported as NaN rather than as concentrated."""
    totals = link_probabilities.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(totals > 0, link_probabilities.max(axis=1) / totals, np.nan)


def symptom_profiles_of_run(run_directory: Path) -> tuple[np.ndarray, np.ndarray]:
    """(profiles [mechanism coordinate, symptom], link probabilities [module, symptom] summed over relations). The
    profile stacks the relations, so a symptom is described by what induces it and what relieves it together."""
    results = json.loads((run_directory / "results.json").read_text())
    links = np.asarray(results["module_symptom_links"], dtype=float)
    if links.ndim == 2:
        links = links[None, :, :]
    return np.concatenate(list(links), axis=0), links.sum(axis=0)


def label_co_occurrence_profiles(evidence_directory: Path, label_selection: Path | None, symptoms: list[str]) -> np.ndarray:
    """The symptom columns of the raw kept outcome matrix, as profiles over perturbations: the baseline of 26e."""
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet", columns=["perturbation_id", "symptom", "relation"])
    records = records[records.relation == "induces"]
    if label_selection is not None:
        selection = pd.read_parquet(label_selection)
        kept_pairs = set(zip(selection[selection.keep].perturbation_id, selection[selection.keep].symptom))
        records = records[[pair in kept_pairs for pair in zip(records.perturbation_id, records.symptom)]]
    perturbations = sorted(set(records.perturbation_id))
    row_of = {perturbation: row for row, perturbation in enumerate(perturbations)}
    column_of = {symptom: column for column, symptom in enumerate(symptoms)}
    outcomes = np.zeros((len(perturbations), len(symptoms)))
    for record in records.itertuples():
        if record.symptom in column_of:
            outcomes[row_of[record.perturbation_id], column_of[record.symptom]] = 1.0
    return outcomes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dirs", type=Path, nargs="+", required=True,
                        help="fitted runs; the first is the full-data fit the readout describes, the rest are the perturbation subsamples of 26 d2")
    parser.add_argument("--permuted-run-dirs", type=Path, nargs="*", default=[], help="permuted-label runs: the null for the link concentration")
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--label-selection", type=Path)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    arguments = parser.parse_args()

    full_fit = json.loads((arguments.run_dirs[0] / "results.json").read_text())
    symptoms = list(full_fit["symptoms"]) if "symptoms" in full_fit else [str(i) for i in range(len(full_fit["symptom_leaks"][0]))]
    profiles, link_probabilities = symptom_profiles_of_run(arguments.run_dirs[0])

    num_clusters, stability_by_count = cluster_count_by_stability(profiles, arguments.seed)
    mechanism_partition = cluster_at(correlation_distances(profiles), num_clusters)

    baseline_profiles = label_co_occurrence_profiles(arguments.evidence_dir, arguments.label_selection, symptoms)
    baseline_count, _ = cluster_count_by_stability(baseline_profiles, arguments.seed)
    baseline_partition = cluster_at(correlation_distances(baseline_profiles), baseline_count)

    across_subsamples = [partition_agreement(mechanism_partition,
                                             cluster_at(correlation_distances(symptom_profiles_of_run(directory)[0]), num_clusters))
                         for directory in arguments.run_dirs[1:]]
    concentration = link_concentration(link_probabilities)
    permuted_concentration = np.concatenate([link_concentration(symptom_profiles_of_run(d)[1]) for d in arguments.permuted_run_dirs]) \
        if arguments.permuted_run_dirs else np.array([])

    record = {
        "reading": "descriptive; no held-out score, in no hypothesis of the decision rule (preregistration item 26)",
        "num_clusters": int(num_clusters), "stability_by_cluster_count": stability_by_count,
        "clusters": {str(label): [symptoms[i] for i in np.flatnonzero(mechanism_partition == label)] for label in np.unique(mechanism_partition)},
        "label_co_occurrence_clusters": {str(label): [symptoms[i] for i in np.flatnonzero(baseline_partition == label)] for label in np.unique(baseline_partition)},
        "agreement_with_the_label_co_occurrence_baseline": partition_agreement(mechanism_partition, baseline_partition),
        "agreement_across_perturbation_subsamples": across_subsamples,
        "link_concentration": [None if np.isnan(value) else float(value) for value in concentration],
        "link_concentration_permuted_null": {"runs": len(arguments.permuted_run_dirs),
                                             "mean": float(np.nanmean(permuted_concentration)) if permuted_concentration.size else None,
                                             "max": float(np.nanmax(permuted_concentration)) if permuted_concentration.size else None},
    }

    lines = ["# Nosology readout (generated by experiments/nosology_readout.py)", "",
             "**Descriptive.** Fitted on all the evidence, so no reading here is a held-out score and none enters H1, H2 or "
             "P1 (preregistration item 26b). It carries no claim that the grouping predicts anything, none about a symptom "
             "outside the set and none that it would recur in another cohort.", "",
             f"Clusters: **{num_clusters}**, taken by the stability criterion over candidate counts "
             + ", ".join(f"{count} ({agreement:.2f})" for count, agreement in stability_by_count.items()) + ".", ""]
    for label, members in record["clusters"].items():
        lines.append(f"- **cluster {label}**: {', '.join(members)}")
    lines += ["", f"Agreement with the label co-occurrence baseline: **{record['agreement_with_the_label_co_occurrence_baseline']:.3f}** "
              f"(the baseline takes {baseline_count} clusters). A cluster the baseline also returns is a grouping the labels "
              "already carry, and is not reported as mechanism-derived; what the graph adds is the difference between the "
              "two partitions, and that difference is this readout's only positive claim.", ""]
    if across_subsamples:
        lines.append(f"Agreement across {len(across_subsamples)} perturbation subsamples: "
                     f"{np.mean(across_subsamples):.3f} (lowest {min(across_subsamples):.3f}).")
    else:
        lines.append("Agreement across perturbation subsamples: not computed, no subsample runs were given.")
    lines += ["", "Link concentration per module (max share of one symptom; 1 means the module's links sit on one symptom, "
              f"{1 / max(len(symptoms), 1):.3f} means they are spread evenly): "
              + ", ".join("n/a" if value is None else f"{value:.3f}" for value in record["link_concentration"]) + ".",
              "", "No threshold is drawn on it: this analysis makes no decision, so the concentration is reported as the "
              "continuous quantity it is."]
    if permuted_concentration.size:
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
