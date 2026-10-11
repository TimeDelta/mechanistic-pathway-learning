"""How the limit on a drug's mechanism targets changes the dataset: drugs, kept pairs, leakage groups and input size.

Part of the label source review (docs/label_source_review.md, sections 4 and 12). Counting only: no label, table or
configuration of the study is changed.

Two things this script reports that `summarise_label_variants.py` does not, and that the review got wrong without them:

  the group that holds the most drugs   the largest leakage group by perturbations is not the group with the most drugs.
                                        Under the three-target rule the largest group (241 perturbations) holds 71
                                        drugs while another group (167 perturbations) holds 114.
  the total input of a drug             each mechanism target enters at magnitude 1, so the total follows the number of
                                        targets.

Leakage groups are those of the registered rule (`disease_cluster_and_targets`, which chains: perturbations that share a
node are joined) on the confirmatory graph.

Usage:
  PYTHONPATH=. python experiments/label_review/target_limit_comparison.py
  PYTHONPATH=. python experiments/label_review/target_limit_comparison.py --variants approved approved_no_cap --added-over approved approved_no_cap
"""
from __future__ import annotations

import argparse
import sys
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from mechanistic_pathway_learning.evaluation.experiment_data import ExperimentData, load_experiment_data  # noqa: E402

REVIEW_DIRECTORY = Path("data/processed/label_review")
CONFIRMATORY_GRAPH = Path("data/processed/graph_full_neuronal_split_binders")
DEFAULT_VARIANTS = ("n_only_cap1", "n_only_cap3", "n_only_cap0", "any_atc_cap1", "any_atc_cap3", "any_atc_cap0",
                    "approved", "approved_no_cap", "approved_with_list", "approved_with_list_no_cap")


def load_variant(variant_name: str) -> tuple[ExperimentData, pd.DataFrame, pd.DataFrame]:
    """The loaded data, the evidence records and the kept drug pairs of one counting build."""
    evidence_directory = REVIEW_DIRECTORY / f"evidence_{variant_name}"
    selection_path = REVIEW_DIRECTORY / f"selection_{variant_name}.parquet"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        data = load_experiment_data(CONFIRMATORY_GRAPH, evidence_directory, label_selection=selection_path, group_by="disease_cluster_and_targets")
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet")
    selection = pd.read_parquet(selection_path)
    kept_drug_pairs = selection[(selection.perturbation_type == "drug") & selection.keep.astype(bool)]
    return data, records, kept_drug_pairs


def total_seed_magnitude_by_drug(data: ExperimentData) -> dict[str, float]:
    """Sum of the absolute seed magnitudes of each drug: 1 per mechanism target."""
    is_drug = np.array(data.perturbation_types) == "drug"
    return {data.perturbation_ids[position]: float(np.abs(np.asarray(data.perturbation_magnitudes[position], dtype=float)).sum()) for position in np.flatnonzero(is_drug)}


def summarise_variant(variant_name: str) -> dict:
    data, records, kept_drug_pairs = load_variant(variant_name)
    perturbation_ids = np.array(data.perturbation_ids)
    group_ids = np.array(data.group_ids)
    is_drug = np.array(data.perturbation_types) == "drug"
    drug_count = int(is_drug.sum())
    kept_pairs_by_drug = Counter(kept_drug_pairs.perturbation_id)
    drugs_by_group = Counter(group_ids[is_drug])
    groups_by_drug_count = drugs_by_group.most_common()
    row = {"variant": variant_name, "perturbations": len(perturbation_ids), "drugs": drug_count, "kept_drug_pairs": len(kept_drug_pairs),
           "drugs_with_a_kept_pair": int(kept_drug_pairs.perturbation_id.nunique()), "groups_holding_a_drug": len(drugs_by_group)}
    for rank, name in ((0, "group_with_most_drugs"), (1, "group_with_second_most_drugs")):
        if rank >= len(groups_by_drug_count):
            continue
        group_id, drugs_in_group = groups_by_drug_count[rank]
        member_drug_ids = perturbation_ids[is_drug & (group_ids == group_id)]
        kept_pairs_in_group = sum(kept_pairs_by_drug[drug_id] for drug_id in member_drug_ids)
        row |= {f"{name}_drugs": drugs_in_group, f"{name}_share_of_drugs": round(drugs_in_group / drug_count, 3),
                f"{name}_share_of_kept_drug_pairs": round(kept_pairs_in_group / max(len(kept_drug_pairs), 1), 3),
                f"{name}_perturbations": int((group_ids == group_id).sum())}
    largest_group_id, largest_group_perturbations = Counter(group_ids).most_common(1)[0]
    row |= {"largest_group_perturbations": largest_group_perturbations, "largest_group_drugs": int((group_ids[is_drug] == largest_group_id).sum())}
    totals = total_seed_magnitude_by_drug(data)
    total_distribution = Counter(int(round(total)) for total in totals.values())
    row |= {"total_seed_magnitude_distribution": dict(sorted(total_distribution.items())), "largest_total_seed_magnitude": max(totals.values()),
            "drugs_with_seeds_of_both_signs": int(sum(len(set(np.sign(np.asarray(data.perturbation_signs[position], dtype=float)).tolist())) > 1 for position in np.flatnonzero(is_drug)))}
    return row


def drugs_added(base_variant: str, wider_variant: str) -> pd.DataFrame:
    """The drugs of the wider build that the base build lacks, with their total input and kept pairs."""
    _, base_records, _ = load_variant(base_variant)
    wider_data, wider_records, wider_kept_drug_pairs = load_variant(wider_variant)
    base_drug_ids = set(base_records[base_records.perturbation_type == "drug"].perturbation_id)
    wider_drug_records = wider_records[wider_records.perturbation_type == "drug"]
    label_by_id = dict(zip(wider_drug_records.perturbation_id, wider_drug_records.perturbation_label))
    totals = total_seed_magnitude_by_drug(wider_data)
    kept_pairs_by_drug = Counter(wider_kept_drug_pairs.perturbation_id)
    rows = [{"drug": label_by_id[drug_id], "mechanism_targets": int(round(totals.get(drug_id, 0.0))), "kept_pairs": kept_pairs_by_drug.get(drug_id, 0),
             "kept_symptoms": ", ".join(sorted(wider_kept_drug_pairs[wider_kept_drug_pairs.perturbation_id == drug_id].symptom.unique()))}
            for drug_id in sorted(set(label_by_id) - base_drug_ids, key=lambda drug_id: -kept_pairs_by_drug.get(drug_id, 0))]
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--variants", nargs="+", default=list(DEFAULT_VARIANTS), help="names of the counting builds under data/processed/label_review")
    parser.add_argument("--added-over", nargs=2, metavar=("BASE", "WIDER"), default=["approved", "approved_no_cap"], help="list the drugs the wider build adds to the base build")
    parser.add_argument("--output", type=Path, default=REVIEW_DIRECTORY / "target_limit_comparison.csv")
    arguments = parser.parse_args()
    rows = []
    for variant_name in arguments.variants:
        if not (REVIEW_DIRECTORY / f"selection_{variant_name}.parquet").exists():
            print(f"{variant_name}: not built, skipped")
            continue
        rows.append(summarise_variant(variant_name))
    summary = pd.DataFrame(rows)
    summary.to_csv(arguments.output, index=False)
    pd.set_option("display.width", 260, "display.max_columns", 40, "display.max_colwidth", 70)
    print(summary.to_string(index=False))
    base_variant, wider_variant = arguments.added_over
    if all((REVIEW_DIRECTORY / f"selection_{name}.parquet").exists() for name in (base_variant, wider_variant)):
        added = drugs_added(base_variant, wider_variant)
        added.to_csv(arguments.output.with_name(f"drugs_added_{wider_variant}_over_{base_variant}.csv"), index=False)
        print(f"\ndrugs in {wider_variant} and not in {base_variant}: {len(added)}, with {int(added.kept_pairs.sum())} kept pairs")
        print(added.to_string(index=False))
    print(f"\nwrote {arguments.output}")


if __name__ == "__main__":
    main()
