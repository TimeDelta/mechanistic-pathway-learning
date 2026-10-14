"""How much the Open Targets class changes the Dawid-Skene fit, and whether option (a)'s missing guard matters.

Decision 8.5 and amendment item 12d (the owner, 13 October 2026: "add it to the Dawid-Skene fit estimates"). The
class reaches the reliability fit and not the grade, so the whole of its effect is on the fitted sensitivities and
on the posterior each item carries. This measures that effect before the wiring is chosen.

The second question is the one that decides between the owner's options (a) and (b). `reliability_reports` skips a
pair the selection has relabelled an evidenced negative, but that set lives downstream of `assemble()`. Option (a)
wires it with the guard inactive. What that costs is exactly the number of Open Targets reports that land on a
pair the selection flags `masks_a_negative`, so that number is computed here rather than assumed small.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evidence.evidence_reliability_model import fit_weighted_dawid_skene
from mechanistic_pathway_learning.evidence.load_open_targets_statements import (
    human_worded_statements,
    reliability_reports,
)

ADVERSE_EFFECTS_PATH = Path("data/raw/opentargets_curation/adverse_effects.tsv")


def build_statements(scoped_rows_path: Path) -> list:
    """The scoping pass already wrote one row per (target symptom, target, curated row); reuse it."""
    rows = pd.read_parquet(scoped_rows_path)
    return human_worded_statements(rows)


def report_matrices(reports: pd.DataFrame, item_keys: list[tuple[str, str, str]],
                    source_names: list[str], implicit_negative_weight: float) -> tuple[np.ndarray, np.ndarray]:
    index_of_item = {key: index for index, key in enumerate(item_keys)}
    index_of_source = {name: index for index, name in enumerate(source_names)}
    positive = np.zeros((len(item_keys), len(source_names)))
    negative = np.zeros((len(item_keys), len(source_names)))
    for row in reports.itertuples():
        item_index = index_of_item.get((row.perturbation_id, row.symptom, row.relation))
        source_index = index_of_source.get(row.source)
        if item_index is None or source_index is None:
            continue
        weight = float(getattr(row, "rubric_weight", 1.0) or 1.0)
        if row.report_value == 1:
            positive[item_index, source_index] += weight
        else:
            negative[item_index, source_index] += weight
    covered = (positive + negative) > 0
    negative += implicit_negative_weight * (~covered)
    return positive, negative


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--scoped-rows", type=Path,
                        default=Path("data/processed/off_target_scoping/open_targets_symptom_rows.parquet"))
    parser.add_argument("--selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
    parser.add_argument("--implicit-negative-weight", type=float, default=0.1)
    parser.add_argument("--output", type=Path, default=Path("docs/open_targets_reliability_effect.json"))
    arguments = parser.parse_args()

    existing_reports = pd.read_parquet(arguments.evidence_dir / "evidence_reports.parquet")
    observations = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet").drop_duplicates(
        ["perturbation_id", "symptom", "relation"])
    item_keys = list(zip(observations.perturbation_id, observations.symptom, observations.relation))
    item_groups = [str(value) for value in observations.evidence_class]

    items_with_a_positive_report = {
        (row.perturbation_id, row.symptom, row.relation)
        for row in existing_reports[existing_reports.report_value == 1].itertuples()}
    gene_perturbation_ids = set(existing_reports[existing_reports.perturbation_type == "gene"].perturbation_id)
    perturbation_of_target = {gene: gene for gene in gene_perturbation_ids}

    statements = build_statements(arguments.scoped_rows)
    open_targets_rows = reliability_reports(statements, items_with_a_positive_report, perturbation_of_target)
    open_targets_frame = pd.DataFrame(open_targets_rows)

    effect: dict = {"curated_statements_read": len(statements),
                    "open_targets_reports_emitted": len(open_targets_frame)}
    if len(open_targets_frame):
        touched_items = set(zip(open_targets_frame.perturbation_id, open_targets_frame.symptom,
                                open_targets_frame.relation))
        effect["distinct_items_touched"] = len(touched_items)
        effect["distinct_targets"] = int(open_targets_frame.perturbation_id.nunique())
        effect["distinct_symptoms"] = int(open_targets_frame.symptom.nunique())
        effect["report_values"] = {str(value): int(count) for value, count
                                   in open_targets_frame.report_value.value_counts().items()}
        selection = pd.read_parquet(arguments.selection)
        masked_pairs = {(row.perturbation_id, row.symptom) for row in selection.itertuples() if row.masks_a_negative}
        landing_on_a_masked_negative = {item for item in touched_items if (item[0], item[1]) in masked_pairs}
        effect["pairs_the_inactive_guard_would_not_skip"] = len(landing_on_a_masked_negative)
        effect["examples_of_those_pairs"] = sorted(f"{p}|{s}|{r}" for p, s, r in landing_on_a_masked_negative)[:10]

        source_names_before = sorted(set(existing_reports.source))
        positive_before, negative_before = report_matrices(existing_reports, item_keys, source_names_before,
                                                           arguments.implicit_negative_weight)
        fit_before = fit_weighted_dawid_skene(positive_before, negative_before, source_names_before, item_keys,
                                              item_groups=item_groups)

        open_targets_frame["rubric_weight"] = 1.0
        combined = pd.concat([existing_reports, open_targets_frame], ignore_index=True)
        source_names_after = sorted(set(combined.source))
        positive_after, negative_after = report_matrices(combined, item_keys, source_names_after,
                                                         arguments.implicit_negative_weight)
        fit_after = fit_weighted_dawid_skene(positive_after, negative_after, source_names_after, item_keys,
                                            item_groups=item_groups)

        sensitivity_before = dict(zip(fit_before.source_names, (float(v) for v in fit_before.sensitivity)))
        sensitivity_after = dict(zip(fit_after.source_names, (float(v) for v in fit_after.sensitivity)))
        effect["sensitivity_before"] = sensitivity_before
        effect["sensitivity_after"] = sensitivity_after
        effect["sensitivity_change"] = {name: round(sensitivity_after[name] - sensitivity_before[name], 5)
                                        for name in sensitivity_before}
        posterior_shift = np.abs(np.asarray(fit_after.posterior) - np.asarray(fit_before.posterior))
        effect["posterior_shift"] = {"items": int(posterior_shift.size),
                                     "mean_absolute": float(posterior_shift.mean()),
                                     "max_absolute": float(posterior_shift.max()),
                                     "items_shifted_more_than_0.01": int((posterior_shift > 0.01).sum()),
                                     "items_shifted_more_than_0.05": int((posterior_shift > 0.05).sum())}
        effect["new_source_name"] = sorted(set(source_names_after) - set(source_names_before))

    arguments.output.write_text(json.dumps(effect, indent=1))
    print(json.dumps(effect, indent=1))


if __name__ == "__main__":
    main()
