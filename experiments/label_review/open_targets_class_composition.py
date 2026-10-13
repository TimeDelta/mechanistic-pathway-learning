"""What the Open Targets target-level class would contribute as training labels, on a named evidence table.

Amendment item 12 fixes the class's loss weight at 0.45 and registers 0.0 as a declared secondary, on the reasoning
that the class is a small part of training. The counts behind that reasoning were measured on the table of the day
(docs/label_source_review.md 9.3 and 27.3), and the table is being rebuilt, so this script recomputes them rather
than repeating them. It reads tables and writes counts: no label, evidence table or configuration is changed.

The class is the human-worded curated triples (target, direction, symptom) that **no study drug label confirms** -
by construction the statements that would add something - restricted to targets the graph holds. A triple's direction
is its own perturbation: "inhibition of X" matches the study's knockout of X, while "activation of X" is a
perturbation the study does not otherwise have, so an activation statement adds a perturbation row rather than a
label on an existing row. The script reports the two apart, because only the second changes the perturbation set and
so the leakage groups the held-out set is drawn from.

It also prints each statement's curated wording beside the symptom it was mapped to. The mapping is the keyword
crosswalk of experiments/scope_open_targets_symptom_labels.py, which that script calls "a scoping estimate, not the
crosswalk a label build would use", so the wording column is what a reviewer reads to accept or reject each mapping
before the class is built.

Usage:
  python experiments/label_review/open_targets_class_composition.py \\
      --evidence-dir data/processed/evidence_full_v4 \\
      --label-selection data/processed/label_selection/better_v2_full_v4.parquet \\
      --graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v4
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.assign_evidence_grades import (
    DEFAULT_GRADE_WEIGHTS,
    DEFAULT_OPEN_TARGETS_WEIGHT,
    OPEN_TARGETS_GRADE,
)

WEIGHT_SCALE_STEP = 0.05  # amendment item 12: the derived weight is rounded down to the nearest step of the weight scale


def rounded_down_to_the_weight_step(value: float) -> float:
    """Round down to the nearest WEIGHT_SCALE_STEP, without the float floor-division error.

    `0.45 // 0.05` is 8.0 in exact arithmetic and 8.999999999999998 in binary floating point, so the naive
    expression reports a weight one step too low at every exact multiple: 0.45 printed as 0.40, 0.30 as 0.25,
    0.15 as 0.10. 0.45 is the registered weight, so the one agreement that would have confirmed the registration
    was the one the script would have contradicted."""
    steps = round(value / WEIGHT_SCALE_STEP, 9)
    import math

    return math.floor(steps) * WEIGHT_SCALE_STEP
CURATED_SYMPTOM_ROWS = Path("data/processed/off_target_scoping/open_targets_symptom_rows.parquet")
SIGN_OF_CURATED_DIRECTION = {"inhibition": -1.0, "activation": 1.0}
GENE_NODE_PREFIX = "GENE:"


def human_worded_triples(curated_rows: pd.DataFrame) -> pd.DataFrame:
    """One row per (target, direction, symptom), keeping only those whose every curated row is worded for people:
    a rodent readout or wording the rodent assay and the clinic share is not a human statement."""
    triples = curated_rows.groupby(["target", "direction", "target_symptom"]).agg(
        every_row_is_an_animal_readout=("animal_readout", "all"),
        every_row_is_species_ambiguous=("species_ambiguous_readout", "all"),
        curated_wording=("curated_symptom", lambda values: " | ".join(sorted(set(values)))),
        references=("reference", lambda values: "; ".join(sorted(set(str(value) for value in values))))).reset_index()
    return triples[~triples.every_row_is_an_animal_readout & ~triples.every_row_is_species_ambiguous]


def drugs_acting_on_each_target(drug_records: pd.DataFrame) -> dict[tuple[str, float], list[str]]:
    acting = {}
    for drug in drug_records.drop_duplicates("perturbation_id").itertuples():
        for node, sign, _occupancy in json.loads(drug.perturbation_nodes):
            acting.setdefault((node, float(sign)), []).append(drug.perturbation_id)
    return acting


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--label-selection", type=Path, required=True)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--curated-rows", type=Path, default=CURATED_SYMPTOM_ROWS)
    parser.add_argument("--markdown-output", type=Path)
    arguments = parser.parse_args()

    triples = human_worded_triples(pd.read_parquet(arguments.curated_rows))
    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet")
    selection = pd.read_parquet(arguments.label_selection)
    kept = selection[selection.keep]
    kept_symptoms_of = kept.groupby("perturbation_id").symptom.apply(set).to_dict()
    graph_genes = {node[len(GENE_NODE_PREFIX):] for node in pd.read_parquet(arguments.graph_dir / "nodes.parquet").node_id
                   if node.startswith(GENE_NODE_PREFIX)}
    gene_perturbations = set(records[records.perturbation_type != "drug"].perturbation_id)
    acting = drugs_acting_on_each_target(records[records.perturbation_type == "drug"])

    rows = []
    for triple in triples.itertuples():
        drugs_in_that_direction = acting.get((f"{GENE_NODE_PREFIX}{triple.target}", SIGN_OF_CURATED_DIRECTION[triple.direction]), [])
        confirmed_by = [drug for drug in drugs_in_that_direction if triple.target_symptom in kept_symptoms_of.get(drug, set())]
        rows.append({"target": triple.target, "direction": triple.direction, "symptom": triple.target_symptom,
                     "curated_wording": triple.curated_wording, "target_is_a_graph_gene": triple.target in graph_genes,
                     "a_study_drug_acts_in_that_direction": bool(drugs_in_that_direction),
                     "a_study_drug_label_confirms_it": bool(confirmed_by),
                     "the_target_is_already_a_perturbation": triple.target in gene_perturbations})
    statements = pd.DataFrame(rows)
    would_train = statements[~statements.a_study_drug_label_confirms_it & statements.target_is_a_graph_gene]
    adds_a_perturbation_row = would_train[would_train.direction == "activation"]
    kept_pairs_on_training_perturbations = int(len(kept))
    testable = statements[statements.a_study_drug_acts_in_that_direction]
    agreement = len(testable[testable.a_study_drug_label_confirms_it]) / len(testable) if len(testable) else float("nan")
    weight_the_agreement_gives = rounded_down_to_the_weight_step(DEFAULT_GRADE_WEIGHTS[OPEN_TARGETS_GRADE] * agreement)

    lines = [f"# Open Targets class composition on {arguments.evidence_dir.name} (generated by experiments/label_review/open_targets_class_composition.py)", "",
             f"| | statements |", "|---|---|",
             f"| human-worded curated triples | {len(statements)} |",
             f"| a study drug acts on the target in that direction | {int(statements.a_study_drug_acts_in_that_direction.sum())} |",
             f"| of those, a kept drug label confirms the symptom | {int(statements.a_study_drug_label_confirms_it.sum())} |",
             f"| **the class: no kept drug label confirms it, and the target is a graph gene** | **{len(would_train)}** |",
             f"| of the class, inhibition (the study's own knockout) | {int((would_train.direction == 'inhibition').sum())} |",
             f"| of the class, activation (a perturbation row the study does not have) | {len(adds_a_perturbation_row)} |",
             f"| of the class, the target is already a gene perturbation | {int(would_train.the_target_is_already_a_perturbation.sum())} |",
             f"| kept pairs in {arguments.label_selection.name} | {kept_pairs_on_training_perturbations} |",
             f"| the class as a share of those | {len(would_train) / kept_pairs_on_training_perturbations:.2%} |"
             if kept_pairs_on_training_perturbations else "| the class as a share of those | no kept pairs in the selection |",
             f"| **agreement of the class with the study's own labels** | **{agreement:.4f}** ({len(testable[testable.a_study_drug_label_confirms_it])} of {len(testable)}) |",
             f"| the weight that agreement gives (grade {OPEN_TARGETS_GRADE} x agreement, rounded down to {WEIGHT_SCALE_STEP}) | **{weight_the_agreement_gives:.2f}** |",
             f"| the weight registered in amendment item 12 | {DEFAULT_OPEN_TARGETS_WEIGHT:.2f} |", "",
             "The agreement is the number amendment item 12 derives its weight from, and this script recomputes it "
             "rather than quoting it. A disagreement with the registered 0.7606 is expected while the table carries "
             "the three-target limit of item 1: a drug's `perturbation_nodes` then names at most three targets, so a "
             "drug that acts on a statement's target through a fourth is neither testable nor confirming here. The "
             "number that settles the weight is this one, computed on the rebuilt table.", "",
             "Every statement of the class, with the curated wording the keyword crosswalk read. The wording column is "
             "what a reviewer checks before the class is built: the crosswalk is a keyword match and "
             "experiments/scope_open_targets_symptom_labels.py calls it a scoping estimate rather than the crosswalk a "
             "label build would use.", "",
             "| target | direction | mapped symptom | curated wording | already a perturbation |", "|---|---|---|---|---|"]
    for statement in would_train.sort_values(["target", "symptom"]).itertuples():
        lines.append(f"| {statement.target} | {statement.direction} | {statement.symptom} | {statement.curated_wording} | "
                     f"{'yes' if statement.the_target_is_already_a_perturbation else 'no'} |")
    table = "\n".join(lines)
    print(table)
    if arguments.markdown_output:
        arguments.markdown_output.write_text(table + "\n")


if __name__ == "__main__":
    main()
