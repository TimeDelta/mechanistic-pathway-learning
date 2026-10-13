"""Evidence class "open_targets": the target-level safety statements of the Open Targets curation, as labels on
perturbations the study already has (preregistration items 10, 12 and 12a).

A curated statement is (target, direction, symptom): "inhibition of X causes Y", "activation of X causes Y". Only
the human-worded ones count, so a rodent readout, or wording the rodent assay and the clinic share, is left out
upstream by experiments/scope_open_targets_symptom_labels.py, as are the wordings the owner ruled out of the
crosswalk on 13 October 2026.

**No Open Targets statement lands on a scored row.** That is the rule this module enforces, and it follows from
what the statements are (the owner, 13 October 2026, of the inhibition statements: "oh ... they're weak evidence",
and then "ok. drop them"). They are review statements about a target class, with no frequency and no recorded
absence, and the weight of 0.45 that prices them in training does nothing where it would matter most: AUPRC reads
a label as a label, so a statement in the lockbox would score exactly like a grade-A monogenic observation. Worse
than its one point, it could tip a symptom over the five-kept-positives line and so decide **which** symptoms the
macro average takes. The two directions resolve differently because the study already represents one of them.

- **Inhibition is dropped.** The study's knockout of that gene is a perturbation of the same node in the same
  direction, and it is scored. Putting the statement there would put weak evidence into the primary score; giving
  it a row of its own would duplicate a perturbation the table already holds. Neither is worth ten labels, so the
  direction carries nothing. `docs/label_source_review.md` 9.3 recommended as much when the class looked five
  times larger: "In the lockbox they would be scored against the weakest labels of the study".
- **Activation becomes its own perturbation row, never scored** (the owner: "what about just adding the receptor
  activations as a new perturbation type?", "what about limiting that perturbation type to the training data since
  it's rows are low and it doesn't have a full mechanism pathway?"). Here there is no duplication, because nothing
  else in the table perturbs these receptors upward: the activation direction is the one thing this curation
  uniquely carries. The rows are `RECEPTOR_ACTIVATION_PERTURBATION_TYPE`, seeded +1 at the gene node, and
  run_main_model.withhold_never_scored_rows keeps them out of every test set and out of training whenever their
  leakage group is held out.

A statement whose pair is already an **evidenced negative** - a never-above-placebo pair relabelled as a negative
under item 3 - does not become a positive. The two are a direct contradiction rather than independent evidence,
and the placebo comparison is the better-controlled of the two.

The output is evidence records plus the perturbation rows the activation statements need; nothing here is a node
or an edge, and assign_evidence_grades.loss_weight_for_record gives the records the registered weight of the class.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.assign_evidence_grades import OPEN_TARGETS_EVIDENCE_CLASS, EvidenceRecord
from mechanistic_pathway_learning.evidence.combined_perturbations import RECEPTOR_ACTIVATION_PERTURBATION_TYPE

CURATED_SYMPTOM_ROWS = Path("data/processed/off_target_scoping/open_targets_symptom_rows.parquet")
SOURCE_NAME = "Open Targets target safety"
GENE_NODE_PREFIX = "GENE:"
SIGN_OF_CURATED_DIRECTION = {"inhibition": -1.0, "activation": 1.0}
DROPPED_DIRECTION = "inhibition"  # the gene knockout already represents it, and it is scored
ROW_OF_ITS_OWN_DIRECTION = "activation"
ACTIVATION_PERTURBATION_PREFIX = "ACTIVATION:"


def activation_perturbation_id(target: str) -> str:
    """The never-scored row that carries agonism at one receptor, distinct from the gene's own knockout row."""
    return f"{ACTIVATION_PERTURBATION_PREFIX}{target}"


@dataclass(frozen=True)
class CuratedStatement:
    target: str
    direction: str
    symptom: str
    curated_wording: str
    references: str


@dataclass(frozen=True)
class StatementPlacement:
    """Where one statement landed, or why it did not land anywhere."""
    statement: CuratedStatement
    perturbation_identifier: str | None
    reason: str


def human_worded_statements(curated_rows: pd.DataFrame) -> list[CuratedStatement]:
    """One statement per (target, direction, symptom) whose every curated row is worded for people."""
    grouped = curated_rows.groupby(["target", "direction", "target_symptom"]).agg(
        every_row_is_an_animal_readout=("animal_readout", "all"),
        every_row_is_species_ambiguous=("species_ambiguous_readout", "all"),
        curated_wording=("curated_symptom", lambda values: " | ".join(sorted(set(values)))),
        references=("reference", lambda values: "; ".join(sorted(set(str(value) for value in values))))).reset_index()
    human = grouped[~grouped.every_row_is_an_animal_readout & ~grouped.every_row_is_species_ambiguous]
    return [CuratedStatement(row.target, row.direction, row.target_symptom, row.curated_wording, row.references)
            for row in human.itertuples()]


def drugs_whose_only_target_is(targets_of_drug: dict[str, set[str]], signs_of_drug: dict[str, dict[str, float]],
                              target: str, direction: str) -> list[str]:
    """Drugs that act on the target in the curated direction and have no other target. A drug with any second target
    is excluded whatever its affinities: the statement names the receptor, and a symptom observed for the drug could
    have come from the other target."""
    node = f"{GENE_NODE_PREFIX}{target}"
    return sorted(drug for drug, targets in targets_of_drug.items()
                  if targets == {node} and signs_of_drug[drug].get(node) == SIGN_OF_CURATED_DIRECTION[direction])


def place_statements(statements: list[CuratedStatement], graph_genes: set[str], gene_perturbations: set[str],
                     targets_of_drug: dict[str, set[str]], signs_of_drug: dict[str, dict[str, float]],
                     symptoms_a_kept_drug_label_confirms: dict[str, set[str]],
                     evidenced_negative_pairs: set[tuple[str, str]] | None = None) -> list[StatementPlacement]:
    """Where each statement goes. A statement a kept drug label already confirms adds nothing and is dropped before
    anything else, which is what makes this class the statements the study's own labels do not make."""
    evidenced_negative_pairs = evidenced_negative_pairs or set()
    placements = []
    for statement in statements:
        node = f"{GENE_NODE_PREFIX}{statement.target}"
        confirming = [drug for drug, signs in signs_of_drug.items()
                      if signs.get(node) == SIGN_OF_CURATED_DIRECTION[statement.direction]
                      and statement.symptom in symptoms_a_kept_drug_label_confirms.get(drug, set())]
        if confirming:
            placements.append(StatementPlacement(statement, None, "a kept drug label already confirms it"))
        elif statement.target not in graph_genes:
            placements.append(StatementPlacement(statement, None, "the target is not a graph gene"))
        elif statement.direction == DROPPED_DIRECTION:
            placements.append(StatementPlacement(statement, None, "inhibition: the gene knockout already represents it and is scored"))
        elif (activation_perturbation_id(statement.target), statement.symptom) in evidenced_negative_pairs:
            placements.append(StatementPlacement(statement, None, "the pair is an evidenced negative"))
        else:
            placements.append(StatementPlacement(statement, activation_perturbation_id(statement.target),
                                                 "a receptor activation row of its own, never scored"))
    return placements


def activation_perturbation_rows(placements: list[StatementPlacement]) -> list[dict]:
    """The perturbation rows the placed activation statements need: one per target, typed so that
    run_main_model.withhold_never_scored_rows keeps it out of every test set, seeded +1 at the gene node, and given
    the leakage group of that gene so a held-out group takes the activation row out of training with it."""
    rows = {}
    for placement in placements:
        if placement.perturbation_identifier is None:
            continue
        rows.setdefault(placement.perturbation_identifier, {
            "perturbation_id": placement.perturbation_identifier,
            "perturbation_type": RECEPTOR_ACTIVATION_PERTURBATION_TYPE,
            "perturbation_label": f"{placement.statement.target} activation",
            "perturbation_nodes": [[f"{GENE_NODE_PREFIX}{placement.statement.target}", 1.0, 1.0]],
            "groups_with_gene": placement.statement.target,
        })
    return sorted(rows.values(), key=lambda row: row["perturbation_id"])


def evidence_records(placements: list[StatementPlacement]) -> list[EvidenceRecord]:
    """One record per placed statement. A pair placed twice - two single-target drugs carrying one statement - gives
    one record each, because the perturbations differ; the same (perturbation, symptom) is never emitted twice."""
    seen, records = set(), []
    for placement in placements:
        if placement.perturbation_identifier is None:
            continue
        pair = (placement.perturbation_identifier, placement.statement.symptom)
        if pair in seen:
            continue
        seen.add(pair)
        records.append(EvidenceRecord(
            perturbation_identifier=placement.perturbation_identifier,
            symptom_identifier=placement.statement.symptom,
            relation="induces",
            evidence_class=OPEN_TARGETS_EVIDENCE_CLASS,
            source=f"{SOURCE_NAME}; {placement.statement.target} {placement.statement.direction}"
                   f" ({placement.statement.curated_wording})",
            distinct_reference_count=len({reference for reference in placement.statement.references.split(";") if reference.strip()}),
        ))
    return records
