"""Evidence class "open_targets": the target-level safety statements of the Open Targets curation, as labels on
perturbations the study already has (preregistration items 10, 12 and 12a).

A curated statement is (target, direction, symptom): "inhibition of X causes Y", "activation of X causes Y". Only
the human-worded ones count, so a rodent readout, or wording the rodent assay and the clinic share, is left out
upstream by experiments/scope_open_targets_symptom_labels.py, as are the wordings the owner ruled out of the
crosswalk on 13 October 2026.

**The class adds no perturbation row.** That is the rule this module exists to enforce, and the reason is the
lockbox: a new perturbation changes the leakage groups, so it changes the held-out draw, so a class that invented
perturbations could not be added once the lockbox was drawn. Each direction is therefore carried by a perturbation
the table already holds, or it is not carried at all.

- **Inhibition** is carried by the study's own knockout of that gene, which is a perturbation of the same node in
  the same direction. What this conflates is registered and not new: a knockout is constitutive and developmental
  where an inhibitor is acute, dosed and reversible, which is the knockout-versus-inhibitor question the
  preregistration already carries as open.
- **Activation** is carried only by a drug that agonises the target **and has no other target at all** (the owner,
  13 October 2026: "make sure to not duplicate or attach to a drug that might not actually cause that symptom due
  to e.g. some other target or somethign (so it would have to at least be limited to drugs that currently have ONLY
  that receptor as a target)"). A drug with several targets cannot carry the statement, because the symptom it is
  being given might come from any of them. On the table of 13 October this admits **nothing**: all 23 activation
  targets are receptors, all 20 carrier pairs run through drugs with 2 to 13 targets, so every activation statement
  is dropped. The rule is implemented rather than the count, because the rebuild changes which drugs are admitted.

A statement whose pair is already an **evidenced negative** - a never-above-placebo pair relabelled as a negative
under item 3 - does not become a positive. The two are a direct contradiction rather than independent evidence,
and the placebo comparison is the better-controlled of the two.

Nothing here is a perturbation, a node or an edge: the output is evidence records, and
assign_evidence_grades.loss_weight_for_record gives them the registered weight of the class.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.assign_evidence_grades import OPEN_TARGETS_EVIDENCE_CLASS, EvidenceRecord

CURATED_SYMPTOM_ROWS = Path("data/processed/off_target_scoping/open_targets_symptom_rows.parquet")
SOURCE_NAME = "Open Targets target safety"
GENE_NODE_PREFIX = "GENE:"
SIGN_OF_CURATED_DIRECTION = {"inhibition": -1.0, "activation": 1.0}
CARRIED_BY_THE_GENE_KNOCKOUT = "inhibition"


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
        elif statement.direction == CARRIED_BY_THE_GENE_KNOCKOUT:
            if statement.target not in gene_perturbations:
                placements.append(StatementPlacement(statement, None, "the target is a graph gene but not a perturbation"))
            elif (statement.target, statement.symptom) in evidenced_negative_pairs:
                placements.append(StatementPlacement(statement, None, "the pair is an evidenced negative"))
            else:
                placements.append(StatementPlacement(statement, statement.target, "carried by the gene knockout"))
        else:
            carriers = drugs_whose_only_target_is(targets_of_drug, signs_of_drug, statement.target, statement.direction)
            carriers = [drug for drug in carriers if (drug, statement.symptom) not in evidenced_negative_pairs]
            if not carriers:
                placements.append(StatementPlacement(statement, None, "no drug whose only target is this receptor"))
            else:
                placements.extend(StatementPlacement(statement, drug, "carried by a single-target agonist") for drug in carriers)
    return placements


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
