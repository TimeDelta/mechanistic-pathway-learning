"""Preregistration item 12a: the Open Targets class adds no perturbation row, and an activation statement is carried
only by a drug whose sole target is that receptor."""
from __future__ import annotations

import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.assign_evidence_grades import (
    DEFAULT_OPEN_TARGETS_WEIGHT,
    OPEN_TARGETS_EVIDENCE_CLASS,
    loss_weight_for_record,
)
from mechanistic_pathway_learning.evidence.combined_perturbations import RECEPTOR_ACTIVATION_PERTURBATION_TYPE
from mechanistic_pathway_learning.evidence.load_open_targets_statements import (
    CuratedStatement,
    activation_perturbation_id,
    activation_perturbation_rows,
    evidence_records,
    human_worded_statements,
    place_statements,
)

GRAPH_GENES = {"GRIN1", "HTR2C", "CHRM1"}
GENE_PERTURBATIONS = {"GRIN1", "HTR2C", "CHRM1"}
TARGETS_OF_DRUG = {"single_target_agonist": {"GENE:HTR2C"},
                   "two_target_agonist": {"GENE:HTR2C", "GENE:CHRM1"},
                   "antagonist": {"GENE:HTR2C"}}
SIGNS_OF_DRUG = {"single_target_agonist": {"GENE:HTR2C": 1.0},
                 "two_target_agonist": {"GENE:HTR2C": 1.0, "GENE:CHRM1": 1.0},
                 "antagonist": {"GENE:HTR2C": -1.0}}


def statement(target: str, direction: str, symptom: str) -> CuratedStatement:
    return CuratedStatement(target, direction, symptom, "some wording", "PMID:1; PMID:2")


def placements_for(statements, confirms=None, evidenced_negatives=None):
    return place_statements(statements, GRAPH_GENES, GENE_PERTURBATIONS, TARGETS_OF_DRUG, SIGNS_OF_DRUG,
                            confirms or {}, evidenced_negatives or set())


def test_an_inhibition_statement_is_dropped_because_the_knockout_already_represents_it_and_is_scored():
    placed = placements_for([statement("GRIN1", "inhibition", "psychosis")])
    assert [p.perturbation_identifier for p in placed] == [None]
    assert "already represents it and is scored" in placed[0].reason


def test_an_activation_statement_gets_a_row_of_its_own_rather_than_a_drug_or_a_gene():
    placed = placements_for([statement("HTR2C", "activation", "anxiety")])
    assert [p.perturbation_identifier for p in placed] == [activation_perturbation_id("HTR2C")]
    assert placed[0].perturbation_identifier not in GENE_PERTURBATIONS and placed[0].perturbation_identifier not in TARGETS_OF_DRUG


def test_the_activation_row_is_of_the_never_scored_type_and_seeds_the_gene_node_upward():
    rows = activation_perturbation_rows(placements_for([statement("HTR2C", "activation", "anxiety"),
                                                        statement("HTR2C", "activation", "insomnia")]))
    assert len(rows) == 1  # one row per target, however many symptoms it carries
    assert rows[0]["perturbation_type"] == RECEPTOR_ACTIVATION_PERTURBATION_TYPE
    assert rows[0]["perturbation_nodes"] == [["GENE:HTR2C", 1.0, 1.0]]
    assert rows[0]["groups_with_gene"] == "HTR2C"  # so a held-out group takes the activation row with it


def test_a_statement_a_kept_drug_label_already_confirms_adds_nothing():
    placed = placements_for([statement("HTR2C", "activation", "anxiety")],
                            confirms={"two_target_agonist": {"anxiety"}, "single_target_agonist": {"anxiety"}})
    assert [p.reason for p in placed] == ["a kept drug label already confirms it"]


def test_a_statement_never_overrides_an_evidenced_negative():
    """A never-above-placebo pair and a review statement are a contradiction, not independent evidence."""
    placed = placements_for([statement("HTR2C", "activation", "anxiety")],
                            evidenced_negatives={(activation_perturbation_id("HTR2C"), "anxiety")})
    assert [p.perturbation_identifier for p in placed] == [None] and placed[0].reason == "the pair is an evidenced negative"


def test_a_target_the_graph_does_not_hold_carries_nothing():
    placed = placements_for([statement("NOT_IN_THE_GRAPH", "activation", "psychosis")])
    assert placed[0].perturbation_identifier is None and "not a graph gene" in placed[0].reason


def test_no_statement_ever_lands_on_a_scored_row():
    """The rule the module exists for, after the owner dropped the inhibitions as weak evidence: an Open Targets
    label must never reach a perturbation the lockbox scores."""
    every_direction = [statement(target, direction, "psychosis")
                       for target in GRAPH_GENES for direction in ("inhibition", "activation")]
    scored_rows = GENE_PERTURBATIONS | set(TARGETS_OF_DRUG)
    assert all(p.perturbation_identifier is None or p.perturbation_identifier not in scored_rows
               for p in placements_for(every_direction))


def test_the_records_carry_the_class_and_its_registered_weight_and_never_repeat_a_pair():
    placed = placements_for([statement("HTR2C", "activation", "anxiety"), statement("HTR2C", "activation", "anxiety")])
    records = evidence_records(placed)
    assert len(records) == 1
    assert records[0].evidence_class == OPEN_TARGETS_EVIDENCE_CLASS and records[0].relation == "induces"
    assert records[0].perturbation_identifier == activation_perturbation_id("HTR2C")
    assert records[0].distinct_reference_count == 2
    assert loss_weight_for_record(records[0]) == pytest.approx(DEFAULT_OPEN_TARGETS_WEIGHT)


def test_only_human_worded_statements_are_read():
    rows = pd.DataFrame({"target": ["A", "A", "B", "C"], "direction": ["inhibition"] * 4,
                         "target_symptom": ["psychosis", "psychosis", "anxiety", "fatigue"],
                         "animal_readout": [False, True, True, False],
                         "species_ambiguous_readout": [False, False, False, True],
                         "curated_symptom": ["w", "w", "w", "w"], "reference": ["PMID:1"] * 4})
    assert [s.target for s in human_worded_statements(rows)] == ["A"]  # B is rodent-only, C is species-ambiguous
