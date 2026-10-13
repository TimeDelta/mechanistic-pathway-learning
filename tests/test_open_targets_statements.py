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
from mechanistic_pathway_learning.evidence.load_open_targets_statements import (
    CuratedStatement,
    drugs_whose_only_target_is,
    evidence_records,
    human_worded_statements,
    place_statements,
)

GRAPH_GENES = {"GRIN1", "HTR2C", "CHRM1", "NOT_A_PERTURBATION"}
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


def test_an_inhibition_statement_is_carried_by_the_gene_knockout():
    placed = placements_for([statement("GRIN1", "inhibition", "psychosis")])
    assert [(p.perturbation_identifier, p.reason) for p in placed] == [("GRIN1", "carried by the gene knockout")]


def test_an_activation_statement_needs_a_drug_whose_only_target_is_that_receptor():
    only_the_two_target_drug_agonises_chrm1 = placements_for([statement("CHRM1", "activation", "irritability_or_aggression")])
    assert [p.perturbation_identifier for p in only_the_two_target_drug_agonises_chrm1] == [None]
    assert only_the_two_target_drug_agonises_chrm1[0].reason == "no drug whose only target is this receptor"

    carried = placements_for([statement("HTR2C", "activation", "anxiety")])
    assert [p.perturbation_identifier for p in carried] == ["single_target_agonist"]


def test_the_direction_of_the_drug_matters_as_well_as_the_target_count():
    assert drugs_whose_only_target_is(TARGETS_OF_DRUG, SIGNS_OF_DRUG, "HTR2C", "activation") == ["single_target_agonist"]
    assert drugs_whose_only_target_is(TARGETS_OF_DRUG, SIGNS_OF_DRUG, "HTR2C", "inhibition") == ["antagonist"]


def test_a_statement_a_kept_drug_label_already_confirms_adds_nothing():
    placed = placements_for([statement("HTR2C", "activation", "anxiety")],
                            confirms={"two_target_agonist": {"anxiety"}, "single_target_agonist": {"anxiety"}})
    assert [p.reason for p in placed] == ["a kept drug label already confirms it"]


def test_a_statement_never_overrides_an_evidenced_negative():
    """A never-above-placebo pair and a review statement are a contradiction, not independent evidence."""
    gene_side = placements_for([statement("GRIN1", "inhibition", "psychosis")], evidenced_negatives={("GRIN1", "psychosis")})
    assert [p.perturbation_identifier for p in gene_side] == [None] and gene_side[0].reason == "the pair is an evidenced negative"
    drug_side = placements_for([statement("HTR2C", "activation", "anxiety")],
                               evidenced_negatives={("single_target_agonist", "anxiety")})
    assert [p.perturbation_identifier for p in drug_side] == [None]


def test_a_target_that_is_a_graph_gene_but_not_a_perturbation_carries_nothing():
    placed = placements_for([statement("NOT_A_PERTURBATION", "inhibition", "psychosis")])
    assert placed[0].perturbation_identifier is None and "not a perturbation" in placed[0].reason


def test_no_placement_ever_names_a_perturbation_the_table_does_not_hold():
    """The rule the module exists for: the class adds no perturbation row, so the leakage groups do not move."""
    every_direction = [statement(target, direction, "psychosis")
                       for target in GRAPH_GENES for direction in ("inhibition", "activation")]
    known = GENE_PERTURBATIONS | set(TARGETS_OF_DRUG)
    assert all(p.perturbation_identifier is None or p.perturbation_identifier in known for p in placements_for(every_direction))


def test_the_records_carry_the_class_and_its_registered_weight_and_never_repeat_a_pair():
    placed = placements_for([statement("GRIN1", "inhibition", "psychosis"), statement("GRIN1", "inhibition", "psychosis")])
    records = evidence_records(placed)
    assert len(records) == 1
    assert records[0].evidence_class == OPEN_TARGETS_EVIDENCE_CLASS and records[0].relation == "induces"
    assert records[0].distinct_reference_count == 2
    assert loss_weight_for_record(records[0]) == pytest.approx(DEFAULT_OPEN_TARGETS_WEIGHT)


def test_only_human_worded_statements_are_read():
    rows = pd.DataFrame({"target": ["A", "A", "B", "C"], "direction": ["inhibition"] * 4,
                         "target_symptom": ["psychosis", "psychosis", "anxiety", "fatigue"],
                         "animal_readout": [False, True, True, False],
                         "species_ambiguous_readout": [False, False, False, True],
                         "curated_symptom": ["w", "w", "w", "w"], "reference": ["PMID:1"] * 4})
    assert [s.target for s in human_worded_statements(rows)] == ["A"]  # B is rodent-only, C is species-ambiguous
