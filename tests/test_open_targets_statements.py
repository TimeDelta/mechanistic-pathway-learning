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
    reliability_reports,
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


def test_a_statement_never_overrides_an_evidenced_negative_on_any_perturbation_it_could_live_on():
    """A never-above-placebo pair and a review statement are a contradiction, not independent evidence. The guard
    has to be keyed on the perturbations an evidenced negative can actually live on - the gene, and the drugs that
    produce this action - and not only on the ACTIVATION id this module invents, which appears in no prior table
    and so could never match."""
    for perturbation in ("HTR2C", "single_target_agonist", "two_target_agonist", activation_perturbation_id("HTR2C")):
        placed = placements_for([statement("HTR2C", "activation", "anxiety")],
                                evidenced_negatives={(perturbation, "anxiety")})
        assert [p.perturbation_identifier for p in placed] == [None], f"the guard missed an evidenced negative on {perturbation}"
        assert placed[0].reason == "the pair is an evidenced negative"
    assert placements_for([statement("HTR2C", "activation", "anxiety")],
                          evidenced_negatives={("antagonist", "anxiety")})[0].perturbation_identifier is not None  # wrong direction


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


EXISTING_ITEMS = {("GRIN1", "psychosis", "induces")}
PERTURBATION_OF_TARGET = {"GRIN1": "GRIN1", "HTR2C": "HTR2C"}


def test_a_reliability_report_is_emitted_only_for_an_item_another_source_already_reports():
    """The assembler makes one observation per item with an explicit report, so a report on a new item would create
    a scored pair out of a review statement - which is what the owner ruled out."""
    reports = reliability_reports([statement("GRIN1", "inhibition", "psychosis"),      # the item exists
                                   statement("GRIN1", "inhibition", "anxiety"),        # the item does not
                                   statement("HTR2C", "activation", "anxiety")],       # nor this one
                                  EXISTING_ITEMS, PERTURBATION_OF_TARGET)
    assert [(r["perturbation_id"], r["symptom"]) for r in reports] == [("GRIN1", "psychosis")]


def test_the_report_carries_the_class_so_it_reaches_the_fit_and_not_the_grade():
    report = reliability_reports([statement("GRIN1", "inhibition", "psychosis")], EXISTING_ITEMS, PERTURBATION_OF_TARGET)[0]
    assert report["evidence_class"] == OPEN_TARGETS_EVIDENCE_CLASS  # the grade reads the pair's primary class only
    assert report["report_value"] == 1 and report["relation"] == "induces"
    assert report["pubmed_reference_count"] == 2


def test_a_target_with_no_perturbation_reports_nothing():
    assert reliability_reports([statement("CHRM1", "inhibition", "psychosis")], EXISTING_ITEMS, PERTURBATION_OF_TARGET) == []


def test_a_reliability_report_needs_an_item_that_already_carries_a_positive():
    """An item can exist as a recorded absence: evidence_full_v4 holds 15 gene observations with no positive
    report. A report of value 1 on one of those would flip a documented absence into a trained and scored
    positive, which is the outcome the module exists to prevent."""
    only_an_absence = set()  # the item exists in the table but carries no positive report
    assert reliability_reports([statement("GRIN1", "inhibition", "psychosis")], only_an_absence, PERTURBATION_OF_TARGET) == []
    assert len(reliability_reports([statement("GRIN1", "inhibition", "psychosis")], EXISTING_ITEMS, PERTURBATION_OF_TARGET)) == 1


def test_a_reliability_report_is_not_emitted_against_an_evidenced_negative():
    assert reliability_reports([statement("GRIN1", "inhibition", "psychosis")], EXISTING_ITEMS, PERTURBATION_OF_TARGET,
                               evidenced_negative_pairs={("GRIN1", "psychosis")}) == []


def test_the_sole_target_rule_is_available_and_differs_from_any_drug():
    from mechanistic_pathway_learning.evidence.load_open_targets_statements import (
        drugs_acting_on_the_target,
        drugs_whose_only_target_is,
    )

    assert drugs_acting_on_the_target(SIGNS_OF_DRUG, "HTR2C", "activation") == {"single_target_agonist", "two_target_agonist"}
    assert drugs_whose_only_target_is(TARGETS_OF_DRUG, SIGNS_OF_DRUG, "HTR2C", "activation") == {"single_target_agonist"}
