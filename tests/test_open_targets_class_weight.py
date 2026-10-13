"""Amendment item 12: the Open Targets statements carry a fixed weight of 0.45, with 0.0 as the declared secondary."""
from __future__ import annotations

import pytest

from mechanistic_pathway_learning.evidence.assign_evidence_grades import (
    DEFAULT_GRADE_WEIGHTS,
    DEFAULT_OPEN_TARGETS_WEIGHT,
    OPEN_TARGETS_EVIDENCE_CLASS,
    OPEN_TARGETS_GRADE,
    OPEN_TARGETS_WEIGHT_WITHOUT_THE_CLASS,
    EvidenceRecord,
    assign_evidence_grade,
    loss_weight_for_record,
)


TRIPLES_A_STUDY_DRUG_CAN_CHECK = 71
TRIPLES_THE_DRUG_LABEL_CONFIRMS = 54
ROUNDING_STEP_OF_THE_WEIGHT_SCALE = 0.05


def open_targets_record(**overrides) -> EvidenceRecord:
    fields = {"perturbation_identifier": "GRIN1", "symptom_identifier": "psychosis", "relation": "induces",
              "evidence_class": OPEN_TARGETS_EVIDENCE_CLASS, "source": "Open Targets target safety"}
    return EvidenceRecord(**{**fields, **overrides})


def test_the_weight_is_the_registered_0_45_and_no_grade_changes_it():
    assert loss_weight_for_record(open_targets_record()) == pytest.approx(DEFAULT_OPEN_TARGETS_WEIGHT)
    assert DEFAULT_OPEN_TARGETS_WEIGHT == pytest.approx(0.45)


def test_0_45_is_grade_b_times_the_measured_agreement_rounded_down():
    """The derivation item 12 states: 0.6 x 0.7606 = 0.4563, rounded down to the nearest 0.05. The rounding goes down
    because the agreement is measured on the triples a study drug can check and the statements that train are the ones
    no study drug label makes, so the measured number is an extrapolation and the lower weight is the conservative
    reading. Rounding to the nearest 0.01 would give 0.46, which is what the item said before 13 October 2026."""
    agreement = TRIPLES_THE_DRUG_LABEL_CONFIRMS / TRIPLES_A_STUDY_DRUG_CAN_CHECK
    assert agreement == pytest.approx(0.7606, abs=0.0001)
    derived = DEFAULT_GRADE_WEIGHTS[OPEN_TARGETS_GRADE] * agreement
    assert derived == pytest.approx(0.4563, abs=0.0001)
    rounded_down = derived // ROUNDING_STEP_OF_THE_WEIGHT_SCALE * ROUNDING_STEP_OF_THE_WEIGHT_SCALE
    assert rounded_down == pytest.approx(DEFAULT_OPEN_TARGETS_WEIGHT)


def test_the_secondary_arm_keeps_the_class_out_of_training():
    assert loss_weight_for_record(open_targets_record(), open_targets_weight=OPEN_TARGETS_WEIGHT_WITHOUT_THE_CLASS) == 0.0
    assert OPEN_TARGETS_WEIGHT_WITHOUT_THE_CLASS == 0.0


def test_the_frequency_scaling_of_grade_b_does_not_reach_the_class():
    """A grade B pharmacological record is scaled by its label frequency; an Open Targets statement has no frequency,
    so the weight has to be the flat one whatever the record carries."""
    with_a_frequency = open_targets_record(label_event_frequency=0.01)
    assert loss_weight_for_record(with_a_frequency) == pytest.approx(DEFAULT_OPEN_TARGETS_WEIGHT)


def test_the_grade_records_where_the_weight_came_from():
    assert assign_evidence_grade(open_targets_record()) == OPEN_TARGETS_GRADE
