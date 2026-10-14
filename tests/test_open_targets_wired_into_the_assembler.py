"""Decision 8.5, built 14 October 2026 on option (a), after the effect was measured rather than assumed.

8.5 and amendment item 12d are the same item: `reliability_reports` was written for the owner's instruction "add
it to the Dawid-Skene fit estimates" and `assemble_evidence_table.py` never called it, so the decision of 13
October was unbuilt. What these tests pin is the property that makes the wiring safe - the class changes a weight
through the reliability fit and creates no label - and the one guard that option (a) leaves inactive.
"""
import inspect

import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence import assemble_evidence_table
from mechanistic_pathway_learning.evidence.load_open_targets_statements import (
    CuratedStatement,
    reliability_reports,
)

EXISTING_ITEM = ("GRIN1", "psychosis", "induces")
PERTURBATION_OF_TARGET = {"GRIN1": "GRIN1"}


def a_statement(target: str, direction: str, symptom: str) -> CuratedStatement:
    return CuratedStatement(target, direction, symptom, f"{direction} of {target}", "pmid:1")


def test_the_assembler_takes_the_scoped_rows_path_so_the_class_can_reach_the_fit():
    """The wiring 8.5 was about: the loader existed, was tested, and had no caller."""
    assert "open_targets_scoped_rows_path" in inspect.signature(assemble_evidence_table.assemble).parameters


def test_the_assembler_imports_the_loader_it_had_never_called():
    assert hasattr(assemble_evidence_table, "reliability_reports")
    assert hasattr(assemble_evidence_table, "human_worded_statements")


def test_a_report_is_emitted_only_for_an_item_that_already_carries_a_positive_report():
    """The load-bearing guard, and the one that is active. Without it a target-safety statement could turn a
    documented absence into a trained and scored positive."""
    on_an_existing_positive = reliability_reports([a_statement("GRIN1", "inhibition", "psychosis")],
                                                 {EXISTING_ITEM}, PERTURBATION_OF_TARGET)
    assert len(on_an_existing_positive) == 1
    on_an_item_with_no_positive = reliability_reports([a_statement("GRIN1", "inhibition", "psychosis")],
                                                      set(), PERTURBATION_OF_TARGET)
    assert on_an_item_with_no_positive == []


def test_a_statement_about_a_target_that_is_not_a_perturbation_emits_nothing():
    assert reliability_reports([a_statement("CHRM1", "inhibition", "psychosis")],
                               {EXISTING_ITEM}, PERTURBATION_OF_TARGET) == []


def test_the_evidenced_negative_guard_still_works_when_a_caller_can_supply_the_set():
    """Option (a) leaves this inactive inside assemble(), because the set is computed downstream. The guard itself
    is intact, so a later caller that has the set can use it - which is what option (b) would have been."""
    skipped = reliability_reports([a_statement("GRIN1", "inhibition", "psychosis")], {EXISTING_ITEM},
                                  PERTURBATION_OF_TARGET, evidenced_negative_pairs={("GRIN1", "psychosis")})
    assert skipped == []


def test_the_emitted_report_names_the_class_and_carries_a_positive_value():
    report = reliability_reports([a_statement("GRIN1", "inhibition", "psychosis")],
                                 {EXISTING_ITEM}, PERTURBATION_OF_TARGET)[0]
    assert report["report_value"] == 1
    assert "Open Targets" in report["source"]
    assert report["perturbation_id"] == "GRIN1"
    assert report["symptom"] == "psychosis"
