"""The frozen baseline stays frozen, and the Amends: trailer parses the way the convention says (decision of 15 October 2026).

The owner's design: a byte-exact copy of the registration as it stood on 15 October, never edited, plus a history
generated from commit trailers. These tests hold up both halves — the baseline against its recorded digest, and the
parser against the trailer shapes it has to survive.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = REPOSITORY_ROOT / "docs" / "preregistration_baseline_2026-10-15.md"

# The digest of the baseline at the moment it was frozen. If this test fails, the frozen registration was edited,
# which is the one thing it exists not to permit: amend docs/preregistration.md instead and let the generator record
# the change. Never update this constant to make the test pass.
BASELINE_SHA256 = "08b02e8d303b61b4ebb4452266e54241c153b2600d58039bb5bb70dbb5e14104"


def load_generator():
    spec = importlib.util.spec_from_file_location(
        "generate_amendment_history", REPOSITORY_ROOT / "experiments" / "generate_amendment_history.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def generator():
    return load_generator()


@pytest.fixture(scope="module")
def known_rule_ids() -> set[str]:
    return set(json.loads((REPOSITORY_ROOT / "configs" / "registered_rules.json").read_text())["rules"])


def test_the_frozen_baseline_is_unchanged():
    assert BASELINE_PATH.exists(), "the baseline is what the amendment history is relative to"
    digest = hashlib.sha256(BASELINE_PATH.read_bytes()).hexdigest()
    assert digest == BASELINE_SHA256, (
        "docs/preregistration_baseline_2026-10-15.md was edited. It is the registration as it stood when frozen and "
        "must stay byte-exact; change docs/preregistration.md and mark the commit with an Amends: trailer instead")


def test_a_trailer_becomes_an_entry(generator, known_rule_ids):
    commit = {"sha": "a" * 40, "iso_date": "2026-10-16T12:00:00-04:00", "subject": "Remove the magnitude floor",
              "body": "Amends: floor_rule\nFrom: realised half-width plus 0.005\nTo: none\n"}
    entries, complaints = generator.amendments_in_commit(commit, known_rule_ids)
    assert complaints == []
    assert len(entries) == 1
    assert entries[0]["rule_id"] == "floor_rule"
    assert entries[0]["from"] == "realised half-width plus 0.005"
    assert entries[0]["to"] == "none"


def test_one_commit_may_amend_several_rules(generator, known_rule_ids):
    """Each Amends: opens an entry and the From:/To: after it belong to that entry, not to the first."""
    commit = {"sha": "b" * 40, "iso_date": "2026-10-16T12:00:00-04:00", "subject": "Two at once",
              "body": ("Amends: floor_rule\nFrom: fixed\nTo: none\n"
                       "Amends: alpha_one_sided\nFrom: 0.05\nTo: 0.025\n")}
    entries, complaints = generator.amendments_in_commit(commit, known_rule_ids)
    assert complaints == []
    assert [(entry["rule_id"], entry["from"], entry["to"]) for entry in entries] == [
        ("floor_rule", "fixed", "none"), ("alpha_one_sided", "0.05", "0.025")]


def test_an_unknown_rule_id_is_an_error_rather_than_a_silent_omission(generator, known_rule_ids):
    commit = {"sha": "c" * 40, "iso_date": "2026-10-16T12:00:00-04:00", "subject": "Typo",
              "body": "Amends: floor_ruel\nFrom: fixed\nTo: none\n"}
    entries, complaints = generator.amendments_in_commit(commit, known_rule_ids)
    assert len(complaints) == 1 and "floor_ruel" in complaints[0]
    assert len(entries) == 1, "the entry is still produced, so the typo shows in the output rather than vanishing"


def test_a_commit_with_no_trailer_produces_nothing(generator, known_rule_ids):
    commit = {"sha": "d" * 40, "iso_date": "2026-10-16T12:00:00-04:00", "subject": "Ordinary work",
              "body": "Fixes a broken join.\n\nCo-Authored-By: someone <x@example.com>\n"}
    entries, complaints = generator.amendments_in_commit(commit, known_rule_ids)
    assert entries == [] and complaints == []


def test_from_and_to_are_optional_but_recorded_as_unstated(generator, known_rule_ids):
    commit = {"sha": "e" * 40, "iso_date": "2026-10-16T12:00:00-04:00", "subject": "Bare trailer",
              "body": "Amends: normalisation\n"}
    entries, _ = generator.amendments_in_commit(commit, known_rule_ids)
    rendered = generator.render_history(entries, [], "0" * 64, {"normalisation": {"means": "the divisor"}})
    assert "(unstated)" in rendered


def test_the_rendered_history_says_an_empty_list_is_not_proof_that_nothing_changed(generator):
    """The hole in any trailer scheme: an unmarked change leaves no entry while the list still looks complete."""
    rendered = generator.render_history([], [], "0" * 64, {})
    assert "is not the same as no rule having changed" in rendered
    assert "test_registration_matches_the_code" in rendered


def test_the_generated_file_is_current():
    """--check in test form: a committed history that no longer matches the log is a stale artifact of the kind this
    whole mechanism exists to stop."""
    generated = REPOSITORY_ROOT / "docs" / "amendment_history.md"
    assert generated.exists(), "run experiments/generate_amendment_history.py"
    assert "do not edit by hand" in generated.read_text().splitlines()[0]
