"""Every registered rule, asserted against the implementation (the owner's diagnosis, 15 October 2026).

Why this exists. He observed that sessions "keep getting tripped up by stale things" and proposed overwriting the
registration rather than appending amendments, with the history generated from marked commits. Both help, and neither
catches the half of the problem that bit hardest: a rule changed in **code** with no amendment written leaves no trace
in the document at all. On 15 October the registration said the early-stopping validation set was "fold 0" while
`run_main_model.py` had long rotated it as `seed % num_validation_folds`, and it described reading
`configs/head_choice.json` when that file did not exist. No documentation standard finds those; a test does.

So `configs/registered_rules.json` holds each load-bearing rule once, and this module asserts each against the code it
names. A code change without a registration change fails here, and so does the reverse.

What it deliberately does not do: parse the registration's prose. Prose is where the staleness lives, and a parser
over it would be the fifth place the floor rule is stated. The JSON is the canonical statement and the prose refers to
it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
REGISTERED_RULES_PATH = REPOSITORY_ROOT / "configs" / "registered_rules.json"
sys.path.insert(0, str(REPOSITORY_ROOT / "experiments"))


@pytest.fixture(scope="module")
def registered_rules() -> dict:
    return json.loads(REGISTERED_RULES_PATH.read_text())["rules"]


def argument_parser_default(module_name: str, flag: str):
    """The default of one command-line flag, read from the module's source with the AST.

    Not by calling the module's build_argument_parser: `score_confirmatory` imports **run_main_model's** parser under
    that name, to reconstruct the arguments a run should have had, and builds its own inline in main(). Calling it
    returns a real parser that is the wrong one, which this helper did on its first attempt and which is the same
    mistake-shape as every broken join in this project - an object that looks right because it has the right type.
    """
    import ast

    source = (REPOSITORY_ROOT / "experiments" / f"{module_name}.py").read_text()
    for node in ast.walk(ast.parse(source)):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "add_argument"):
            continue
        flags = [argument.value for argument in node.args if isinstance(argument, ast.Constant)]
        if flag not in flags:
            continue
        for keyword in node.keywords:
            if keyword.arg == "default":
                return ast.literal_eval(keyword.value)
        return None  # a flag declared with no default, which for a registered rule is itself a finding
    raise AssertionError(f"{module_name} declares no flag {flag}")


def test_every_rule_names_a_value_a_meaning_and_a_code_location(registered_rules):
    """A rule with no code_location cannot be checked, which is how a rule goes stale unnoticed."""
    for rule_id, rule in registered_rules.items():
        assert "value" in rule, f"{rule_id} states no value"
        assert rule.get("means"), f"{rule_id} states no meaning, so a reader cannot tell what it governs"
        assert rule.get("code_location"), f"{rule_id} names no code location, so nothing can check it"


def test_the_floor_rule_is_none(registered_rules):
    """Registered 15 October 2026: no magnitude floor, the paired interval carries the decision."""
    assert registered_rules["floor_rule"]["value"] == "none"
    assert argument_parser_default("score_confirmatory", "--floor-rule") == "none"


def test_alpha_is_one_sided_025(registered_rules):
    assert registered_rules["alpha_one_sided"]["value"] == 0.025
    assert argument_parser_default("score_confirmatory", "--alpha") == 0.025


def test_the_bootstrap_is_a_cluster_bootstrap_over_leakage_groups(registered_rules):
    assert registered_rules["bootstrap_unit"]["value"] == "group"
    assert argument_parser_default("score_confirmatory", "--bootstrap-unit") == "group"


def test_the_deciding_bootstrap_draws_four_thousand_resamples(registered_rules):
    assert registered_rules["num_bootstrap"]["value"] == 4000
    assert argument_parser_default("score_confirmatory", "--num-bootstrap") == 4000


def test_two_models_are_tested_and_the_head_choice_names_them(registered_rules):
    """Two, not four: the sigmoid head was dropped because only noisy-OR has modules readable as pathways."""
    registered = registered_rules["tested_models"]["value"]
    assert len(registered) == 2, "four-model arithmetic was superseded by the head decision"
    head_choice = json.loads((REPOSITORY_ROOT / "configs" / "head_choice.json").read_text())
    assert sorted(head_choice["chosen"].values()) == sorted(registered)
    assert all("noisy_or" in model for model in registered)


def test_the_held_out_fraction_is_one_seventh(registered_rules):
    from mechanistic_pathway_learning.evaluation.target_held_out_split import REGISTERED_HELD_OUT_FRACTION

    assert registered_rules["held_out_fraction"]["value"] == pytest.approx(REGISTERED_HELD_OUT_FRACTION)
    assert REGISTERED_HELD_OUT_FRACTION == pytest.approx(1 / 7)


def test_the_minimum_positives_to_score_agree_between_the_rule_and_the_scorer(registered_rules):
    import importlib.util

    spec = importlib.util.spec_from_file_location("score_confirmatory", REPOSITORY_ROOT / "experiments" / "score_confirmatory.py")
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    assert registered_rules["minimum_positives_to_score"]["value"] == scorer.MINIMUM_POSITIVES_TO_SCORE


def test_the_normalisation_is_the_total_in_degree_divisor(registered_rules):
    assert registered_rules["normalisation"]["value"] == "total_in_degree"
    source = (REPOSITORY_ROOT / "experiments" / "run_main_model.py").read_text()
    assert '"--normalisation"' in source and 'default="total_in_degree"' in source


def test_the_validation_fold_rotates_with_the_seed(registered_rules):
    """The registration said "fold 0" until 15 October 2026; the code has long rotated it. This is the case that
    motivated the whole module, so it is asserted against the source rather than trusted."""
    assert registered_rules["validation_fold_rule"]["value"] == "seed % num_validation_folds"
    source = (REPOSITORY_ROOT / "experiments" / "run_main_model.py").read_text()
    assert "arguments.seed % num_validation_folds" in source, (
        "the validation fold no longer rotates with the seed; the registered rule and configs/registered_rules.json "
        "both need changing, and folds 5 and 6 may now be reachable")


def test_the_registration_points_at_the_rule_file_rather_than_restating_the_rules():
    """The standard this file enforces: prose refers to the canonical values, it does not duplicate them."""
    registration = (REPOSITORY_ROOT / "docs" / "preregistration.md").read_text()
    assert "configs/registered_rules.json" in registration, (
        "docs/preregistration.md does not reference the machine-checked rule file, so a reader has no way to know "
        "which of its statements is canonical")
