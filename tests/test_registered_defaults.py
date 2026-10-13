"""The trainer's defaults are the registered v4 model (the owner, 13 October 2026: "all of the v4 stuff should be ON by
default"), held equal to the explicit flag lists of experiments/run_main_model_batch.py so that a default cannot drift
from the registration without this test failing (docs/confirmatory_architecture_review.md, row 8)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.run_main_model import VALUES_BEFORE_THE_V4_DEFAULTS, build_argument_parser, resolve_registered_defaults  # noqa: E402
from experiments.run_main_model_batch import V4_LINEAR_RESPONSE_MODEL, V4_MESSAGE_PASSING_MODEL  # noqa: E402

# data paths stand in for the rebuilt v4 graph's files; the parser does not open them
DATA_PATHS = ["--graph-dir", "graph", "--evidence-dir", "evidence", "--node-descriptors", "descriptors.parquet"]
LINEAR_RESPONSE_DATA = [*DATA_PATHS, "--cell-class-weights", "weights.parquet"]


def resolved(argument_list):
    arguments = build_argument_parser().parse_args(argument_list)
    resolve_registered_defaults(arguments)
    return vars(arguments)


@pytest.mark.parametrize("data_paths, encoder, registered", [
    (DATA_PATHS, "message_passing", V4_MESSAGE_PASSING_MODEL),
    (LINEAR_RESPONSE_DATA, "linear_response", V4_LINEAR_RESPONSE_MODEL),
])
def test_the_defaults_are_the_registered_model(data_paths, encoder, registered):
    by_default = resolved([*data_paths, "--encoder", encoder])
    named = resolved([*data_paths, *registered])
    differing = {key: (by_default[key], named[key]) for key in named if by_default.get(key) != named[key]}
    assert not differing, f"a default differs from the registered model: {differing}"


def test_the_registered_choices_of_13_october_are_on_by_default():
    linear_response = resolved([*LINEAR_RESPONSE_DATA, "--encoder", "linear_response"])
    assert linear_response["class_weight_rule"] == "limiting_step"
    assert linear_response["input_weights"] == "positive"
    assert linear_response["descriptor_treatment"] == "seed_masked"
    assert linear_response["fixed_point_tolerance"] == pytest.approx(1e-3)
    assert linear_response["extracellular_coupling"] and linear_response["cofactor_relations"]
    message_passing = resolved([*DATA_PATHS, "--encoder", "message_passing"])
    assert message_passing["descriptor_treatment"] == "seed_masked"
    assert message_passing["message_passing_edge_signs"] == "graph"


def test_without_descriptors_there_is_nothing_to_mask_and_an_explicit_request_still_fails_loudly():
    assert resolved(["--graph-dir", "graph", "--evidence-dir", "evidence"])["descriptor_treatment"] == "plain"
    assert resolved(["--graph-dir", "graph", "--evidence-dir", "evidence", "--encoder", "linear_response"])["extracellular_coupling"] is False


def test_each_v4_default_can_still_be_turned_off_explicitly():
    arguments = resolved([*LINEAR_RESPONSE_DATA, "--encoder", "linear_response", "--class-weight-rule", "product", "--input-weights", "signed",
                          "--fixed-point-tolerance", "none", "--no-extracellular-coupling", "--no-cofactor-relations", "--no-refit-on-validation",
                          "--no-equal-drug-shares", "--no-start-at-weighted-optimum", "--descriptor-treatment", "plain"])
    assert arguments["class_weight_rule"] == "product" and arguments["input_weights"] == "signed"
    assert arguments["fixed_point_tolerance"] is None
    assert not any(arguments[key] for key in ("extracellular_coupling", "cofactor_relations", "refit_on_validation", "equal_drug_shares",
                                              "start_at_weighted_optimum"))


def test_a_run_made_before_the_v4_defaults_is_rebuilt_as_it_was_trained():
    """experiments/evaluate_on_rewired_graph.py fills an argument a results.json does not record from these values, so
    a change of default cannot rebuild an old run as a different model."""
    from experiments.evaluate_on_rewired_graph import stored_arguments
    old_record = {"arguments": {"graph_dir": "graph", "evidence_dir": "evidence", "encoder": "linear_response"}}
    rebuilt = vars(stored_arguments(old_record))
    for key, value in VALUES_BEFORE_THE_V4_DEFAULTS.items():
        assert rebuilt[key] == value, key
    defaults = vars(build_argument_parser().parse_args([]))
    assert set(VALUES_BEFORE_THE_V4_DEFAULTS) <= set(defaults)  # no stale name in the map
