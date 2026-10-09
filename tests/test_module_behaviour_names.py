"""Tests for naming a module from its behaviour: its drivers, their reach into its support and the symptoms it feeds."""
import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.module_behaviour import (
    behaviour_label,
    domain_shares,
    hops_into_support,
    module_drivers,
    module_symptom_side,
    propagation_hops,
    reverse_adjacency,
    symptom_rdoc_domains,
)

CROSSWALK = pd.DataFrame({
    "target_symptom": ["parkinsonism", "catatonia", "insomnia", "fatigue", "anxiety"],
    "rdoc_domain": ["sensorimotor", "arousal;sensorimotor", "arousal_regulatory", "arousal_regulatory", "negative_valence"],
})


def test_a_module_whose_activation_does_not_vary_has_no_drivers() -> None:
    flat = np.full(50, 0.3) + np.linspace(0.0, 0.01, 50)  # standard deviation well under the responding threshold
    reading = module_drivers(flat, [f"GENE{i}" for i in range(50)], ["gene"] * 50)
    assert reading["responds"] is False
    assert "num_drivers" not in reading
    assert "no driver" in behaviour_label({"available": True, "drivers": reading, "symptom_side": {}})


def test_drivers_are_the_top_fifth_and_a_curated_module_names_them() -> None:
    activations = np.concatenate([np.full(40, 0.1), np.linspace(0.6, 0.9, 10)])
    labels = [f"OTHER{i}" for i in range(40)] + ["ALAD", "FECH", "HMBS", "UROD", "CPOX"] + [f"REST{i}" for i in range(5)]
    reading = module_drivers(activations, labels, ["gene"] * 50,
                             curated_module_genes={"heme_porphyrin": {"ALAD", "FECH", "HMBS", "UROD", "CPOX"}})
    assert reading["responds"] is True
    assert reading["num_drivers"] == 10
    assert reading["group"]["name"] == "heme_porphyrin"
    assert reading["group"]["drivers"] == 5


def test_a_group_holding_no_more_than_the_background_does_not_name_the_drivers() -> None:
    activations = np.concatenate([np.full(40, 0.1), np.linspace(0.6, 0.9, 10)])
    labels = [f"DRUG{i}" for i in range(50)]
    every_perturbation_is_a_drug = ["drug"] * 50
    reading = module_drivers(activations, labels, every_perturbation_is_a_drug)
    assert reading["group"] is None  # the drivers are all drugs, but so is every other perturbation


def test_reach_is_read_against_the_perturbations_that_are_not_drivers() -> None:
    activations = np.concatenate([np.full(40, 0.1), np.linspace(0.6, 0.9, 10)])
    labels = [f"GENE{i}" for i in range(50)]
    reaching_drivers_only = {position: 1 for position in range(40, 48)}
    reading = module_drivers(activations, labels, ["gene"] * 50, hops_of_perturbation=reaching_drivers_only, hop_budget=2)
    assert reading["drivers_reaching_the_support"] == 8
    assert reading["others_reaching_the_support"] == 0
    assert reading["reach_is_enriched"] is True
    everything_reaches = {position: 1 for position in range(50)}
    flat_reach = module_drivers(activations, labels, ["gene"] * 50, hops_of_perturbation=everything_reaches, hop_budget=2)
    assert flat_reach["reach_is_enriched"] is False
    assert "does not follow the path" in behaviour_label({"available": True, "drivers": flat_reach, "symptom_side": {}})


def test_hops_into_support_walks_backwards_and_stops_at_the_budget() -> None:
    sources = np.array(["a", "b", "c", "x"])
    targets = np.array(["b", "c", "d", "d"])
    reverse = reverse_adjacency(sources, targets, excluded_nodes={"x"})
    distance, truncated = hops_into_support({"d"}, reverse, max_hops=2)
    assert truncated is False
    assert distance == {"d": 0, "c": 1, "b": 2}  # "a" is three hops away and "x" is currency
    assert hops_into_support({"d"}, reverse, max_hops=0)[0] == {"d": 0}


def test_a_link_profile_matching_the_other_modules_names_no_symptom_family() -> None:
    symptoms = ["parkinsonism", "catatonia", "insomnia", "fatigue", "anxiety"]
    domains = symptom_rdoc_domains(CROSSWALK)
    nearly_identical = np.array([[0.25, 0.24, 0.10, 0.11, 0.30],
                                 [0.26, 0.25, 0.11, 0.10, 0.31],
                                 [0.24, 0.23, 0.10, 0.12, 0.29]])
    reading = module_symptom_side(nearly_identical, 0, symptoms, domains)
    assert reading["distinct_profile"] is False
    assert reading["domain"] is None
    assert "matches the other modules" in behaviour_label({"available": True, "drivers": {"responds": False, "activation_standard_deviation": 0.0}, "symptom_side": reading})


def test_a_distinct_profile_takes_the_domain_of_the_symptoms_it_feeds() -> None:
    symptoms = ["parkinsonism", "catatonia", "insomnia", "fatigue", "anxiety"]
    domains = symptom_rdoc_domains(CROSSWALK)
    links = np.array([[0.60, 0.50, 0.05, 0.05, 0.05],
                      [0.10, 0.10, 0.40, 0.40, 0.40],
                      [0.10, 0.10, 0.40, 0.40, 0.40]])
    reading = module_symptom_side(links, 0, symptoms, domains)
    assert reading["distinct_profile"] is True
    assert reading["symptoms"] == ["parkinsonism", "catatonia"]
    assert reading["domain"]["name"] == "sensorimotor"
    assert "feeds sensorimotor symptoms" in behaviour_label({"available": True, "drivers": {"responds": False, "activation_standard_deviation": 0.0}, "symptom_side": reading})


def test_one_symptom_alone_does_not_name_a_domain() -> None:
    symptoms = ["parkinsonism", "catatonia", "insomnia", "fatigue", "anxiety"]
    domains = symptom_rdoc_domains(CROSSWALK)
    links = np.array([[0.60, 0.10, 0.10, 0.10, 0.10],
                      [0.10, 0.30, 0.30, 0.30, 0.30],
                      [0.10, 0.30, 0.30, 0.30, 0.30]])
    reading = module_symptom_side(links, 0, symptoms, domains)
    assert reading["symptoms"] == ["parkinsonism"]
    assert reading["domain"] is None


def test_domain_shares_split_a_symptom_between_its_domains() -> None:
    domains = symptom_rdoc_domains(CROSSWALK)
    assert domains["catatonia"] == ["arousal", "sensorimotor"]
    shares = domain_shares(["parkinsonism", "catatonia"], domains)
    assert shares["sensorimotor"] == 0.75
    assert shares["arousal"] == 0.25


def test_the_hop_budget_follows_the_encoder_the_run_used() -> None:
    assert propagation_hops({"encoder": "linear_response", "propagation_steps": 8, "num_layers": 2}) == 8
    assert propagation_hops({"encoder": "message_passing", "propagation_steps": 8, "num_layers": 2}) == 2


def test_a_curated_module_that_holds_as_many_non_drivers_does_not_name_the_drivers() -> None:
    activations = np.concatenate([np.full(40, 0.1), np.linspace(0.6, 0.9, 10)])
    # the curated module holds three of the ten drivers and twelve of the forty perturbations that are not drivers,
    # so its share of the drivers is no larger than its share of the scored perturbations
    everywhere = {f"OTHER{i}" for i in range(12)} | {"ALAD", "FECH", "HMBS"}
    labels = [f"OTHER{i}" for i in range(40)] + ["ALAD", "FECH", "HMBS"] + [f"REST{i}" for i in range(7)]
    reading = module_drivers(activations, labels, ["gene"] * 50, curated_module_genes={"heme_porphyrin": everywhere})
    assert reading["group"] is None
