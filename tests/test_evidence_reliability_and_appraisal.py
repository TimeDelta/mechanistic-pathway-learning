"""Tests: the reliability model recovers source quality from agreement; the appraisal module parses, caches and scores."""
import json
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evidence.evidence_reliability_model import (
    anchor_pseudo_counts_from_posterior,
    fit_dawid_skene,
    fit_feature_dependent_sensitivity,
    observation_weights_from_posterior,
    per_item_reliability,
    posterior_truth_probability,
)
from mechanistic_pathway_learning.evidence.llm_evidence_appraisal import (
    EvidenceRubric,
    appraise_document,
    build_appraisal_prompt,
    gold_set_agreement,
    parse_rubric_response,
)


def simulate_reports(num_items: int, true_sensitivity: list[float], true_specificity: list[float], prevalence: float, seed: int = 0):
    generator = np.random.default_rng(seed)
    truth = generator.random(num_items) < prevalence
    reports = np.zeros((num_items, len(true_sensitivity)))
    for source_index, (sensitivity, specificity) in enumerate(zip(true_sensitivity, true_specificity)):
        positive_draw = generator.random(num_items) < sensitivity
        negative_draw = generator.random(num_items) < (1.0 - specificity)
        reports[:, source_index] = np.where(truth, positive_draw, negative_draw)
    coverage = np.ones_like(reports)
    return truth, reports, coverage


def test_dawid_skene_recovers_reliable_versus_noisy_sources() -> None:
    # three sources: a two-class latent model is not identified with fewer
    truth, reports, coverage = simulate_reports(3000, [0.95, 0.80, 0.55], [0.95, 0.85, 0.55], prevalence=0.3)
    parameters, posterior = fit_dawid_skene(reports, coverage, ["monogenic", "label", "comention"])
    assert parameters.sensitivity[0] > parameters.sensitivity[2] + 0.2
    assert parameters.specificity[0] > parameters.specificity[2] + 0.2
    assert abs(parameters.prevalence - 0.3) < 0.08
    recovered = (posterior > 0.5) == truth
    assert recovered.mean() > 0.9
    assert parameters.log_likelihood_trace[-1] >= parameters.log_likelihood_trace[0] - 1e-6


def test_posterior_respects_coverage() -> None:
    reports = np.array([[1.0, 0.0], [1.0, 0.0]])
    coverage = np.array([[1.0, 1.0], [1.0, 0.0]])  # second item not covered by the second source
    posterior = posterior_truth_probability(reports, coverage, 0.3, np.array([0.9, 0.9]), np.array([0.9, 0.9]))
    assert posterior[1] > posterior[0]  # a silent uncovered source is not a 0 report


def test_feature_dependent_sensitivity_learns_design_effect() -> None:
    generator = np.random.default_rng(1)
    num_items = 2000
    good_design = generator.random(num_items) < 0.5
    truth = generator.random(num_items) < 0.4
    sensitivity = np.where(good_design, 0.95, 0.5)
    reports = np.where(truth, generator.random(num_items) < sensitivity, generator.random(num_items) < 0.1).astype(float)
    features = good_design[:, None].astype(float)
    sensitivity_model, specificity_model = fit_feature_dependent_sensitivity(features, reports, truth.astype(float))
    per_item_sensitivity, _ = per_item_reliability(sensitivity_model, specificity_model, features, 0.8, 0.8)
    assert per_item_sensitivity[good_design].mean() > per_item_sensitivity[~good_design].mean() + 0.3


def test_weights_and_anchor_counts() -> None:
    posterior = np.array([0.5, 0.99, 0.999999])
    assert np.allclose(observation_weights_from_posterior(posterior, global_scale=2.0), [1.0, 1.98, 1.999998])
    counts = anchor_pseudo_counts_from_posterior(posterior, maximum_pseudo_count=99.0)
    assert abs(counts[0] - 1.0) < 1e-9 and abs(counts[1] - 99.0) < 1e-6 and counts[2] == 99.0


def test_appraisal_parses_caches_and_scores(tmp_path: Path) -> None:
    calls = []

    def fake_client(prompt: str) -> str:
        calls.append(prompt)
        assert "Do not use outside knowledge" in prompt
        return json.dumps({
            "study_design": "human_monogenic_case_series", "species": "human", "sample_size": 12,
            "measurement_level": "diagnosis", "effect_direction": "induces", "cohort_identifier": "not_reported",
            "mechanistic_specificity": "single_target", "supporting_span": "psychotic symptoms were present in 12 patients",
            "extraction_confidence": 0.8,
        })

    rubric = appraise_document(fake_client, "fake-model", "Twelve patients with HMBS mutations had psychotic symptoms.", "HMBS loss of function", "psychosis", cache_directory=tmp_path)
    assert rubric.study_design == "human_monogenic_case_series" and rubric.sample_size == 12
    assert len(rubric.feature_vector()) == len(rubric.feature_vector())  # stable length
    again = appraise_document(fake_client, "fake-model", "Twelve patients with HMBS mutations had psychotic symptoms.", "HMBS loss of function", "psychosis", cache_directory=tmp_path)
    assert again == rubric and len(calls) == 1  # second call served from cache
    assert "psychosis" in build_appraisal_prompt("text", "HMBS", "psychosis")
    fenced = parse_rubric_response("```json\n{\"study_design\": \"made_up\", \"effect_direction\": \"relieves\"}\n```")
    assert fenced.study_design == "not_reported" and fenced.effect_direction == "relieves"
    human = EvidenceRubric("human_monogenic_case_series", "human", 12, "diagnosis", "induces", "not_reported", "single_target", "", 1.0)
    agreement = gold_set_agreement([rubric, fenced], [human, human])
    assert 0.0 <= agreement["sample_size_exact_match"] <= 1.0 and "study_design_kappa" in agreement
