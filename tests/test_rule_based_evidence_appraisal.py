from mechanistic_pathway_learning.evidence.rule_based_evidence_appraisal import appraise_from_metadata, rubric_features, rubric_weight_from_metadata, study_design_from_metadata


def test_study_design_precedence_and_species() -> None:
    assert study_design_from_metadata(["Journal Article", "Randomized Controlled Trial"], ["9606"], "drug") == "human_randomized_trial"
    assert study_design_from_metadata(["Review", "Randomized Controlled Trial"], ["9606"], "drug") == "review_or_secondary"  # a review of trials is secondary
    assert study_design_from_metadata(["Case Reports"], ["9606"], "gene") == "human_case_report"
    assert study_design_from_metadata(["Journal Article"], ["10090"], "gene") == "animal_genetic_perturbation"
    assert study_design_from_metadata(["Journal Article"], ["10090"], "drug") == "animal_pharmacological"
    assert study_design_from_metadata(["Journal Article"], ["10090", "9606"], "gene") == "not_reported"  # mixed species, no design type
    assert study_design_from_metadata([], [], "gene") == "not_reported"
    assert study_design_from_metadata(["Comparative Study", "Journal Article"], ["10116"], "drug") == "animal_pharmacological"  # species before publication types
    assert study_design_from_metadata(["Case Reports"], ["9615"], "drug") == "animal_pharmacological"  # a veterinary case report is an animal study
    assert study_design_from_metadata(["Comparative Study"], ["9606"], "drug") == "not_reported"  # Comparative Study is scope, not a design
    assert study_design_from_metadata(["Review"], ["10090"], "gene") == "review_or_secondary"  # the review check still comes first


def test_weights_and_limitations() -> None:
    trial = appraise_from_metadata(["Randomized Controlled Trial"], ["9606"], ["Homo sapiens"], "drug", "symptom", True)
    assert rubric_weight_from_metadata(trial) == 1.0 and rubric_features(trial)["rubric_randomized_trial"] == 1.0
    diagnosis_review = appraise_from_metadata(["Review"], ["9606"], ["Homo sapiens"], "drug", "diagnosis", False)
    assert abs(rubric_weight_from_metadata(diagnosis_review) - 0.4 * 0.7) < 1e-9
    unfetched = appraise_from_metadata([], [], [], "gene", "symptom", True, metadata_available=False)
    assert unfetched.study_design == "metadata_not_fetched" and rubric_weight_from_metadata(unfetched) == 0.3 and rubric_features(unfetched)["rubric_metadata_available"] == 0.0
    assert "document metadata not fetched" in unfetched.limitations and "publication type not reported" not in unfetched.limitations
    fetched_uninformative = appraise_from_metadata(["Journal Article"], [], [], "gene", "symptom", True)
    assert fetched_uninformative.study_design == "not_reported" and rubric_weight_from_metadata(fetched_uninformative) == 0.35 and rubric_features(fetched_uninformative)["rubric_metadata_available"] == 1.0
    mixed_polarity = appraise_from_metadata([], ["9606"], ["Homo sapiens"], "drug", "mixed_polarity_diagnosis", True)
    assert mixed_polarity.measurement_level == "diagnosis"
    assert "secondary literature" in diagnosis_review.limitations and "descriptor is a diagnosis" in diagnosis_review.limitations and "publication year unknown" in diagnosis_review.limitations
    mouse = appraise_from_metadata([], ["10090"], ["Mus musculus"], "gene", "symptom", True)
    assert mouse.animal_only and "publication type not reported" in mouse.limitations and rubric_features(mouse)["rubric_animal_only"] == 1.0
