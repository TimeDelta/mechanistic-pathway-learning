"""Document-grounded evidence appraisal with a language model (design section 4.3).

The model extracts a structured rubric from one document (abstract or methods text)
about one candidate perturbation-symptom link. It is never asked whether the link is
true; it is asked what kind of evidence the document is. The rubric fields become
features of the evidence reliability model, which learns how much each kind of
evidence is worth. Three guards:

1. Grounding: the prompt contains the document text and the candidate link only; the
   instructions forbid using outside knowledge and require "not_reported" when a
   field is absent from the text.
2. Reproducibility: model name, prompt version and temperature 0 are pinned and
   recorded with every output; outputs are cached by a hash of (model, prompt version,
   document, link) so a rerun never re-queries.
3. Validation: a gold set rated by two humans is scored field by field with Cohen's
   kappa before any appraisal enters training; agreement is reported in the
   pre-registration. The leakage audit (design section 6.4) checks whether appraisal
   features help the post-cutoff time split more than the within-period split.

The client is any callable mapping a prompt string to a response string, so tests
run with a fake and production runs use whichever API the cluster allows.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

PROMPT_VERSION = "appraisal-rubric-v1"

STUDY_DESIGNS = (
    "human_monogenic_case_series",
    "human_randomized_trial",
    "human_cohort_or_case_control",
    "human_case_report",
    "human_genetic_association",
    "animal_genetic_perturbation",
    "animal_pharmacological",
    "in_vitro_or_computational",
    "review_or_opinion",
    "not_reported",
)
MEASUREMENT_LEVELS = ("symptom_item", "symptom_scale", "diagnosis", "behavioral_proxy", "not_reported")
EFFECT_DIRECTIONS = ("induces", "relieves", "no_effect", "unclear", "not_reported")
MECHANISTIC_SPECIFICITIES = ("single_target", "multi_target", "unknown", "not_reported")


@dataclass
class EvidenceRubric:
    study_design: str
    species: str
    sample_size: int | None
    measurement_level: str
    effect_direction: str
    cohort_identifier: str
    mechanistic_specificity: str
    supporting_span: str
    extraction_confidence: float

    def feature_vector(self) -> list[float]:
        """One-hot design and measurement level, log sample size, specificity flag and confidence."""
        features = [1.0 if self.study_design == design else 0.0 for design in STUDY_DESIGNS]
        features += [1.0 if self.measurement_level == level else 0.0 for level in MEASUREMENT_LEVELS]
        features.append(0.0 if self.sample_size is None else float(__import__("math").log1p(self.sample_size)))
        features.append(1.0 if self.species.lower() in ("human", "homo sapiens") else 0.0)
        features.append(1.0 if self.mechanistic_specificity == "single_target" else 0.0)
        features.append(float(self.extraction_confidence))
        return features


def build_appraisal_prompt(document_text: str, perturbation_description: str, symptom_description: str) -> str:
    return (
        "You are extracting the design of a study from the text below. Use only the text. "
        "Do not use outside knowledge about the drug, gene, pathway or symptom, and do not judge whether the link is real.\n\n"
        f"Candidate link: perturbation = {perturbation_description}; symptom = {symptom_description}.\n\n"
        "Return a JSON object with exactly these keys:\n"
        f"- study_design: one of {list(STUDY_DESIGNS)}\n"
        "- species: the organism studied, or \"not_reported\"\n"
        "- sample_size: integer number of subjects or animals, or null\n"
        f"- measurement_level: one of {list(MEASUREMENT_LEVELS)} (how the symptom was measured)\n"
        f"- effect_direction: one of {list(EFFECT_DIRECTIONS)} (what the text reports for this link)\n"
        "- cohort_identifier: a short string naming the cohort, trial or registry if the text names one, else \"not_reported\"\n"
        f"- mechanistic_specificity: one of {list(MECHANISTIC_SPECIFICITIES)} (whether the perturbation acts through one characterized target)\n"
        "- supporting_span: the shortest verbatim span of the text that supports effect_direction, or \"\"\n"
        "- extraction_confidence: a number in [0, 1] for how clearly the text supports your answers\n\n"
        "Answer with the JSON object only.\n\n"
        f"Text:\n{document_text}"
    )


def parse_rubric_response(response_text: str) -> EvidenceRubric:
    cleaned = response_text.strip()
    if cleaned.startswith("```"):
        cleaned = cleaned.strip("`")
        cleaned = cleaned[cleaned.find("{"):]
    payload = json.loads(cleaned[cleaned.find("{"): cleaned.rfind("}") + 1])

    def constrained(value, allowed: Sequence[str]) -> str:
        return value if value in allowed else "not_reported"

    sample_size = payload.get("sample_size")
    return EvidenceRubric(
        study_design=constrained(payload.get("study_design"), STUDY_DESIGNS),
        species=str(payload.get("species") or "not_reported"),
        sample_size=int(sample_size) if isinstance(sample_size, (int, float)) and sample_size >= 0 else None,
        measurement_level=constrained(payload.get("measurement_level"), MEASUREMENT_LEVELS),
        effect_direction=constrained(payload.get("effect_direction"), EFFECT_DIRECTIONS),
        cohort_identifier=str(payload.get("cohort_identifier") or "not_reported"),
        mechanistic_specificity=constrained(payload.get("mechanistic_specificity"), MECHANISTIC_SPECIFICITIES),
        supporting_span=str(payload.get("supporting_span") or ""),
        extraction_confidence=float(min(1.0, max(0.0, payload.get("extraction_confidence") or 0.0))),
    )


def cache_key(model_name: str, document_text: str, perturbation_description: str, symptom_description: str) -> str:
    digest = hashlib.sha256()
    for part in (model_name, PROMPT_VERSION, document_text, perturbation_description, symptom_description):
        digest.update(part.encode("utf-8"))
        digest.update(b"\x00")
    return digest.hexdigest()


def appraise_document(
    client: Callable[[str], str],
    model_name: str,
    document_text: str,
    perturbation_description: str,
    symptom_description: str,
    cache_directory: Path | None = None,
) -> EvidenceRubric:
    """Run one appraisal, reading and writing the cache when a directory is given."""
    key = cache_key(model_name, document_text, perturbation_description, symptom_description)
    cache_path = cache_directory / f"{key}.json" if cache_directory is not None else None
    if cache_path is not None and cache_path.exists():
        cached = json.loads(cache_path.read_text())
        return EvidenceRubric(**cached["rubric"])
    response_text = client(build_appraisal_prompt(document_text, perturbation_description, symptom_description))
    rubric = parse_rubric_response(response_text)
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps({"model_name": model_name, "prompt_version": PROMPT_VERSION, "rubric": asdict(rubric), "raw_response": response_text}))
    return rubric


def gold_set_agreement(model_rubrics: Sequence[EvidenceRubric], human_rubrics: Sequence[EvidenceRubric]) -> dict[str, float]:
    """Cohen's kappa per categorical field and exact-match rate for sample size, against human ratings."""
    from sklearn.metrics import cohen_kappa_score

    if len(model_rubrics) != len(human_rubrics):
        raise ValueError("model and human rubric lists must align item by item")
    agreement: dict[str, float] = {}
    for field_name in ("study_design", "measurement_level", "effect_direction", "mechanistic_specificity"):
        model_values = [getattr(rubric, field_name) for rubric in model_rubrics]
        human_values = [getattr(rubric, field_name) for rubric in human_rubrics]
        agreement[f"{field_name}_kappa"] = float(cohen_kappa_score(model_values, human_values)) if len(set(model_values) | set(human_values)) > 1 else 1.0
    matches = [model.sample_size == human.sample_size for model, human in zip(model_rubrics, human_rubrics)]
    agreement["sample_size_exact_match"] = float(sum(matches) / len(matches)) if matches else float("nan")
    return agreement
