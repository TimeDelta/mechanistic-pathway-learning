"""Evidence class E2, later label slice: OnSIDES v3.1.1 label events as report-level rows (design section 4.2).

OnSIDES (Tanaka et al. 2025, reference [19]) extracts adverse-event terms from current drug labels with a
PubMedBERT model: from the Adverse Reactions (AR), Boxed Warning (BW) and Warnings and Precautions (WP)
sections of US labels, and from the undivided undesirable-effects text of UK, EU and Japanese labels, whose
events carry label_section "NA" (every non-US event in v3.1.1 and no US event). Its v3.1.1 data release
(2026-04-22, CC BY 4.0) has 6.9 million label-event mentions over 41,119 labels and 1,866 RxNorm
ingredients. Generic products repeat the same label text, so the unit of evidence here is the deduplicated
label statement: one report per (ingredient, target symptom, MedDRA term, label section, label country)
with the number of distinct labels as a feature. A non-US statement gets evidence_code non_us_label_section
and counts as adverse-reactions text for rubric_adverse_reactions_section (the extraction covers that
section), so the section rubric features are not a second copy of rubric_us_label; rubric_boxed_warning can
only be set on US labels. OnSIDES carries no frequencies and no label dates; the data release date is an
upper bound on when a statement appeared, which the limitations text says.

Events are coded at PT or LLT level. A PT is matched to a target symptom by its crosswalk name; an LLT is
matched through its PT, using the LLT-to-PT pairing SIDER's meddra_all_se.tsv.gz carries per label term
(sider_lower_level_to_preferred_terms), and the statement keeps the source term type. Without that pairing
only PT-coded events match, since no LLT name equals a PT name in the OnSIDES vocabulary.

Drugs are keyed by RxNorm ingredient; the identifier bridge (onsides_identifier_bridge.py) takes each
ingredient to its ChEMBL parent, its mechanism targets and, when the parent or the name is a SIDER drug's,
to that drug's PubChem identifier so SIDER and OnSIDES report on one perturbation. The qualification rule
is the one SIDER rows pass (assemble_evidence_table.py): a single mechanism target and the ATC N proxy for
central action, so an ingredient without any ATC code does not qualify. A statement carried only by
combination-product labels is disqualified (combination_label_only): the label of carbidopa/levodopa
says nothing about carbidopa alone.

Rows are grade B label evidence like SIDER rows, never literature. MedDRA term names stay inside
data/processed and are blanked by the data release.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS
from mechanistic_pathway_learning.evidence.load_drug_label_events import is_nervous_system_atc
from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import IngredientBridge, mechanism_targets_for_parent, passes_single_target_rule
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import GraphNodeLookup, graph_compound_targets
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import DrugTarget

ONSIDES_RELEASE = "v3.1.1"
ONSIDES_RELEASE_DATE = "2026-04-22"
ONSIDES_SOURCE = "OnSIDES-label"
NON_US_LABEL_SECTION = "NA"  # OnSIDES' section value for the undivided non-US label (UK, EU, JP), not an unsectioned mention
NON_US_SECTION_EVIDENCE_CODE = "non_us_label_section"
SECTION_EVIDENCE_CODES = {"AR": "adverse_reactions_section", "BW": "boxed_warning_section", "WP": "warnings_precautions_section", NON_US_LABEL_SECTION: NON_US_SECTION_EVIDENCE_CODE}
OTHER_SECTION_EVIDENCE_CODE = "other_section"
ADVERSE_REACTIONS_RUBRIC_SECTIONS = ("AR", NON_US_LABEL_SECTION)  # sections whose text is the adverse-reactions text of the label
COMBINATION_LABEL_ONLY_REASON = "combination_label_only"
TARGET_RECORD_MISSING_REASON = "target_record_missing"
ONSIDES_EXTRA_COLUMNS: tuple[str, ...] = (
    "ingredient_id", "meddra_term_type", "label_section", "source_country", "label_count", "single_ingredient_label_count", "labels_of_ingredient_in_country",
    "max_pred1", "mean_pred1", "match_methods", "chembl_parent", "unified_with_sider", "qualifies", "disqualified_reason",
)


def load_onsides_tables(parquet_directory: Path) -> dict[str, pd.DataFrame]:
    names = ("product_adverse_effect", "product_label", "product_to_rxnorm", "vocab_rxnorm_ingredient_to_product", "vocab_rxnorm_ingredient", "vocab_meddra_adverse_effect")
    return {name: pd.read_parquet(parquet_directory / f"{name}.parquet") for name in names}


def sider_lower_level_to_preferred_terms(sider_directory: Path) -> dict[str, set[str]]:
    """lowercase LLT name -> lowercase PT names it falls under, from the LLT and PT rows SIDER lists per label term.

    meddra_all_se.tsv.gz carries, for every (drug, label term) pair, one LLT row and one PT row that share the
    UMLS label concept (column 3); pairing them gives a MedDRA LLT-to-PT table restricted to the terms SIDER
    saw, without redistributing MedDRA. Names are used in memory only and never written out.
    """
    lower_level_names_by_concept: dict[str, set[str]] = {}
    preferred_names_by_concept: dict[str, set[str]] = {}
    with gzip.open(sider_directory / "meddra_all_se.tsv.gz", "rt", encoding="utf-8") as side_effect_file:
        for row in csv.reader(side_effect_file, delimiter="\t"):
            label_concept, concept_type, term_name = row[2], row[3], row[5].strip().lower()
            (lower_level_names_by_concept if concept_type == "LLT" else preferred_names_by_concept).setdefault(label_concept, set()).add(term_name)
    pairing: dict[str, set[str]] = {}
    for label_concept, lower_level_names in lower_level_names_by_concept.items():
        for lower_level_name in lower_level_names:
            pairing.setdefault(lower_level_name, set()).update(preferred_names_by_concept.get(label_concept, set()))
    return {lower_level_name: preferred_names for lower_level_name, preferred_names in pairing.items() if preferred_names}


def meddra_terms_to_target_symptoms(vocabulary: pd.DataFrame, preferred_term_to_symptom: Mapping[str, str], lower_level_to_preferred_terms: Mapping[str, set[str]] | None = None) -> pd.DataFrame:
    """MedDRA ids matched to a target symptom: a PT row by its crosswalk name, an LLT row through its PT names.

    An LLT row is matched when lower_level_to_preferred_terms takes its name to a crosswalk preferred term (one
    row per target symptom reached); without the pairing LLT rows match only on an exact name, which in the
    v3.1.1 vocabulary never happens. meddra_term_type keeps the source term type of the row.
    """
    names = vocabulary.meddra_name.astype(str).str.strip().str.lower()
    is_lower_level = vocabulary.meddra_term_type.astype(str).str.upper() == "LLT"
    symptoms = names.map(preferred_term_to_symptom).where(~is_lower_level).map(lambda symptom: {symptom} if isinstance(symptom, str) else set())
    if lower_level_to_preferred_terms:
        symptoms = symptoms.where(~is_lower_level, names.map(lambda name: {preferred_term_to_symptom[preferred] for preferred in lower_level_to_preferred_terms.get(name, set()) if preferred in preferred_term_to_symptom}))
    else:
        symptoms = symptoms.where(~is_lower_level, names.map(lambda name: {preferred_term_to_symptom[name]} if name in preferred_term_to_symptom else set()))
    matched = vocabulary.assign(target_symptom=symptoms.map(sorted))[symptoms.map(len) > 0].explode("target_symptom")
    return matched[["meddra_id", "meddra_name", "meddra_term_type", "target_symptom"]].reset_index(drop=True)


def aggregate_label_statements(tables: Mapping[str, pd.DataFrame], preferred_term_to_symptom: Mapping[str, str], lower_level_to_preferred_terms: Mapping[str, set[str]] | None = None) -> pd.DataFrame:
    """One row per (ingredient, target symptom, MedDRA id, label section, country): the deduplicated label statement.

    Joins product_adverse_effect to the label's country, the label's RxNorm products and their ingredients;
    a combination product contributes its statement to each ingredient. label_count is the number of distinct
    labels carrying the statement, single_ingredient_label_count the number of those labels whose products
    hold exactly one ingredient (0 when every label behind the statement is a combination product) and
    labels_of_ingredient_in_country the ingredient's labels in that country.
    """
    terms = meddra_terms_to_target_symptoms(tables["vocab_meddra_adverse_effect"], preferred_term_to_symptom, lower_level_to_preferred_terms)
    events = tables["product_adverse_effect"][["product_label_id", "label_section", "effect_meddra_id", "match_method", "pred1"]]
    events = events.merge(terms, left_on="effect_meddra_id", right_on="meddra_id", how="inner")
    labels = tables["product_label"][["label_id", "source"]].rename(columns={"source": "source_country"})
    events = events.merge(labels, left_on="product_label_id", right_on="label_id", how="inner")
    label_products = tables["product_to_rxnorm"].rename(columns={"rxnorm_product_id": "product_id"})
    product_ingredients = tables["vocab_rxnorm_ingredient_to_product"]
    label_ingredients = label_products.merge(product_ingredients, on="product_id", how="inner")[["label_id", "ingredient_id"]].drop_duplicates()
    ingredients_per_label = label_ingredients.groupby("label_id").ingredient_id.nunique()
    label_ingredients["single_ingredient_label_id"] = label_ingredients.label_id.where(label_ingredients.label_id.map(ingredients_per_label) == 1)
    events = events.merge(label_ingredients, on="label_id", how="inner")
    ingredient_names = tables["vocab_rxnorm_ingredient"].set_index("rxnorm_id")["rxnorm_name"]
    labels_per_ingredient_country = label_ingredients.merge(labels, on="label_id").groupby(["ingredient_id", "source_country"]).label_id.nunique().rename("labels_of_ingredient_in_country")
    grouped = events.groupby(["ingredient_id", "target_symptom", "meddra_id", "meddra_name", "meddra_term_type", "label_section", "source_country"], as_index=False).agg(
        label_count=("label_id", "nunique"), single_ingredient_label_count=("single_ingredient_label_id", "nunique"), max_pred1=("pred1", "max"), mean_pred1=("pred1", "mean"),
        match_methods=("match_method", lambda values: ";".join(sorted(set(str(value) for value in values)))),
    )
    grouped = grouped.merge(labels_per_ingredient_country, on=["ingredient_id", "source_country"], how="left")
    grouped["ingredient_name"] = grouped.ingredient_id.map(ingredient_names).fillna("")
    return grouped.sort_values(["ingredient_id", "target_symptom", "meddra_id", "label_section", "source_country"]).reset_index(drop=True)


@dataclass
class OnsidesQualification:
    qualifies: bool
    reason: str
    drug_targets: list[DrugTarget]
    atc_known: bool
    missing_target_ids: list[str] = field(default_factory=list)  # mechanism targets without a record in targets.json


def qualify_ingredient(bridge: IngredientBridge, mechanisms_by_molecule: Mapping[str, list[dict]], targets: Mapping[str, dict],
                       compounds_by_chembl_id: Mapping[str, list[str]] | None = None) -> OnsidesQualification:
    """The E2 qualification of section 4.2 for one ingredient, the rule SIDER rows pass: a single mechanism target and an ATC N code.

    An ingredient with no mechanism target at all whose parent is listed in configs/drugs_acting_as_graph_compounds.csv
    (compounds_by_chembl_id) takes that graph compound as its one target, as a SIDER drug does.
    An ingredient with no ATC code at all fails with no_atc_code, as a SIDER drug without one does; atc_known
    is kept on the result for the rubric_atc_known feature.
    """
    missing_target_ids: set[str] = set()
    drug_targets = mechanism_targets_for_parent(bridge.chembl_parent, mechanisms_by_molecule, targets, missing_target_ids=missing_target_ids)
    drug_targets = graph_compound_targets(drug_targets, None, [bridge.chembl_parent] if bridge.chembl_parent else [], {}, compounds_by_chembl_id or {})
    atc_known = bool(bridge.atc_codes)
    missing = sorted(missing_target_ids)
    if bridge.chembl_parent is None:
        return OnsidesQualification(False, "no_chembl_parent", drug_targets, atc_known, missing)
    if not passes_single_target_rule(drug_targets):
        return OnsidesQualification(False, "no_single_mechanism_target", drug_targets, atc_known, missing)
    if not atc_known:
        return OnsidesQualification(False, "no_atc_code", drug_targets, atc_known, missing)
    if not is_nervous_system_atc(list(bridge.atc_codes)):
        return OnsidesQualification(False, "atc_not_nervous_system", drug_targets, atc_known, missing)
    return OnsidesQualification(True, "", drug_targets, atc_known, missing)


def section_limitation_clause(label_section: str) -> str:
    if label_section == "BW":
        return "boxed warning"
    if label_section == "WP":
        return "warnings and precautions section"
    if label_section == "AR":
        return "adverse reactions section"
    if label_section == NON_US_LABEL_SECTION:
        return "non-US label, section not recorded by OnSIDES"
    return "unsectioned mention"


def limitations_for_statement(statement: pd.Series, atc_known: bool, scaled_score: float) -> str:
    clauses = ["label statement without frequency", f"{int(statement.label_count)} label(s) in {statement.source_country}"]
    clauses.append(section_limitation_clause(str(statement.label_section)))
    if int(statement.single_ingredient_label_count) == 0:
        clauses.append("combination-product labels only")
    clauses.append("string match without model score" if scaled_score == 0.0 and "PMB" not in str(statement.match_methods) else f"extraction score {scaled_score:.2f} of 1")
    if not atc_known:
        clauses.append("no ATC code known, central action not tested")
    clauses.append(f"dated by the OnSIDES {ONSIDES_RELEASE} release ({ONSIDES_RELEASE_DATE}), an upper bound")
    return "; ".join(clauses)


def onsides_reports(
    statements: pd.DataFrame,
    bridges: Mapping[str, IngredientBridge],
    mechanisms_by_molecule: Mapping[str, list[dict]],
    targets: Mapping[str, dict],
    node_lookup: GraphNodeLookup | Mapping[str, str],
    compounds_by_chembl_id: Mapping[str, list[str]] | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Report rows (REPORT_COLUMNS first, then ONSIDES_EXTRA_COLUMNS) for every statement; qualifies marks the rows that enter the evidence table.

    node_lookup places each target on the graph (GraphNodeLookup); a gene symbol -> node id mapping places protein
    targets only."""
    if not isinstance(node_lookup, GraphNodeLookup):
        node_lookup = GraphNodeLookup(dict(node_lookup), {}, {})
    scored = statements.max_pred1.where(statements.match_methods.astype(str).str.contains("PMB"))
    score_low, score_high = float(scored.min()) if scored.notna().any() else 0.0, float(scored.max()) if scored.notna().any() else 1.0
    counts = {"statements": int(len(statements)), "ingredients_without_bridge": 0, "disqualified": {}, "qualifying_reports": 0}
    qualifications: dict[str, OnsidesQualification] = {}
    rows: list[dict] = []
    for statement in statements.itertuples(index=False):
        bridge = bridges.get(str(statement.ingredient_id))
        if bridge is None:
            counts["ingredients_without_bridge"] += 1
            continue
        qualification = qualifications.get(bridge.ingredient_id) or qualify_ingredient(bridge, mechanisms_by_molecule, targets, compounds_by_chembl_id)
        qualifications[bridge.ingredient_id] = qualification
        perturbation_nodes = node_lookup.perturbation_nodes(qualification.drug_targets)
        combination_only = int(statement.single_ingredient_label_count) == 0
        qualifies = qualification.qualifies and not combination_only and bool(perturbation_nodes)
        if qualification.reason:
            reason = qualification.reason
        elif combination_only:
            reason = COMBINATION_LABEL_ONLY_REASON
        elif perturbation_nodes:
            reason = ""
        else:
            reason = TARGET_RECORD_MISSING_REASON if qualification.missing_target_ids else "no_graph_node"
        if not qualifies:
            counts["disqualified"][reason] = counts["disqualified"].get(reason, 0) + 1
        is_model_scored = "PMB" in str(statement.match_methods) and not math.isnan(float(statement.max_pred1))
        scaled_score = (float(statement.max_pred1) - score_low) / (score_high - score_low) if is_model_scored and score_high > score_low else 0.0
        target_description = ";".join(sorted(f"{target.target_chembl_id}:{target.action_type}" for target in qualification.drug_targets))
        row = {column: "" for column in REPORT_COLUMNS}
        row.update({
            "report_id": f"{ONSIDES_SOURCE}|{bridge.perturbation_id}|{statement.ingredient_id}|{statement.meddra_id}|{statement.target_symptom}|{statement.label_section}|{statement.source_country}",
            "perturbation_id": bridge.perturbation_id, "perturbation_type": "drug", "perturbation_label": bridge.sider_drug_name or statement.ingredient_name,
            "symptom": statement.target_symptom, "relation": "induces", "evidence_class": "pharmacological", "source": ONSIDES_SOURCE, "report_value": 1,
            "source_record_id": str(statement.ingredient_id), "source_record_label": statement.ingredient_name,
            "source_term_id": str(statement.meddra_id), "source_term_label": statement.meddra_name,
            "evidence_code": SECTION_EVIDENCE_CODES.get(statement.label_section, OTHER_SECTION_EVIDENCE_CODE),
            "model_description": f"human; drug label (OnSIDES {ONSIDES_RELEASE}, {statement.source_country}); ChEMBL targets {target_description}",
            "frequency": None, "frequency_denominator": None, "placebo_flag": False, "onset": None, "sex": None, "references": "", "pubmed_reference_count": 0,
            "evidence_date": ONSIDES_RELEASE_DATE,
            "rubric_log_sample_size": 0.0, "rubric_frequency_known": 0.0, "rubric_evidence_code_pcs": 0.0, "rubric_evidence_code_tas": 0.0, "rubric_evidence_code_iea": 0.0,
            "rubric_placebo_controlled": 0.0, "rubric_curated_synopsis": 0.0, "association_type": "", "rubric_causal_association": 0.0,
            "limitations": limitations_for_statement(statement, qualification.atc_known, scaled_score),
            "perturbation_nodes": json.dumps(perturbation_nodes),
            "rubric_log_label_count": float(math.log1p(int(statement.label_count))),
            "rubric_boxed_warning": 1.0 if statement.label_section == "BW" else 0.0,
            "rubric_adverse_reactions_section": 1.0 if statement.label_section in ADVERSE_REACTIONS_RUBRIC_SECTIONS else 0.0,
            "rubric_nlp_score": scaled_score, "rubric_us_label": 1.0 if statement.source_country == "US" else 0.0,
            "rubric_atc_known": 1.0 if qualification.atc_known else 0.0,
            "ingredient_id": str(statement.ingredient_id), "meddra_term_type": statement.meddra_term_type, "label_section": statement.label_section,
            "source_country": statement.source_country, "label_count": int(statement.label_count), "single_ingredient_label_count": int(statement.single_ingredient_label_count),
            "labels_of_ingredient_in_country": int(statement.labels_of_ingredient_in_country) if not math.isnan(float(statement.labels_of_ingredient_in_country)) else 0,
            "max_pred1": float(statement.max_pred1), "mean_pred1": float(statement.mean_pred1), "match_methods": str(statement.match_methods),
            "chembl_parent": bridge.chembl_parent or "", "unified_with_sider": bool(bridge.unified_with_sider), "qualifies": qualifies, "disqualified_reason": reason,
        })
        rows.append(row)
        if qualifies:
            counts["qualifying_reports"] += 1
    columns = list(REPORT_COLUMNS) + ["rubric_log_label_count", "rubric_boxed_warning", "rubric_adverse_reactions_section", "rubric_nlp_score", "rubric_us_label", "rubric_atc_known"] + list(ONSIDES_EXTRA_COLUMNS)
    table = pd.DataFrame(rows, columns=columns)
    counts["qualifying_ingredients"] = int(table[table.qualifies].ingredient_id.nunique()) if len(table) else 0
    counts["qualifying_perturbations"] = int(table[table.qualifies].perturbation_id.nunique()) if len(table) else 0
    return table, counts
