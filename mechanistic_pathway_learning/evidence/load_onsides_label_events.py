"""Evidence class E2, later label slice: OnSIDES v3.1.1 label events as report-level rows (design section 4.2).

OnSIDES (Tanaka et al. 2025, reference [29]) extracts adverse-event terms from the Adverse Reactions, Boxed
Warning and Warnings and Precautions sections of current drug labels (US, UK, EU and Japan) with a
PubMedBERT model; its v3.1.1 data release (2026-04-22, CC BY 4.0) has 6.9 million label-event mentions
over 41,119 labels and 1,866 RxNorm ingredients. Generic products repeat the same label text, so the
unit of evidence here is the deduplicated label statement: one report per (ingredient, target symptom,
MedDRA term, label section, label country) with the number of distinct labels as a feature. OnSIDES
carries no frequencies and no label dates; the data release date is an upper bound on when a statement
appeared, which the limitations text says.

Drugs are keyed by RxNorm ingredient; the identifier bridge (onsides_identifier_bridge.py) takes each
ingredient to its ChEMBL parent, its mechanism targets and, when the parent is also a SIDER drug's,
to that drug's PubChem identifier so SIDER and OnSIDES report on one perturbation. The qualification
rule is the one SIDER rows pass: a single mechanism target and the ATC N proxy for central action,
with the ATC test waived (and flagged) when no ATC code is known for the ingredient.

Rows are grade B label evidence like SIDER rows, never literature. MedDRA term names stay inside
data/processed and are blanked by the data release.
"""
from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS
from mechanistic_pathway_learning.evidence.load_drug_label_events import is_nervous_system_atc
from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import IngredientBridge, mechanism_targets_for_parent, passes_single_target_rule
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import DrugTarget

ONSIDES_RELEASE = "v3.1.1"
ONSIDES_RELEASE_DATE = "2026-04-22"
ONSIDES_SOURCE = "OnSIDES-label"
SECTION_EVIDENCE_CODES = {"AR": "adverse_reactions_section", "BW": "boxed_warning_section", "WP": "warnings_precautions_section"}
OTHER_SECTION_EVIDENCE_CODE = "other_section"
ONSIDES_EXTRA_COLUMNS: tuple[str, ...] = (
    "ingredient_id", "meddra_term_type", "label_section", "source_country", "label_count", "labels_of_ingredient_in_country",
    "max_pred1", "mean_pred1", "match_methods", "chembl_parent", "unified_with_sider", "qualifies", "disqualified_reason",
)


def load_onsides_tables(parquet_directory: Path) -> dict[str, pd.DataFrame]:
    names = ("product_adverse_effect", "product_label", "product_to_rxnorm", "vocab_rxnorm_ingredient_to_product", "vocab_rxnorm_ingredient", "vocab_meddra_adverse_effect")
    return {name: pd.read_parquet(parquet_directory / f"{name}.parquet") for name in names}


def meddra_terms_to_target_symptoms(vocabulary: pd.DataFrame, preferred_term_to_symptom: Mapping[str, str]) -> pd.DataFrame:
    """MedDRA ids whose name is a crosswalk preferred term; PT rows always, LLT rows only on an exact name match."""
    names = vocabulary.meddra_name.astype(str).str.strip().str.lower()
    symptoms = names.map(preferred_term_to_symptom)
    matched = vocabulary.assign(target_symptom=symptoms)[symptoms.notna()]
    return matched[["meddra_id", "meddra_name", "meddra_term_type", "target_symptom"]]


def aggregate_label_statements(tables: Mapping[str, pd.DataFrame], preferred_term_to_symptom: Mapping[str, str]) -> pd.DataFrame:
    """One row per (ingredient, target symptom, MedDRA id, label section, country): the deduplicated label statement.

    Joins product_adverse_effect to the label's country, the label's RxNorm products and their ingredients;
    a combination product contributes its statement to each ingredient. label_count is the number of distinct
    labels carrying the statement, labels_of_ingredient_in_country the ingredient's labels in that country.
    """
    terms = meddra_terms_to_target_symptoms(tables["vocab_meddra_adverse_effect"], preferred_term_to_symptom)
    events = tables["product_adverse_effect"][["product_label_id", "label_section", "effect_meddra_id", "match_method", "pred1"]]
    events = events.merge(terms, left_on="effect_meddra_id", right_on="meddra_id", how="inner")
    labels = tables["product_label"][["label_id", "source"]].rename(columns={"source": "source_country"})
    events = events.merge(labels, left_on="product_label_id", right_on="label_id", how="inner")
    label_products = tables["product_to_rxnorm"].rename(columns={"rxnorm_product_id": "product_id"})
    product_ingredients = tables["vocab_rxnorm_ingredient_to_product"]
    label_ingredients = label_products.merge(product_ingredients, on="product_id", how="inner")[["label_id", "ingredient_id"]].drop_duplicates()
    events = events.merge(label_ingredients, on="label_id", how="inner")
    ingredient_names = tables["vocab_rxnorm_ingredient"].set_index("rxnorm_id")["rxnorm_name"]
    labels_per_ingredient_country = label_ingredients.merge(labels, on="label_id").groupby(["ingredient_id", "source_country"]).label_id.nunique().rename("labels_of_ingredient_in_country")
    grouped = events.groupby(["ingredient_id", "target_symptom", "meddra_id", "meddra_name", "meddra_term_type", "label_section", "source_country"], as_index=False).agg(
        label_count=("label_id", "nunique"), max_pred1=("pred1", "max"), mean_pred1=("pred1", "mean"),
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


def qualify_ingredient(bridge: IngredientBridge, mechanisms_by_molecule: Mapping[str, list[dict]], targets: Mapping[str, dict]) -> OnsidesQualification:
    """The E2 qualification of section 4.2 for one ingredient: a single mechanism target and the ATC N proxy (waived when no ATC code is known)."""
    drug_targets = mechanism_targets_for_parent(bridge.chembl_parent, mechanisms_by_molecule, targets)
    if bridge.chembl_parent is None:
        return OnsidesQualification(False, "no_chembl_parent", drug_targets, bool(bridge.atc_codes))
    if not passes_single_target_rule(drug_targets):
        return OnsidesQualification(False, "no_single_mechanism_target", drug_targets, bool(bridge.atc_codes))
    atc_known = bool(bridge.atc_codes)
    if atc_known and not is_nervous_system_atc(list(bridge.atc_codes)):
        return OnsidesQualification(False, "atc_not_nervous_system", drug_targets, atc_known)
    return OnsidesQualification(True, "", drug_targets, atc_known)


def limitations_for_statement(statement: pd.Series, atc_known: bool, scaled_score: float) -> str:
    clauses = ["label statement without frequency", f"{int(statement.label_count)} label(s) in {statement.source_country}"]
    clauses.append("boxed warning" if statement.label_section == "BW" else "warnings and precautions section" if statement.label_section == "WP" else "adverse reactions section" if statement.label_section == "AR" else "unsectioned mention")
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
    node_by_symbol: Mapping[str, str],
) -> tuple[pd.DataFrame, dict]:
    """Report rows (REPORT_COLUMNS first, then ONSIDES_EXTRA_COLUMNS) for every statement; qualifies marks the rows that enter the evidence table."""
    scored = statements.max_pred1.where(statements.match_methods.astype(str).str.contains("PMB"))
    score_low, score_high = float(scored.min()) if scored.notna().any() else 0.0, float(scored.max()) if scored.notna().any() else 1.0
    counts = {"statements": int(len(statements)), "ingredients_without_bridge": 0, "disqualified": {}, "no_graph_node": 0, "qualifying_reports": 0}
    qualifications: dict[str, OnsidesQualification] = {}
    rows: list[dict] = []
    for statement in statements.itertuples(index=False):
        bridge = bridges.get(str(statement.ingredient_id))
        if bridge is None:
            counts["ingredients_without_bridge"] += 1
            continue
        qualification = qualifications.get(bridge.ingredient_id) or qualify_ingredient(bridge, mechanisms_by_molecule, targets)
        qualifications[bridge.ingredient_id] = qualification
        perturbation_nodes = []
        for target in qualification.drug_targets:
            mapped = [node_by_symbol[symbol] for symbol in target.gene_symbols if symbol in node_by_symbol]
            perturbation_nodes.extend([node, target.sign, 1.0 / len(mapped)] for node in mapped)
        qualifies = qualification.qualifies and bool(perturbation_nodes)
        reason = qualification.reason or ("no_graph_node" if not perturbation_nodes else "")
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
            "rubric_adverse_reactions_section": 1.0 if statement.label_section == "AR" else 0.0,
            "rubric_nlp_score": scaled_score, "rubric_us_label": 1.0 if statement.source_country == "US" else 0.0,
            "rubric_atc_known": 1.0 if qualification.atc_known else 0.0,
            "ingredient_id": str(statement.ingredient_id), "meddra_term_type": statement.meddra_term_type, "label_section": statement.label_section,
            "source_country": statement.source_country, "label_count": int(statement.label_count),
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
