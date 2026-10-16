"""Trial outcome measures as a therapeutic-direction source (decision 9.60, the owner, 16 October 2026).

SIDER states an indication as a disorder, so the table's `relieves` pairs pass through the existing nosology. A
posted results record on ClinicalTrials.gov states, per arm, the change on a named scale or subscale, which is a
drug-to-construct measurement: a placebo-controlled difference on the ADHD-RS hyperactivity/impulsivity subscale
says what the drug did to hyperactivity, whatever disorder the trial enrolled for.

What a report is. One outcome measure of one placebo-controlled study, for one study drug, mapped by
configs/trial_outcome_scale_crosswalk.csv (a regular expression on the outcome title, the construct it measures and
whether a higher score is worse). The drug arms and placebo arms are chosen by the same title rule as the adverse-event
arms (clinical_trial_adverse_events.select_arms: a group naming the drug and not placebo, against a group naming
placebo and not the drug; a group naming another study drug is left out). The comparison is the sponsor's own
between-group analysis where it names exactly those two groups and gives a mean difference with an interval; else the
difference of the two arms' means with a Welch interval from their dispersions and denominators. The difference is
signed so that a negative value is improvement on the drug, and standardised by the pooled standard deviation where
the dispersions allow it.

What the report says. relation `relieves`, report_value 1 when the interval excludes 0 on the improving side;
relation `relieves`, report_value 0 (an absence of relief) when the standardised interval lies within
NO_EFFECT_MARGIN of 0 on both sides, which needs the dispersions; relation `induces`, report_value 1 when the interval
excludes 0 on the worsening side. A comparison that decides none of these is undecided: it carries relation
`relieves` and report_value -1 in the comparisons table the fetch writes beside the report table, and never enters a
report table, since every reader of report_value takes 0 as an absence. Class human_association,
source ClinicalTrials-outcomes, grade C: a soft prior, never a label (9.49, 9.50).

What is not read. Trials without a placebo arm (an active comparator says nothing about the construct); disorder-level
composites (CGI, PANSS total, a total ADHD score that fuses inattention with hyperactivity) unless the crosswalk names
them; measures reported as medians or counts; a trial whose arms name two study drugs (9.48's combination rule).
"""
from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from mechanistic_pathway_learning.evidence.clinical_trial_adverse_events import name_pattern, study_descriptors
from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS

EVIDENCE_CLASS = "human_association"
SOURCE_NAME = "ClinicalTrials-outcomes"
EVIDENCE_CODE = "trial_outcome_comparison"
GRADE = "C"
NO_EFFECT_MARGIN = 0.2  # standardised mean difference inside which a drug is read as not moving the construct (Cohen's small effect)
Z_95 = 1.959964
DEFAULT_SCALE_CROSSWALK = Path("configs/trial_outcome_scale_crosswalk.csv")
OUTCOME_COLUMNS = ("nct_id", "outcome_title", "scale", "construct", "higher_is_worse", "param_type", "dispersion_type", "drug_group", "placebo_group",
                   "drug_n", "placebo_n", "drug_value", "placebo_value", "difference", "ci_lower", "ci_upper", "standardised_difference",
                   "standardised_ci_lower", "standardised_ci_upper", "comparison_source", "verdict")


@dataclass(frozen=True)
class ScaleRule:
    pattern: re.Pattern
    scale: str
    construct: str
    higher_is_worse: bool
    note: str


def read_scale_crosswalk(path: Path = DEFAULT_SCALE_CROSSWALK) -> list[ScaleRule]:
    """Rows: pattern (a regular expression matched against the outcome title, case-insensitive), scale, construct,
    higher_is_worse (true/false), note. The first matching row decides, so a subscale row precedes its total."""
    rules = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if not row.get("pattern", "").strip():
                continue
            rules.append(ScaleRule(re.compile(row["pattern"], re.IGNORECASE), row["scale"].strip(), row["construct"].strip(),
                                   row["higher_is_worse"].strip().lower() in ("true", "1", "yes"), row.get("note", "")))
    return rules


def rule_for(title: str, rules: list[ScaleRule]) -> ScaleRule | None:
    for rule in rules:
        if rule.pattern.search(title):
            return rule
    return None


@dataclass(frozen=True)
class GroupSelection:
    drug_groups: tuple[str, ...]
    placebo_groups: tuple[str, ...]
    excluded: tuple[str, ...]


def select_outcome_groups(outcome: dict, drug_name: str, other_study_drugs: Iterable[str] = ()) -> GroupSelection:
    """The adverse-event arm rule applied to an outcome measure's own groups (their ids differ from the event groups')."""
    pattern = name_pattern(drug_name)
    others = [name_pattern(other) for other in other_study_drugs if other.strip().lower() != drug_name.strip().lower()]
    drug_groups, placebo_groups, excluded = [], [], []
    for group in outcome.get("groups", []):
        text = (str(group.get("title", "")) + " " + str(group.get("description", ""))).lower()
        names_drug, names_placebo = bool(pattern.search(text)), "placebo" in text
        if any(other.search(text) for other in others) and (names_drug or names_placebo):
            excluded.append(group.get("id", ""))
            continue
        if names_drug and not names_placebo:
            drug_groups.append(group["id"])
        elif names_placebo and not names_drug:
            placebo_groups.append(group["id"])
    return GroupSelection(tuple(drug_groups), tuple(placebo_groups), tuple(excluded))


def to_float(value) -> float | None:
    try:
        number = float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def group_counts(outcome: dict) -> dict[str, int]:
    """Participants analysed per group, from the measure's denominators (the first 'Participants' denominator)."""
    for denominator in outcome.get("denoms", []):
        if str(denominator.get("units", "")).lower().startswith("participant"):
            return {count.get("groupId"): int(float(count.get("value"))) for count in denominator.get("counts", []) if to_float(count.get("value")) is not None}
    return {}


def first_measurements(outcome: dict) -> dict[str, dict]:
    """groupId -> the first class's first category's measurement (a single-row measure; a measure with several rows
    reports its first, which the title describes)."""
    for measure_class in outcome.get("classes", []):
        for category in measure_class.get("categories", []):
            measurements = {m.get("groupId"): m for m in category.get("measurements", []) if m.get("groupId")}
            if measurements:
                return measurements
    return {}


def standard_deviation(measurement: dict, dispersion_type: str, count: int | None) -> float | None:
    """A standard deviation from the measure's dispersion: itself, a standard error times root n, or a 95 percent
    interval's half-width times root n over z."""
    dispersion = (dispersion_type or "").lower()
    spread = to_float(measurement.get("spread"))
    if "standard deviation" in dispersion:
        return spread
    if "standard error" in dispersion:
        return spread * math.sqrt(count) if spread is not None and count else None
    lower, upper = to_float(measurement.get("lowerLimit")), to_float(measurement.get("upperLimit"))
    if "95" in dispersion and lower is not None and upper is not None and count:
        return (upper - lower) / (2 * Z_95) * math.sqrt(count)
    return None


def welch_interval(drug_mean: float, drug_sd: float, drug_n: int, placebo_mean: float, placebo_sd: float, placebo_n: int) -> tuple[float, float, float]:
    """(difference, lower, upper) of drug minus placebo with a normal 95 percent interval on Welch's standard error."""
    difference = drug_mean - placebo_mean
    standard_error = math.sqrt(drug_sd ** 2 / drug_n + placebo_sd ** 2 / placebo_n)
    return difference, difference - Z_95 * standard_error, difference + Z_95 * standard_error


def sponsor_comparison(outcome: dict, drug_group: str, placebo_group: str) -> tuple[float, float, float] | None:
    """The sponsor's mean-difference analysis over exactly these two groups, as drug minus placebo, when it has an interval."""
    for analysis in outcome.get("analyses", []):
        group_ids = list(analysis.get("groupIds", []))
        if sorted(group_ids) != sorted([drug_group, placebo_group]):
            continue
        parameter = str(analysis.get("paramType", "")).lower()
        if "difference" not in parameter or "median" in parameter:
            continue
        value, lower, upper = (to_float(analysis.get(key)) for key in ("paramValue", "ciLowerLimit", "ciUpperLimit"))
        if value is None or lower is None or upper is None:
            continue
        if str(analysis.get("ciPctValue", "95")).strip() not in ("95", "95.0"):
            continue
        sign = 1.0 if group_ids[0] == drug_group else -1.0  # the API states the difference in the order the groups are listed
        if sign < 0:
            value, lower, upper = -value, -upper, -lower
        return value, lower, upper
    return None


def verdict_of(standardised: tuple[float, float, float] | None, raw: tuple[float, float, float], higher_is_worse: bool) -> str:
    """relief, worsening, no_effect or undecided, with the interval signed so that negative is improvement."""
    orientation = 1.0 if higher_is_worse else -1.0
    lower, upper = sorted((orientation * raw[1], orientation * raw[2]))
    if upper < 0:
        return "relief"
    if lower > 0:
        return "worsening"
    if standardised is not None:
        s_lower, s_upper = sorted((orientation * standardised[1], orientation * standardised[2]))
        if -NO_EFFECT_MARGIN <= s_lower and s_upper <= NO_EFFECT_MARGIN:
            return "no_effect"
    return "undecided"


def compare_outcome(outcome: dict, selection: GroupSelection, rule: ScaleRule) -> dict | None:
    """One comparison row for an outcome measure, or None when no placebo-controlled mean comparison can be read."""
    if len(selection.drug_groups) != 1 or len(selection.placebo_groups) != 1:
        return None  # several dose arms or several placebo arms: which pair to read is a choice this module does not make
    drug_group, placebo_group = selection.drug_groups[0], selection.placebo_groups[0]
    parameter = str(outcome.get("paramType", "")).upper()
    if "MEAN" not in parameter:
        return None
    measurements, counts = first_measurements(outcome), group_counts(outcome)
    drug, placebo = measurements.get(drug_group), measurements.get(placebo_group)
    if drug is None or placebo is None:
        return None
    drug_value, placebo_value = to_float(drug.get("value")), to_float(placebo.get("value"))
    if drug_value is None or placebo_value is None:
        return None
    drug_n, placebo_n = counts.get(drug_group), counts.get(placebo_group)
    dispersion = str(outcome.get("dispersionType", ""))
    drug_sd, placebo_sd = standard_deviation(drug, dispersion, drug_n), standard_deviation(placebo, dispersion, placebo_n)
    pooled_sd = None
    if drug_sd is not None and placebo_sd is not None and drug_n and placebo_n and drug_n + placebo_n > 2 and (drug_sd > 0 or placebo_sd > 0):
        pooled_sd = math.sqrt(((drug_n - 1) * drug_sd ** 2 + (placebo_n - 1) * placebo_sd ** 2) / (drug_n + placebo_n - 2))
    sponsor = sponsor_comparison(outcome, drug_group, placebo_group)
    if sponsor is not None:
        raw, source = sponsor, "sponsor analysis"
    elif pooled_sd is not None:
        raw, source = welch_interval(drug_value, drug_sd, drug_n, placebo_value, placebo_sd, placebo_n), "arm means, Welch interval"
    else:
        return None
    standardised = tuple(value / pooled_sd for value in raw) if pooled_sd else None
    verdict = verdict_of(standardised, raw, rule.higher_is_worse)
    return {"outcome_title": outcome.get("title", ""), "scale": rule.scale, "construct": rule.construct, "higher_is_worse": rule.higher_is_worse,
            "param_type": parameter, "dispersion_type": dispersion, "drug_group": drug_group, "placebo_group": placebo_group,
            "drug_n": drug_n, "placebo_n": placebo_n, "drug_value": drug_value, "placebo_value": placebo_value,
            "difference": raw[0], "ci_lower": raw[1], "ci_upper": raw[2],
            "standardised_difference": standardised[0] if standardised else None, "standardised_ci_lower": standardised[1] if standardised else None,
            "standardised_ci_upper": standardised[2] if standardised else None, "comparison_source": source, "verdict": verdict}


def build_reports(studies: Iterable[dict], drug_name: str, perturbation: dict, rules: list[ScaleRule], other_study_drugs: Iterable[str] = ()) -> tuple[list[dict], dict]:
    """Every readable placebo-controlled comparison of the drug's studies on a crosswalked scale, as report rows, with counts."""
    rows, counts = [], {"studies": 0, "outcome_measures": 0, "crosswalked": 0, "placebo_controlled": 0, "compared": 0,
                        "relief": 0, "worsening": 0, "no_effect": 0, "undecided": 0}
    others = list(other_study_drugs)
    for study in studies:
        counts["studies"] += 1
        nct_id = study.get("protocolSection", {}).get("identificationModule", {}).get("nctId", "")
        title = study.get("protocolSection", {}).get("identificationModule", {}).get("briefTitle", "")
        for outcome in study.get("resultsSection", {}).get("outcomeMeasuresModule", {}).get("outcomeMeasures", []):
            counts["outcome_measures"] += 1
            rule = rule_for(str(outcome.get("title", "")), rules)
            if rule is None:
                continue
            counts["crosswalked"] += 1
            selection = select_outcome_groups(outcome, drug_name, others)
            if not selection.drug_groups or not selection.placebo_groups:
                continue
            counts["placebo_controlled"] += 1
            comparison = compare_outcome(outcome, selection, rule)
            if comparison is None:
                continue
            counts["compared"] += 1
            counts[comparison["verdict"]] += 1
            relation, value = {"relief": ("relieves", 1), "worsening": ("induces", 1), "no_effect": ("relieves", 0), "undecided": ("relieves", -1)}[comparison["verdict"]]
            descriptors = study_descriptors(study)
            row = {column: "" for column in REPORT_COLUMNS}
            row.update({
                "report_id": f"{SOURCE_NAME}|{perturbation['perturbation_id']}|{nct_id}|{comparison['scale']}|{comparison['construct']}|{relation}|{value}",
                "perturbation_id": perturbation["perturbation_id"], "perturbation_type": "drug", "perturbation_label": perturbation["perturbation_label"],
                "symptom": comparison["construct"], "relation": relation, "evidence_class": EVIDENCE_CLASS, "source": SOURCE_NAME, "report_value": value,
                "source_record_id": nct_id, "source_record_label": title, "source_term_id": comparison["scale"], "source_term_label": comparison["outcome_title"],
                "evidence_code": EVIDENCE_CODE,
                "model_description": (f"human; placebo-controlled trial {nct_id}; {comparison['scale']} ({'higher is worse' if comparison['higher_is_worse'] else 'higher is better'}); "
                                      f"drug minus placebo {comparison['difference']:.4g} [{comparison['ci_lower']:.4g}, {comparison['ci_upper']:.4g}] by {comparison['comparison_source']}; "
                                      f"n {comparison['drug_n']} against {comparison['placebo_n']}; verdict {comparison['verdict']}; ChEMBL targets {perturbation.get('targets', '')}"),
                "frequency": None, "frequency_denominator": comparison["drug_n"], "placebo_flag": True, "onset": None, "sex": None,
                "references": f"NCT:{nct_id}", "pubmed_reference_count": 0,
                "evidence_date": descriptors["results_first_posted"] or None, "evidence_date_source": "release" if descriptors["results_first_posted"] else "",
                "rubric_log_sample_size": math.log1p(comparison["drug_n"] or 0), "rubric_frequency_known": 0.0,
                "rubric_evidence_code_pcs": 0.0, "rubric_evidence_code_tas": 0.0, "rubric_evidence_code_iea": 0.0,
                "rubric_placebo_controlled": 1.0, "rubric_curated_synopsis": 0.0, "association_type": "", "rubric_causal_association": 0.0,
                "limitations": (f"the outcome measure's first reported row is read; allocation {descriptors['allocation'] or 'unstated'}, masking {descriptors['masking'] or 'unstated'}, "
                                f"status {descriptors['overall_status'] or 'unstated'}; comparison by {comparison['comparison_source']}"),
                "perturbation_nodes": perturbation.get("perturbation_nodes", ""),
                "nct_id": nct_id, **{key: comparison[key] for key in OUTCOME_COLUMNS if key in comparison},
            })
            rows.append(row)
    return rows, counts


def rows_to_dataframe(rows: list[dict]):
    """Report rows as a table with the report columns, then the outcome columns, with stable dtypes; empty-safe."""
    import pandas as pd

    table = pd.DataFrame(rows, columns=[*REPORT_COLUMNS, *OUTCOME_COLUMNS])
    table["report_value"] = pd.to_numeric(table.report_value, errors="coerce").fillna(0).astype("int64")
    table["pubmed_reference_count"] = pd.to_numeric(table.pubmed_reference_count, errors="coerce").fillna(0).astype("int64")
    table["frequency"] = pd.to_numeric(table.frequency, errors="coerce").astype("float64")
    for column in ("frequency_denominator", "drug_n", "placebo_n"):
        table[column] = pd.array([None if value in ("", None) or (isinstance(value, float) and math.isnan(value)) else int(value) for value in table[column]], dtype="Int64")
    table["placebo_flag"] = table.placebo_flag.map(lambda value: bool(value) if value != "" else False).astype(bool)
    table["higher_is_worse"] = table.higher_is_worse.map(lambda value: bool(value) if value != "" else False).astype(bool)
    for column in ("drug_value", "placebo_value", "difference", "ci_lower", "ci_upper", "standardised_difference", "standardised_ci_lower", "standardised_ci_upper"):
        table[column] = pd.to_numeric(table[column], errors="coerce").astype("float64")
    for column in ("rubric_log_sample_size", "rubric_frequency_known", "rubric_evidence_code_pcs", "rubric_evidence_code_tas", "rubric_evidence_code_iea",
                   "rubric_placebo_controlled", "rubric_curated_synopsis", "rubric_causal_association"):
        table[column] = pd.to_numeric(table[column], errors="coerce").fillna(0.0).astype("float64")
    for column in ("onset", "sex", "evidence_date"):
        table[column] = table[column].map(lambda value: None if value in ("", None) else value).astype(object)
    return table
