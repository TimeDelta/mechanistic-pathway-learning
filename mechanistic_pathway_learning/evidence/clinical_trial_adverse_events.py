"""Adverse-event tables of placebo-controlled trials on ClinicalTrials.gov as evidence reports, positive and negative.

Why this source. Every other label-derived source here records presences only: a drug label lists an event, SIDER
counts it, OnSIDES extracts it, and silence is silence. A trial's results section is different in kind: it reports,
for each event term above the frequency threshold, the number affected and the number at risk in every arm. A drug
arm whose rate is no higher than the placebo arm's is a measured absence with a known exposure, which is the
negative evidence the reliability fit cannot estimate specificity without (docs/evidence_reports_spec.md section 3.4;
the owner, 15 October 2026: "it doesn't have to be positive only if we find a way to incorporate negatives").

What a report is. One (drug, symptom, trial, term, seriousness) row: the drug arms of the trial pooled against its
placebo arms, the term mapped to a target symptom through the crosswalk's MedDRA preferred terms (and SIDER's
lowest-level-term to preferred-term table, since trials code at either level), report_value 1 when the drug-arm rate
exceeds the placebo-arm rate and 0 when it does not and the drug arm's exposure reaches the floor. A term the trial
does not list is not a report of anything: the results section lists only terms above its frequency threshold
(`frequencyThreshold`, usually 5 percent), so an unlisted term is below threshold in every arm, not absent.

Arms and the leakage groups. The held-out split groups drugs by mechanism target, so an arm that names two study
drugs (a combination, or one drug on top of another as background) would carry one arm's events into two drugs'
records, one of which may be held out while the other trains. Such an arm is excluded from both sides
(`other_study_drugs` in select_arms), and the exclusion is counted. Matching is by the drug's SIDER label name as a
whole word; a trial that names the drug by a synonym is not found, which under-fetches and never contaminates.

What this is not. Not a causal claim per trial: arms are pooled across doses, terms are the sponsor's MedDRA coding,
the threshold hides rare events, and a terminated or open-label trial is weaker than a completed double-blind one.
Those facts travel in the limitations text and the rubric features; the tier and the reliability fit decide the
weight.

Licence: ClinicalTrials.gov is run by the U.S. National Library of Medicine; "Works produced by the U.S. government
are not subject to copyright protection in the United States" and "may be freely used or reproduced without
permission in the U.S." (NLM web policies, https://www.nlm.nih.gov/web_policies.html, read 15 October 2026). The
site's own terms page is rendered by script and could not be read from the shell; docs/data_sources.md records that.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from mechanistic_pathway_learning.evidence.comparison_statistics import STATISTIC_COLUMNS, empty_statistic, log_risk_ratio
from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS

SOURCE_NAME = "ClinicalTrials-arms"
# the arm counts and the risk difference with its 95 percent Newcombe interval, beside the REPORT_COLUMNS
ARM_COLUMNS = ("drug_affected", "drug_at_risk", "placebo_affected", "placebo_at_risk", "risk_difference", "risk_difference_lower_95", "risk_difference_upper_95")
EVIDENCE_CODE = "trial_arm_comparison"
# the design's grade-C class (assign_evidence_grades): never a label by default, enters the reliability fit and the priors
EVIDENCE_CLASS = "human_association"
DEFAULT_MINIMUM_EXPOSURE_FOR_ABSENCE = 100  # drug-arm participants at risk before a non-excess counts as a measured absence
DEFAULT_FREQUENCY_THRESHOLD = 5.0


@dataclass(frozen=True)
class ArmSelection:
    nct_id: str
    drug_group_ids: tuple[str, ...]
    placebo_group_ids: tuple[str, ...]
    drug_group_titles: tuple[str, ...]
    placebo_group_titles: tuple[str, ...]
    excluded_for_another_study_drug: tuple[str, ...] = ()  # titles of arms naming a second study drug, left out of both sides

    @property
    def usable(self) -> bool:
        return bool(self.drug_group_ids) and bool(self.placebo_group_ids)


@dataclass(frozen=True)
class TermComparison:
    nct_id: str
    term: str
    organ_system: str
    seriousness: str  # "other" or "serious"
    drug_affected: int
    drug_at_risk: int
    placebo_affected: int
    placebo_at_risk: int
    source_vocabulary: str

    @property
    def drug_rate(self) -> float | None:
        return self.drug_affected / self.drug_at_risk if self.drug_at_risk else None

    @property
    def placebo_rate(self) -> float | None:
        return self.placebo_affected / self.placebo_at_risk if self.placebo_at_risk else None

    @property
    def risk_difference_interval(self) -> tuple[float, float, float] | None:
        """(difference, lower, upper): the drug-arm rate minus the placebo-arm rate with its 95 percent Newcombe interval."""
        if not self.drug_at_risk or not self.placebo_at_risk:
            return None
        return newcombe_risk_difference_interval(self.drug_affected, self.drug_at_risk, self.placebo_affected, self.placebo_at_risk)


def wilson_interval(affected: int, at_risk: int, z: float = 1.959964) -> tuple[float, float]:
    """The Wilson score interval of a proportion; unlike the Wald interval it has width at 0 of n."""
    proportion = affected / at_risk
    centre = (proportion + z * z / (2 * at_risk)) / (1 + z * z / at_risk)
    half_width = z * math.sqrt(proportion * (1 - proportion) / at_risk + z * z / (4 * at_risk * at_risk)) / (1 + z * z / at_risk)
    return max(0.0, centre - half_width), min(1.0, centre + half_width)


def newcombe_risk_difference_interval(drug_affected: int, drug_at_risk: int, placebo_affected: int, placebo_at_risk: int) -> tuple[float, float, float]:
    """The difference of two proportions with Newcombe's hybrid score interval (Newcombe 1998, *Stat Med* 17:873, method
    10): the Wilson limits of each proportion combined, which keeps a width at zero events where the Wald interval has
    none. Returns (difference, lower, upper)."""
    drug_rate, placebo_rate = drug_affected / drug_at_risk, placebo_affected / placebo_at_risk
    drug_lower, drug_upper = wilson_interval(drug_affected, drug_at_risk)
    placebo_lower, placebo_upper = wilson_interval(placebo_affected, placebo_at_risk)
    difference = drug_rate - placebo_rate
    lower = difference - math.sqrt((drug_rate - drug_lower) ** 2 + (placebo_upper - placebo_rate) ** 2)
    upper = difference + math.sqrt((drug_upper - drug_rate) ** 2 + (placebo_rate - placebo_lower) ** 2)
    return difference, lower, upper


def name_pattern(drug_name: str) -> re.Pattern:
    """The drug name as a whole word, so "doxepin" does not match "doxepinone" and "lu" does not match "flu"."""
    return re.compile(r"(?<![a-z])" + re.escape(drug_name.strip().lower()) + r"(?![a-z])")


def select_arms(study: dict, drug_name: str, other_study_drugs: Iterable[str] = ()) -> ArmSelection:
    """Event groups whose title names the drug and not placebo, against groups that name placebo and not the drug.

    A group naming both (placebo added to the drug as background in every arm) belongs to neither side, which is
    what keeps an adjunctive trial of another drug from being read as a trial of the background drug. A group naming
    any other study drug is left out of both sides, so that one arm's events never reach two drugs' records across
    the held-out split.
    """
    nct_id = study.get("protocolSection", {}).get("identificationModule", {}).get("nctId", "")
    groups = study.get("resultsSection", {}).get("adverseEventsModule", {}).get("eventGroups", [])
    pattern = name_pattern(drug_name)
    other_patterns = [name_pattern(other) for other in other_study_drugs if other.strip().lower() != drug_name.strip().lower()]
    drug_groups, placebo_groups, excluded = [], [], []
    for group in groups:
        title = str(group.get("title", "")) + " " + str(group.get("description", ""))
        lowered = title.lower()
        names_drug = bool(pattern.search(lowered))
        names_placebo = "placebo" in lowered
        names_another = any(other.search(lowered) for other in other_patterns)
        if names_another and (names_drug or names_placebo):
            excluded.append(group.get("title", ""))
            continue
        if names_drug and not names_placebo:
            drug_groups.append((group["id"], group.get("title", "")))
        elif names_placebo and not names_drug:
            placebo_groups.append((group["id"], group.get("title", "")))
    return ArmSelection(nct_id, tuple(identifier for identifier, _ in drug_groups), tuple(identifier for identifier, _ in placebo_groups),
                        tuple(title for _, title in drug_groups), tuple(title for _, title in placebo_groups), tuple(excluded))


def term_comparisons(study: dict, selection: ArmSelection) -> list[TermComparison]:
    """Every listed term, with the drug arms pooled and the placebo arms pooled."""
    module = study.get("resultsSection", {}).get("adverseEventsModule", {})
    comparisons = []
    for seriousness, key in (("other", "otherEvents"), ("serious", "seriousEvents")):
        for event in module.get(key, []) or []:
            drug_affected = drug_at_risk = placebo_affected = placebo_at_risk = 0
            seen_drug = seen_placebo = False
            for stat in event.get("stats", []) or []:
                group_id = stat.get("groupId")
                affected = int(stat.get("numAffected") or 0)
                at_risk = int(stat.get("numAtRisk") or 0)
                if group_id in selection.drug_group_ids:
                    drug_affected += affected
                    drug_at_risk += at_risk
                    seen_drug = True
                elif group_id in selection.placebo_group_ids:
                    placebo_affected += affected
                    placebo_at_risk += at_risk
                    seen_placebo = True
            if seen_drug and seen_placebo and drug_at_risk > 0 and placebo_at_risk > 0:
                comparisons.append(TermComparison(selection.nct_id, str(event.get("term", "")).strip(), str(event.get("organSystem", "")),
                                                  seriousness, drug_affected, drug_at_risk, placebo_affected, placebo_at_risk,
                                                  str(event.get("sourceVocabulary", ""))))
    return comparisons


def read_lowest_level_to_preferred_term(sider_directory: Path) -> dict[str, str]:
    """SIDER's meddra_all_se.tsv names each term at the level it was found (LLT or PT) with the preferred-term concept
    beside it, so a lowest-level term name maps to its preferred term's name."""
    preferred_name_by_concept: dict[str, str] = {}
    preferred_concept_by_lowest_name: dict[str, str] = {}
    path = Path(sider_directory) / "meddra_all_se.tsv.gz"
    if not path.exists():
        return {}
    with gzip.open(path, "rt", encoding="utf-8") as se_file:
        for row in csv.reader(se_file, delimiter="\t"):
            if len(row) < 6:
                continue
            concept_type, preferred_concept, name = row[3], row[4], row[5].strip().lower()
            if concept_type == "PT":
                preferred_name_by_concept[preferred_concept] = name
            elif concept_type == "LLT":
                preferred_concept_by_lowest_name[name] = preferred_concept
    return {lowest: preferred_name_by_concept[concept] for lowest, concept in preferred_concept_by_lowest_name.items()
            if concept in preferred_name_by_concept}


def symptom_of_term(term: str, preferred_term_to_symptom: dict[str, str], lowest_to_preferred: dict[str, str]) -> str | None:
    lowered = term.strip().lower()
    if lowered in preferred_term_to_symptom:
        return preferred_term_to_symptom[lowered]
    preferred = lowest_to_preferred.get(lowered)
    return preferred_term_to_symptom.get(preferred) if preferred else None


def study_descriptors(study: dict) -> dict:
    protocol = study.get("protocolSection", {})
    design = protocol.get("designModule", {}).get("designInfo", {})
    status = protocol.get("statusModule", {})
    module = study.get("resultsSection", {}).get("adverseEventsModule", {})
    return {
        "overall_status": status.get("overallStatus", ""),
        "results_first_posted": (status.get("resultsFirstPostDateStruct") or {}).get("date", ""),
        "allocation": design.get("allocation", ""),
        "masking": (design.get("maskingInfo") or {}).get("masking", ""),
        "frequency_threshold": float(module.get("frequencyThreshold") or DEFAULT_FREQUENCY_THRESHOLD),
        "time_frame": module.get("timeFrame", ""),
    }


def build_reports(studies: Iterable[dict], drug_name: str, perturbation: dict, preferred_term_to_symptom: dict[str, str],
                  lowest_to_preferred: dict[str, str], minimum_exposure_for_absence: int = DEFAULT_MINIMUM_EXPOSURE_FOR_ABSENCE,
                  other_study_drugs: Iterable[str] = ()) -> tuple[list[dict], dict]:
    """Report rows (REPORT_COLUMNS) for one drug across its studies, and counts of what was kept and why not.

    `perturbation` carries the identity the assembled table already uses for the drug: perturbation_id (the STITCH
    flat id), perturbation_label, perturbation_nodes (JSON) and the ChEMBL target description, so these reports join
    the drug's existing reports rather than admitting it afresh.
    """
    other_study_drugs = tuple(other_study_drugs)
    rows: list[dict] = []
    counts = {"studies": 0, "studies_with_usable_arms": 0, "arms_excluded_for_another_study_drug": 0, "terms": 0, "terms_on_a_symptom": 0,
              "reports_positive": 0, "reports_negative": 0, "non_excess_below_exposure_floor": 0}
    for study in studies:
        counts["studies"] += 1
        selection = select_arms(study, drug_name, other_study_drugs)
        counts["arms_excluded_for_another_study_drug"] += len(selection.excluded_for_another_study_drug)
        if not selection.usable:
            continue
        counts["studies_with_usable_arms"] += 1
        descriptors = study_descriptors(study)
        for ordinal, comparison in enumerate(term_comparisons(study, selection)):
            counts["terms"] += 1
            symptom = symptom_of_term(comparison.term, preferred_term_to_symptom, lowest_to_preferred)
            if symptom is None:
                continue
            counts["terms_on_a_symptom"] += 1
            excess = comparison.drug_rate > comparison.placebo_rate
            if not excess and comparison.drug_at_risk < minimum_exposure_for_absence:
                counts["non_excess_below_exposure_floor"] += 1
                continue
            report_value = 1 if excess else 0
            counts["reports_positive" if excess else "reports_negative"] += 1
            limitations = (f"frequency threshold {descriptors['frequency_threshold']:g} percent, so unlisted terms are below it in every arm; "
                           f"arms pooled across doses: drug {', '.join(selection.drug_group_titles)} against placebo {', '.join(selection.placebo_group_titles)}; "
                           f"allocation {descriptors['allocation'] or 'unstated'}, masking {descriptors['masking'] or 'unstated'}, status {descriptors['overall_status'] or 'unstated'}; "
                           f"{comparison.seriousness} events; sponsor coding ({comparison.source_vocabulary or 'vocabulary unstated'})")
            row = {column: "" for column in REPORT_COLUMNS}
            row.update({
                "report_id": f"{SOURCE_NAME}|{perturbation['perturbation_id']}|{comparison.nct_id}|{symptom}|{comparison.seriousness}|{ordinal}",
                "perturbation_id": perturbation["perturbation_id"], "perturbation_type": "drug", "perturbation_label": perturbation["perturbation_label"],
                "symptom": symptom, "relation": "induces", "evidence_class": EVIDENCE_CLASS, "source": SOURCE_NAME, "report_value": report_value,
                "source_record_id": comparison.nct_id, "source_record_label": study.get("protocolSection", {}).get("identificationModule", {}).get("briefTitle", ""),
                "source_term_id": "", "source_term_label": comparison.term, "evidence_code": EVIDENCE_CODE,
                "model_description": (f"human; placebo-controlled trial {comparison.nct_id}; drug arm {comparison.drug_affected}/{comparison.drug_at_risk} "
                                      f"({comparison.drug_rate:.4f}) against placebo {comparison.placebo_affected}/{comparison.placebo_at_risk} ({comparison.placebo_rate:.4f}); "
                                      f"ChEMBL targets {perturbation.get('targets', '')}"),
                "frequency": comparison.drug_rate, "frequency_denominator": comparison.drug_at_risk, "placebo_flag": True, "onset": None, "sex": None,
                "references": f"NCT:{comparison.nct_id}", "pubmed_reference_count": 0,
                "evidence_date": descriptors["results_first_posted"] or None, "evidence_date_source": "release" if descriptors["results_first_posted"] else "",
                "rubric_log_sample_size": math.log1p(comparison.drug_at_risk), "rubric_frequency_known": 1.0,
                "rubric_evidence_code_pcs": 0.0, "rubric_evidence_code_tas": 0.0, "rubric_evidence_code_iea": 0.0,
                "rubric_placebo_controlled": 1.0, "rubric_curated_synopsis": 0.0, "association_type": "", "rubric_causal_association": 0.0,
                "limitations": limitations, "perturbation_nodes": perturbation.get("perturbation_nodes", ""),
                # the arm counts and the risk difference with its interval, structured, so a reader of the table (the soft
                # priors) can set its own bar for an excess or an absence instead of taking any excess at its word
                "drug_affected": comparison.drug_affected, "drug_at_risk": comparison.drug_at_risk,
                "placebo_affected": comparison.placebo_affected, "placebo_at_risk": comparison.placebo_at_risk,
            })
            interval = comparison.risk_difference_interval
            if interval is not None:
                row["risk_difference"], row["risk_difference_lower_95"], row["risk_difference_upper_95"] = interval
            row.update(comparison_statistic(comparison))
            rows.append(row)
    return rows, counts


def comparison_statistic(comparison: TermComparison) -> dict:
    """The arm comparison's own statistic (decision 9.48, option 3): the log risk ratio, or empty when it is undefined."""
    return log_risk_ratio(comparison.drug_affected, comparison.drug_at_risk, comparison.placebo_affected, comparison.placebo_at_risk) or empty_statistic()


COMPARISON_COLUMNS = ("perturbation_id", "perturbation_type", "perturbation_label", "symptom", "relation", "evidence_class", "source", "source_record_id",
                      "source_term_label", "seriousness", "drug_affected", "drug_at_risk", "placebo_affected", "placebo_at_risk", "verdict", *STATISTIC_COLUMNS)


def build_comparisons(studies: Iterable[dict], drug_name: str, perturbation: dict, preferred_term_to_symptom: dict[str, str],
                      lowest_to_preferred: dict[str, str], minimum_exposure_for_absence: int = DEFAULT_MINIMUM_EXPOSURE_FOR_ABSENCE,
                      other_study_drugs: Iterable[str] = ()) -> list[dict]:
    """Every arm comparison of one drug on a target symptom (COMPARISON_COLUMNS), the ones build_reports leaves out for
    lack of exposure included: the counts, the verdict of the rule in use (1 an excess, 0 an absence at the exposure
    floor, None a non-excess below it) and the comparison's own statistic (decision 9.48, option 3)."""
    other_study_drugs = tuple(other_study_drugs)
    rows = []
    for study in studies:
        selection = select_arms(study, drug_name, other_study_drugs)
        if not selection.usable:
            continue
        for comparison in term_comparisons(study, selection):
            symptom = symptom_of_term(comparison.term, preferred_term_to_symptom, lowest_to_preferred)
            if symptom is None:
                continue
            excess = comparison.drug_rate > comparison.placebo_rate
            verdict = 1 if excess else (0 if comparison.drug_at_risk >= minimum_exposure_for_absence else None)
            rows.append({"perturbation_id": perturbation["perturbation_id"], "perturbation_type": "drug", "perturbation_label": perturbation["perturbation_label"],
                         "symptom": symptom, "relation": "induces", "evidence_class": EVIDENCE_CLASS, "source": SOURCE_NAME, "source_record_id": comparison.nct_id,
                         "source_term_label": comparison.term, "seriousness": comparison.seriousness, "drug_affected": comparison.drug_affected,
                         "drug_at_risk": comparison.drug_at_risk, "placebo_affected": comparison.placebo_affected, "placebo_at_risk": comparison.placebo_at_risk,
                         "verdict": verdict, **comparison_statistic(comparison)})
    return rows


def comparisons_to_dataframe(rows: list[dict], columns: tuple[str, ...]):
    """Comparison rows as a table with the given columns, the statistic columns as floats and the verdict as a nullable
    integer, empty-safe."""
    import pandas as pd

    table = pd.DataFrame(rows, columns=list(columns))
    table["verdict"] = pd.array([None if value is None or (isinstance(value, float) and math.isnan(value)) else int(value) for value in table.verdict], dtype="Int64")
    for column in STATISTIC_COLUMNS:
        if column in ("statistic_name", "statistic_censored"):
            table[column] = table[column].fillna("").astype(str)
        else:
            table[column] = pd.to_numeric(table[column], errors="coerce").astype("float64")
    return table


def perturbations_from_reports(reports_path: Path) -> dict[str, dict]:
    """drug label (lowercase) -> the identity the assembled table uses, read from its report table."""
    import pandas as pd

    reports = pd.read_parquet(reports_path, columns=["perturbation_id", "perturbation_type", "perturbation_label", "perturbation_nodes", "model_description"])
    drugs = reports[reports.perturbation_type == "drug"].drop_duplicates("perturbation_id")
    identities = {}
    for row in drugs.itertuples(index=False):
        match = re.search(r"ChEMBL targets ([^;]*)", str(row.model_description))
        identities[str(row.perturbation_label).lower()] = {"perturbation_id": row.perturbation_id, "perturbation_label": row.perturbation_label,
                                                            "perturbation_nodes": row.perturbation_nodes, "targets": match.group(1) if match else ""}
    return identities


def rows_to_dataframe(rows: list[dict]):
    """Report rows as a table with the columns of REPORT_COLUMNS and the dtypes of evidence_reports.parquet, empty-safe."""
    import pandas as pd

    table = pd.DataFrame(rows, columns=[*REPORT_COLUMNS, *ARM_COLUMNS, *STATISTIC_COLUMNS])
    table["report_value"] = pd.to_numeric(table.report_value, errors="coerce").fillna(0).astype("int64")
    for column in STATISTIC_COLUMNS:  # each comparison's own statistic (decision 9.48, option 3); empty for a row that states none
        if column in ("statistic_name", "statistic_censored"):
            table[column] = table[column].fillna("").astype(str)
        else:
            table[column] = pd.to_numeric(table[column], errors="coerce").astype("float64")
    for column in ARM_COLUMNS[:4]:
        table[column] = pd.array([None if value in ("", None) or (isinstance(value, float) and math.isnan(value)) else int(value) for value in table[column]], dtype="Int64")
    for column in ARM_COLUMNS[4:]:
        table[column] = pd.to_numeric(table[column], errors="coerce").astype("float64")
    table["pubmed_reference_count"] = pd.to_numeric(table.pubmed_reference_count, errors="coerce").fillna(0).astype("int64")
    table["frequency"] = pd.to_numeric(table.frequency, errors="coerce").astype("float64")
    table["frequency_denominator"] = pd.array([None if value in ("", None) or (isinstance(value, float) and math.isnan(value)) else int(value)
                                               for value in table.frequency_denominator], dtype="Int64")
    table["placebo_flag"] = table.placebo_flag.map(lambda value: bool(value) if value != "" else False).astype(bool)
    for column in ("rubric_log_sample_size", "rubric_frequency_known", "rubric_evidence_code_pcs", "rubric_evidence_code_tas", "rubric_evidence_code_iea",
                   "rubric_placebo_controlled", "rubric_curated_synopsis", "rubric_causal_association"):
        table[column] = pd.to_numeric(table[column], errors="coerce").fillna(0.0).astype("float64")
    for column in ("onset", "sex", "evidence_date"):
        table[column] = table[column].map(lambda value: None if value in ("", None) else value).astype(object)
    return table


def write_summary(path: Path, per_drug: dict[str, dict], minimum_exposure_for_absence: int) -> None:
    totals: dict[str, int] = {}
    for counts in per_drug.values():
        for key, value in counts.items():
            if isinstance(value, int):
                totals[key] = totals.get(key, 0) + value
    path.write_text(json.dumps({"source": SOURCE_NAME, "minimum_exposure_for_absence": minimum_exposure_for_absence, "drugs_queried": len(per_drug),
                                "totals": totals, "per_drug": per_drug}, indent=1) + "\n", encoding="utf-8")
