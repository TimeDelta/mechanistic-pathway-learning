"""FAERS spontaneous reports (openFDA) as disproportionality evidence, positive and negative, per drug and symptom.

What this measures and what it cannot. A spontaneous-report system has no denominator of exposed patients, so the
only comparison available is disproportionality: whether a reaction is reported with the drug more often than with
everything else. The proportional reporting ratio (PRR) and the reporting odds ratio (ROR) are the standard
statistics (Evans, Waller and Davis 2001; Rothman, Lanes and Sacks 2004); the usual signal rule is PRR >= 2 with a
chi-square of at least 4 and at least 3 reports. A reaction reported with the drug *less* often than background,
with a confidence bound that excludes 1 and a large number of reports behind the drug, is as much a measurement as
a signal is: the system was exposed to the drug many thousands of times and the reaction did not surface in
proportion. That is the negative this source supplies. Its biases are the ones every reader of FAERS knows and the
limitations text names them on each report: stimulated reporting after publicity, confounding by indication (an
antidepressant's "depression" reports), duplicates, co-medication (one report names several drugs, so counts are
correlated across drugs that are prescribed together, a structure the target-based leakage groups do not see), and
the absence of any exposure count. For those reasons the rows enter at the lowest evidence tier.

Licence: "Unless otherwise noted, the content, data, documentation, code, and related materials on openFDA is public
domain and made available with a Creative Commons CC0 1.0 Universal dedication" (https://open.fda.gov/terms/, read
15 October 2026).
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS

SOURCE_NAME = "FAERS-disproportionality"
EVIDENCE_CODE = "faers_disproportionality"
DEFAULT_SIGNAL_PRR = 2.0
DEFAULT_SIGNAL_CHI_SQUARE = 4.0
DEFAULT_SIGNAL_MINIMUM_REPORTS = 3
DEFAULT_NEGATIVE_MINIMUM_DRUG_REPORTS = 2000  # reports naming the drug before an under-reported reaction counts as an absence


@dataclass(frozen=True)
class Disproportionality:
    reaction: str
    with_drug: int          # a: reports with the drug and the reaction
    drug_total: int         # a + b: reports with the drug
    reaction_total: int     # a + c: reports with the reaction
    all_reports: int        # N
    bounded: bool = False   # with_drug is an upper bound: the reaction is not among the drug's listed top reactions, so its count is at most the smallest listed

    @property
    def table(self) -> tuple[float, float, float, float]:
        a = float(self.with_drug)
        b = float(self.drug_total - self.with_drug)
        c = float(self.reaction_total - self.with_drug)
        d = float(self.all_reports - a - b - c)
        return a, b, c, d

    @property
    def prr(self) -> float | None:
        """None when a margin is empty; infinite when the reaction is reported with the drug and never without it."""
        a, b, c, d = self.table
        if a + b <= 0 or c + d <= 0:
            return None
        if c <= 0:
            return math.inf if a > 0 else None
        return (a / (a + b)) / (c / (c + d))

    @property
    def ror(self) -> float | None:
        a, b, c, d = (value + 0.5 for value in self.table)  # Haldane–Anscombe correction keeps a zero cell finite
        return (a * d) / (b * c) if b * c > 0 else None

    @property
    def ror_interval(self) -> tuple[float, float] | None:
        a, b, c, d = (value + 0.5 for value in self.table)
        if min(a, b, c, d) <= 0:
            return None
        log_ror = math.log((a * d) / (b * c))
        half_width = 1.96 * math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        return math.exp(log_ror - half_width), math.exp(log_ror + half_width)

    @property
    def chi_square(self) -> float:
        """Pearson chi-square with Yates's correction on the 2x2 table."""
        a, b, c, d = self.table
        n = a + b + c + d
        row1, row2, col1, col2 = a + b, c + d, a + c, b + d
        if min(row1, row2, col1, col2) <= 0:
            return 0.0
        return n * (abs(a * d - b * c) - n / 2) ** 2 / (row1 * row2 * col1 * col2)

    def classify(self, signal_prr: float = DEFAULT_SIGNAL_PRR, signal_chi_square: float = DEFAULT_SIGNAL_CHI_SQUARE,
                 signal_minimum_reports: int = DEFAULT_SIGNAL_MINIMUM_REPORTS,
                 negative_minimum_drug_reports: int = DEFAULT_NEGATIVE_MINIMUM_DRUG_REPORTS) -> int | None:
        """1 for a signal, 0 for a measured under-reporting, None for neither. A bounded count can never be a signal."""
        prr = self.prr
        if not self.bounded and prr is not None and prr >= signal_prr and self.chi_square >= signal_chi_square and self.with_drug >= signal_minimum_reports:
            return 1
        interval = self.ror_interval
        if interval is not None and interval[1] < 1.0 and self.drug_total >= negative_minimum_drug_reports:
            return 0
        return None


def disproportionalities(reaction_counts: dict[str, int], drug_total: int, background_counts: dict[str, int], all_reports: int,
                         unlisted_upper_bound: int | None = None) -> list[Disproportionality]:
    """One Disproportionality per reaction that is either reported with the drug or in the background set.

    openFDA's count endpoint lists a drug's most reported reactions (100 without a key), so a background reaction the
    list leaves out has an unknown count no larger than the smallest listed one; with `unlisted_upper_bound` it enters
    at that bound, marked, which can only make an under-reporting more conservative and can never make a signal.
    """
    reactions = set(reaction_counts) | set(background_counts)
    items = []
    for reaction in sorted(reactions):
        if reaction in reaction_counts:
            items.append(Disproportionality(reaction, int(reaction_counts[reaction]), int(drug_total), int(background_counts.get(reaction, 0)), int(all_reports)))
        elif unlisted_upper_bound is not None:
            items.append(Disproportionality(reaction, int(unlisted_upper_bound), int(drug_total), int(background_counts.get(reaction, 0)), int(all_reports), bounded=True))
        else:
            items.append(Disproportionality(reaction, 0, int(drug_total), int(background_counts.get(reaction, 0)), int(all_reports)))
    return items


def build_reports(drug_label: str, perturbation: dict, reaction_counts: dict[str, int], drug_total: int, background_counts: dict[str, int],
                  all_reports: int, preferred_term_to_symptom: dict[str, str], indicated_symptoms: set[str], retrieved: str,
                  thresholds: dict | None = None, unlisted_upper_bound: int | None = None) -> tuple[list[dict], dict]:
    """Report rows (REPORT_COLUMNS) for one drug: a signal row per reaction on a target symptom, an absence row per
    under-reported reaction on a target symptom, nothing for the symptoms the drug is indicated for.

    Reaction names are FAERS preferred terms in upper case; `preferred_term_to_symptom` is keyed by lower case."""
    thresholds = thresholds or {}
    rows: list[dict] = []
    counts = {"reactions_with_drug": len(reaction_counts), "reactions_on_a_symptom": 0, "signals": 0, "absences": 0, "neither": 0, "indication_confounded": 0}
    for item in disproportionalities(reaction_counts, drug_total, background_counts, all_reports, unlisted_upper_bound):
        symptom = preferred_term_to_symptom.get(item.reaction.strip().lower())
        if symptom is None:
            continue
        counts["reactions_on_a_symptom"] += 1
        if symptom in indicated_symptoms:
            counts["indication_confounded"] += 1
            continue
        report_value = item.classify(**thresholds)
        if report_value is None:
            counts["neither"] += 1
            continue
        counts["signals" if report_value == 1 else "absences"] += 1
        interval = item.ror_interval
        row = {column: "" for column in REPORT_COLUMNS}
        row.update({
            "report_id": f"{SOURCE_NAME}|{perturbation['perturbation_id']}|{item.reaction.lower().replace(' ', '_')}|{symptom}",
            "perturbation_id": perturbation["perturbation_id"], "perturbation_type": "drug", "perturbation_label": perturbation["perturbation_label"],
            "symptom": symptom, "relation": "induces", "evidence_class": "pharmacological", "source": SOURCE_NAME, "report_value": report_value,
            "source_record_id": f"openFDA:{drug_label}", "source_record_label": drug_label, "source_term_id": "", "source_term_label": item.reaction,
            "evidence_code": EVIDENCE_CODE,
            "model_description": (f"human; FAERS spontaneous reports via openFDA (retrieved {retrieved}); {'at most ' if item.bounded else ''}{item.with_drug} of {item.drug_total} reports naming the drug "
                                  f"carry the reaction against {item.reaction_total} of {item.all_reports} overall; PRR {item.prr if item.prr is None else round(item.prr, 3)}, "
                                  f"ROR {None if item.ror is None else round(item.ror, 3)} [{'-' if interval is None else f'{interval[0]:.3f}, {interval[1]:.3f}'}], "
                                  f"chi-square {item.chi_square:.1f}; ChEMBL targets {perturbation.get('targets', '')}"),
            "frequency": None, "frequency_denominator": None, "placebo_flag": False, "onset": None, "sex": None, "references": f"openFDA drug/event, retrieved {retrieved}",
            "pubmed_reference_count": 0, "evidence_date": retrieved[:10] if retrieved else None, "evidence_date_source": "release" if retrieved else "",
            "rubric_log_sample_size": math.log1p(item.drug_total), "rubric_frequency_known": 0.0,
            "rubric_evidence_code_pcs": 0.0, "rubric_evidence_code_tas": 0.0, "rubric_evidence_code_iea": 0.0, "rubric_placebo_controlled": 0.0,
            "rubric_curated_synopsis": 0.0, "association_type": "", "rubric_causal_association": 0.0,
            "limitations": ("spontaneous reports: no exposure denominator, stimulated reporting after publicity, duplicates, co-medication (one report names several drugs), "
                            "reactions coded by the reporter; a signal is disproportion against all other drugs, not causation; an absence is under-reporting relative to background "
                            f"with at least {thresholds.get('negative_minimum_drug_reports', DEFAULT_NEGATIVE_MINIMUM_DRUG_REPORTS)} reports naming the drug, not a measured rate"
                            + ("; the reaction is not among the drug's listed top reactions, so its count is an upper bound" if item.bounded else "")),
            "perturbation_nodes": perturbation.get("perturbation_nodes", ""),
        })
        rows.append(row)
    return rows, counts


def indicated_symptoms_by_drug(reports_path: Path) -> dict[str, set[str]]:
    """perturbation_id -> symptoms the assembled table holds a 'relieves' report for (the drug's indications)."""
    import pandas as pd

    reports = pd.read_parquet(reports_path, columns=["perturbation_id", "symptom", "relation"])
    relieves = reports[reports.relation == "relieves"]
    return {perturbation: set(group.symptom) for perturbation, group in relieves.groupby("perturbation_id")}


def write_summary(path: Path, per_drug: dict[str, dict], background: dict, thresholds: dict) -> None:
    totals: dict[str, int] = {}
    for counts in per_drug.values():
        for key, value in counts.items():
            if isinstance(value, int):
                totals[key] = totals.get(key, 0) + value
    path.write_text(json.dumps({"source": SOURCE_NAME, "thresholds": thresholds, "background": background, "drugs_queried": len(per_drug),
                                "totals": totals, "per_drug": per_drug}, indent=1) + "\n", encoding="utf-8")
