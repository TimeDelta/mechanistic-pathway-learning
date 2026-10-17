"""Soft priors: the evidence below the label bar, turned into a target and a weight for the loss (design section 4.2).

A label (grade A or B) says a pair is positive and is trained toward the positive target at its evidence weight. The
evidence below that bar - an HPO annotation that missed the grade-A criteria, a placebo-controlled trial arm, a FAERS
disproportionality, an IMPC knockout, a text-mined relation - says less, and until now said nothing to the loss: grade
C pairs were set aside and grades D and E never left the report table. Here each such pair becomes a soft prior,

    target     = the tier's prior probability (positive evidence) or 0 (an absence),
    weight     = global strength * tier unit strength * min(n, cap) ** tempering,

where n is the number of independent records behind the pair (distinct PMIDs where the reports cite one, else distinct
source records). The loss that takes them is the weighted binary cross-entropy already in use: with a fractional
target q the per-pair loss is the expectation of the hard loss over a Bernoulli(q) truth, which is the proper way to
say "probably positive". The weight is a power-prior exponent times a tempered count (Ibrahim et al. 2015, "The power
prior: theory and applications"): the step from one independent publication to two carries the signal (Cusick et al.
2009), and beyond a small cap a count measures how much a pair was studied, not whether it is true (Schaefer et al.
2015; Gillis and Pavlidis 2012). The tiers are discrete, capped and chosen by evidence type, after ClinGen's
gene-disease validity scheme (Strande et al. 2017), which caps points so that "a preponderance of weak evidence"
cannot inflate a classification.

Rules the arrays enforce:
- A soft prior never touches a labelled positive, a labelled negative, or a pair the selection set aside for a
  controlled comparison; it applies to unlabelled pairs and to pairs set aside only for lacking a label grade.
- A pair with evidence in both directions (a trial arm above placebo and a FAERS absence, say) is contested and takes
  no prior; it is counted and listed.
- One class decides a pair's prior: the highest-ranked class present (assemble_evidence_table.EVIDENCE_CLASS_RANK), so
  a human source is never averaged with a mouse or a sentence.
- The prior is applied to the rows a fit trains on and to nothing scored; the trainer does that.
- The degree-preserving permutation (permute_soft_priors) keeps every perturbation's and every symptom's number of
  soft pairs and moves the pairs among eligible cells: a prior that is only popularity survives it, one that carries
  pair-specific information does not (Zietz et al. 2024, the degree-only edge prior).
"""
from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evidence.assemble_evidence_table import EVIDENCE_CLASS_RANK

DIRECTIONS = ("positive", "negative")
DEFAULT_TIERS_PATH = Path("configs/soft_prior_tiers.csv")
DEFAULT_TEMPERING_EXPONENT = 0.5
DEFAULT_COUNT_CAP = 2
DEFAULT_RECORD_GRADES = ("C",)
# the literature relations admitted as evidence of induction; co-mention ("associated_with") is not a direction
LITERATURE_SOURCES_ADMITTED = ("PubTator3-cause", "PubTator3-stimulate", "PubTator3-inhibit", "CTD-curated")
# A trial-arm comparison is positive in its report table on any excess over placebo (3 of 100 against 2 of 100 counts),
# which is the right record of the measurement and too weak a bar for a prior. Here an excess enters only when the
# lower 95 percent bound of the risk difference (Newcombe's score interval, clinical_trial_adverse_events) is above 0,
# and an absence only when the upper bound is below this margin: the drug adds at most five in a hundred.
TRIAL_ARM_SOURCE = "ClinicalTrials-arms"
ABSENCE_RISK_DIFFERENCE_MARGIN = 0.05
# Decision 9.48, option 3 (the owner, 16 October 2026: "Carrying stat to calibrated proper is the way"): a report that
# carries its comparison's own statistic (evidence/comparison_statistics.py) passes it on, so the calibrated reliability
# (evidence/calibrated_reliability.py) reads the statistic and not only the verdict cut from it. Per pair and class the
# statistic of the most informative comparison is kept, with the number of comparisons behind it: the comparisons of
# one pair (several MedDRA terms of one symptom in one source, several arms of one trial) share reports or patients,
# so their information is not added.
STATISTIC_EVIDENCE_COLUMNS = ("statistic_probability_of_excess", "statistic_information", "statistic_comparisons")
EVIDENCE_KEYS = ["perturbation_id", "symptom", "evidence_class", "direction"]


@dataclass(frozen=True)
class Tier:
    evidence_class: str
    direction: str
    prior_probability: float
    unit_strength: float
    note: str = ""


def read_tiers(path: Path = DEFAULT_TIERS_PATH) -> dict[tuple[str, str], Tier]:
    """configs/soft_prior_tiers.csv as {(evidence_class, direction): Tier}; a class without a row takes no prior."""
    tiers = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            direction = row["direction"].strip()
            if direction not in DIRECTIONS:
                raise ValueError(f"{path}: direction {direction!r} is not one of {DIRECTIONS}")
            probability, strength = float(row["prior_probability"]), float(row["unit_strength"])
            if not 0.0 <= probability <= 1.0 or strength < 0.0:
                raise ValueError(f"{path}: {row['evidence_class']} {direction} needs a probability in [0, 1] and a strength >= 0")
            if direction == "negative" and probability != 0.0:
                raise ValueError(f"{path}: an absence is trained toward 0, as the registered placebo rule is; {row['evidence_class']} negative has {probability}")
            tiers[(row["evidence_class"].strip(), direction)] = Tier(row["evidence_class"].strip(), direction, probability, strength, row.get("note", ""))
    return tiers


def tempered_count(count, cap: int = DEFAULT_COUNT_CAP, exponent: float = DEFAULT_TEMPERING_EXPONENT):
    """min(count, cap) ** exponent, so the second independent record adds and the tenth does not."""
    return np.minimum(np.asarray(count, dtype=float), float(cap)) ** float(exponent)


def record_key(row) -> str:
    """What makes two reports of one pair independent: the cited PMID where there is one, else the source record."""
    pmid = getattr(row, "pmid", None)
    if pmid is not None and not pd.isna(pmid) and str(pmid).strip():
        return f"PMID:{str(pmid).strip()}"
    references = getattr(row, "references", None)
    if references is not None and not pd.isna(references) and str(references).strip():
        first = str(references).split(";")[0].strip()
        if first.startswith("PMID:"):
            return first
    record = getattr(row, "source_record_id", None)
    if record is not None and not pd.isna(record) and str(record).strip():
        return f"{row.source}|{record}"
    return str(row.report_id)


def admit_trial_arm_reports(table: pd.DataFrame, counts: dict | None = None, absence_margin: float = ABSENCE_RISK_DIFFERENCE_MARGIN) -> pd.DataFrame:
    """The trial-arm rows whose interval supports them: an excess with its lower bound above 0, an absence with its upper
    bound below the margin. Rows of other sources pass; a trial-arm table without the interval columns passes as it is
    and is counted as such."""
    if "source" not in table.columns or not (table.source == TRIAL_ARM_SOURCE).any():
        return table
    is_trial = table.source == TRIAL_ARM_SOURCE
    if "risk_difference_lower_95" not in table.columns or "risk_difference_upper_95" not in table.columns:
        if counts is not None:
            counts["trial_arm_rows_without_an_interval"] = int(is_trial.sum())
        return table
    lower, upper = pd.to_numeric(table.risk_difference_lower_95, errors="coerce"), pd.to_numeric(table.risk_difference_upper_95, errors="coerce")
    positive = table.report_value.astype(int) == 1
    admitted = ~is_trial | (positive & (lower > 0)) | (~positive & (upper < absence_margin))
    if counts is not None:
        counts["trial_arm_excesses_admitted"] = int((is_trial & positive & admitted).sum())
        counts["trial_arm_excesses_not_supported_by_the_interval"] = int((is_trial & positive & ~admitted).sum())
        counts["trial_arm_absences_admitted"] = int((is_trial & ~positive & admitted).sum())
        counts["trial_arm_absences_with_an_excess_left_open"] = int((is_trial & ~positive & ~admitted).sum())
    return table[admitted]


def soft_evidence_from_reports(reports: pd.DataFrame, relation: str = "induces", counts: dict | None = None) -> pd.DataFrame:
    """Per (perturbation, symptom, evidence_class, direction): the number of independent records and the sources.

    Takes a report table of the evidence_reports schema (one row per report, report_value 1 or 0). Literature rows
    enter only from the admitted sources (LITERATURE_SOURCES_ADMITTED) and only on the relation asked for; the
    co-mention relation "associated_with" is not evidence of induction. Trial-arm rows enter only where their interval
    supports them (admit_trial_arm_reports). Reports with the symptom-level crosswalk flagged as indicated (an
    indication, not an adverse event) are the builders' business and are not re-filtered here. `counts`, when given,
    receives the admission counts.
    """
    empty = pd.DataFrame(columns=[*EVIDENCE_KEYS, "n_distinct", "sources", *STATISTIC_EVIDENCE_COLUMNS])
    if reports.empty:
        return empty
    table = reports[reports.relation == relation]
    if "evidence_class" in table.columns:
        is_literature = table.evidence_class == "literature"
        table = table[~is_literature | table.source.isin(LITERATURE_SOURCES_ADMITTED)]
    table = admit_trial_arm_reports(table, counts)
    if table.empty:
        return empty
    has_statistic = "statistic_probability_of_excess" in table.columns and "statistic_information" in table.columns
    rows = []
    for row in table.itertuples(index=False):
        probability = float(row.statistic_probability_of_excess) if has_statistic and not pd.isna(row.statistic_probability_of_excess) else np.nan
        information = float(row.statistic_information) if has_statistic and not pd.isna(row.statistic_information) else np.nan
        rows.append((str(row.perturbation_id), str(row.symptom), str(row.evidence_class), "positive" if int(row.report_value) == 1 else "negative", record_key(row),
                     str(row.source), probability, information))
    keyed = pd.DataFrame(rows, columns=[*EVIDENCE_KEYS, "record", "source", "statistic_probability_of_excess", "statistic_information"])
    grouped = keyed.groupby(EVIDENCE_KEYS, sort=True)
    evidence = grouped.agg(n_distinct=("record", "nunique"), sources=("source", lambda values: ";".join(sorted(set(values)))),
                           statistic_comparisons=("statistic_information", "count")).reset_index()
    return with_most_informative_statistic(evidence, keyed)


def with_most_informative_statistic(evidence: pd.DataFrame, rows: pd.DataFrame) -> pd.DataFrame:
    """evidence (one row per pair, class and direction) with the probability of excess and the information of the most
    informative comparison among rows of the same keys; NaN where none carries a statistic."""
    columns = ["statistic_probability_of_excess", "statistic_information"]
    evidence = evidence.drop(columns=[column for column in columns if column in evidence.columns])
    carrying = rows.dropna(subset=["statistic_information"]) if "statistic_information" in rows.columns else rows.iloc[0:0]
    if carrying.empty:
        return evidence.assign(statistic_probability_of_excess=np.nan, statistic_information=np.nan)[[*evidence.columns, *columns]]
    best = carrying.sort_values("statistic_information", ascending=False, kind="stable").drop_duplicates(EVIDENCE_KEYS)[[*EVIDENCE_KEYS, *columns]]
    return evidence.merge(best, on=EVIDENCE_KEYS, how="left")


def soft_evidence_from_records(records: pd.DataFrame, grades=DEFAULT_RECORD_GRADES, relation: str = "induces") -> pd.DataFrame:
    """The evidence_records rows of the soft grades (grade C by default), in the same shape as the report evidence.

    A record is positive when it holds positive reports and no negative one, negative when the reverse; both is
    contested and is kept as two rows so resolve_pairs can list it. The count is the distinct PubMed references where
    the record cites any, else its reports.
    """
    columns = ["perturbation_id", "symptom", "evidence_class", "direction", "n_distinct", "sources"]
    if records.empty or "grade" not in records.columns:
        return pd.DataFrame(columns=columns)
    table = records[(records.relation == relation) & records.grade.isin(grades)]
    rows = []
    for row in table.itertuples(index=False):
        positive = int(getattr(row, "positive_report_count", 1) or 0)
        negative = int(getattr(row, "negative_report_count", 0) or 0)
        references = getattr(row, "distinct_pubmed_reference_count", None)
        references = 0 if references is None or pd.isna(references) else int(references)
        if positive > 0:
            rows.append((str(row.perturbation_id), str(row.symptom), str(row.evidence_class), "positive", max(references, positive, 1), str(row.source)))
        if negative > 0:
            rows.append((str(row.perturbation_id), str(row.symptom), str(row.evidence_class), "negative", max(negative, 1), str(row.source)))
    return pd.DataFrame(rows, columns=columns)


@dataclass
class ResolvedPriors:
    pairs: pd.DataFrame  # perturbation_id, symptom, evidence_class, direction, n_distinct, prior_probability, strength, sources
    contested: pd.DataFrame  # the pairs with evidence in both directions, which take no prior
    without_a_tier: dict[str, int]  # evidence classes present in the evidence that the tier table does not name


def resolve_pairs(evidence: pd.DataFrame, tiers: dict[tuple[str, str], Tier], tempering_exponent: float = DEFAULT_TEMPERING_EXPONENT,
                  count_cap: int = DEFAULT_COUNT_CAP) -> ResolvedPriors:
    """One prior per (perturbation, symptom): the highest-ranked class with a tier decides; both directions is contested.

    The direction is read over every class with a tier, so a trial-arm absence (human_association, negative) against
    an IMPC hit (animal, positive) is contested rather than resolved by rank: the two kinds of evidence do not share a
    species, and the rule for that disagreement is a decision, not a sort order.
    """
    columns = ["perturbation_id", "symptom", "evidence_class", "direction", "n_distinct", "prior_probability", "strength", "sources", *STATISTIC_EVIDENCE_COLUMNS]
    if evidence.empty:
        return ResolvedPriors(pd.DataFrame(columns=columns), pd.DataFrame(columns=["perturbation_id", "symptom", "classes", "directions"]), {})
    for column in STATISTIC_EVIDENCE_COLUMNS:  # tables built before 9.48 (3), and grade-C records, carry no statistic
        if column not in evidence.columns:
            evidence = evidence.assign(**{column: 0 if column == "statistic_comparisons" else np.nan})
    statistics = evidence[[*EVIDENCE_KEYS, "statistic_probability_of_excess", "statistic_information"]]
    evidence = evidence.groupby(EVIDENCE_KEYS, sort=True).agg(
        n_distinct=("n_distinct", "sum"), sources=("sources", lambda values: ";".join(sorted(set(";".join(values).split(";"))))),
        statistic_comparisons=("statistic_comparisons", lambda values: int(pd.to_numeric(values, errors="coerce").fillna(0).sum()))).reset_index()
    evidence = with_most_informative_statistic(evidence, statistics)
    has_tier = [(row.evidence_class, row.direction) in tiers for row in evidence.itertuples(index=False)]
    without_a_tier = evidence[[not value for value in has_tier]].groupby("evidence_class").size().to_dict()
    tiered = evidence[has_tier]
    resolved, contested = [], []
    for (perturbation_id, symptom), pair in tiered.groupby(["perturbation_id", "symptom"], sort=True):
        directions = set(pair.direction)
        if len(directions) > 1:
            contested.append((perturbation_id, symptom, ";".join(sorted(set(pair.evidence_class))), ";".join(sorted(directions))))
            continue
        ranks = pair.evidence_class.map(lambda name: EVIDENCE_CLASS_RANK.get(name, len(EVIDENCE_CLASS_RANK)))
        first = pair.loc[ranks.idxmin()]
        tier = tiers[(first.evidence_class, first.direction)]
        strength = float(tier.unit_strength * tempered_count(first.n_distinct, count_cap, tempering_exponent))
        resolved.append((perturbation_id, symptom, first.evidence_class, first.direction, int(first.n_distinct), tier.prior_probability, strength, first.sources,
                         first.statistic_probability_of_excess, first.statistic_information, int(first.statistic_comparisons)))
    return ResolvedPriors(pd.DataFrame(resolved, columns=columns), pd.DataFrame(contested, columns=["perturbation_id", "symptom", "classes", "directions"]),
                          {str(key): int(value) for key, value in without_a_tier.items()})


def soft_prior_matrices(resolved: ResolvedPriors, perturbation_ids: list[str], symptoms: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """[num_perturbations, num_symptoms] prior probability (NaN where no prior) and strength (0 where none), aligned to
    the perturbation and symptom order of an ExperimentData; pairs outside that universe are dropped."""
    row_of = {perturbation_id: index for index, perturbation_id in enumerate(perturbation_ids)}
    column_of = {symptom: index for index, symptom in enumerate(symptoms)}
    probability = np.full((len(perturbation_ids), len(symptoms)), np.nan)
    strength = np.zeros_like(probability)
    for row in resolved.pairs.itertuples(index=False):
        position, column = row_of.get(row.perturbation_id), column_of.get(row.symptom)
        if position is None or column is None:
            continue
        probability[position, column] = row.prior_probability
        strength[position, column] = row.strength
    return probability, strength


def eligible_cells(outcomes: np.ndarray, label_mask: np.ndarray | None, evidenced_negative: np.ndarray | None, masked_for_soft_evidence: np.ndarray | None) -> np.ndarray:
    """Where a soft prior may apply: not a labelled positive, not an evidenced negative, and either labelled (so an
    unlabelled negative) or set aside only for lacking a label grade."""
    eligible = np.asarray(outcomes) == 0
    if evidenced_negative is not None:
        eligible &= ~np.asarray(evidenced_negative, dtype=bool)
    if label_mask is not None:
        allowed = np.asarray(label_mask, dtype=bool)
        if masked_for_soft_evidence is not None:
            allowed = allowed | np.asarray(masked_for_soft_evidence, dtype=bool)
        eligible &= allowed
    return eligible


def apply_soft_priors(outcomes: np.ndarray, weights: np.ndarray, probability: np.ndarray, strength: np.ndarray, eligible: np.ndarray, fit_rows: np.ndarray,
                      global_strength: float) -> tuple[np.ndarray, np.ndarray, dict]:
    """Copies of the outcome and weight matrices with the priors written into the eligible cells of fit_rows.

    A positive prior sets the target to its probability and the weight to the scaled strength; an absence keeps the
    target at 0 and raises the weight to the scaled strength where that is more than the pair's weight, as the
    registered placebo rule does. Rows outside fit_rows are untouched, so validation and test rows never see a prior.
    """
    new_outcomes, new_weights = np.array(outcomes, dtype=float, copy=True), np.array(weights, dtype=float, copy=True)
    in_fit = np.zeros(len(outcomes), dtype=bool)
    in_fit[np.asarray(fit_rows, dtype=int)] = True
    if global_strength <= 0.0:  # off: a target moved at weight 0 would still change the recorded matrices
        return new_outcomes, new_weights, {"positive_cells": 0, "absence_cells": 0, "fit_rows": int(in_fit.sum()), "priors_outside_eligible_cells": 0,
                                           "mean_positive_weight": None, "mean_positive_target": None}
    applicable = eligible & in_fit[:, None] & (strength > 0)
    positive = applicable & np.isfinite(probability) & (probability > 0)
    negative = applicable & ~positive
    new_outcomes[positive] = probability[positive]
    new_weights[positive] = global_strength * strength[positive]
    new_weights[negative] = np.maximum(new_weights[negative], global_strength * strength[negative])
    summary = {"positive_cells": int(positive.sum()), "absence_cells": int(negative.sum()), "fit_rows": int(in_fit.sum()),
               "priors_outside_eligible_cells": int(((strength > 0) & in_fit[:, None] & ~eligible).sum()),
               "mean_positive_weight": float(new_weights[positive].mean()) if positive.any() else None,
               "mean_positive_target": float(new_outcomes[positive].mean()) if positive.any() else None}
    return new_outcomes, new_weights, summary


def permute_soft_priors(probability: np.ndarray, strength: np.ndarray, eligible: np.ndarray, rows: np.ndarray, random_seed: int = 0,
                        swaps_per_pair: int = 10) -> tuple[np.ndarray, np.ndarray, dict]:
    """The priors of the given rows moved among those rows' eligible cells, each perturbation keeping its number of
    soft pairs and each symptom keeping its number, by column swaps between pairs of prior cells (a degree-preserving
    rewiring of the bipartite prior graph); a prior's probability and strength travel with it. Rows outside `rows`
    keep their priors. Returns the new matrices and a record with the SHA-256 of the moved layout."""
    probability, strength = np.array(probability, copy=True), np.array(strength, copy=True)
    in_rows = np.zeros(len(strength), dtype=bool)
    in_rows[np.asarray(rows, dtype=int)] = True
    cells = [tuple(cell) for cell in np.argwhere((strength > 0) & in_rows[:, None])]
    occupied = set(cells)
    generator = np.random.default_rng(random_seed)
    attempted, swapped = 0, 0
    for _ in range(swaps_per_pair * len(cells)):
        if len(cells) < 2:
            break
        first, second = generator.choice(len(cells), size=2, replace=False)
        (row_a, column_a), (row_b, column_b) = cells[first], cells[second]
        attempted += 1
        if column_a == column_b or row_a == row_b:
            continue
        target_a, target_b = (row_a, column_b), (row_b, column_a)
        if target_a in occupied or target_b in occupied or not eligible[target_a] or not eligible[target_b]:
            continue
        probability[target_a], probability[target_b] = probability[row_a, column_a], probability[row_b, column_b]
        strength[target_a], strength[target_b] = strength[row_a, column_a], strength[row_b, column_b]
        probability[row_a, column_a], probability[row_b, column_b] = np.nan, np.nan
        strength[row_a, column_a], strength[row_b, column_b] = 0.0, 0.0
        occupied -= {(row_a, column_a), (row_b, column_b)}
        occupied |= {target_a, target_b}
        cells[first], cells[second] = target_a, target_b
        swapped += 1
    layout = np.argwhere(strength > 0)
    record = {"seed": int(random_seed), "pairs": len(cells), "swaps_attempted": attempted, "swaps_made": swapped,
              "layout_sha256": hashlib.sha256(layout.astype(np.int64).tobytes()).hexdigest()}
    return probability, strength, record


def literature_popularity(perturbation_ids: list[str], literature_reports: pd.DataFrame) -> np.ndarray:
    """Distinct PMIDs per perturbation in the literature report table, any symptom and relation: the popularity axis
    the degree bins do not measure. 0 for a perturbation the table never mentions."""
    if literature_reports.empty or "pmid" not in literature_reports.columns:
        return np.zeros(len(perturbation_ids))
    counts = literature_reports.groupby("perturbation_id").pmid.nunique()
    return np.array([float(counts.get(perturbation_id, 0)) for perturbation_id in perturbation_ids])


def summarise(resolved: ResolvedPriors, probability: np.ndarray, strength: np.ndarray, eligible: np.ndarray) -> dict:
    """Counts by class and direction, the contested pairs, and how many resolved pairs fall on eligible cells."""
    pairs = resolved.pairs
    by_class = {}
    if len(pairs):
        for (evidence_class, direction), group in pairs.groupby(["evidence_class", "direction"]):
            by_class[f"{evidence_class}:{direction}"] = {"pairs": int(len(group)), "mean_strength": float(group.strength.mean()),
                                                         "n_distinct_at_least_two": int((group.n_distinct >= 2).sum())}
    placed = int(((strength > 0)).sum())
    return {"resolved_pairs": int(len(pairs)), "pairs_in_the_perturbation_universe": placed, "pairs_on_eligible_cells": int(((strength > 0) & eligible).sum()),
            "pairs_on_labelled_or_set_aside_cells": int(((strength > 0) & ~eligible).sum()), "contested_pairs": int(len(resolved.contested)),
            "contested_examples": resolved.contested.head(10).to_dict("records"), "classes_without_a_tier": resolved.without_a_tier, "by_class_and_direction": by_class,
            "pairs_with_a_statistic": int(pd.to_numeric(pairs.get("statistic_information", pd.Series(dtype=float)), errors="coerce").notna().sum())}
