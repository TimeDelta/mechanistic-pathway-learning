"""Write the "better (not more)" label selection: which positive pairs train and score, and which are set aside.

The rule was fixed on 7 October 2026, before any full-graph model was scored, from evidence properties only (never
from a model's predictions), and is recorded in docs/preregistration.md:

  - gene pairs (grade A, HPO): kept when the frequency is at least 0.30, the lower bound of the HPO qualifier
    "frequent" (30 to 79 percent), so "frequent", "very frequent", "obligate" and counted frequencies of 30 percent or
    more are kept; "occasional", "very rare", counted frequencies below 30 percent and pairs with no frequency are set
    aside;
  - drug pairs (grade B, SIDER and OnSIDES labels): kept when the label frequency is at least 1 percent (the CIOMS
    "common" band) or when both SIDER and OnSIDES list the pair (two independent label extractions); the rest,
    mostly label statements with no frequency, are set aside.

A pair set aside is masked, not relabelled: it stays a positive in the outcome matrix and load_experiment_data's
label_mask removes it from the loss and from every metric, so it is neither a positive nor a negative. A
perturbation whose every positive is set aside still contributes its unlabelled pairs, as any perturbation does.

--mask-grades C (better_v2, the user's decision of 8 October 2026) also masks the pairs whose only evidence is of those
grades (grade C: a human association, not a label). Without it such a pair is a negative in scoring and a weak
negative in training. Each gets a row with keep False and masks_a_negative True, which load_experiment_data reads as
"neither positive nor negative"; a pair that also has grade A or B evidence is a positive and gets no such row.

Usage:
  python experiments/build_label_selection.py --evidence-dir data/processed/evidence_full \
      --output data/processed/label_selection/better_v1_full.parquet
  python experiments/build_label_selection.py --evidence-dir data/processed/evidence_full_v2 --mask-grades C \
      --output data/processed/label_selection/better_v2_full_v2.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import DEFAULT_LABEL_GRADES

GENE_MINIMUM_FREQUENCY = 0.30
DRUG_MINIMUM_FREQUENCY = 0.01


def select_pairs(records: pd.DataFrame) -> pd.DataFrame:
    """One row per positive (perturbation_id, symptom) with keep and the reason."""
    positives = records[(records.relation == "induces") & records.grade.isin(DEFAULT_LABEL_GRADES)]
    if "positive_report_count" in positives.columns:
        positives = positives[positives.positive_report_count > 0]
    rows = []
    for (perturbation_id, symptom), group in positives.groupby(["perturbation_id", "symptom"], sort=True):
        perturbation_type = group.perturbation_type.iloc[0]
        frequency = pd.to_numeric(group.label_frequency, errors="coerce").max()
        sources = ";".join(group.source.astype(str))
        if perturbation_type == "gene":
            if pd.isna(frequency):
                keep, reason = False, "gene pair with no frequency"
            elif frequency >= GENE_MINIMUM_FREQUENCY:
                keep, reason = True, f"gene pair, frequency {frequency:.2f} at least {GENE_MINIMUM_FREQUENCY}"
            else:
                keep, reason = False, f"gene pair, frequency {frequency:.2f} below {GENE_MINIMUM_FREQUENCY}"
        else:
            both_sources = "OnSIDES-label" in sources and "SIDER-label" in sources
            if not pd.isna(frequency) and frequency >= DRUG_MINIMUM_FREQUENCY:
                keep, reason = True, f"drug pair, frequency {frequency:.3f} at least {DRUG_MINIMUM_FREQUENCY}"
            elif both_sources:
                keep, reason = True, "drug pair listed by both SIDER and OnSIDES"
            elif pd.isna(frequency):
                keep, reason = False, "drug pair with no frequency, one label source"
            else:
                keep, reason = False, f"drug pair, frequency {frequency:.3f} below {DRUG_MINIMUM_FREQUENCY}, one label source"
        rows.append({"perturbation_id": perturbation_id, "symptom": symptom, "perturbation_type": perturbation_type,
                     "label_frequency": frequency, "keep": keep, "reason": reason})
    return pd.DataFrame(rows)


def masked_negative_pairs(records: pd.DataFrame, mask_grades: list[str]) -> pd.DataFrame:
    """One row (keep False, masks_a_negative True) per (perturbation_id, symptom) with evidence of mask_grades and none of
    the label grades, so that the pair is neither a positive nor a negative."""
    claims = records[records.relation == "induces"]
    if "positive_report_count" in claims.columns:
        claims = claims[claims.positive_report_count > 0]
    labelled = set(zip(*[claims[claims.grade.isin(DEFAULT_LABEL_GRADES)][column] for column in ("perturbation_id", "symptom")]))
    rows = []
    for (perturbation_id, symptom), group in claims[claims.grade.isin(mask_grades)].groupby(["perturbation_id", "symptom"], sort=True):
        if (perturbation_id, symptom) in labelled:
            continue
        rows.append({"perturbation_id": perturbation_id, "symptom": symptom, "perturbation_type": group.perturbation_type.iloc[0],
                     "label_frequency": float("nan"), "keep": False, "reason": f"grade {'/'.join(sorted(group.grade.unique()))} pair only, masked",
                     "masks_a_negative": True})
    return pd.DataFrame(rows, columns=["perturbation_id", "symptom", "perturbation_type", "label_frequency", "keep", "reason", "masks_a_negative"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_selection/better_v1_full.parquet"))
    parser.add_argument("--mask-grades", nargs="*", default=[], help="also mask the pairs whose only evidence is of these grades (C for better_v2)")
    arguments = parser.parse_args()
    if set(arguments.mask_grades) & set(DEFAULT_LABEL_GRADES):
        raise SystemExit(f"--mask-grades {arguments.mask_grades} names a label grade")
    if arguments.output.exists():
        raise SystemExit(f"{arguments.output} exists; a selection is written once (jobs and the lockbox record its hash)")
    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet")
    positive_selection = select_pairs(records)
    masked_negatives = masked_negative_pairs(records, arguments.mask_grades) if arguments.mask_grades else None
    selection = positive_selection if masked_negatives is None else pd.concat([positive_selection.assign(masks_a_negative=False), masked_negatives], ignore_index=True)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    selection.to_parquet(arguments.output, index=False)
    kept = positive_selection[positive_selection.keep]
    summary = {
        "output": str(arguments.output), "evidence_dir": str(arguments.evidence_dir),
        "gene_minimum_frequency": GENE_MINIMUM_FREQUENCY, "drug_minimum_frequency": DRUG_MINIMUM_FREQUENCY,
        "positive_pairs": len(positive_selection), "kept_pairs": len(kept), "set_aside_pairs": int((~positive_selection.keep).sum()),
        "kept_by_type": kept.groupby("perturbation_type").size().to_dict(),
        "perturbations_with_a_kept_positive_by_type": kept.groupby("perturbation_type").perturbation_id.nunique().to_dict(),
        "kept_by_symptom": kept.groupby("symptom").size().to_dict(),
        "all_by_symptom": positive_selection.groupby("symptom").size().to_dict(),
        "set_aside_reasons": positive_selection[~positive_selection.keep].reason.str.replace(r"[0-9.]+ below", "x below", regex=True).value_counts().to_dict(),
    }
    if masked_negatives is not None:
        summary["mask_grades"] = sorted(arguments.mask_grades)
        summary["masked_negative_pairs"] = len(masked_negatives)
        summary["masked_negative_pairs_by_type"] = masked_negatives.groupby("perturbation_type").size().to_dict()
        summary["masked_negative_pairs_by_symptom"] = masked_negatives.groupby("symptom").size().to_dict()
    arguments.output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
