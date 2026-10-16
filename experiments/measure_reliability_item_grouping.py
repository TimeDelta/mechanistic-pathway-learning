"""Measure what the reliability fit's item grouping does on a joined evidence table (decision 9.55).

Decision 9.55 found that joining the three absence sources (trial arms, FAERS, IMPC) to the label sources makes the
weighted Dawid-Skene fit identified and badly behaved: EM hits its iteration cap, IMPC comes out degenerate and the
animal and human_association item groups take a prevalence near zero, because an item group is an evidence class and
those classes are overwhelmingly absences. The question it left open is whether the items should be grouped another
way. This script fits the same report matrices under four groupings and prints what each gives, so the owner rules on
a measurement rather than on the argument:

  evidence_class    the assembler's grouping (one prevalence per primary class of the item)
  perturbation_type one prevalence for genes and one for drugs, the original reason for grouping at all (no source
                    covers both, so a pooled prevalence would hand the HPO agreement to every drug pair)
  coverage_family   one prevalence per combination of source families that cover the item (genetics, label, trial,
                    faers, impc): an item reported only by IMPC is compared with items reported only by IMPC, and an
                    item the human genetics sources also cover with those
  pooled            one prevalence for every item

Nothing here changes a table or a rule; the output is a JSON of the four fits and a table on standard output.

Usage:
  python experiments/measure_reliability_item_grouping.py --evidence-dir <joined table directory> --output <json>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evidence.assemble_evidence_table import explicit_report_cells
from mechanistic_pathway_learning.evidence.evidence_reliability_model import fit_weighted_dawid_skene, report_count_matrices
from mechanistic_pathway_learning.evidence.evidence_reports import SOURCE_NAMES, RubricWeightDefaults

SOURCE_FAMILIES = {
    "HPO-OMIM": "genetics", "HPO-Orphanet": "genetics",
    "SIDER-label": "label", "SIDER-indication": "label", "OnSIDES-label": "label",
    "ClinicalTrials-arms": "trial", "FAERS-disproportionality": "faers", "IMPC-knockout": "impc",
}
GROUPINGS = ("evidence_class", "perturbation_type", "coverage_family", "pooled")


def coverage_family_labels(coverage: np.ndarray, source_names: list[str]) -> list[str]:
    """One label per item: the sorted, '+'-joined families of the sources that cover it ('none' when no source does)."""
    families = [SOURCE_FAMILIES.get(name, name) for name in source_names]
    labels = []
    for row in np.asarray(coverage) > 0:
        covering = sorted({family for family, covered in zip(families, row) if covered})
        labels.append("+".join(covering) if covering else "none")
    return labels


def item_groups_for(grouping: str, observations: pd.DataFrame, coverage: np.ndarray, source_names: list[str]) -> list[str]:
    if grouping == "evidence_class":
        return [str(value) for value in observations.evidence_class]
    if grouping == "perturbation_type":
        return [str(value) for value in observations.perturbation_type]
    if grouping == "coverage_family":
        return coverage_family_labels(coverage, source_names)
    if grouping == "pooled":
        return ["all"] * len(observations)
    raise ValueError(f"unknown grouping {grouping}")


def posterior_by_grade(posterior: np.ndarray, observations: pd.DataFrame) -> dict:
    """Median and 5-95 percent band of the posterior per grade, and how far apart the label grades and the rest sit."""
    result = {}
    for grade, rows in observations.groupby("grade").indices.items():
        values = posterior[rows]
        result[str(grade)] = {"items": int(len(values)), "median": float(np.median(values)), "p5": float(np.percentile(values, 5)), "p95": float(np.percentile(values, 95))}
    return result


def fit_under_grouping(grouping: str, observations: pd.DataFrame, reports: pd.DataFrame, defaults: RubricWeightDefaults) -> dict:
    """The assembler's fit (assemble_evidence_table.fit_report_reliability) with the item grouping swapped."""
    item_keys = list(zip(observations.perturbation_id, observations.symptom, observations.relation))
    present = set(reports.source)
    source_names = [name for name in SOURCE_NAMES if name in present] + sorted(present - set(SOURCE_NAMES))
    positive_weights, negative_weights, coverage = report_count_matrices(reports, item_keys, source_names, defaults.implicit_negative_weight)
    explicit = explicit_report_cells(reports, item_keys, source_names)
    explicit_negative_weight = (negative_weights - defaults.implicit_negative_weight * coverage * (~explicit)).sum(axis=0)
    groups = item_groups_for(grouping, observations, coverage, source_names)
    fit = fit_weighted_dawid_skene(positive_weights, negative_weights, source_names, item_keys, item_groups=groups, estimate_specificity=explicit_negative_weight > 0)
    group_sizes = pd.Series(groups).value_counts().to_dict()
    return {
        "grouping": grouping,
        "groups": len(group_sizes),
        "converged": bool(fit.converged),
        "iterations": int(fit.num_iterations),
        "degenerate_sources": list(fit.degenerate_sources),
        "prevalence_pooled": float(fit.prevalence),
        "prevalence_by_group": {group: {"items": int(group_sizes[group]), "prevalence": float(value)} for group, value in sorted(fit.prevalence_by_group.items())},
        "per_source": {name: {"sensitivity": float(fit.sensitivity[column]), "specificity": float(fit.specificity[column]),
                              "specificity_held": fit.specificity_hold_reasons.get(name)} for column, name in enumerate(source_names)},
        "posterior_by_grade": posterior_by_grade(fit.posterior, observations),
    }


def render(results: list[dict]) -> str:
    lines = [f"{'grouping':18s} {'groups':>6s} {'converged':>9s} {'iterations':>10s}  degenerate", ""]
    for result in results:
        lines.append(f"{result['grouping']:18s} {result['groups']:6d} {str(result['converged']):>9s} {result['iterations']:10d}  {', '.join(result['degenerate_sources']) or '-'}")
    lines.append("")
    sources = list(results[0]["per_source"])
    lines.append(f"{'source':26s} " + " ".join(f"{result['grouping'][:17]:>19s}" for result in results))
    for source in sources:
        cells = []
        for result in results:
            entry = result["per_source"][source]
            held = "" if entry["specificity_held"] is None else "*"
            cells.append(f"{entry['sensitivity']:.3f} / {entry['specificity']:.3f}{held:1s}".rjust(19))
        lines.append(f"{source:26s} " + " ".join(cells))
    lines.append("(sensitivity / specificity; * specificity held at the prior mean)")
    lines.append("")
    for result in results:
        lines.append(f"{result['grouping']}: posterior median by grade " + ", ".join(f"{grade} {entry['median']:.3f} (n={entry['items']})" for grade, entry in sorted(result["posterior_by_grade"].items())))
        if result["groups"] <= 12:
            lines.append("  prevalence by group: " + ", ".join(f"{group} {entry['prevalence']:.4f} (n={entry['items']})" for group, entry in result["prevalence_by_group"].items()))
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, required=True, help="a directory holding evidence_records.parquet and evidence_reports.parquet, assembled with --extra-reports")
    parser.add_argument("--output", type=Path, default=None, help="JSON of the four fits")
    parser.add_argument("--groupings", nargs="*", default=list(GROUPINGS), choices=GROUPINGS)
    arguments = parser.parse_args()
    observations = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet")
    reports = pd.read_parquet(arguments.evidence_dir / "evidence_reports.parquet")
    defaults = RubricWeightDefaults()
    results = [fit_under_grouping(grouping, observations, reports, defaults) for grouping in arguments.groupings]
    print(render(results))
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps({"evidence_dir": str(arguments.evidence_dir), "observations": len(observations), "reports": len(reports), "fits": results}, indent=1) + "\n", encoding="utf-8")
        print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
