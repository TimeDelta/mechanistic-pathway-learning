#!/usr/bin/env python
"""Which drugs act centrally, by the rule the owner approved on 12 October 2026, measured rather than by ATC group.

The rule. A drug is admitted when both hold:

  1. a mechanism target of it is a node of the confirmatory graph that GTEx calls brain-expressed, and
  2. the drug has a measured unbound brain-to-plasma ratio, Kp,uu,brain, at or above --minimum-kpuu-brain.

The threshold is not decoration. Without it the rule tests that a measurement exists rather than what it says, and the
first run of this script admitted atenolol at 0.030, the drug the rule was asked to exclude. 0.1 is the default because
the beta blockers, the class the owner raised, are measured either side of a gap: metoprolol 0.64, propranolol 0.95,
pindolol 0.50 and oxprenolol 0.20 pass, atenolol 0.030, nadolol 0.037 and labetalol 0.054 fail, and nothing measured
sits between 0.054 and 0.20. So the cut is read off the measured distribution of the class rather than chosen, and it
is a property of the drug measured before any outcome is looked at.

Why this rule and not the ATC letter. The list it replaces (configs/admitted_atc_prefixes_draft.csv) admitted a drug
when its ATC group name said "centrally acting", which excluded every beta blocker although ADRB1 is brain-expressed at
12.7 median TPM and ADRB2 at 2.5, and although measured Kp,uu,brain separates metoprolol (0.43 to 0.93) from atenolol
(0.026 to 0.036) by a factor of twenty. The owner, 12 October 2026: "Can't exclude all beta blockers … just ones that
don't act centrally because beta receptors are all throughout the brain", then "The revised rule seems good for the beta
blockers". Both halves of the rule are properties of the drug's pharmacology and of the graph, read before any outcome
is looked at, so the rule cannot select the drug set on the symptoms the study predicts, which is the objection that the
ATC list had to be guarded against by hand (docs/label_source_review.md 9.4 and 28).

Both halves are also already read by the registered occupancy anchor, which turns fraction unbound and Kp,uu,brain into
an unbound brain concentration at a target: a drug the rule admits has a concentration at a real target, and a drug it
turns away has no occupancy to compute. So the rule adds no threshold of its own.

What a verdict means. "admitted" is both halves. "below the threshold" is a measured value under
--minimum-kpuu-brain, a statement about the drug. "target only" is a brain-expressed target with no measured
Kp,uu,brain at all, which is a gap in the exposure tables and not a finding about the drug, and is reported apart so the
owner can decide those one at a time. "no brain-expressed target" is the first half failing, which is a statement about
the drug's pharmacology.

What this rule must not be used for. Only 89 of the 222 drugs already in the table have a measured Kp,uu,brain, so
applied as a filter on the drugs already in it would turn away 148 of them for want of a measurement. Whether a drug
has been measured tracks how well it has been studied, which is the confound the measurement-coverage strata exist for
(evaluation/measurement_coverage_strata.py). So this is an admission rule for drugs the nervous-system filter left out,
and the asymmetry it creates has to be disclosed: the drug set then holds drugs admitted by an ATC rule beside drugs
admitted by measured penetration.

Measurement and a table; no evidence table, graph or model is changed. Idempotent.

Usage:
  python experiments/apply_central_action_rule.py \
      --candidate-evidence-dir data/processed/label_review/evidence_any_atc_cap3 \
      --admitted-evidence-dir data/processed/evidence_full_v4 \
      --graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_drug_exposure_table import brain_exposure_summary, study_drugs

sys.path.insert(0, str(Path(__file__).resolve().parent / "label_review"))
from effect_site_coverage import normalise

OUTPUT_TABLE = Path("configs/centrally_acting_drugs.csv")
OUTPUT_DOCUMENT = Path("docs/central_action_rule.md")
OUTPUT_SUMMARY = Path("data/processed/label_review/central_action_rule.json")
DEFAULT_MINIMUM_KPUU_BRAIN = 0.1
VERDICTS = ("admitted", "below the Kp,uu,brain threshold", "target only, no measured Kp,uu,brain", "no brain-expressed target")


def name_keys(drugs: pd.DataFrame) -> tuple[dict[str, str], list[str]]:
    """(normalised name -> perturbation_id, the names two candidates share).

    effect_site_coverage.study_keys refuses a shared name, which is right for the study table and wrong here: the
    candidate pool is every drug of the label sources, and it holds pairs such as candesartan and its prodrug that the
    nervous-system filter had removed before that check ever ran. A shared name goes to the first candidate and is
    reported, because an exposure value matched through it could belong to either.
    """
    key_to_drug: dict[str, str] = {}
    shared: list[str] = []
    for drug in drugs.itertuples(index=False):
        for name in [drug.perturbation_label, *[name for name in str(drug.rxnorm_names).split(";") if name]]:
            key = normalise(name)
            if not key:
                continue
            if key_to_drug.setdefault(key, drug.perturbation_id) != drug.perturbation_id:
                shared.append(key)
    return key_to_drug, sorted(set(shared))


def brain_expressed_nodes(graph_directory: Path) -> set[str]:
    """The node ids the graph's GTEx annotation calls brain-expressed."""
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    if "brain_expressed" not in nodes.columns:
        raise ValueError(f"{graph_directory} has no brain_expressed column; build it with experiments/build_brain_expression_descriptors.py")
    return set(nodes.node_id[nodes.brain_expressed.fillna(False).astype(bool)])


def target_nodes_of_each_drug(evidence_directory: Path) -> dict[str, list[str]]:
    """perturbation_id -> the graph node ids its mechanism targets seed, from the evidence table's own triples."""
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet",
                              columns=["perturbation_id", "perturbation_type", "perturbation_nodes"])
    drugs = records[records.perturbation_type == "drug"].drop_duplicates("perturbation_id")
    return {row.perturbation_id: [node_id for node_id, _, _ in json.loads(row.perturbation_nodes)] for row in drugs.itertuples(index=False)}


def verdicts(candidate_evidence_dir: Path, admitted_evidence_dir: Path, graph_directory: Path,
             minimum_kpuu_brain: float = DEFAULT_MINIMUM_KPUU_BRAIN) -> pd.DataFrame:
    candidates = study_drugs(candidate_evidence_dir)
    key_to_drug, shared_names = name_keys(candidates)
    brain_exposure = brain_exposure_summary(key_to_drug)
    expressed = brain_expressed_nodes(graph_directory)
    targets_of = target_nodes_of_each_drug(candidate_evidence_dir)
    already_admitted = set(target_nodes_of_each_drug(admitted_evidence_dir))
    protein_of_gene = {}
    gene_to_protein = graph_directory / "gene_to_protein.parquet"
    if gene_to_protein.exists():
        table = pd.read_parquet(gene_to_protein)
        protein_of_gene = table.groupby("gene_node_id", sort=False).protein_node_id.agg(list).to_dict()

    rows = []
    for drug in candidates.itertuples(index=False):
        target_nodes = targets_of.get(drug.perturbation_id, [])
        # a target names a gene node; on a split graph the drug acts on that gene's protein, so either may carry the flag
        passing = sorted({node for node in target_nodes
                          if node in expressed or any(protein in expressed for protein in protein_of_gene.get(node, []))})
        exposure = brain_exposure.get(drug.perturbation_id, {})
        kpuu = exposure.get("kpuu_brain")
        if not passing:
            verdict = VERDICTS[3]
        elif kpuu is None:
            verdict = VERDICTS[2]
        else:
            verdict = VERDICTS[0] if kpuu >= minimum_kpuu_brain else VERDICTS[1]
        rows.append({
            "perturbation_id": drug.perturbation_id, "drug": drug.perturbation_label,
            "already_admitted": drug.perturbation_id in already_admitted, "verdict": verdict,
            "mechanism_targets": ";".join(target_nodes), "brain_expressed_targets": ";".join(passing),
            "kpuu_brain": kpuu, "kpuu_brain_species": exposure.get("kpuu_brain_species", ""),
            "kpuu_brain_values": exposure.get("kpuu_brain_values", 0),
            "kpuu_brain_min": exposure.get("kpuu_brain_min"), "kpuu_brain_max": exposure.get("kpuu_brain_max"),
            "kpuu_brain_files": exposure.get("kpuu_brain_files", ""),
        })
    table = pd.DataFrame(rows).sort_values(["verdict", "drug"], ignore_index=True)
    table.attrs["names_shared_by_two_candidates"] = shared_names
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidate-evidence-dir", type=Path, required=True, help="a table built without the nervous-system filter, the pool the rule chooses from")
    parser.add_argument("--admitted-evidence-dir", type=Path, required=True, help="the table in use, whose drugs are marked already admitted")
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--minimum-kpuu-brain", type=float, default=DEFAULT_MINIMUM_KPUU_BRAIN,
                        help="the measured unbound brain-to-plasma ratio a drug must reach; see the module docstring for where the default comes from")
    parser.add_argument("--output-table", type=Path, default=OUTPUT_TABLE)
    parser.add_argument("--output-document", type=Path, default=OUTPUT_DOCUMENT)
    parser.add_argument("--output-summary", type=Path, default=OUTPUT_SUMMARY)
    arguments = parser.parse_args()

    table = verdicts(arguments.candidate_evidence_dir, arguments.admitted_evidence_dir, arguments.graph_dir, arguments.minimum_kpuu_brain)
    new = table[~table.already_admitted]
    summary = {
        "candidate_evidence_dir": str(arguments.candidate_evidence_dir), "admitted_evidence_dir": str(arguments.admitted_evidence_dir),
        "graph_dir": str(arguments.graph_dir), "minimum_kpuu_brain": arguments.minimum_kpuu_brain,
        "candidates": len(table), "already_admitted": int(table.already_admitted.sum()),
        "already_admitted_with_a_measured_kpuu_brain": int(table.already_admitted.to_numpy().__and__(table.kpuu_brain.notna().to_numpy()).sum()),
        "by_verdict": {verdict: int((table.verdict == verdict).sum()) for verdict in VERDICTS},
        "new_by_verdict": {verdict: int((new.verdict == verdict).sum()) for verdict in VERDICTS},
        "already_admitted_that_the_rule_would_turn_away": sorted(table.drug[table.already_admitted & (table.verdict != VERDICTS[0])]),
        "names_shared_by_two_candidates": table.attrs.get("names_shared_by_two_candidates", []),
    }
    arguments.output_table.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(arguments.output_table, index=False)
    arguments.output_summary.parent.mkdir(parents=True, exist_ok=True)
    arguments.output_summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items()
                      if key not in ("already_admitted_that_the_rule_would_turn_away", "names_shared_by_two_candidates")}, indent=2))
    print(f"already admitted, and the rule would turn away: {len(summary['already_admitted_that_the_rule_would_turn_away'])}")
    write_document(arguments.output_document, table, summary)


def write_document(path: Path, table: pd.DataFrame, summary: dict) -> None:
    new = table[~table.already_admitted]
    provenance = (f"Written by `experiments/apply_central_action_rule.py` over `{summary['candidate_evidence_dir']}` against "
                  f"`{summary['graph_dir']}`. The rule, approved by the owner on 12 October 2026 for the beta blockers: a drug is "
                  f"admitted when a mechanism target of it is a brain-expressed node of the graph and the drug has a measured "
                  f"Kp,uu,brain of at least {summary['minimum_kpuu_brain']}. Nothing is applied by this script; "
                  "`configs/centrally_acting_drugs.csv` carries the evidence per drug.")
    counted = f"{summary['candidates']} candidate drugs, {summary['already_admitted']} of them already in the table in use."
    lines = ["# Which drugs act centrally, by the measured rule", "", provenance, "", counted, "",
             "| verdict | all candidates | not already admitted |", "| --- | --- | --- |"]
    for verdict in VERDICTS:
        lines.append(f"| {verdict} | {summary['by_verdict'][verdict]} | {summary['new_by_verdict'][verdict]} |")
    turned_away = summary["already_admitted_that_the_rule_would_turn_away"]
    lines += ["", "## Drugs already in the table that the rule would turn away", "",
              f"{len(turned_away)} of the {summary['already_admitted']} admitted drugs. The rule is for admitting drugs the "
              "nervous-system filter left out; it is not proposed as a filter on the drugs already in, and this list is "
              "here so the owner can see what applying it both ways would cost." if turned_away else
              "None: every drug already in the table passes the rule.", ""]
    if turned_away:
        lines += ["| drug |", "| --- |"] + [f"| {drug} |" for drug in turned_away] + [""]
    admitted = new[new.verdict == VERDICTS[0]].sort_values("kpuu_brain", ascending=False)
    lines += ["## The drugs the rule would admit", "",
              "| drug | brain-expressed targets | Kp,uu,brain | species | values | source files |", "| --- | --- | --- | --- | --- | --- |"]
    for row in admitted.itertuples(index=False):
        targets = row.brain_expressed_targets.replace("GENE:", "").replace(";", ", ")
        lines.append(f"| {row.drug} | {targets} | {row.kpuu_brain:.3g} ({row.kpuu_brain_min:.3g} to {row.kpuu_brain_max:.3g}) | "
                     f"{row.kpuu_brain_species} | {int(row.kpuu_brain_values)} | {row.kpuu_brain_files} |")
    below = new[new.verdict == VERDICTS[1]].sort_values("kpuu_brain", ascending=False)
    lines += ["", "## Measured and below the threshold", "",
              f"{len(below)} drugs, at a threshold of {summary['minimum_kpuu_brain']}.", "",
              "| drug | brain-expressed targets | Kp,uu,brain |", "| --- | --- | --- |"]
    for row in below.itertuples(index=False):
        lines.append(f"| {row.drug} | {row.brain_expressed_targets.replace('GENE:', '').replace(';', ', ')} | {row.kpuu_brain:.3g} |")
    target_only = new[new.verdict == VERDICTS[2]]
    unmeasured = (f"{len(target_only)} drugs. The second half of the rule is unmet for want of a measurement, which is a gap in "
                  "the exposure tables rather than a finding about the drug, so these are listed and not admitted.")
    lines += ["", "## A brain-expressed target and no measured Kp,uu,brain", "", unmeasured, "",
              "| drug | brain-expressed targets |", "| --- | --- |"]
    for row in target_only.sort_values("drug").itertuples(index=False):
        lines.append(f"| {row.drug} | {row.brain_expressed_targets.replace('GENE:', '').replace(';', ', ')} |")
    lines.append("")
    path.write_text("\n".join(lines))
    print(f"wrote {path} and {OUTPUT_TABLE}")


if __name__ == "__main__":
    main()
