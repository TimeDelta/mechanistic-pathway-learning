"""Interpretability deliverables (design section 6.6): one mechanism card per symptom.

A card gathers, across the splits of a noisy-OR run analysed by experiments/analyze_pathway_modules.py,
the modules linked to the symptom (link above the threshold), each with the name that script reads off its
support (module_label, which also prints the share the name rests on), its support size, node types,
top subsystems, top nodes, curated-module overlap and the equifinality indices of the symptom in that
split; the sufficiency-test rows for the symptom; and how often a comparable module appeared across
splits (consensus is judged by the cross-split Jaccard matches of the analysis). Cards are markdown
files under <run>/mechanism_cards/. When the report table (evidence_reports.parquet) is passed, a card
also lists the supporting evidence behind the symptom with each pair's limitations; without it the card
renders as before.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def evidence_limitations_for_symptom(reports: pd.DataFrame, symptom: str, max_pairs: int | None = None) -> list[str]:
    """One line per perturbation with a report on the symptom: its reports' values and the limitations of the evidence.

    A pair that rests on one disease entry says so, because its several reports are descendant terms of one
    annotation rather than independent studies.
    """
    lines: list[str] = []
    symptom_reports = reports[reports.symptom == symptom]
    for (perturbation_id, relation), pair_reports in symptom_reports.groupby(["perturbation_id", "relation"], sort=False):
        positive = int((pair_reports.report_value == 1).sum())
        negative = int(len(pair_reports) - positive)
        clauses: list[str] = []
        for limitations in pair_reports.limitations:
            for clause in str(limitations).split("; "):
                if clause and clause not in clauses:
                    clauses.append(clause)
        if pair_reports.evidence_class.iloc[0] == "monogenic" and pair_reports.source_record_id.nunique() == 1:
            clauses.append(f"rests on one disease entry ({len(pair_reports)} annotation rows)")
        label = pair_reports.perturbation_label.iloc[0]
        lines.append(f"{label} ({perturbation_id}, {relation}): {positive} positive, {negative} negative report(s) from {', '.join(sorted(set(pair_reports.source)))}; limitations: {'; '.join(clauses) or 'none recorded'}")
        if max_pairs is not None and len(lines) >= max_pairs:
            break
    return lines


def build_mechanism_card(symptom_identifier: str, analysis: dict, top_nodes: int = 10, reports: pd.DataFrame | None = None) -> str:
    lines = [f"# Mechanism card: {symptom_identifier}", "", f"Run: {analysis.get('run_dir', '')}", ""]
    if reports is not None:
        lines += ["## Supporting evidence", "", *(f"- {line}" for line in evidence_limitations_for_symptom(reports, symptom_identifier)), ""]
    stability = analysis.get("stability", {})
    lines += [f"Splits analysed: {len(analysis.get('splits', {}))}; non-empty supports: {stability.get('num_nonempty_supports', 'n/a')}; mean best cross-split Jaccard of supports: {stability.get('mean_best_match_jaccard', float('nan')):.2f}", ""]
    for split_name, split_entry in analysis.get("splits", {}).items():
        symptom_entry = split_entry.get("symptoms", {}).get(symptom_identifier)
        if symptom_entry is None:
            continue
        active = ", ".join(f"{name} ({split_entry['modules'].get(name, {}).get('label', '')})"
                           for name in symptom_entry["active_modules"]) or "none above the link threshold"
        lines += [f"## {split_name}", "",
                  f"Active modules: {active}; links: {symptom_entry['links']}; "
                  f"independence index: {symptom_entry['independence_index']:.2f}; convergence index: {symptom_entry['convergence_index']:.2f}", ""]
        for module_name in symptom_entry["active_modules"]:
            module = split_entry["modules"][module_name]
            lines += [f"### {module_name}, {module.get('label', '')} ({module.get('support_rule', '')}; support {module['size']} nodes)", "",
                      f"Node types: {module.get('by_node_type', {})}; compartments: {module.get('by_compartment', {})}; top subsystems: {module.get('top_subsystems', [])}; curated overlap: {module.get('curated_module_overlap', {})}", "",
                      "| node | type | name | compartment | subsystem | gate |", "|---|---|---|---|---|---|"]
            for row in module.get("top_nodes", [])[:top_nodes]:
                lines.append(f"| {row['node_id']} | {row['node_type']} | {str(row['name'])[:50]} | {row['compartment']} | {str(row['subsystem'])[:40]} | {row['gate']:.2f} |")
            lines.append("")
        sufficiency_rows = [row for row in split_entry.get("sufficiency_test", []) if row["symptom"] == symptom_identifier]
        if sufficiency_rows:
            lines += ["Sufficiency test:", "", "| dominant module | group size | positives | AUPRC full | ablate own | ablate others | own-ablation drop [95% CI] |", "|---|---|---|---|---|---|---|"]
            for row in sufficiency_rows:
                lines.append(f"| {row['module']} | {row['group_size']} | {row['positives']} | {row['auprc_full']:.3f} | {row['auprc_ablate_own']:.3f} | {row['auprc_ablate_others_mean']:.3f} | {row['own_ablation_drop']:+.3f} [{row['own_ablation_drop_ci'][0]:+.3f}, {row['own_ablation_drop_ci'][1]:+.3f}] |")
            lines.append("")
    return "\n".join(lines) + "\n"


def write_mechanism_cards(analysis_path: Path, output_directory: Path | None = None) -> list[Path]:
    """Render one card per symptom from <run>/module_analysis.json; returns the written paths."""
    analysis = json.loads(Path(analysis_path).read_text())
    output_directory = output_directory or Path(analysis_path).parent / "mechanism_cards"
    output_directory.mkdir(parents=True, exist_ok=True)
    symptoms = sorted({symptom for split in analysis.get("splits", {}).values() for symptom in split.get("symptoms", {})})
    written = []
    for symptom in symptoms:
        path = output_directory / f"{symptom}.md"
        path.write_text(build_mechanism_card(symptom, analysis))
        written.append(path)
    return written


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--analysis", type=Path, required=True, help="<run>/module_analysis.json written by experiments/analyze_pathway_modules.py")
    parser.add_argument("--output-dir", type=Path, default=None)
    arguments = parser.parse_args()
    for path in write_mechanism_cards(arguments.analysis, arguments.output_dir):
        print(path)
