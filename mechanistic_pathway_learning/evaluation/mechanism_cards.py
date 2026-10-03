"""Interpretability deliverables (design section 6.6): one mechanism card per symptom.

A card gathers, across the splits of a noisy-OR run analysed by experiments/analyze_pathway_modules.py,
the modules linked to the symptom (link above the threshold), each with its support size, node types,
top subsystems, top nodes, curated-module overlap and the equifinality indices of the symptom in that
split; the sufficiency-test rows for the symptom; and how often a comparable module appeared across
splits (consensus is judged by the cross-split Jaccard matches of the analysis). Cards are markdown
files under <run>/mechanism_cards/.
"""
from __future__ import annotations

import json
from pathlib import Path


def build_mechanism_card(symptom_identifier: str, analysis: dict, top_nodes: int = 10) -> str:
    lines = [f"# Mechanism card: {symptom_identifier}", "", f"Run: {analysis.get('run_dir', '')}", ""]
    stability = analysis.get("stability", {})
    lines += [f"Splits analysed: {len(analysis.get('splits', {}))}; non-empty supports: {stability.get('num_nonempty_supports', 'n/a')}; mean best cross-split Jaccard of supports: {stability.get('mean_best_match_jaccard', float('nan')):.2f}", ""]
    for split_name, split_entry in analysis.get("splits", {}).items():
        symptom_entry = split_entry.get("symptoms", {}).get(symptom_identifier)
        if symptom_entry is None:
            continue
        lines += [f"## {split_name}", "",
                  f"Active modules: {symptom_entry['active_modules'] or 'none above the link threshold'}; links: {symptom_entry['links']}; "
                  f"independence index: {symptom_entry['independence_index']:.2f}; convergence index: {symptom_entry['convergence_index']:.2f}", ""]
        for module_name in symptom_entry["active_modules"]:
            module = split_entry["modules"][module_name]
            lines += [f"### {module_name} ({module.get('support_rule', '')}; support {module['size']} nodes)", "",
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
