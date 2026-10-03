"""Interpret the learned pathway modules of a noisy-OR run (design sections 6.5 and 6.6).

For one run directory of experiments/run_main_model.py (noisy-OR head), per split:
  - module supports: nodes whose deterministic gate exceeds --support-threshold, annotated with node type,
    display name, compartment and subsystem; top nodes per module by gate value;
  - module-to-symptom links above --link-threshold and the active module set M_s per symptom;
  - equifinality indices per symptom: independence (1 - max pairwise Jaccard of supports) and convergence
    (mean pairwise Jaccard of downstream node sets within --downstream-hops), from
    equifinality_independence_and_convergence.py;
  - overlap of each module support with the curated modules (gene nodes of docs/curated_pathway_modules.csv)
    and the Human-GEM subsystems of the reactions in the support;
  - stability across splits: for every pair of modules from different splits, Jaccard of supports; modules are
    matched greedily and the selection frequency of each consensus support is reported.

Writes <run-dir>/module_analysis.json and docs/<name>_modules.md (one card per symptom).
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.equifinality_independence_and_convergence import (
    convergence_index,
    independence_index,
    jaccard_index,
)
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import read_curated_modules


def load_split(split_directory: Path) -> tuple[dict, np.ndarray]:
    results = json.loads((split_directory / "results.json").read_text())
    support = np.load(split_directory / "module_support.npy")
    return results, support


def directed_adjacency(edges: pd.DataFrame, exclude_nodes: set[str]) -> dict[str, list[str]]:
    adjacency: dict[str, list[str]] = defaultdict(list)
    for source, target in zip(edges.source_id, edges.target_id):
        if source in exclude_nodes or target in exclude_nodes:
            continue
        adjacency[source].append(target)
    return adjacency


def describe_support(node_ids: list[str], nodes: pd.DataFrame, gate_values: np.ndarray, top_k: int) -> dict:
    table = nodes.set_index("node_id")
    rows = []
    for node_id in node_ids:
        record = table.loc[node_id]
        rows.append({"node_id": node_id, "node_type": record.node_type, "name": str(record.display_name) if pd.notna(record.display_name) else node_id,
                     "compartment": str(record.compartment) if pd.notna(record.get("compartment")) else "", "subsystem": str(record.subsystem) if pd.notna(record.get("subsystem")) else "",
                     "gate": float(gate_values[node_id]) if isinstance(gate_values, dict) else float("nan")})
    rows.sort(key=lambda row: -row["gate"])
    return {
        "size": len(rows),
        "by_node_type": dict(Counter(row["node_type"] for row in rows)),
        "by_compartment": dict(Counter(row["compartment"] for row in rows if row["compartment"])),
        "top_subsystems": Counter(row["subsystem"] for row in rows if row["subsystem"]).most_common(5),
        "top_nodes": rows[:top_k],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--curated-modules", type=Path, default=Path("docs/curated_pathway_modules.csv"))
    parser.add_argument("--support-threshold", type=float, default=0.5)
    parser.add_argument("--link-threshold", type=float, default=0.5)
    parser.add_argument("--downstream-hops", type=int, default=3)
    parser.add_argument("--top-k", type=int, default=15)
    parser.add_argument("--minimum-support-nodes", type=int, default=25, help="when fewer gates exceed the threshold, the top-ranked gates up to this many form the support (reported as such)")
    parser.add_argument("--markdown-output", type=Path, default=None)
    arguments = parser.parse_args()
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet")
    node_ids = list(nodes.node_id)
    currency = set(nodes.loc[nodes.is_currency.fillna(False).astype(bool), "node_id"])
    adjacency = directed_adjacency(edges, currency)
    curated = {module: {f"GENE:{symbol}" for symbol in symbols} for module, symbols in read_curated_modules(arguments.curated_modules).items()}

    split_directories = sorted(path for path in arguments.run_dir.iterdir() if path.is_dir() and (path / "module_support.npy").exists() and (path / "results.json").exists())
    analysis = {"run_dir": str(arguments.run_dir), "splits": {}, "stability": {}}
    all_supports: list[tuple[str, int, set[str]]] = []
    for split_directory in split_directories:
        results, support = load_split(split_directory)
        symptoms = results["symptoms"]
        links = np.array(results["module_symptom_links"])  # [K, S]
        supports, gate_values, support_rule = {}, {}, {}
        for k in range(support.shape[0]):
            above = np.where(support[k] > arguments.support_threshold)[0]
            if len(above) < arguments.minimum_support_nodes:
                above = np.argsort(-support[k])[: arguments.minimum_support_nodes]
                support_rule[f"module_{k}"] = f"top {arguments.minimum_support_nodes} gates (max gate {support[k].max():.2f})"
            else:
                support_rule[f"module_{k}"] = f"gate > {arguments.support_threshold}"
            supports[f"module_{k}"] = {node_ids[i] for i in above}
            gate_values[f"module_{k}"] = {node_ids[i]: float(support[k, i]) for i in above}
        split_entry = {"modules": {}, "symptoms": {}, "gate_summary": {f"module_{k}": {"max": float(support[k].max()), "mean": float(support[k].mean()), "above_threshold": int((support[k] > arguments.support_threshold).sum())} for k in range(support.shape[0])}}
        for module_name, module_support in supports.items():
            k = int(module_name.split("_")[1])
            description = describe_support(sorted(module_support), nodes, gate_values[module_name], arguments.top_k) if module_support else {"size": 0}
            description["support_rule"] = support_rule[module_name]
            description["links_above_threshold"] = {symptoms[s]: round(float(links[k, s]), 3) for s in range(links.shape[1]) if links[k, s] > arguments.link_threshold}
            description["links"] = {symptoms[s]: round(float(links[k, s]), 3) for s in range(links.shape[1])}
            description["curated_module_overlap"] = {module: len(module_support & genes) for module, genes in curated.items() if module_support & genes}
            split_entry["modules"][module_name] = description
            all_supports.append((split_directory.name, k, module_support))
        for s, symptom in enumerate(symptoms):
            active = [f"module_{k}" for k in range(links.shape[0]) if links[k, s] > arguments.link_threshold and supports[f"module_{k}"]]
            split_entry["symptoms"][symptom] = {
                "active_modules": active, "links": {module: round(float(links[int(module.split('_')[1]), s]), 3) for module in active},
                "independence_index": independence_index(supports, active), "convergence_index": convergence_index(supports, active, adjacency, arguments.downstream_hops) if len(active) >= 2 else float("nan"),
            }
        analysis["splits"][split_directory.name] = split_entry
    # stability: match supports across splits by Jaccard
    pairwise = []
    for (split_a, k_a, support_a), (split_b, k_b, support_b) in combinations(all_supports, 2):
        if split_a == split_b or not support_a or not support_b:
            continue
        pairwise.append({"a": f"{split_a}/module_{k_a}", "b": f"{split_b}/module_{k_b}", "jaccard": jaccard_index(support_a, support_b)})
    pairwise.sort(key=lambda entry: -entry["jaccard"])
    analysis["stability"] = {"num_nonempty_supports": sum(1 for _, _, s in all_supports if s), "best_cross_split_matches": pairwise[:20],
                             "mean_best_match_jaccard": float(np.mean([max([p["jaccard"] for p in pairwise if p["a"].startswith(split) or p["b"].startswith(split)] or [0.0]) for split in {s for s, _, _ in all_supports}])) if pairwise else float("nan")}
    (arguments.run_dir / "module_analysis.json").write_text(json.dumps(analysis, indent=1, default=str))

    markdown_output = arguments.markdown_output or Path("docs") / f"{arguments.run_dir.name}_modules.md"
    lines = [f"# Learned pathway modules: {arguments.run_dir.name} (generated by experiments/analyze_pathway_modules.py)", "",
             f"Support threshold {arguments.support_threshold} on the deterministic gate; link threshold {arguments.link_threshold}; downstream hops {arguments.downstream_hops}; currency metabolites removed from the downstream walk.", ""]
    for split_name, split_entry in analysis["splits"].items():
        lines += [f"## {split_name}", "", "| module | support rule | gates above threshold | node types | top subsystems | links above threshold | curated overlap |", "|---|---|---|---|---|---|---|"]
        for module_name, description in split_entry["modules"].items():
            gate_summary = split_entry["gate_summary"][module_name]
            lines.append(f"| {module_name} | {description['support_rule']} | {gate_summary['above_threshold']} | {description.get('by_node_type', {})} | {[name for name, _ in description.get('top_subsystems', [])][:3]} | {description['links_above_threshold']} | {description['curated_module_overlap']} |")
        lines += ["", "| symptom | active modules | independence index | convergence index |", "|---|---|---|---|"]
        for symptom, entry in split_entry["symptoms"].items():
            lines.append(f"| {symptom} | {entry['active_modules']} | {entry['independence_index']:.2f} | {entry['convergence_index']:.2f} |")
        lines.append("")
        for module_name, description in split_entry["modules"].items():
            if description["size"] == 0:
                continue
            lines += [f"### {split_name} {module_name}: top nodes", "", "| node | type | name | compartment | subsystem | gate |", "|---|---|---|---|---|---|"]
            for row in description["top_nodes"]:
                lines.append(f"| {row['node_id']} | {row['node_type']} | {row['name'][:60]} | {row['compartment']} | {row['subsystem'][:40]} | {row['gate']:.2f} |")
            lines.append("")
    lines += ["## Stability across splits", "", f"Non-empty supports: {analysis['stability']['num_nonempty_supports']}; mean best cross-split Jaccard: {analysis['stability']['mean_best_match_jaccard']:.2f}", "",
              "| support A | support B | Jaccard |", "|---|---|---|"] + [f"| {p['a']} | {p['b']} | {p['jaccard']:.2f} |" for p in analysis["stability"]["best_cross_split_matches"]]
    markdown_output.write_text("\n".join(lines) + "\n")
    print(f"wrote {markdown_output}; splits analysed: {len(analysis['splits'])}; non-empty supports: {analysis['stability']['num_nonempty_supports']}")


if __name__ == "__main__":
    main()
