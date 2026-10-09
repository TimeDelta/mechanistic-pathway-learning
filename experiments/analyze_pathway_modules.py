"""Interpret the learned pathway modules of a noisy-OR run (design sections 6.5 and 6.6).

For one run directory of experiments/run_main_model.py (noisy-OR head), per split:
  - module supports: nodes whose evaluation gate exceeds --support-threshold, annotated with node type,
    display name, compartment and subsystem; top nodes per module by gate value;
  - module-to-symptom links above --link-threshold and the active module set M_s per symptom;
  - equifinality indices per symptom: independence (1 - max pairwise Jaccard of supports) and convergence
    (mean pairwise Jaccard of downstream node sets within --downstream-hops), from
    equifinality_independence_and_convergence.py;
  - overlap of each module support with the curated modules (gene nodes of docs/curated_pathway_modules.csv)
    and the Human-GEM subsystems of the reactions in the support;
  - a name per module, from its support (module_label) and from its behaviour (module_behaviour.py: the held-out
    perturbations that switch it on, how far they are from its support and the symptoms it feeds more than the
    other modules do);
  - stability across splits: for every pair of modules from different splits, Jaccard of supports;
  - the sufficiency test of design section 6.5, computed post hoc from the saved test module activations: for
    symptom s, P_k is the set of held-out perturbations whose largest contribution link_{k,s} * a_k comes from
    module k; the symptom probability is recomputed with one module removed from the noisy-OR; ablating k should
    lower AUPRC on P_k and ablating any other module should not (paired bootstrap over perturbations).

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

from sklearn.metrics import average_precision_score

from mechanistic_pathway_learning.evaluation.equifinality_independence_and_convergence import (
    convergence_index,
    independence_index,
    jaccard_index,
)
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.module_behaviour import (
    behaviour_label,
    hops_into_support,
    module_drivers,
    module_symptom_side,
    propagation_hops,
    reverse_adjacency,
    symptom_rdoc_domains,
)
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import read_curated_modules


def noisy_or_probabilities(activations: np.ndarray, links: np.ndarray, leaks: np.ndarray, ablated_module: int | None = None) -> np.ndarray:
    """P(symptom) from activations [n, K], links [K, S] and leaks [S] or [n, S] (degree-dependent), with one module removed when given."""
    contributions = activations[:, :, None] * links[None, :, :]  # [n, K, S]
    if ablated_module is not None:
        contributions = np.delete(contributions, ablated_module, axis=1)
    leak_rows = leaks if leaks.ndim == 2 else leaks[None, :]
    return 1.0 - (1.0 - leak_rows) * np.prod(1.0 - np.clip(contributions, 0.0, 1.0 - 1e-6), axis=1)


def sufficiency_test(activations: np.ndarray, links: np.ndarray, leaks: np.ndarray, outcomes: np.ndarray, symptoms: list[str], min_group_size: int = 8, min_positives: int = 3, num_bootstrap: int = 500, seed: int = 0) -> list[dict]:
    """Design section 6.5: per (symptom, dominant module) group, AUPRC with the full model, with the own module ablated and with each other module ablated."""
    generator = np.random.default_rng(seed)
    num_modules = links.shape[0]
    full = noisy_or_probabilities(activations, links, leaks)
    ablated = [noisy_or_probabilities(activations, links, leaks, k) for k in range(num_modules)]
    contributions = activations[:, :, None] * links[None, :, :]
    dominant = contributions.argmax(axis=1)  # [n, S]
    rows_out = []
    for s, symptom in enumerate(symptoms):
        for k in range(num_modules):
            group = np.where(dominant[:, s] == k)[0]
            positives = outcomes[group, s].sum() if len(group) else 0
            if len(group) < min_group_size or positives < min_positives or positives == len(group):
                continue
            auprc_full = average_precision_score(outcomes[group, s], full[group, s])
            auprc_own = average_precision_score(outcomes[group, s], ablated[k][group, s])
            others = [average_precision_score(outcomes[group, s], ablated[other][group, s]) for other in range(num_modules) if other != k]
            differences = []
            for _ in range(num_bootstrap):
                sample = generator.choice(group, size=len(group), replace=True)
                if outcomes[sample, s].sum() in (0, len(sample)):
                    continue
                differences.append(average_precision_score(outcomes[sample, s], full[sample, s]) - average_precision_score(outcomes[sample, s], ablated[k][sample, s]))
            lower, upper = (float(np.quantile(differences, 0.025)), float(np.quantile(differences, 0.975))) if differences else (float("nan"), float("nan"))
            rows_out.append({"symptom": symptom, "module": f"module_{k}", "group_size": int(len(group)), "positives": int(positives), "auprc_full": float(auprc_full),
                             "auprc_ablate_own": float(auprc_own), "auprc_ablate_others_mean": float(np.mean(others)) if others else float("nan"),
                             "own_ablation_drop": float(auprc_full - auprc_own), "own_ablation_drop_ci": [lower, upper],
                             "mean_contribution": float(contributions[group, k, s].mean())})
    return rows_out


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


# A module is a learned gate, so a name for it has to be read off the support it ends up holding rather than assigned.
# The rules below take the first name the support actually carries, and every name states the share it rests on, so the
# name cannot claim more than the content behind it. The dictionary key stays module_<k>: the index is what identifies
# a gate inside a split and what the sufficiency test, the links and the stability matrix join on. Two splits whose
# modules take the same name are not thereby the same module, because the gates are fitted per split; whether they are
# is what the cross-split Jaccard in the stability section measures (the user, 9 October 2026: "can you please
# identify better names for the modules than a simple indexing?").
CURATED_GENES_FOR_A_NAME = 3  # below this a curated-module overlap is a coincidence of one or two genes
SUBSYSTEM_SHARE_FOR_A_NAME = 0.25
# A raw share is not enough on its own: Transport reactions is the largest subsystem of the graph, so a module holding
# a quarter of its support in it may hold no more than the graph does. The subsystem also has to be this many times
# over-represented against its share of the graph, when that background is known.
SUBSYSTEM_ENRICHMENT_FOR_A_NAME = 2.0
COMPARTMENT_SHARE_FOR_A_NAME = 0.5
NODE_TYPE_SHARE_FOR_A_NAME = 0.75
# Human-GEM's own compartment names (data/raw/Human-GEM/model/Human-GEM.xml, listOfCompartments), plus the vesicle
# compartment the Reactome import adds (mechanistic_pathway_learning/graph/reactome_import.py: "vesicle lumens become v")
COMPARTMENT_NAMES = {"c": "cytosol", "e": "extracellular", "l": "lysosome", "r": "endoplasmic reticulum",
                     "m": "mitochondria", "x": "peroxisome", "n": "nucleus", "g": "Golgi apparatus",
                     "i": "inner mitochondria", "v": "vesicle"}


def subsystem_shares_of_graph(nodes: pd.DataFrame) -> dict[str, float]:
    """Each subsystem's share of the nodes that carry one, the background a module's share is judged against."""
    subsystems = nodes.subsystem.dropna().astype(str)
    subsystems = subsystems[subsystems != ""]
    return (subsystems.value_counts() / len(subsystems)).to_dict() if len(subsystems) else {}


def module_label(description: dict, subsystem_background: dict[str, float] | None = None) -> str:
    """A readable name for one module, derived from its support and carrying the share it rests on.

    Without `subsystem_background` the subsystem rule judges a raw share, which the largest subsystem of the graph can
    reach without being characteristic of the module; with it the subsystem also has to be over-represented.
    """
    size = description.get("size", 0)
    if not size:
        return "empty support"
    curated = description.get("curated_module_overlap", {}) or {}
    if curated:
        module, genes = max(curated.items(), key=lambda item: item[1])
        if genes >= CURATED_GENES_FOR_A_NAME:
            return f"{module} ({genes} of its genes, {size} nodes)"
    subsystems = description.get("top_subsystems", []) or []
    if subsystems and subsystems[0][1] >= SUBSYSTEM_SHARE_FOR_A_NAME * size:
        subsystem, count = subsystems[0]
        background = (subsystem_background or {}).get(subsystem)
        enriched = background is None or count / size >= SUBSYSTEM_ENRICHMENT_FOR_A_NAME * background
        if not enriched:
            return (f"mixed support ({size} nodes; its largest subsystem, {subsystem.lower()}, holds "
                    f"{count / size:.0%} against {background:.0%} of the graph)")
        name = f"{subsystem.lower()} ({count} of {size} nodes)"
        compartments = description.get("by_compartment", {}) or {}
        if compartments:
            compartment, there = max(compartments.items(), key=lambda item: item[1])
            if there >= COMPARTMENT_SHARE_FOR_A_NAME * size:
                name = f"{subsystem.lower()} in the {COMPARTMENT_NAMES.get(compartment, compartment)} ({count} of {size} nodes)"
        return name
    node_types = description.get("by_node_type", {}) or {}
    if node_types:
        node_type, count = max(node_types.items(), key=lambda item: item[1])
        if count >= NODE_TYPE_SHARE_FOR_A_NAME * size:
            return f"{node_type} nodes across subsystems ({count} of {size} nodes)"
    return f"mixed support ({size} nodes, no subsystem above {SUBSYSTEM_SHARE_FOR_A_NAME:.0%})"


def behaviour_of_module(module_index: int, module_support: set[str], activations: np.ndarray | None, test_rows: np.ndarray | None,
                        data, links: np.ndarray, symptoms: list[str], rdoc_domains: dict[str, list[str]],
                        curated_driver_genes: dict[str, set[str]], reverse_edges: dict[str, list[str]], hop_budget: int) -> dict:
    """One module's behaviour reading: who switches it on, how far they are from its support, what it feeds."""
    if activations is None or test_rows is None:
        return {"available": False, "reason": "no saved module activations, so no driver reading"}
    labels = [data.perturbation_labels[row] for row in test_rows]
    types = [data.perturbation_types[row] for row in test_rows]
    hops_of_perturbation, truncated = None, False
    if module_support and hop_budget:
        distance, truncated = hops_into_support(module_support, reverse_edges, hop_budget)
        hops_of_perturbation = {}
        for position, row in enumerate(test_rows):
            reached = [distance[data.node_ids[index]] for index in data.perturbation_seeds[row]
                       if data.node_ids[index] in distance]
            if reached:
                hops_of_perturbation[position] = min(reached)
    return {"available": True, "reach_truncated": truncated,
            "drivers": module_drivers(activations[:, module_index], labels, types, curated_driver_genes,
                                      hops_of_perturbation, hop_budget),
            "symptom_side": module_symptom_side(links, module_index, symptoms, rdoc_domains)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--graph-dir", type=Path, default=None, help="default: the graph directory the run was trained on (its results.json arguments)")
    parser.add_argument("--curated-modules", type=Path, default=Path("docs/curated_pathway_modules.csv"))
    parser.add_argument("--support-threshold", type=float, default=0.5)
    parser.add_argument("--link-threshold", type=float, default=0.5)
    parser.add_argument("--downstream-hops", type=int, default=3)
    parser.add_argument("--top-k", type=int, default=15)
    parser.add_argument("--minimum-support-nodes", type=int, default=25, help="when fewer gates exceed the threshold, the top-ranked gates up to this many form the support (reported as such)")
    parser.add_argument("--evidence-dir", type=Path, default=None, help="default: the evidence directory the run was trained on")
    parser.add_argument("--markdown-output", type=Path, default=None)
    parser.add_argument("--symptom-crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"),
                        help="the RDoC domain per symptom, which groups the symptoms a module feeds into a family")
    arguments = parser.parse_args()
    split_directories = sorted(path for path in arguments.run_dir.iterdir() if path.is_dir() and (path / "module_support.npy").exists() and (path / "results.json").exists())
    # the module supports index the nodes of the graph the run was trained on; reading another graph's node list
    # (the old fixed default, data/processed/graph, for a graph_neuronal run) misnames nodes or fails
    recorded_arguments = [json.loads((path / "results.json").read_text()).get("arguments", {}) for path in split_directories]
    trained_on = {(recorded.get("graph_dir"), recorded.get("evidence_dir")) for recorded in recorded_arguments}
    if len(trained_on) > 1:
        raise SystemExit(f"the splits of {arguments.run_dir} were trained on different graph or evidence directories: {sorted(trained_on, key=str)}")
    recorded_graph_dir, recorded_evidence_dir = next(iter(trained_on)) if trained_on else (None, None)
    for name, recorded in (("graph_dir", recorded_graph_dir), ("evidence_dir", recorded_evidence_dir)):
        if getattr(arguments, name) is None:
            setattr(arguments, name, Path(recorded or ("data/processed/graph" if name == "graph_dir" else "data/processed/evidence")))
        elif recorded is not None and Path(recorded) != getattr(arguments, name):
            print(f"warning: --{name.replace('_', '-')} {getattr(arguments, name)} differs from the run's {recorded}")
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir)
    position_of = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet")
    node_ids = list(nodes.node_id)
    currency = set(nodes.loc[nodes.is_currency.fillna(False).astype(bool), "node_id"])
    adjacency = directed_adjacency(edges, currency)
    curated_symbols = read_curated_modules(arguments.curated_modules)
    curated = {module: {f"GENE:{symbol}" for symbol in symbols} for module, symbols in curated_symbols.items()}
    curated_driver_genes = {module: set(symbols) for module, symbols in curated_symbols.items()}
    subsystem_background = subsystem_shares_of_graph(nodes)
    reverse_edges = reverse_adjacency(edges.source_id.to_numpy(), edges.target_id.to_numpy(), currency)
    rdoc_domains = symptom_rdoc_domains(pd.read_csv(arguments.symptom_crosswalk))

    analysis = {"run_dir": str(arguments.run_dir), "graph_dir": str(arguments.graph_dir), "splits": {}, "stability": {}}
    all_supports: list[tuple[str, int, set[str]]] = []
    for split_directory in split_directories:
        results, support = load_split(split_directory)
        if support.shape[1] != len(node_ids):
            raise SystemExit(f"{split_directory}: module supports cover {support.shape[1]} nodes but {arguments.graph_dir} has {len(node_ids)}")
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
        activations_path = split_directory / "test_module_activations.npy"
        activations, test_rows = None, None
        if activations_path.exists():
            saved = np.load(activations_path)
            rows = np.array([position_of[p] for p in results["test_perturbation_ids"] if p in position_of])
            if len(rows) == saved.shape[0]:
                activations, test_rows = saved, rows
        hop_budget = propagation_hops(results.get("arguments", {}))
        sufficiency_rows = []
        if activations is not None and "symptom_leaks" in results:
            leaks_path = split_directory / "test_leaks.npy"  # per-perturbation leaks of a degree-offset run
            leaks = np.load(leaks_path) if leaks_path.exists() else np.array(results["symptom_leaks"])
            sufficiency_rows = sufficiency_test(activations, links, leaks, data.outcomes[test_rows], symptoms)
        split_entry = {"sufficiency_test": sufficiency_rows, "modules": {}, "symptoms": {}, "gate_summary": {f"module_{k}": {"max": float(support[k].max()), "mean": float(support[k].mean()), "above_threshold": int((support[k] > arguments.support_threshold).sum())} for k in range(support.shape[0])}}
        for module_name, module_support in supports.items():
            k = int(module_name.split("_")[1])
            description = describe_support(sorted(module_support), nodes, gate_values[module_name], arguments.top_k) if module_support else {"size": 0}
            description["support_rule"] = support_rule[module_name]
            description["links_above_threshold"] = {symptoms[s]: round(float(links[k, s]), 3) for s in range(links.shape[1]) if links[k, s] > arguments.link_threshold}
            description["links"] = {symptoms[s]: round(float(links[k, s]), 3) for s in range(links.shape[1])}
            description["curated_module_overlap"] = {module: len(module_support & genes) for module, genes in curated.items() if module_support & genes}
            description["label"] = module_label(description, subsystem_background)
            description["behaviour"] = behaviour_of_module(k, module_support, activations, test_rows, data, links, symptoms,
                                                           rdoc_domains, curated_driver_genes, reverse_edges, hop_budget)
            description["behaviour_label"] = behaviour_label(description["behaviour"])
            split_entry["modules"][module_name] = description
            all_supports.append((split_directory.name, k, module_support))
        driver_sets = {name: set(entry.get("behaviour", {}).get("drivers", {}).get("drivers", []))
                       for name, entry in split_entry["modules"].items()
                       if entry.get("behaviour", {}).get("available") and entry["behaviour"]["drivers"].get("responds")}
        driver_overlaps = [jaccard_index(left, right) for (left, right) in combinations(driver_sets.values(), 2) if left and right]
        split_entry["driver_overlap"] = {"responding_modules": len(driver_sets),
                                         "median_pairwise_jaccard": float(np.median(driver_overlaps)) if driver_overlaps else float("nan")}
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
             f"Support threshold {arguments.support_threshold} on the evaluation gate saved by the run; link threshold {arguments.link_threshold}; downstream hops {arguments.downstream_hops}; currency metabolites removed from the downstream walk.", "",
             "Each module keeps its index, which is what the tables below join on, and carries two names. The first is "
             "read off its own support: the curated module it overlaps by at least "
             f"{CURATED_GENES_FOR_A_NAME} genes, else the subsystem holding at least {SUBSYSTEM_SHARE_FOR_A_NAME:.0%} of "
             f"the support and at least {SUBSYSTEM_ENRICHMENT_FOR_A_NAME:g} times its share of the graph (with the "
             "compartment when one holds half), else the node type holding at least "
             f"{NODE_TYPE_SHARE_FOR_A_NAME:.0%}, else \"mixed support\". The second is read off its behaviour "
             "(mechanistic_pathway_learning/evaluation/module_behaviour.py): the held-out perturbations in the top "
             "fifth of its activation, the group they share, how many of them have a directed path into its support "
             "no longer than the propagation the encoder runs, and the symptoms it links to more strongly than the "
             "other modules do, grouped by the crosswalk's RDoC domain. The share is printed with each name so the "
             "name claims no more than its support, and a reading that does not concentrate gives the number behind "
             "its refusal instead of a name. Two splits whose modules take the same name need not be the same module: "
             "the gates are fitted per split, and the stability section is where that is judged.", ""]
    for split_name, split_entry in analysis["splits"].items():
        lines += [f"## {split_name}", "", "| module | what its support holds | how it behaves | support rule | gates above threshold | node types | top subsystems | links above threshold | curated overlap |", "|---|---|---|---|---|---|---|---|---|"]
        for module_name, description in split_entry["modules"].items():
            gate_summary = split_entry["gate_summary"][module_name]
            lines.append(f"| {module_name} | {description.get('label', '')} | {description.get('behaviour_label', '')} | {description['support_rule']} | {gate_summary['above_threshold']} | {description.get('by_node_type', {})} | {[name for name, _ in description.get('top_subsystems', [])][:3]} | {description['links_above_threshold']} | {description['curated_module_overlap']} |")
        overlap = split_entry["driver_overlap"]
        lines += ["", f"Modules whose activation varies across the held-out perturbations: {overlap['responding_modules']} of "
                  f"{len(split_entry['modules'])}; median pairwise Jaccard of their driver sets: "
                  f"{overlap['median_pairwise_jaccard']:.2f} (near 1 the same perturbations drive every module, so a "
                  "driver name does not separate them).", ""]
        lines += ["| symptom | active modules | independence index | convergence index |", "|---|---|---|---|"]
        for symptom, entry in split_entry["symptoms"].items():
            lines.append(f"| {symptom} | {entry['active_modules']} | {entry['independence_index']:.2f} | {entry['convergence_index']:.2f} |")
        lines.append("")
        if split_entry["sufficiency_test"]:
            lines += ["Sufficiency test (design section 6.5): held-out perturbations grouped by the module with the largest contribution to the symptom; AUPRC with the full noisy-OR, with that module ablated and with each other module ablated (mean); 95 percent paired bootstrap interval of the own-ablation drop.", "",
                      "| symptom | dominant module | group size | positives | AUPRC full | ablate own | ablate others (mean) | own-ablation drop [95% CI] | mean contribution |", "|---|---|---|---|---|---|---|---|---|"]
            for row in split_entry["sufficiency_test"]:
                lines.append(f"| {row['symptom']} | {row['module']} | {row['group_size']} | {row['positives']} | {row['auprc_full']:.3f} | {row['auprc_ablate_own']:.3f} | {row['auprc_ablate_others_mean']:.3f} | {row['own_ablation_drop']:+.3f} [{row['own_ablation_drop_ci'][0]:+.3f}, {row['own_ablation_drop_ci'][1]:+.3f}] | {row['mean_contribution']:.3f} |")
            lines.append("")
        for module_name, description in split_entry["modules"].items():
            if description["size"] == 0:
                continue
            lines += [f"### {split_name} {module_name}, {description.get('label', '')}: top nodes", "",
                      f"How it behaves: {description.get('behaviour_label', '')}", "",
                      "| node | type | name | compartment | subsystem | gate |", "|---|---|---|---|---|---|"]
            for row in description["top_nodes"]:
                lines.append(f"| {row['node_id']} | {row['node_type']} | {row['name'][:60]} | {row['compartment']} | {row['subsystem'][:40]} | {row['gate']:.2f} |")
            lines.append("")
    lines += ["## Stability across splits", "", f"Non-empty supports: {analysis['stability']['num_nonempty_supports']}; mean best cross-split Jaccard: {analysis['stability']['mean_best_match_jaccard']:.2f}", "",
              "| support A | support B | Jaccard |", "|---|---|---|"] + [f"| {p['a']} | {p['b']} | {p['jaccard']:.2f} |" for p in analysis["stability"]["best_cross_split_matches"]]
    markdown_output.write_text("\n".join(lines) + "\n")
    print(f"wrote {markdown_output}; splits analysed: {len(analysis['splits'])}; non-empty supports: {analysis['stability']['num_nonempty_supports']}")


if __name__ == "__main__":
    main()
