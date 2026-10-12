"""What a measured concentration costs and what it admits: coverage, and how far the leakage groups merge.

Two questions the owner put on 12 October 2026.

First: "Why does the concentration per target cost 67 drugs? We can still keep the drugs that don't specify and train
them the old way." The cost is coverage, and this script prints the arithmetic: the placeholder of
docs/off_target_scoping.md needs only an affinity at a mechanism gene, while a measured concentration needs a
therapeutic range and a fraction unbound as well.

Second: "Or are you saying because this allows the off targets to include additional mechanisms, which grows the
leakage group?" That is a different cost, and this script measures it. Each scheme admits a different set of
off-targets at each occupancy threshold, and each admitted off-target that is a labelled gene perturbation joins the
drug's leakage group to that gene's. The three joining rules of experiments/scope_off_target_binding.py are reused
unchanged, so the numbers here sit beside the ones that script reports under the placeholder.

Three concentration schemes, over the same affinities and the same thresholds:
  placeholder  C = 9 * K of the drug's best-bound mechanism gene (occupancy of that gene fixed at 0.9)
  measured     C = the unbound concentration at each target's own site (brain when the graph records the gene as
               expressed in a brain tissue and the drug has a brain concentration, plasma otherwise)
  mixed        measured where the drug has one, placeholder for every other drug

Counting only: no label, table, graph or model input is changed.

Usage:
  PYTHONPATH=. python experiments/label_review/measured_occupancy_leakage.py
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, merge_drugs_with_their_targets

EXPERIMENTS = Path(__file__).resolve().parent.parent


def load_script(file_name: str):
    specification = importlib.util.spec_from_file_location(Path(file_name).stem, EXPERIMENTS / file_name)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


scoping = load_script("scope_off_target_binding.py")
occupancy_measurement = load_script("measure_target_occupancy.py")

SCHEMES = ("placeholder", "measured", "mixed")


def unbound_concentrations(exposure: pd.DataFrame) -> dict[str, dict[str, float]]:
    """perturbation id -> the high end of the unbound range in plasma and in brain, in micromolar, where it exists."""
    concentrations = {}
    for perturbation_id, row in exposure.iterrows():
        plasma, brain = row.get("unbound_plasma_high_micromolar"), row.get("unbound_brain_high_micromolar")
        entry = {}
        if not pd.isna(plasma):
            entry["plasma"] = float(plasma)
        if not pd.isna(brain):
            entry["brain"] = float(brain)
        if entry:
            concentrations[perturbation_id] = entry
    return concentrations


def occupancy_under(scheme: str, gene: str, target_p_affinity: float, primary_p_affinity: float | None,
                    concentration: dict[str, float] | None, brain_genes: set[str]) -> float:
    """The occupancy this scheme gives the drug at this gene, or NaN where the scheme has no concentration for it."""
    measured_is_available = concentration is not None
    use_measured = scheme == "measured" or (scheme == "mixed" and measured_is_available)
    if use_measured:
        if not measured_is_available:
            return float("nan")
        site = "brain" if (gene in brain_genes and "brain" in concentration) else "plasma"
        if site not in concentration:
            return float("nan")
        return occupancy_measurement.occupancy_at(concentration[site], target_p_affinity)
    if primary_p_affinity is None:
        return float("nan")
    return scoping.occupancy(primary_p_affinity, target_p_affinity)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", default="v4")
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_review/measured_occupancy_leakage.json"))
    arguments = parser.parse_args()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster",
                                    label_selection=arguments.label_selection)
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet", columns=["node_id", "node_type", "gene_symbol"])
    gene_node_of_symbol = dict(zip(nodes.gene_symbol[nodes.node_type == "gene"], nodes.node_id[nodes.node_type == "gene"]))
    node_index = {node_id: index for index, node_id in enumerate(data.node_ids)}
    brain_genes = occupancy_measurement.brain_expressed_genes(arguments.graph_dir)

    exposure = pd.read_csv(occupancy_measurement.EXPOSURE_DIRECTORY / f"drug_exposure_{arguments.name}.csv").set_index("perturbation_id")
    concentrations = unbound_concentrations(exposure)
    affinities = occupancy_measurement.affinities_by_drug(data, arguments.graph_dir)
    position_of_perturbation = {perturbation_id: position for position, perturbation_id in enumerate(data.perturbation_ids)}

    per_drug = []
    for perturbation_id, entry in affinities.items():
        mechanism_affinities = entry["mechanism"]
        per_drug.append({
            "position": position_of_perturbation[perturbation_id],
            "perturbation_id": perturbation_id,
            "label": entry["label"],
            "primary_p_affinity": max(mechanism_affinities.values()) if mechanism_affinities else None,
            "has_measured_concentration": perturbation_id in concentrations,
            "mechanism_affinities": mechanism_affinities,
            "off_target_affinities": entry["other"],
        })

    coverage = {
        "drugs": len(per_drug),
        "with_an_affinity_at_a_mechanism_gene": sum(entry["primary_p_affinity"] is not None for entry in per_drug),
        "with_a_measured_concentration": sum(entry["has_measured_concentration"] for entry in per_drug),
        "with_both": sum(entry["primary_p_affinity"] is not None and entry["has_measured_concentration"] for entry in per_drug),
    }
    coverage["lost_by_requiring_a_measured_concentration"] = (
        coverage["with_an_affinity_at_a_mechanism_gene"] - coverage["with_both"])
    print("Coverage of a mechanism-target occupancy, over the drugs of the evidence table:")
    for name, value in coverage.items():
        print(f"  {name}: {value}")

    # The owner's proposal: "keep the drugs that don't specify and train them the old way", that is, a measured
    # occupancy for the drugs that have one and a unit input per target for the rest. The objection is that the two
    # groups then carry inputs on different scales, and which group a drug falls in follows how well it was studied.
    # The comparable quantity is a drug's total input magnitude over the same seeds, so both are computed that way:
    # under the unit rule it is the seed count, and under an occupancy rule the sum of each seed's gene's occupancy.
    # A seed whose gene has no affinity has no occupancy, which is a coverage hole inside the drug, counted here too.
    graph_symbol_of_node = dict(zip(nodes.node_id, nodes.gene_symbol))
    total_input_occupancy_rule, total_input_unit_rule, seeds_without_an_affinity = [], [], []
    for entry in per_drug:
        concentration = concentrations.get(entry["perturbation_id"])
        seeds = data.perturbation_seeds[entry["position"]]
        seed_genes = [graph_symbol_of_node.get(data.node_ids[seed]) for seed in seeds]
        per_seed_occupancy = [
            occupancy_under("measured", gene, entry["mechanism_affinities"][gene], entry["primary_p_affinity"],
                            concentration, brain_genes)
            if gene in entry["mechanism_affinities"] else float("nan")
            for gene in seed_genes]
        usable = [value for value in per_seed_occupancy if not np.isnan(value)]
        if usable:
            total_input_occupancy_rule.append(float(sum(usable)))
            seeds_without_an_affinity.append(len(seeds) - len(usable))
        else:
            total_input_unit_rule.append(len(seeds))

    def spread(values):
        return {"drugs": len(values), "median": float(np.median(values)),
                "quartiles": [float(value) for value in np.percentile(values, [25, 75])],
                "max": float(max(values))}

    two_scales = {
        "occupancy_rule_where_there_is_a_concentration": spread(total_input_occupancy_rule),
        "unit_rule_for_the_rest": spread(total_input_unit_rule),
        "seeds_of_a_covered_drug_with_no_affinity_and_so_no_occupancy": spread(seeds_without_an_affinity),
    }
    print("\nThe two scales of a mixed input: a drug's total input magnitude under each rule, over the same seeds")
    for name, value in two_scales.items():
        print(f"  {name}: {value}")

    base_groups = merge_drugs_with_their_targets(data.perturbation_ids, data.perturbation_types, data.group_ids,
                                                 data.perturbation_seeds, data.node_ids)
    kept_positive = ((data.outcomes > 0) & (data.label_mask > 0 if data.label_mask is not None else True)).sum(axis=1)
    gene_perturbation_of_node = {data.node_ids[seeds[0]]: position
                                 for position, (kind, seeds) in enumerate(zip(data.perturbation_types, data.perturbation_seeds))
                                 if kind == "gene" and len(seeds)}

    summary = {"evidence_dir": str(arguments.evidence_dir), "coverage": coverage, "two_scales": two_scales,
               "without_off_targets": scoping.largest_group(base_groups, kept_positive, data.perturbation_types),
               "by_scheme": {}}
    print("\nWithout any off-target, the registered chained rule:",
          json.dumps(summary["without_off_targets"]["largest_by_perturbations"]))

    for scheme in SCHEMES:
        summary["by_scheme"][scheme] = {}
        print(f"\n=== scheme: {scheme}")
        for threshold in scoping.OCCUPANCY_THRESHOLDS:
            admitted_pairs, drugs_with_any, in_graph = 0, 0, 0
            seeds_with_off_targets = list(data.perturbation_seeds)
            extra_nodes_of_drug = {}
            for entry in per_drug:
                concentration = concentrations.get(entry["perturbation_id"])
                kept_genes = [gene for gene, p_affinity in entry["off_target_affinities"].items()
                              if occupancy_under(scheme, gene, p_affinity, entry["primary_p_affinity"],
                                                 concentration, brain_genes) >= threshold]
                admitted_pairs += len(kept_genes)
                drugs_with_any += len(kept_genes) > 0
                kept_nodes = [gene_node_of_symbol[gene] for gene in kept_genes
                              if gene in gene_node_of_symbol and gene_node_of_symbol[gene] in node_index]
                in_graph += len(kept_nodes)
                extra_nodes_of_drug[entry["position"]] = kept_nodes
                extra_positions = [node_index[node] for node in kept_nodes]
                seeds_with_off_targets[entry["position"]] = np.concatenate([
                    np.asarray(data.perturbation_seeds[entry["position"]], dtype=int),
                    np.asarray(extra_positions, dtype=int)])

            every_off_target_joins = merge_drugs_with_their_targets(data.perturbation_ids, data.perturbation_types,
                                                                    data.group_ids, seeds_with_off_targets, data.node_ids)
            genes_only = scoping.join_off_target_genes_only(base_groups, data.perturbation_ids, extra_nodes_of_drug,
                                                            gene_perturbation_of_node)
            drugs_joined_to_a_gene = sum(any(node in gene_perturbation_of_node for node in kept)
                                         for kept in extra_nodes_of_drug.values())
            record = {
                "off_target_pairs": admitted_pairs, "drugs_with_any": drugs_with_any, "gene_node_in_graph": in_graph,
                "drugs_with_an_off_target_that_is_a_labelled_gene": drugs_joined_to_a_gene,
                "every_off_target_joins": scoping.largest_group(every_off_target_joins, kept_positive, data.perturbation_types),
                "off_target_joins_only_its_gene_perturbation": scoping.largest_group(genes_only, kept_positive, data.perturbation_types),
            }
            summary["by_scheme"][scheme][str(threshold)] = record
            largest_chained = record["every_off_target_joins"]["largest_by_perturbations"]
            largest_one_step = record["off_target_joins_only_its_gene_perturbation"]["largest_by_perturbations"]
            print(f"  threshold {threshold}: {admitted_pairs} off-target pairs over {drugs_with_any} drugs, "
                  f"{in_graph} of them graph genes, {drugs_joined_to_a_gene} drugs touch a labelled gene")
            print(f"    largest leakage group, every off-target joins: {largest_chained['perturbations']} perturbations "
                  f"({largest_chained['share_of_perturbations']:.1%}), {largest_chained['drugs']} drugs")
            print(f"    largest leakage group, off-target joins only its own gene: {largest_one_step['perturbations']} "
                  f"perturbations ({largest_one_step['share_of_perturbations']:.1%}), {largest_one_step['drugs']} drugs")

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(summary, indent=1, sort_keys=True, default=str))
    print(f"\nwrote {arguments.output}")


if __name__ == "__main__":
    main()
