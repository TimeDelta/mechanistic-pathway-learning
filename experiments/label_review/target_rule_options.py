"""Alternatives to the limit on a drug's number of targets, with the leakage groups each one leaves.

The owner, 12 October 2026: "I really don't like having a target limit. Can you present another option if leakage
changes are needed?"

What the current rule is, which matters for what can replace it. `has_dominant_target` is an *admission* test, not a
way of choosing targets: a drug with more than --max-drug-targets distinct ChEMBL mechanism targets is dropped from
the study, not truncated. Its own docstring says "Affinity margins are not yet used", so the rule is named after a
dominance it does not measure. That is the opening: a rule that measures dominance can admit the drugs the count
excludes and still seed each drug on the targets that carry its action.

Four rules are compared, all over the same table, the approved changes with no limit (230 drugs against 222):

  no_limit          every drug with a mechanism target, seeded on all of them
  three_targets     the current rule: a drug with at most three ChEMBL mechanism targets, seeded on all of them
  affinity_margin   every drug with a measured affinity, seeded only on the genes within <margin> log units of its
                    best-bound gene. A target 100-fold weaker than the best is not what the drug is doing at a
                    therapeutic dose, and this says so with a measurement instead of a count
  curator_primary   every drug with a target that DrugCentral flags MOA = 1 or Guide to PHARMACOLOGY flags as a
                    primary target, seeded only on those genes

For each the script reports how many drugs are admitted, how many genes each seeds and the leakage groups that
result, under the registered chained rule and under the one-step rule. Counting only: no label, table, graph or model
input is changed, and no evidence table is rebuilt, so the admission counts are over the drugs the no-limit table
already holds.

Usage:
  PYTHONPATH=. python experiments/label_review/target_rule_options.py
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
specification = importlib.util.spec_from_file_location("scope_off_target_binding", EXPERIMENTS / "scope_off_target_binding.py")
scoping = importlib.util.module_from_spec(specification)
specification.loader.exec_module(scoping)

AFFINITY_MARGINS = (0.5, 1.0, 1.5, 2.0)  # log units below the drug's best-bound gene that a target may fall and stay


def measurements_by_drug(data, symbol_of_node) -> dict[int, dict]:
    """position -> {"genes": seeded gene symbols, "affinity": {gene: pK}, "primary": {genes a curator flags}}."""
    drugcentral = scoping.drugcentral_measurements(scoping.DRUGCENTRAL_INTERACTIONS)
    gtopdb = scoping.gtopdb_measurements(scoping.GTOPDB_INTERACTIONS, scoping.GTOPDB_LIGAND_MAPPING)
    drugcentral_id_of_name = drugcentral.groupby("drug_name").drugcentral_id.first().to_dict()
    by_position = {}
    for position, (perturbation_id, kind, label, seeds) in enumerate(
            zip(data.perturbation_ids, data.perturbation_types, data.perturbation_labels, data.perturbation_seeds)):
        if kind != "drug":
            continue
        pubchem_cid, name = scoping.drug_keys(perturbation_id, label)
        from_gtopdb = gtopdb[((gtopdb.pubchem_cid == pubchem_cid) if pubchem_cid else False) | (gtopdb.drug_name == name) | (gtopdb.inn == name)]
        drugcentral_ids = set(from_gtopdb.drugcentral_id.dropna().astype(int)) | (
            {drugcentral_id_of_name[name]} if name in drugcentral_id_of_name else set())
        measurements = pd.concat([drugcentral[drugcentral.drugcentral_id.isin(drugcentral_ids)], from_gtopdb], ignore_index=True)
        seeded_genes = {symbol_of_node.get(data.node_ids[seed]) for seed in seeds} - {None}
        affinity, primary = {}, set()
        if len(measurements):
            best = measurements.dropna(subset=["p_affinity"]).groupby("gene").p_affinity.max().to_dict()
            affinity = {gene: value for gene, value in best.items() if gene in seeded_genes}
            primary = {gene for gene in measurements.gene[measurements.mechanism.fillna(False)] if gene in seeded_genes}
        by_position[position] = {"label": label, "genes": sorted(seeded_genes), "affinity": affinity, "primary": primary}
    return by_position


def seeded_genes_under(rule: str, entry: dict, three_target_drugs: set[str], margin: float) -> list[str] | None:
    """The genes this rule seeds for this drug, or None where the rule does not admit the drug."""
    if rule == "no_limit":
        return entry["genes"]
    if rule == "three_targets":
        return entry["genes"] if entry["label"] in three_target_drugs else None
    if rule == "affinity_margin":
        if not entry["affinity"]:
            return None
        best = max(entry["affinity"].values())
        return sorted(gene for gene, value in entry["affinity"].items() if value >= best - margin)
    if rule == "curator_primary":
        return sorted(entry["primary"]) or None
    if rule == "affinity_margin_or_all":
        # No admission test at all, which is what "no target limit" means: every drug enters, and a drug whose
        # affinities are known is seeded only on the targets within <margin> of its best. A drug with no measured
        # affinity keeps all of its targets, because there is nothing to prune it by.
        if not entry["affinity"]:
            return entry["genes"]
        best = max(entry["affinity"].values())
        pruned = sorted(gene for gene, value in entry["affinity"].items() if value >= best - margin)
        unmeasured = sorted(set(entry["genes"]) - set(entry["affinity"]))
        return pruned + unmeasured
    if rule == "affinity_margin_or_three":
        # The same pruning, with the old count as the fallback for a drug that has no affinity to prune by.
        if not entry["affinity"]:
            return entry["genes"] if entry["label"] in three_target_drugs else None
        best = max(entry["affinity"].values())
        pruned = sorted(gene for gene, value in entry["affinity"].items() if value >= best - margin)
        unmeasured = sorted(set(entry["genes"]) - set(entry["affinity"]))
        return pruned + unmeasured
    raise ValueError(rule)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--no-limit-dir", type=Path, default=Path("data/processed/label_review/evidence_approved_no_cap"))
    parser.add_argument("--three-target-dir", type=Path, default=Path("data/processed/label_review/evidence_approved"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_review/target_rule_options.json"))
    arguments = parser.parse_args()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        data = load_experiment_data(arguments.graph_dir, arguments.no_limit_dir, group_by="disease_cluster")
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet", columns=["node_id", "gene_symbol"])
    symbol_of_node = dict(zip(nodes.node_id, nodes.gene_symbol))
    node_index = {node_id: index for index, node_id in enumerate(data.node_ids)}
    gene_node_of_symbol = {symbol: node for node, symbol in symbol_of_node.items() if isinstance(symbol, str)}

    three_target_records = pd.read_parquet(arguments.three_target_dir / "evidence_records.parquet",
                                           columns=["perturbation_type", "perturbation_label"])
    three_target_drugs = set(three_target_records.perturbation_label[three_target_records.perturbation_type == "drug"])

    by_position = measurements_by_drug(data, symbol_of_node)
    kept_positive = (data.outcomes > 0).sum(axis=1)
    gene_perturbation_of_node = {data.node_ids[seeds[0]]: position
                                 for position, (kind, seeds) in enumerate(zip(data.perturbation_types, data.perturbation_seeds))
                                 if kind == "gene" and len(seeds)}

    print(f"the no-limit table holds {len(by_position)} drugs; the three-target table holds {len(three_target_drugs)}")
    print(f"drugs the count excludes: {sorted(entry['label'] for entry in by_position.values() if entry['label'] not in three_target_drugs)}")

    summary = {}
    rules = ([("no_limit", None), ("three_targets", None)]
             + [("affinity_margin", margin) for margin in AFFINITY_MARGINS]
             + [("affinity_margin_or_all", margin) for margin in AFFINITY_MARGINS]
             + [("affinity_margin_or_three", margin) for margin in AFFINITY_MARGINS]
             + [("curator_primary", None)])
    for rule, margin in rules:
        name = rule if margin is None else f"{rule}_{margin}"
        admitted, seeded_counts, base_groups_seeds = 0, [], list(data.perturbation_seeds)
        dropped_genes = 0
        for position, entry in by_position.items():
            genes = seeded_genes_under(rule, entry, three_target_drugs, margin or 1.0)
            if genes is None:
                base_groups_seeds[position] = np.asarray([], dtype=int)  # not admitted: seeds nothing, joins nothing
                continue
            admitted += 1
            seeded_counts.append(len(genes))
            dropped_genes += len(entry["genes"]) - len(genes)
            positions = [node_index[gene_node_of_symbol[gene]] for gene in genes
                         if gene in gene_node_of_symbol and gene_node_of_symbol[gene] in node_index]
            base_groups_seeds[position] = np.asarray(sorted(set(positions)), dtype=int)

        groups = merge_drugs_with_their_targets(data.perturbation_ids, data.perturbation_types, data.group_ids,
                                                base_groups_seeds, data.node_ids)
        largest = scoping.largest_group(groups, kept_positive, data.perturbation_types)
        summary[name] = {"drugs_admitted": admitted, "genes_seeded_median": float(np.median(seeded_counts)) if seeded_counts else 0.0,
                         "genes_seeded_max": int(max(seeded_counts)) if seeded_counts else 0,
                         "gene_seats_dropped": dropped_genes, "groups": largest}
        by_perturbations = largest["largest_by_perturbations"]
        print(f"\n{name}: {admitted} drugs admitted, median {summary[name]['genes_seeded_median']:.0f} genes seeded "
              f"(max {summary[name]['genes_seeded_max']}), {dropped_genes} gene seats dropped")
        print(f"  leakage groups: {largest['leakage_groups']}; largest holds {by_perturbations['perturbations']} perturbations "
              f"({by_perturbations['share_of_perturbations']:.1%}) and {by_perturbations['drugs']} drugs "
              f"({by_perturbations['share_of_drugs']:.1%} of drugs)")

    # Which drugs hold the largest group together. If a handful of drugs bridge two blocks, then "no limit, with
    # these drugs named and set aside" is an option, and it is a rule that can be written down in advance rather
    # than a count chosen to make the split come out. Measured by removing one drug at a time and asking how large
    # the biggest group is without it.
    for rule, margin in (("no_limit", None), ("affinity_margin_or_all", 0.5)):
        name = rule if margin is None else f"{rule}_{margin}"
        seeds_of_position = {}
        for position, entry in by_position.items():
            genes = seeded_genes_under(rule, entry, three_target_drugs, margin or 1.0)
            seeds_of_position[position] = ([node_index[gene_node_of_symbol[gene]] for gene in (genes or [])
                                            if gene in gene_node_of_symbol and gene_node_of_symbol[gene] in node_index])

        def largest_drug_count(positions_left_out: set[int]) -> int:
            trial_seeds = list(data.perturbation_seeds)
            for position in by_position:
                trial_seeds[position] = np.asarray([] if position in positions_left_out else sorted(set(seeds_of_position[position])), dtype=int)
            trial_groups = merge_drugs_with_their_targets(data.perturbation_ids, data.perturbation_types,
                                                          data.group_ids, trial_seeds, data.node_ids)
            return scoping.largest_group(trial_groups, kept_positive, data.perturbation_types)["largest_by_perturbations"]["drugs"]

        whole = largest_drug_count(set())
        effects = sorted(((whole - largest_drug_count({position}), by_position[position]["label"])
                          for position in by_position), reverse=True)
        print(f"\n{name}: removing one drug at a time, the ten whose removal shrinks the largest group most "
              f"(whole group holds {whole} drugs)")
        for shrink, label in effects[:10]:
            print(f"    {label}: -{shrink} drugs")
        summary.setdefault(name, {})["bridging_drugs"] = [{"drug": label, "drugs_removed_from_largest_group": shrink}
                                                          for shrink, label in effects[:10]]

        # Then the menu: no limit at all, with the named bridges set aside one after another.
        label_of_position = {position: entry["label"] for position, entry in by_position.items()}
        position_of_label = {label: position for position, label in label_of_position.items()}
        set_aside, menu = [], []
        for _shrink, label in effects:
            if len(set_aside) >= 4:
                break
            set_aside.append(label)
            remaining = largest_drug_count({position_of_label[name_] for name_ in set_aside})
            menu.append({"set_aside": list(set_aside), "largest_group_drugs": remaining})
            print(f"    setting aside {', '.join(set_aside)}: largest group holds {remaining} drugs")
        summary[name]["set_aside_menu"] = menu

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(summary, indent=1, sort_keys=True, default=str))
    print(f"\nwrote {arguments.output}")


if __name__ == "__main__":
    main()
