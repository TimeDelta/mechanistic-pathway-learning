"""Leakage group rules for the drugs: how far each rule merges the drugs, and what leakage each one leaves behind.

The registered rule (group_by "disease_cluster_and_targets", merge_drugs_with_their_targets) is union-find over shared
target nodes: a drug joins every gene perturbation whose node it targets and every drug that shares a target node with
it. That closes the drug-gene leakage of docs/drug_target_leakage.md, but it is all-or-nothing, and once measured
off-targets join as well it collapses 94 percent of the drugs into one group (docs/off_target_scoping.md).

Target families are the alternative: a drug's group is the GtoPdb family of its mechanism target
(data/raw/gtopdb_2026/targets_and_families.csv). 40 of the 222 drugs have mechanism targets in more than one family, so
the rule needs a tie-break. The user's decision of 9 October 2026 ("I don't like the alphabetical thing. Try the other
version of family tie breaking and also make sure to include a measurement of the remaining leakage confound from this
decision") is the highest-affinity tie-break: the family of the mechanism target the drug binds best.

Family grouping bounds leakage rather than removing it, so every rule is measured on the same five grouped folds
(assign_grouped_folds, the study's own fold builder) and reported with what crosses a fold boundary under it:
  - test drugs sharing a mechanism target with a training drug, and the kept positives they hold;
  - test drugs sharing any measured target (mechanism or off-target above the occupancy threshold) with a training drug;
  - test drugs whose mechanism target is itself a training gene perturbation, and the kept positives they hold;
  - test gene perturbations whose gene is a mechanism target of a training drug;
  - test drugs sharing any mechanism family with a training drug, which the family rules do not bound because a drug is
    assigned one family and may bind targets in others.

    PYTHONPATH=. python experiments/measure_leakage_group_rules.py
"""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from experiments.scope_off_target_binding import (GTOPDB_INTERACTIONS, GTOPDB_LIGAND_MAPPING, DRUGCENTRAL_INTERACTIONS,
                                                  drug_keys, drugcentral_measurements, gtopdb_measurements, occupancy)
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, merge_drugs_with_their_targets
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds

GTOPDB_FAMILIES = Path("data/raw/gtopdb_2026/targets_and_families.csv")
MARKUP = re.compile(r"<[^>]+>")  # GtoPdb family names carry HTML (5-HT<sub>3</sub> receptors)
OCCUPANCY_THRESHOLD = 0.25
NUM_FOLDS = 5
UNASSIGNED_FAMILY = "unassigned"
# The measurement of experiments/measure_family_in_graph.py, which bounds what any grouping of drugs can remove.
FAMILY_IS_IN_THE_GRAPH_CAVEAT = (
    "the family is recoverable from what the model reads, so no grouping of drugs removes the family confound: docs/family_in_graph.md recovers a held-out mechanism target's family 0.90 of the time from the 64 protein components, 0.91 from its graph neighbours and 0.97 from its Reactome protein-entity memberships, against a chance rate of 0.09. Grouping changes which drugs share a fold, not what the graph tells the model about a held-out target.")


def target_families(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    """HGNC symbol -> GtoPdb family name, and HGNC symbol -> target type."""
    table = pd.read_csv(path, skiprows=1, low_memory=False)
    table = table[table["HGNC symbol"].notna()]
    family_names = [MARKUP.sub("", str(name)) for name in table["Family name"]]
    return dict(zip(table["HGNC symbol"], family_names)), dict(zip(table["HGNC symbol"], table["Type"]))


def measurements_per_drug(data, symbol_of_node: dict[str, str]) -> dict[int, dict]:
    """Position of each drug perturbation -> its mechanism genes, every gene with a measured affinity, and the
    off-target genes it occupies above OCCUPANCY_THRESHOLD. Affinity per pair is the maximum over the measurements
    (the user's decision of 9 October 2026, as in experiments/scope_off_target_binding.py)."""
    drugcentral = drugcentral_measurements(DRUGCENTRAL_INTERACTIONS)
    gtopdb = gtopdb_measurements(GTOPDB_INTERACTIONS, GTOPDB_LIGAND_MAPPING)
    drugcentral_ids_of_name = drugcentral.groupby("drug_name").drugcentral_id.first().to_dict()
    per_drug = {}
    for position, (perturbation_id, kind, label, seeds) in enumerate(zip(data.perturbation_ids, data.perturbation_types,
                                                                        data.perturbation_labels, data.perturbation_seeds)):
        if kind != "drug":
            continue
        mechanism_genes = sorted({symbol_of_node.get(data.node_ids[seed]) for seed in seeds} - {None})
        pubchem_cid, name = drug_keys(perturbation_id, label)
        from_gtopdb = gtopdb[((gtopdb.pubchem_cid == pubchem_cid) if pubchem_cid else False) | (gtopdb.drug_name == name) | (gtopdb.inn == name)]
        drugcentral_ids = set(from_gtopdb.drugcentral_id.dropna().astype(int)) | ({drugcentral_ids_of_name[name]} if name in drugcentral_ids_of_name else set())
        measurements = pd.concat([drugcentral[drugcentral.drugcentral_id.isin(drugcentral_ids)], from_gtopdb], ignore_index=True)
        entry = {"mechanism_genes": mechanism_genes, "p_affinity_of_gene": {}, "num_measurements_of_gene": {}, "kept_off_target_genes": []}
        if not measurements.empty:
            by_gene = measurements.groupby("gene").agg(p_affinity=("p_affinity", "max"), num_measurements=("p_affinity", "count"))
            with_affinity = by_gene[by_gene.p_affinity.notna()]
            entry["p_affinity_of_gene"] = {gene: float(value) for gene, value in with_affinity.p_affinity.items()}
            entry["num_measurements_of_gene"] = {gene: int(value) for gene, value in with_affinity.num_measurements.items()}
            on_mechanism = [value for gene, value in entry["p_affinity_of_gene"].items() if gene in set(mechanism_genes)]
            primary = max(on_mechanism) if on_mechanism else None
            if primary is not None:
                entry["kept_off_target_genes"] = sorted(gene for gene, value in entry["p_affinity_of_gene"].items()
                                                        if gene not in set(mechanism_genes) and occupancy(primary, value) >= OCCUPANCY_THRESHOLD)
        per_drug[position] = entry
    return per_drug


def family_of_drug(entry: dict, family_of_gene: dict[str, str], tie_break: str) -> tuple[str, str]:
    """(family name, how it was chosen). "alphabetical" takes the first family name of the drug's mechanism targets.
    "highest_affinity" takes the family of the mechanism target with the largest measured affinity, ties broken by the
    number of measurements and then by gene symbol; a drug with no measured affinity on a mechanism target that carries
    a family falls back to the alphabetical rule."""
    with_a_family = [gene for gene in entry["mechanism_genes"] if gene in family_of_gene]
    if not with_a_family:
        return UNASSIGNED_FAMILY, "no mechanism target carries a family"
    alphabetical = sorted({family_of_gene[gene] for gene in with_a_family})[0]
    if tie_break == "alphabetical":
        return alphabetical, "first family name"
    measured = [gene for gene in with_a_family if gene in entry["p_affinity_of_gene"]]
    if not measured:
        return alphabetical, "no measured affinity on a mechanism target; first family name"
    best = max(measured, key=lambda gene: (entry["p_affinity_of_gene"][gene], entry["num_measurements_of_gene"][gene], gene))
    return family_of_gene[best], f"best-bound mechanism target {best} (p{entry['p_affinity_of_gene'][best]:.2f})"


def groups_under_rule(rule: str, data, per_drug: dict[int, dict], family_of_gene: dict[str, str],
                      gene_perturbation_of_symbol: dict[str, str], gene_node_of_symbol: dict[str, str]) -> list[str]:
    """One group label per perturbation under the named rule."""
    if rule == "disease_cluster_and_targets":
        return merge_drugs_with_their_targets(data.perturbation_ids, data.perturbation_types, data.group_ids, data.perturbation_seeds, data.node_ids)
    if rule == "disease_cluster_and_targets_with_off_targets":
        seeds_with_off_targets = list(data.perturbation_seeds)
        node_index = {node_id: index for index, node_id in enumerate(data.node_ids)}
        for position, entry in per_drug.items():
            extra = [node_index[gene_node_of_symbol[gene]] for gene in entry["kept_off_target_genes"]
                     if gene in gene_node_of_symbol and gene_node_of_symbol[gene] in node_index]
            seeds_with_off_targets[position] = list(data.perturbation_seeds[position]) + extra
        return merge_drugs_with_their_targets(data.perturbation_ids, data.perturbation_types, data.group_ids, seeds_with_off_targets, data.node_ids)
    tie_break = "alphabetical" if rule.endswith("alphabetical") else "highest_affinity"
    groups = ["group:" + group_id for group_id in data.group_ids]
    for position, entry in per_drug.items():
        groups[position] = "family:" + family_of_drug(entry, family_of_gene, tie_break)[0]
    if rule.endswith("_with_target_genes"):  # keep the drug-gene join of the registered rule, union-find over it only
        parent: dict[str, str] = {}

        def find(item: str) -> str:
            parent.setdefault(item, item)
            while parent[item] != item:
                parent[item] = parent[parent[item]]
                item = parent[item]
            return item

        def union(first: str, second: str) -> None:
            parent[find(first)] = find(second)

        for perturbation_id, group in zip(data.perturbation_ids, groups):
            union("perturbation:" + perturbation_id, group)
        for position, entry in per_drug.items():
            for gene in entry["mechanism_genes"]:
                if gene in gene_perturbation_of_symbol:
                    union("perturbation:" + data.perturbation_ids[position], "perturbation:" + gene_perturbation_of_symbol[gene])
        names: dict[str, str] = {}
        for perturbation_id, group in zip(data.perturbation_ids, groups):
            root = find("perturbation:" + perturbation_id)
            names[root] = min(names.get(root, group), group)
        return [names[find("perturbation:" + perturbation_id)] for perturbation_id in data.perturbation_ids]
    return groups


def group_sizes(groups: list[str], data, kept_positive: np.ndarray) -> dict:
    """Group counts over the drugs and the largest group's share of drugs and of kept drug positives."""
    drug_positions = [position for position, kind in enumerate(data.perturbation_types) if kind == "drug"]
    drugs_per_group, positives_per_group = Counter(), Counter()
    for position in drug_positions:
        drugs_per_group[groups[position]] += 1
        positives_per_group[groups[position]] += int(kept_positive[position])
    total_positives = max(sum(positives_per_group.values()), 1)
    return {"num_groups_holding_a_drug": len(drugs_per_group), "num_groups_over_all_perturbations": len(set(groups)),
            "largest_group_drugs": max(drugs_per_group.values()), "largest_group_share_of_drugs": round(max(drugs_per_group.values()) / len(drug_positions), 4),
            "largest_group_share_of_kept_drug_positives": round(max(positives_per_group.values()) / total_positives, 4),
            "groups_with_five_or_more_drugs": sum(1 for size in drugs_per_group.values() if size >= 5),
            "drugs_in_groups_with_five_or_more": sum(size for size in drugs_per_group.values() if size >= 5)}


def residual_leakage(groups: list[str], data, per_drug: dict[int, dict], family_of_gene: dict[str, str],
                     kept_positive: np.ndarray, gene_perturbation_of_symbol: dict[str, str]) -> dict:
    """What crosses a fold boundary under this rule, as the mean over the five grouped folds."""
    fold_of_perturbation = assign_grouped_folds(list(data.perturbation_ids), groups, num_folds=NUM_FOLDS, random_seed=0)
    folds = np.array([fold_of_perturbation[perturbation_id] for perturbation_id in data.perturbation_ids])
    symbol_of_gene_perturbation = {perturbation_id: symbol for symbol, perturbation_id in gene_perturbation_of_symbol.items()}
    measured_genes_of_drug = {position: set(entry["mechanism_genes"]) | set(entry["kept_off_target_genes"]) for position, entry in per_drug.items()}
    families_of_drug = {position: {family_of_gene[gene] for gene in entry["mechanism_genes"] if gene in family_of_gene} for position, entry in per_drug.items()}
    readings = defaultdict(list)
    for fold in range(NUM_FOLDS):
        test_drugs = [position for position in per_drug if folds[position] == fold]
        training_drugs = [position for position in per_drug if folds[position] != fold]
        if not test_drugs:
            continue
        training_mechanism_genes = {gene for position in training_drugs for gene in per_drug[position]["mechanism_genes"]}
        training_measured_genes = {gene for position in training_drugs for gene in measured_genes_of_drug[position]}
        training_families = {family for position in training_drugs for family in families_of_drug[position]}
        training_gene_perturbation_symbols = {symbol_of_gene_perturbation[data.perturbation_ids[position]]
                                              for position, kind in enumerate(data.perturbation_types)
                                              if kind == "gene" and folds[position] != fold and data.perturbation_ids[position] in symbol_of_gene_perturbation}
        shares_mechanism = [position for position in test_drugs if set(per_drug[position]["mechanism_genes"]) & training_mechanism_genes]
        shares_measured = [position for position in test_drugs if measured_genes_of_drug[position] & training_measured_genes]
        shares_family = [position for position in test_drugs if families_of_drug[position] & training_families]
        targets_a_training_gene = [position for position in test_drugs if set(per_drug[position]["mechanism_genes"]) & training_gene_perturbation_symbols]
        test_positives = max(sum(int(kept_positive[position]) for position in test_drugs), 1)
        test_genes = [position for position, kind in enumerate(data.perturbation_types) if kind == "gene" and folds[position] == fold]
        training_drug_mechanism_genes = {gene for position in training_drugs for gene in per_drug[position]["mechanism_genes"]}
        genes_targeted_by_a_training_drug = [position for position in test_genes
                                             if symbol_of_gene_perturbation.get(data.perturbation_ids[position]) in training_drug_mechanism_genes]
        readings["test_drugs"].append(len(test_drugs))
        readings["test_drugs_sharing_a_mechanism_target_with_a_training_drug"].append(len(shares_mechanism) / len(test_drugs))
        readings["test_drugs_sharing_any_measured_target_with_a_training_drug"].append(len(shares_measured) / len(test_drugs))
        readings["test_drugs_sharing_a_mechanism_family_with_a_training_drug"].append(len(shares_family) / len(test_drugs))
        readings["test_drugs_whose_mechanism_target_is_a_training_gene_perturbation"].append(len(targets_a_training_gene) / len(test_drugs))
        readings["share_of_test_drug_positives_on_drugs_sharing_a_mechanism_target"].append(sum(int(kept_positive[position]) for position in shares_mechanism) / test_positives)
        readings["share_of_test_drug_positives_on_drugs_targeting_a_training_gene"].append(sum(int(kept_positive[position]) for position in targets_a_training_gene) / test_positives)
        readings["test_gene_perturbations_targeted_by_a_training_drug"].append(len(genes_targeted_by_a_training_drug) / max(len(test_genes), 1))
    summary = {name: round(float(np.mean(values)), 4) for name, values in readings.items() if name != "test_drugs"}
    summary["test_drugs_per_fold"] = [int(value) for value in readings["test_drugs"]]
    return summary


def markdown_report(result: dict) -> str:
    lines = ["# Leakage group rules for the drugs", "",
             f"Measured on {result['evidence_dir']} with {result['num_drugs']} drug and "
             f"{result['num_gene_perturbations']} gene perturbations, off-targets kept at occupancy "
             f"{result['occupancy_threshold']}, {result['num_folds']} grouped folds (seed 0). "
             "Generated by experiments/measure_leakage_group_rules.py.", "",
             "## How far each rule merges the drugs", "",
             "| rule | groups holding a drug | largest group: drugs | share of drugs | share of kept drug positives | groups with 5+ drugs |",
             "|---|---|---|---|---|---|"]
    for rule, entry in result["rules"].items():
        sizes = entry["group_sizes"]
        lines.append(f"| `{rule}` | {sizes['num_groups_holding_a_drug']} | {sizes['largest_group_drugs']} | "
                     f"{sizes['largest_group_share_of_drugs']:.2f} | {sizes['largest_group_share_of_kept_drug_positives']:.2f} | "
                     f"{sizes['groups_with_five_or_more_drugs']} |")
    lines += ["", "## Residual leakage under each rule", "",
              "Mean over the five folds. A drug or gene counted here sits in the test fold while something it shares a "
              "target with sits in training; the rule bounds the leakage it does not remove.", "",
              "| rule | test drugs sharing a mechanism target | their share of test drug positives | sharing any measured target | "
              "targeting a training gene perturbation | their share of test drug positives | sharing a mechanism family | "
              "test genes targeted by a training drug |",
              "|---|---|---|---|---|---|---|---|"]
    for rule, entry in result["rules"].items():
        leakage = entry["residual_leakage"]
        lines.append(f"| `{rule}` | {leakage['test_drugs_sharing_a_mechanism_target_with_a_training_drug']:.3f} | "
                     f"{leakage['share_of_test_drug_positives_on_drugs_sharing_a_mechanism_target']:.3f} | "
                     f"{leakage['test_drugs_sharing_any_measured_target_with_a_training_drug']:.3f} | "
                     f"{leakage['test_drugs_whose_mechanism_target_is_a_training_gene_perturbation']:.3f} | "
                     f"{leakage['share_of_test_drug_positives_on_drugs_targeting_a_training_gene']:.3f} | "
                     f"{leakage['test_drugs_sharing_a_mechanism_family_with_a_training_drug']:.3f} | "
                     f"{leakage['test_gene_perturbations_targeted_by_a_training_drug']:.3f} |")
    tie = result["tie_break"]
    lines += ["", "## The tie-break", "",
              f"{tie['drugs_with_a_family']} of {result['num_drugs']} drugs have a mechanism target in a GtoPdb family; "
              f"{tie['drugs_spanning_more_than_one_family']} have mechanism targets in more than one family. The "
              f"highest-affinity rule and the alphabetical rule put {tie['drugs_the_two_rules_place_differently']} drugs "
              f"in different families. {tie['drugs_without_a_measured_mechanism_affinity']} drugs have no measured affinity on a "
              "mechanism target that carries a family, so the highest-affinity rule falls back to the first family name for them.", "",
              "Why the choice is made per drug: how it was chosen, for the drugs the two rules place differently:", ""]
    for label, chosen in tie["examples_of_disagreement"]:
        lines.append(f"- {label}: {chosen}")
    registered = result["rules"]["disease_cluster_and_targets"]
    alphabetical = result["rules"]["target_family_alphabetical"]
    best_bound = result["rules"]["target_family_highest_affinity"]
    with_genes = result["rules"]["target_family_highest_affinity_with_target_genes"]
    moved = result["tie_break"]["drugs_the_two_rules_place_differently"]
    tie_break_reading = (
        f"No drug in this evidence table has mechanism targets in more than one family, so the two tie-breaks give the "
        f"same {best_bound['group_sizes']['num_groups_holding_a_drug']} groups and the choice between them does not "
        "arise here."
        if moved == 0 else
        f"The tie-break by best-bound mechanism target takes the largest group to "
        f"{best_bound['group_sizes']['largest_group_share_of_drugs']:.2f} of the drugs from "
        f"{alphabetical['group_sizes']['largest_group_share_of_drugs']:.2f} alphabetically and moves {moved} drugs, each "
        "for a reason a pharmacologist would give; it leaves the leakage where it was "
        f"({best_bound['residual_leakage']['test_drugs_sharing_a_mechanism_target_with_a_training_drug']:.3f} of test "
        f"drugs sharing a mechanism target with a training drug against "
        f"{alphabetical['residual_leakage']['test_drugs_sharing_a_mechanism_target_with_a_training_drug']:.3f}).")
    lines += ["", "## Reading", "", tie_break_reading, "",
              "The confound the family rules leave is the join they drop. Under the registered rule two drugs that share "
              "a target node are one group and a drug is held out with the gene perturbation of its target, so both "
              f"readings are 0.000 by construction. Under the best-bound family rule "
              f"{best_bound['residual_leakage']['test_drugs_sharing_a_mechanism_target_with_a_training_drug']:.3f} of test "
              "drugs share a mechanism target with a training drug (carrying "
              f"{best_bound['residual_leakage']['share_of_test_drug_positives_on_drugs_sharing_a_mechanism_target']:.3f} of "
              "the test fold's kept drug positives) and "
              f"{best_bound['residual_leakage']['test_drugs_whose_mechanism_target_is_a_training_gene_perturbation']:.3f} "
              "target a gene whose own perturbation is in training (carrying "
              f"{best_bound['residual_leakage']['share_of_test_drug_positives_on_drugs_targeting_a_training_gene']:.3f}), "
              "which is the leakage docs/drug_target_leakage.md exists to close: a model learns \"perturbing this node "
              "causes this symptom\" from one perturbation and is credited for it on another.", "",
              "What the family rule buys is fold balance: its largest group holds "
              f"{best_bound['group_sizes']['largest_group_share_of_drugs']:.2f} of the drugs against "
              f"{registered['group_sizes']['largest_group_share_of_drugs']:.2f} under the registered rule, whose largest "
              f"group also holds {registered['group_sizes']['largest_group_share_of_kept_drug_positives']:.2f} of the kept "
              "drug positives and makes the group bootstrap coarse. Putting the registered rule's drug-to-gene join back "
              "inside the family rule closes that leak "
              f"({with_genes['residual_leakage']['test_drugs_whose_mechanism_target_is_a_training_gene_perturbation']:.3f}) "
              f"and takes the largest group to {with_genes['group_sizes']['largest_group_share_of_drugs']:.2f} of the drugs. "
              "So the balance is bought with leakage rather than with a better partition, and the registered grouping "
              "stays until the user decides otherwise.", "",
              "Caveats on the family rule itself:", "",
              f"- {result['tie_break']['drugs_without_a_measured_mechanism_affinity']} drugs have no measured affinity on a "
              "mechanism target that carries a family, so the best-bound rule falls back to the first family name for them; "
              "the fallback is the arbitrary rule the user rejected, applied to the drugs the data cannot place.",
              "- the grouping inherits GtoPdb's family coverage, and a mechanism gene outside it cannot be grouped at all.",
              f"- a drug is assigned one family while it may bind targets in others, so "
              f"{best_bound['residual_leakage']['test_drugs_sharing_a_mechanism_family_with_a_training_drug']:.3f} of test "
              "drugs share a mechanism family with a training drug under the family rule (0 where every drug's targets sit "
              "in one family).",
              "- the numbers are for the evidence table named at the top of this document, not necessarily the one a "
              "confirmatory run reads.",
              FAMILY_IS_IN_THE_GRAPH_CAVEAT, ""]
    lines += ["", "## Unmatched genes", "",
              f"{len(result['mechanism_genes_without_a_family'])} mechanism genes carry no GtoPdb family: "
              f"{', '.join(result['mechanism_genes_without_a_family'])}.", ""]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v3_max3_targets"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v3_max3_targets.parquet"))
    parser.add_argument("--json-output", type=Path, default=Path("data/processed/off_target_scoping/leakage_group_rules.json"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/leakage_group_rules.md"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster", label_selection=arguments.label_selection)
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet", columns=["node_id", "node_type", "gene_symbol"])
    symbol_of_node = dict(zip(nodes.node_id, nodes.gene_symbol))
    gene_node_of_symbol = dict(zip(nodes.gene_symbol[nodes.node_type == "gene"], nodes.node_id[nodes.node_type == "gene"]))
    family_of_gene, _ = target_families(GTOPDB_FAMILIES)
    kept_positive = ((data.outcomes > 0) & (data.label_mask > 0 if data.label_mask is not None else True)).sum(axis=1)
    gene_perturbation_of_symbol = {symbol_of_node.get(data.node_ids[seeds[0]]): perturbation_id
                                   for perturbation_id, kind, seeds in zip(data.perturbation_ids, data.perturbation_types, data.perturbation_seeds)
                                   if kind == "gene" and len(seeds)}
    gene_perturbation_of_symbol.pop(None, None)
    per_drug = measurements_per_drug(data, symbol_of_node)

    rules = ["disease_cluster_and_targets", "disease_cluster_and_targets_with_off_targets",
             "target_family_alphabetical", "target_family_highest_affinity", "target_family_highest_affinity_with_target_genes"]
    result = {"evidence_dir": str(arguments.evidence_dir), "graph_dir": str(arguments.graph_dir),
              "num_drugs": len(per_drug), "num_gene_perturbations": int(sum(kind == "gene" for kind in data.perturbation_types)),
              "occupancy_threshold": OCCUPANCY_THRESHOLD, "num_folds": NUM_FOLDS, "rules": {}}
    for rule in rules:
        groups = groups_under_rule(rule, data, per_drug, family_of_gene, gene_perturbation_of_symbol, gene_node_of_symbol)
        result["rules"][rule] = {"group_sizes": group_sizes(groups, data, kept_positive),
                                 "residual_leakage": residual_leakage(groups, data, per_drug, family_of_gene, kept_positive, gene_perturbation_of_symbol)}
        print(f"{rule}: {json.dumps(result['rules'][rule]['group_sizes'])}")

    alphabetical = {position: family_of_drug(entry, family_of_gene, "alphabetical") for position, entry in per_drug.items()}
    highest = {position: family_of_drug(entry, family_of_gene, "highest_affinity") for position, entry in per_drug.items()}
    disagreement = [position for position in per_drug if alphabetical[position][0] != highest[position][0]]
    result["tie_break"] = {
        "drugs_with_a_family": sum(1 for position, entry in per_drug.items() if any(gene in family_of_gene for gene in entry["mechanism_genes"])),
        "drugs_spanning_more_than_one_family": sum(1 for entry in per_drug.values()
                                                   if len({family_of_gene[gene] for gene in entry["mechanism_genes"] if gene in family_of_gene}) > 1),
        "drugs_the_two_rules_place_differently": len(disagreement),
        "drugs_without_a_measured_mechanism_affinity": sum(1 for position, entry in per_drug.items()
                                                           if not [gene for gene in entry["mechanism_genes"] if gene in family_of_gene and gene in entry["p_affinity_of_gene"]]),
        "examples_of_disagreement": [[data.perturbation_labels[position],
                                      f"{highest[position][0]} by {highest[position][1]}, against {alphabetical[position][0]} alphabetically"]
                                     for position in disagreement[:12]]}
    result["mechanism_genes_without_a_family"] = sorted({gene for entry in per_drug.values() for gene in entry["mechanism_genes"] if gene not in family_of_gene})

    arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
    arguments.json_output.write_text(json.dumps(result, indent=1) + "\n")
    arguments.markdown_output.write_text(markdown_report(result))
    print(f"wrote {arguments.json_output} and {arguments.markdown_output}")


if __name__ == "__main__":
    main()
