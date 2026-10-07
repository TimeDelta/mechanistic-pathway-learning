"""How much mechanism a held-out gene shares with the training genes of the disease-cluster split.

The disease-cluster grouping keeps the genes of one disease in one fold, which blocks leakage inside a disease. It does
not block leakage across diseases: two genes whose products form one protein complex, or catalyse the same reaction,
can cause different diseases and land on opposite sides of the split. Network propagation benchmarks have found that
standard cross-validation is over-optimistic for exactly this reason (Picart-Armada et al. 2018, protein complexes).
The opposite pressure is a held-out gene whose pathway has no training gene at all, which the subsystem hold-out
measures as a split of its own.

For every gene perturbation this script records, against the genes in the other folds:

- a complex partner: a gene joined to it by "and" in one clause of a reaction's gene-reaction rule (a subunit of the
  same enzyme complex), after the rule is expanded to disjunctive normal form
- a shared catalysed reaction: any gene catalysing a reaction it catalyses (complex partners and isoenzymes)
- a shared pathway subsystem: any Human-GEM subsystem in common, catch-all subsystems excluded
- the same primary subsystem (the subsystem most of its reactions belong to, as the subsystem hold-out assigns it)
- a primary subsystem with no training gene at all

Sharing is leakage only if it carries label information, so the script also compares the symptom labels of each
held-out gene with those of its training partners under each kind of sharing, against all held-out x training pairs:
the Jaccard similarity of positive symptom sets, and the rate at which a partner shares at least one positive symptom.

Usage: python experiments/measure_fold_mechanism_sharing.py [--graph-dir ...] [--evidence-dir ...] [--output ...]
"""
from __future__ import annotations

import argparse
import itertools
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import (
    CATCH_ALL_SUBSYSTEMS,
    assign_grouped_folds,
    primary_subsystem_by_gene_node,
)

GENE_RULE_TOKEN = re.compile(r"\(|\)|\band\b|\bor\b|[A-Za-z0-9_.:-]+")
MAXIMUM_CLAUSES_PER_RULE = 4096


def gene_rule_clauses(rule: str) -> tuple[list[frozenset[str]], bool]:
    """Expand a gene-reaction rule to disjunctive normal form: a list of gene sets, each a complex that suffices.

    Returns the clauses and whether the expansion was cut at MAXIMUM_CLAUSES_PER_RULE (then the cut clauses are
    still valid complexes, only not all of them)."""
    tokens = GENE_RULE_TOKEN.findall(rule)
    position = 0
    was_truncated = False

    def parse_disjunction() -> list[frozenset[str]]:
        nonlocal position
        clauses = parse_conjunction()
        while position < len(tokens) and tokens[position] == "or":
            position += 1
            clauses = clauses + parse_conjunction()
        return clauses

    def parse_conjunction() -> list[frozenset[str]]:
        nonlocal position, was_truncated
        clauses = parse_factor()
        while position < len(tokens) and tokens[position] == "and":
            position += 1
            right_clauses = parse_factor()
            combined = []
            for left, right in itertools.product(clauses, right_clauses):
                if len(combined) >= MAXIMUM_CLAUSES_PER_RULE:
                    was_truncated = True
                    break
                combined.append(left | right)
            clauses = combined
        return clauses

    def parse_factor() -> list[frozenset[str]]:
        nonlocal position
        token = tokens[position]
        position += 1
        if token == "(":
            clauses = parse_disjunction()
            position += 1  # the closing parenthesis
            return clauses
        return [frozenset([token])]

    if not tokens:
        return [], False
    return parse_disjunction(), was_truncated


def jaccard_of_positive_sets(first_positive: np.ndarray, second_positive: np.ndarray) -> float:
    union = np.logical_or(first_positive, second_positive).sum()
    return float(np.logical_and(first_positive, second_positive).sum() / union) if union else 0.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=Path("docs/fold_mechanism_sharing.md"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster")
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.seed)
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")

    catalysed_by_index = data.relation_types.index("catalyzed_by")
    reactions_by_gene_node: dict[int, set[int]] = defaultdict(set)
    genes_by_reaction_node: dict[int, set[int]] = defaultdict(set)
    for source, target, relation in zip(data.edge_source, data.edge_target, data.edge_relation):
        if relation != catalysed_by_index:
            continue
        gene_node, reaction_node = (source, target) if data.node_types[source] == "gene" else (target, source)
        reactions_by_gene_node[int(gene_node)].add(int(reaction_node))
        genes_by_reaction_node[int(reaction_node)].add(int(gene_node))

    gene_node_by_ensembl = {ensembl: data.node_index[node_id] for node_id, ensembl in zip(nodes.node_id, nodes.ensembl_gene_id)
                            if isinstance(ensembl, str) and ensembl and node_id in data.node_index}
    complex_partners_by_gene_node: dict[int, set[int]] = defaultdict(set)
    rules_truncated = 0
    for rule in nodes.loc[nodes.node_type == "reaction", "gene_reaction_rule"].dropna():
        if " and " not in rule:
            continue
        clauses, was_truncated = gene_rule_clauses(rule)
        rules_truncated += was_truncated
        for clause in clauses:
            members = [gene_node_by_ensembl[ensembl] for ensembl in clause if ensembl in gene_node_by_ensembl]
            for first, second in itertools.permutations(members, 2):
                complex_partners_by_gene_node[first].add(second)

    pathway_subsystems_by_gene_node = {
        gene_node: {data.node_subsystem[reaction] for reaction in reactions if data.node_subsystem[reaction]
                    and data.node_subsystem[reaction] not in CATCH_ALL_SUBSYSTEMS}
        for gene_node, reactions in reactions_by_gene_node.items()}
    primary_subsystem = primary_subsystem_by_gene_node(data.node_subsystem, data.edge_source, data.edge_target,
                                                       data.edge_relation, catalysed_by_index)

    gene_rows = [index for index, kind in enumerate(data.perturbation_types) if kind == "gene" and len(data.perturbation_seeds[index]) == 1]
    gene_node_of_row = {index: int(data.perturbation_seeds[index][0]) for index in gene_rows}
    positive = data.outcomes >= 0.5
    fold_of_row = {index: fold_by_perturbation[data.perturbation_ids[index]] for index in gene_rows}

    sharing_kinds = ("complex partner", "shared catalysed reaction", "shared pathway subsystem", "same primary subsystem")

    def shares(kind: str, first_node: int, second_node: int) -> bool:
        if kind == "complex partner":
            return second_node in complex_partners_by_gene_node.get(first_node, ())
        if kind == "shared catalysed reaction":
            return bool(reactions_by_gene_node.get(first_node, set()) & reactions_by_gene_node.get(second_node, set()))
        if kind == "shared pathway subsystem":
            return bool(pathway_subsystems_by_gene_node.get(first_node, set()) & pathway_subsystems_by_gene_node.get(second_node, set()))
        first_primary, second_primary = primary_subsystem.get(first_node), primary_subsystem.get(second_node)
        # a catch-all primary subsystem (a transporter with no pathway) is not a pathway, as in the subsystem hold-out
        return first_primary is not None and first_primary not in CATCH_ALL_SUBSYSTEMS and first_primary == second_primary

    held_out_with_training_partner = {kind: defaultdict(int) for kind in sharing_kinds}
    held_out_with_same_fold_partner_only = {kind: defaultdict(int) for kind in sharing_kinds}
    pathway_absent_from_training = defaultdict(int)
    held_out_with_a_primary_subsystem = defaultdict(int)
    held_out_count = defaultdict(int)
    pair_jaccard = {kind: [] for kind in sharing_kinds + ("any held-out x training pair",)}
    # per held-out gene: mean Jaccard to its training partners of a kind minus its mean Jaccard to every training gene,
    # so each gene is compared with its own baseline and pairs from one gene are not treated as independent
    paired_jaccard_excess = {kind: [] for kind in sharing_kinds}
    pair_shares_a_positive = {kind: [] for kind in sharing_kinds + ("any held-out x training pair",)}

    for held_out_row in gene_rows:
        fold = fold_of_row[held_out_row]
        held_out_node = gene_node_of_row[held_out_row]
        held_out_count[fold] += 1
        training_rows = [row for row in gene_rows if fold_of_row[row] != fold]
        same_fold_rows = [row for row in gene_rows if fold_of_row[row] == fold and row != held_out_row]
        for kind in sharing_kinds:
            training_partners = [row for row in training_rows if shares(kind, held_out_node, gene_node_of_row[row])]
            if training_partners:
                held_out_with_training_partner[kind][fold] += 1
            elif any(shares(kind, held_out_node, gene_node_of_row[row]) for row in same_fold_rows):
                held_out_with_same_fold_partner_only[kind][fold] += 1
            if positive[held_out_row].any():
                # the same population as the baseline row: both genes of the pair have a positive symptom
                for row in (row for row in training_partners if positive[row].any()):
                    pair_jaccard[kind].append(jaccard_of_positive_sets(positive[held_out_row], positive[row]))
                    pair_shares_a_positive[kind].append(bool(np.logical_and(positive[held_out_row], positive[row]).any()))
        if held_out_node in primary_subsystem and primary_subsystem[held_out_node] not in CATCH_ALL_SUBSYSTEMS:
            held_out_with_a_primary_subsystem[fold] += 1
            if not any(primary_subsystem.get(gene_node_of_row[row]) == primary_subsystem[held_out_node] for row in training_rows):
                pathway_absent_from_training[fold] += 1
        if positive[held_out_row].any():
            jaccard_to_every_training_gene = np.mean([jaccard_of_positive_sets(positive[held_out_row], positive[row])
                                                      for row in training_rows if positive[row].any()])
            for kind in sharing_kinds:
                partner_jaccards = [jaccard_of_positive_sets(positive[held_out_row], positive[row]) for row in training_rows
                                    if positive[row].any() and shares(kind, held_out_node, gene_node_of_row[row])]
                if partner_jaccards:
                    paired_jaccard_excess[kind].append(float(np.mean(partner_jaccards)) - float(jaccard_to_every_training_gene))
            for row in training_rows:
                if positive[row].any():
                    pair_jaccard["any held-out x training pair"].append(jaccard_of_positive_sets(positive[held_out_row], positive[row]))
                    pair_shares_a_positive["any held-out x training pair"].append(bool(np.logical_and(positive[held_out_row], positive[row]).any()))

    folds = sorted(held_out_count)
    total_held_out = sum(held_out_count.values())
    lines = [
        "# Mechanism shared across the disease-cluster folds (generated by experiments/measure_fold_mechanism_sharing.py)",
        "",
        f"{total_held_out} gene perturbations, {arguments.num_folds} disease-cluster folds (seed {arguments.seed}), graph "
        f"{arguments.graph_dir}, evidence {arguments.evidence_dir}. Each gene is counted in the fold that holds it out, "
        "against the genes of the other folds. Complex partners come from the gene-reaction rules expanded to disjunctive "
        f"normal form ({rules_truncated} rules cut at {MAXIMUM_CLAUSES_PER_RULE} clauses); catch-all subsystems "
        f"({', '.join(sorted(CATCH_ALL_SUBSYSTEMS))}) are excluded from the subsystem comparisons.",
        "",
        "## Held-out genes sharing mechanism with a training gene",
        "",
        "| fold | held-out genes | " + " | ".join(sharing_kinds) + " | primary subsystem absent from training |",
        "|---|---|" + "---|" * len(sharing_kinds) + "---|",
    ]
    for fold in folds:
        cells = [f"{held_out_with_training_partner[kind][fold]} ({held_out_with_training_partner[kind][fold] / held_out_count[fold]:.0%})"
                 for kind in sharing_kinds]
        absent = f"{pathway_absent_from_training[fold]} of {held_out_with_a_primary_subsystem[fold]}"
        lines.append(f"| {fold} | {held_out_count[fold]} | " + " | ".join(cells) + f" | {absent} |")
    total_cells = [f"{sum(held_out_with_training_partner[kind].values())} ({sum(held_out_with_training_partner[kind].values()) / total_held_out:.0%})"
                   for kind in sharing_kinds]
    total_absent = f"{sum(pathway_absent_from_training.values())} of {sum(held_out_with_a_primary_subsystem.values())}"
    lines.append(f"| all | {total_held_out} | " + " | ".join(total_cells) + f" | {total_absent} |")
    lines += [
        "",
        "Shared only with a gene of the same fold (so the grouping or the fold assignment kept the pair together):",
        "",
        "| " + " | ".join(sharing_kinds) + " |",
        "|" + "---|" * len(sharing_kinds),
        "| " + " | ".join(str(sum(held_out_with_same_fold_partner_only[kind].values())) for kind in sharing_kinds) + " |",
        "",
        "## Do partners across the split share labels?",
        "",
        "Pairs of a held-out gene and a training gene in which both have at least one positive symptom (outcome at least 0.5), "
        "by the kind of mechanism they share; the last row is every such pair whatever they share. If sharing leaks labels, "
        "partner pairs should be more alike than that baseline. The same primary subsystem excludes catch-all subsystems.",
        "",
        "| pair kind | pairs | mean Jaccard of positive symptom sets | share at least one positive symptom |",
        "|---|---|---|---|",
    ]
    for kind, values in pair_jaccard.items():
        if values:
            lines.append(f"| {kind} | {len(values)} | {np.mean(values):.3f} | {np.mean(pair_shares_a_positive[kind]):.3f} |")
        else:
            lines.append(f"| {kind} | 0 | - | - |")
    bootstrap_generator = np.random.default_rng(arguments.seed)
    lines += [
        "",
        "The pair table counts every pair, so a held-out gene with many partners weighs more and the pairs are not independent. "
        "The paired form below takes, for each held-out gene with at least one partner of the kind, its mean Jaccard to those "
        "partners minus its mean Jaccard to every training gene, and bootstraps the mean of that excess over held-out genes "
        "(2,000 resamples, 95 percent percentile interval).",
        "",
        "| pair kind | held-out genes | mean excess Jaccard over the gene's own baseline [95% CI] |",
        "|---|---|---|",
    ]
    for kind in sharing_kinds:
        excess = np.asarray(paired_jaccard_excess[kind])
        if len(excess) < 2:
            lines.append(f"| {kind} | {len(excess)} | - |")
            continue
        resampled_means = [excess[bootstrap_generator.integers(0, len(excess), len(excess))].mean() for _ in range(2000)]
        low, high = np.percentile(resampled_means, [2.5, 97.5])
        lines.append(f"| {kind} | {len(excess)} | {excess.mean():+.3f} [{low:+.3f}, {high:+.3f}] |")
    arguments.output.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
