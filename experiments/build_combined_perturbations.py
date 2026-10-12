#!/usr/bin/env python
"""Write the combination rows of an evidence table: a drug carried with a lowered gene, or with a second drug.

What this builds. Rows 4 to 6 of section 8 of docs/label_source_review.md name three label kinds that no single
perturbation can express, and the statements for all three are already curated with their quotations in
configs/comparative_checks.csv (experiments/build_comparative_checks.py). This script turns those statements into
labelled examples of the combined perturbation type (mechanistic_pathway_learning/evidence/combined_perturbations.py),
which experiments/run_main_model.py reads with --combined-perturbations and which nothing reads by default.

The label rule, and what it does and does not claim. A statement says the drug's effects are at least as large when the
second factor is present. Read as a label that is a weaker claim than the statement makes: the combination is a
positive on the symptoms the statement reads, which are the one it names or, where it names none, the drug's own kept
positive symptoms, and it carries that drug's own weight and grade for each. So the row asserts that the symptom still
occurs under the combination, not that it occurs more. The ordering, which is what the statement actually says, is
read by experiments/run_comparative_checks.py and is not a label here.

Two consequences worth stating before any run uses these rows:

  - a combination of this kind carries labels the drug already carries, so a drug with a statement would take two
    drugs' share of the loss. experiments/run_main_model.py --equal-drug-shares puts a combination in its drug's slot
    for that reason (drug_loss_slots there), and a run without that option gives the drug the larger share;
  - a negative control of the check table (expectation "none") is not written as a label at all. Its statement says
    nothing is expected, and a labelled negative would say the symptom does not occur, which no source here states.

Leakage. A combination inherits the groups of its parts, which is section 8's condition for adding these rows after a
held-out set is drawn, and the loader makes it exact under group_by "disease_cluster_and_targets"
(merge_drugs_with_their_targets). Under the simpler groupings a lowered gene that is itself a perturbation of the table
(ALB is one) would sit in its own group while the combination sat in the drug's, so the summary names those genes and
the script warns.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evaluation.comparative_checks import (
    CHECK_KINDS,
    ComparativeCheck,
    load_checks,
)
from mechanistic_pathway_learning.evaluation.experiment_data import read_protein_of_gene
from mechanistic_pathway_learning.evidence.combined_perturbations import (
    COMBINATION_PERTURBATION_TYPE,
    LOWERED_GENE_SIGN,
    PARTS_COLUMN,
    PerturbationPart,
    combined_perturbation_id,
    nodes_json_of_parts,
    parts_json,
    parts_with_disjoint_nodes,
)

EVIDENCE_CLASS = "combined_perturbation"
OUTPUT_FILE = "combined_perturbation_records.parquet"
SUMMARY_FILE = "combined_perturbation_summary.json"
KINDS_WITH_A_GENE_PART: tuple[str, ...] = ("carrier_level", "enzyme_activity", "displacement")


def gene_part(gene_symbols: tuple[str, ...], protein_of_gene: dict[str, list[str]], node_ids: set[str]) -> list[PerturbationPart]:
    """One part per lowered gene, seeded on the nodes that stand for its product with sign LOWERED_GENE_SIGN.

    The nodes and the shared magnitude are those evaluation/comparative_checks.py seeds for the same gene, so the label
    and the check put the same input to a model: the gene's protein nodes on a split graph, its gene node otherwise,
    with a total magnitude of 1 shared among them.
    """
    parts = []
    for gene_symbol in gene_symbols:
        nodes = [node_id for node_id in (protein_of_gene.get(f"GENE:{gene_symbol}") or [f"GENE:{gene_symbol}"]) if node_id in node_ids]
        if not nodes:
            continue
        parts.append(PerturbationPart(gene_symbol, "gene", tuple((node_id, LOWERED_GENE_SIGN, 1.0 / len(nodes)) for node_id in nodes)))
    return parts


def parts_of_check(check: ComparativeCheck, drug_nodes: dict[str, str], protein_of_gene: dict[str, list[str]], node_ids: set[str]) -> tuple[list[PerturbationPart], str]:
    """(the parts of the combination a statement describes, the reason it has none).

    The order is the one the statement reads in: the drug, then the second drug where there is one, then the lowered
    genes. A drug part's nodes are the triples the evidence table gives that drug, so its input is that drug's input.
    """
    parts: list[PerturbationPart] = []
    for perturbation_id in check.drug_perturbation_ids:
        if perturbation_id not in drug_nodes:
            return [], f"{perturbation_id} is not a drug of the evidence table"
        triples = [(node_id, float(sign), float(magnitude)) for node_id, sign, magnitude in json.loads(drug_nodes[perturbation_id]) if node_id in node_ids]
        if not triples:
            return [], f"{perturbation_id} seeds no node of this graph"
        parts.append(PerturbationPart(perturbation_id, "drug", tuple(triples)))
    if check.kind in KINDS_WITH_A_GENE_PART:
        lowered = gene_part(check.lowered_genes, protein_of_gene, node_ids)
        if len(lowered) != len(check.lowered_genes):
            return [], f"no node for {', '.join(check.lowered_genes)}"
        parts += lowered
    if len(parts) < 2:
        return [], f"a {check.kind} statement gives one part only"
    # Two parts can name one node (two drugs on one transporter): the node then carries the sum of their signed
    # magnitudes and belongs to the first part that names it, so each part's nodes are its own.
    disjoint, _ = parts_with_disjoint_nodes(parts)
    if not any(part.nodes for part in disjoint):
        return [], "every node of the combination cancels to magnitude zero"
    return disjoint, ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--checks", type=Path, default=Path("configs/comparative_checks.csv"))
    parser.add_argument("--evidence-dir", type=Path, required=True, help="the evidence table the drugs and their labels are read from")
    parser.add_argument("--graph-dir", type=Path, required=True, help="the graph the nodes are resolved in")
    parser.add_argument("--kinds", nargs="+", choices=CHECK_KINDS, default=list(CHECK_KINDS))
    parser.add_argument("--label-selection", type=Path, default=None,
                        help="the selection a run uses (experiments/build_label_selection.py); a positive it sets aside is not a label of a combination either")
    parser.add_argument("--relation", default="induces")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, default=None, help="markdown written beside the parquet (default: none)")
    arguments = parser.parse_args()

    evidence = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet")
    evidence = evidence[evidence.relation == arguments.relation]
    if "positive_report_count" in evidence.columns:
        evidence = evidence[evidence.positive_report_count > 0]
    node_ids = set(pd.read_parquet(arguments.graph_dir / "nodes.parquet").node_id)
    protein_of_gene = read_protein_of_gene(arguments.graph_dir)
    # A pair the selection does not keep as a positive is not a positive of a combination either. That covers the pairs
    # it sets aside (keep False) and the never-above-placebo pairs it turns into negatives (counts_as_negative): a
    # placebo comparison that says the drug alone is not above placebo on a symptom cannot leave the combination a
    # positive on it. Rows with masks_a_negative name pairs that are not positives at all, so they never arise here.
    not_a_kept_positive = set()
    if arguments.label_selection is not None:
        selection = pd.read_parquet(arguments.label_selection)
        kept = selection.keep.astype(bool)
        if "counts_as_negative" in selection.columns:
            kept = kept & ~selection.counts_as_negative.fillna(False).astype(bool)
        not_a_kept_positive = set(zip(selection.perturbation_id[~kept], selection.symptom[~kept]))

    drug_rows = evidence[evidence.perturbation_type == "drug"]
    drug_nodes = dict(zip(drug_rows.perturbation_id, drug_rows.perturbation_nodes))
    drug_labels = {perturbation_id: group for perturbation_id, group in drug_rows.groupby("perturbation_id")}
    perturbations_in_the_table = set(evidence.perturbation_id)
    quotes = pd.read_csv(arguments.checks, dtype=str).fillna("").set_index("statement_id")

    rows, skipped = [], []
    for check in load_checks(arguments.checks):
        if check.kind not in arguments.kinds:
            continue
        if check.expectation != "at_least_as_high":
            skipped.append({"statement_id": check.statement_id, "kind": check.kind, "reason": f"expectation {check.expectation!r}: a negative control is not written as a label"})
            continue
        parts, reason = parts_of_check(check, drug_nodes, protein_of_gene, node_ids)
        if reason:
            skipped.append({"statement_id": check.statement_id, "kind": check.kind, "reason": reason})
            continue
        nodes_per_part = [len(part.nodes) for part in parts]
        labels = drug_labels[check.drug_perturbation_id]
        if check.symptom:
            labels = labels[labels.symptom == check.symptom]
        labels = labels[[(check.drug_perturbation_id, symptom) not in not_a_kept_positive for symptom in labels.symptom]]
        if not len(labels):
            skipped.append({"statement_id": check.statement_id, "kind": check.kind,
                            "reason": f"{check.drug_perturbation_id} has no kept positive symptom" + (f" at {check.symptom}" if check.symptom else "")})
            continue
        statement = quotes.loc[check.statement_id]
        perturbation_id = combined_perturbation_id([part.part_id for part in parts])
        label = " + ".join([statement["drug"]] + ([statement["second_drug"]] if statement["second_drug"] else []) +
                           ([f"{gene} lowered" for gene in check.lowered_genes]))
        for row in labels.itertuples(index=False):
            rows.append({"perturbation_id": perturbation_id, "perturbation_type": COMBINATION_PERTURBATION_TYPE, "perturbation_label": label,
                         "group_id": row.group_id, "disease_cluster_id": row.disease_cluster_id, "symptom": row.symptom, "relation": arguments.relation,
                         "evidence_class": EVIDENCE_CLASS, "grade": row.grade, "weight": float(row.weight), "label_frequency": row.label_frequency,
                         "in_metabolic_layer": bool(row.in_metabolic_layer), "perturbation_nodes": nodes_json_of_parts(parts), PARTS_COLUMN: parts_json(parts),
                         "statement_id": check.statement_id, "statement_kind": check.kind, "source": statement["source"], "source_id": statement["source_id"],
                         "source_version": statement["source_version"], "quote": statement["quote"],
                         "nodes_per_part": ";".join(str(count) for count in nodes_per_part),
                         "positive_report_count": 1, "report_count": 1, "negative_report_count": 0})

    table = pd.DataFrame(rows)
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_parquet(arguments.output_dir / OUTPUT_FILE, index=False)
    lowered_also_a_perturbation = sorted({gene for check in load_checks(arguments.checks) for gene in check.lowered_genes} & perturbations_in_the_table)
    summary = {
        "checks": str(arguments.checks), "evidence_dir": str(arguments.evidence_dir), "graph_dir": str(arguments.graph_dir),
        "label_selection": None if arguments.label_selection is None else str(arguments.label_selection),
        "kinds": list(arguments.kinds), "combinations": int(table.perturbation_id.nunique()) if len(table) else 0, "positive_pairs": len(table),
        "by_kind": {} if not len(table) else {kind: {"combinations": int(group.perturbation_id.nunique()), "positive_pairs": len(group)}
                                             for kind, group in table.groupby("statement_kind")},
        "statements_written": 0 if not len(table) else int(table.statement_id.nunique()),
        "statements_skipped": len(skipped), "skipped": skipped,
        "lowered_genes_that_are_themselves_perturbations": lowered_also_a_perturbation,
    }
    (arguments.output_dir / SUMMARY_FILE).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "skipped"}, indent=2))
    if lowered_also_a_perturbation:
        print(f"warning: {lowered_also_a_perturbation} are lowered by a statement and are perturbations of the table; load these rows with "
              f"group_by disease_cluster_and_targets, which holds a combination out with every part of it")
    if arguments.report is not None:
        write_report(arguments.report, table, summary)


def write_report(path: Path, table: pd.DataFrame, summary: dict) -> None:
    provenance = (f"Written by `experiments/build_combined_perturbations.py` from `{summary['checks']}` on the labels of "
                  f"`{summary['evidence_dir']}` and the nodes of `{summary['graph_dir']}`.")
    counted = (f"{summary['combinations']} combinations over {summary['statements_written']} statements, "
               f"{summary['positive_pairs']} positive pairs. {summary['statements_skipped']} statements wrote nothing.")
    lines = ["# Combination rows of the evidence table", "", provenance, "", counted, "",
             "| kind | combinations | positive pairs |", "| --- | --- | --- |"]
    for kind, counts in sorted(summary["by_kind"].items()):
        lines.append(f"| {kind} | {counts['combinations']} | {counts['positive_pairs']} |")
    lines += ["", "## The combinations", "", "| perturbation | kind | symptoms | statement |", "| --- | --- | --- | --- |"]
    if len(table):
        for (perturbation_id, kind, statement_id), group in table.groupby(["perturbation_label", "statement_kind", "statement_id"]):
            lines.append(f"| {perturbation_id} | {kind} | {len(group)} | `{statement_id}` |")
    lines += ["", "## Statements that wrote nothing", "", "| statement | kind | reason |", "| --- | --- | --- |"]
    for row in summary["skipped"]:
        lines.append(f"| `{row['statement_id']}` | {row['kind']} | {row['reason']} |")
    lines.append("")
    path.write_text("\n".join(lines))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
