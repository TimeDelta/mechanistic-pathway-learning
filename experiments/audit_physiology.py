"""The physiology audit: does the graph and the evidence behind it represent the underlying physiology as truly as the
sources allow? (The owner, 15 October 2026: "The top priority should always be representing the underlying physiology
as truly as possible.")

Seven sections, each a function returning a dict, rendered to docs/physiology_audit.md (generated, stamped) and
data/processed/audit/physiology_audit.json:

  1. graph integrity: duplicate edges, self-loops, dangling endpoints, contradictory signs, unsigned fractions,
     degree concentration, whether stoichiometric coefficients survive into the graph the model reads;
  2. compartments: a non-transport reaction whose substrates and products sit in different compartments;
  3. mass and charge balance of the Human-GEM reactions, from the SBML through cobra, boundary reactions excluded;
  4. the cited validation set (configs/physiology_validation_set.csv): established chains the evidence should hold
     and established negatives it should not contradict, with the drug's target and sign checked against the table;
  5. evidence coverage per symptom, by class, grade and negatives;
  6. IntAct's non-interacting pairs against the protein edges (mechanistic_pathway_learning/graph/edge_negatives.py);
  7. the rules applied to the sources that could have entered and did not (imaging, mouse proxies without validity).

Usage:
  python experiments/audit_physiology.py --graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v4 --evidence-dir data/processed/evidence_full_v4
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.documentation import provenance_stamp

DEFAULT_VALIDATION_SET = Path("configs/physiology_validation_set.csv")
EXPECTED_SIGN = {"antagonist": -1.0, "inhibitor": -1.0, "loss_of_function": -1.0, "negative_allosteric_modulator": -1.0,
                 "agonist": 1.0, "positive_allosteric_modulator": 1.0, "gain_of_function": 1.0, "activator": 1.0,
                 # a substrate-type releaser drives the transporter backwards, so the pathway sees more transmitter, not
                 # less: the sign convention is the effect on the pathway and ChEMBL types these RELEASING AGENT at +1
                 "releasing_agent": 1.0}
LABEL_GRADES = ("A", "B")


# ----------------------------------------------------------------------------------------------------------------------
# 1. graph integrity


def node_symbol(node_id: str) -> str:
    return str(node_id).split(":", 1)[1] if ":" in str(node_id) else str(node_id)


def graph_integrity(nodes: pd.DataFrame, edges: pd.DataFrame) -> dict:
    node_ids = set(nodes.node_id)
    key = ["source_id", "target_id", "relation_type"]
    exact_duplicates = edges[edges.duplicated(subset=key + ["sign"], keep=False)]
    signed = edges.groupby(key).sign.agg(lambda values: set(float(value) for value in values))
    contradictory = signed[signed.map(lambda signs: 1.0 in signs and -1.0 in signs)]
    unsigned = edges.groupby("relation_type").sign.apply(lambda values: float((values == 0).mean())).round(4).to_dict()
    self_loops = int((edges.source_id == edges.target_id).sum())
    dangling = int((~edges.source_id.isin(node_ids) | ~edges.target_id.isin(node_ids)).sum())
    degree = pd.concat([edges.source_id, edges.target_id]).value_counts()
    isolated = int((~nodes.node_id.isin(degree.index)).sum())
    top = []
    currency = dict(zip(nodes.node_id, nodes.get("is_currency", pd.Series(False, index=nodes.index))))
    names = dict(zip(nodes.node_id, nodes.display_name)) if "display_name" in nodes.columns else {}
    for node_id, count in degree.head(15).items():
        top.append({"node": str(node_id), "name": str(names.get(node_id, "")), "degree": int(count), "is_currency": bool(currency.get(node_id, False))})
    stoichiometry = {"column_present": "stoichiometric_coefficient" in edges.columns}
    if stoichiometry["column_present"]:
        metabolic = edges[edges.relation_type.isin(["substrate_of", "product_of"])]
        stoichiometry["metabolic_edges"] = int(len(metabolic))
        stoichiometry["edges_with_coefficient_other_than_one"] = int((metabolic.stoichiometric_coefficient.fillna(1.0) != 1.0).sum())
    return {"nodes": int(len(nodes)), "nodes_by_type": nodes.node_type.value_counts().to_dict(), "edges": int(len(edges)),
            "edges_by_relation": edges.relation_type.value_counts().to_dict(),
            "edges_by_source": edges.evidence_source.astype(str).str.split(" ").str[0].value_counts().to_dict(),
            "exact_duplicate_edge_rows": int(len(exact_duplicates)), "exact_duplicate_edge_pairs": int(exact_duplicates.groupby(key).ngroups) if len(exact_duplicates) else 0,
            "contradictory_sign_pairs": int(len(contradictory)),
            "contradictory_examples": [{"source_id": s, "target_id": t, "relation_type": r} for (s, t, r) in list(contradictory.index[:10])],
            "unsigned_fraction_by_relation": unsigned, "self_loops": self_loops, "dangling_endpoint_edges": dangling, "isolated_nodes": isolated,
            "top_degree_nodes": top, "stoichiometry": stoichiometry}


# ----------------------------------------------------------------------------------------------------------------------
# 2. compartments


def compartment_consistency(nodes: pd.DataFrame, edges: pd.DataFrame) -> dict:
    compartment = dict(zip(nodes.node_id, nodes.compartment)) if "compartment" in nodes.columns else {}
    is_transport = dict(zip(nodes.node_id, nodes.is_transport)) if "is_transport" in nodes.columns else {}
    reaction_ids = set(nodes.loc[nodes.node_type == "reaction", "node_id"])
    substrates: dict[str, set[str]] = defaultdict(set)
    products: dict[str, set[str]] = defaultdict(set)
    for edge in edges[edges.relation_type.isin(["substrate_of", "product_of"])].itertuples(index=False):
        if edge.relation_type == "substrate_of" and edge.target_id in reaction_ids:
            substrates[edge.target_id].add(str(compartment.get(edge.source_id, "")))
        elif edge.relation_type == "product_of" and edge.source_id in reaction_ids:
            products[edge.source_id].add(str(compartment.get(edge.target_id, "")))
    crossing_without_transport = []
    for reaction in reaction_ids:
        compartments = substrates.get(reaction, set()) | products.get(reaction, set())
        compartments.discard("")
        if len(compartments) > 1 and not bool(is_transport.get(reaction, False)):
            crossing_without_transport.append(reaction)
    metabolites = nodes[nodes.node_type == "metabolite"]
    return {"reactions": len(reaction_ids), "reactions_spanning_compartments_without_transport_flag": len(crossing_without_transport),
            "examples": sorted(crossing_without_transport)[:10],
            "metabolites_by_compartment": metabolites.compartment.value_counts().to_dict() if "compartment" in metabolites.columns else {},
            "metabolites_without_compartment": int(metabolites.compartment.isna().sum()) if "compartment" in metabolites.columns else None}


# ----------------------------------------------------------------------------------------------------------------------
# 3. mass and charge balance


def mass_charge_balance(sbml_path: Path, limit_examples: int = 15) -> dict:
    """Every reaction's elemental and charge balance from the SBML through cobra; boundary, exchange, demand and sink
    reactions are not balanced by construction and are counted apart."""
    import cobra

    model = cobra.io.read_sbml_model(str(sbml_path))
    boundary = {reaction.id for reaction in model.boundary}
    unbalanced, proton_only, charge_only = [], [], []
    for reaction in model.reactions:
        if reaction.id in boundary:
            continue
        try:
            imbalance = reaction.check_mass_balance()
        except Exception as error:  # a metabolite without a formula: reported, not hidden
            unbalanced.append({"reaction": reaction.id, "imbalance": f"unreadable: {error}"})
            continue
        if not imbalance:
            continue
        elements = {element: value for element, value in imbalance.items() if element != "charge"}
        if set(elements) <= {"H"} and "charge" in imbalance and abs(imbalance.get("H", 0) - imbalance.get("charge", 0)) < 1e-6:
            proton_only.append(reaction.id)
        elif not elements:
            charge_only.append(reaction.id)
        else:
            unbalanced.append({"reaction": reaction.id, "imbalance": {key: round(float(value), 3) for key, value in imbalance.items()}})
    internal = len(model.reactions) - len(boundary)
    return {"reactions": len(model.reactions), "boundary_reactions": len(boundary), "internal_reactions": internal,
            "proton_and_charge_only_imbalances": len(proton_only), "charge_only_imbalances": len(charge_only),
            "element_imbalances": len(unbalanced), "element_imbalance_fraction_of_internal": round(len(unbalanced) / internal, 4) if internal else None,
            "examples": unbalanced[:limit_examples], "proton_examples": proton_only[:10]}


# ----------------------------------------------------------------------------------------------------------------------
# 4. the validation set


def read_validation_set(path: Path = DEFAULT_VALIDATION_SET) -> list[dict]:
    with open(path, encoding="utf-8") as validation_file:
        return list(csv.DictReader(validation_file))


def perturbation_targets(perturbation_nodes_json: str) -> dict[str, float]:
    try:
        entries = json.loads(perturbation_nodes_json) if perturbation_nodes_json else []
    except (TypeError, ValueError):
        return {}
    targets: dict[str, float] = {}
    for entry in entries:
        if isinstance(entry, (list, tuple)) and len(entry) >= 2:
            targets[node_symbol(entry[0])] = float(entry[1])
    return targets


def hpo_parents(obo_path: Path) -> dict[str, set[str]]:
    parents: dict[str, set[str]] = defaultdict(set)
    term = None
    with open(obo_path, encoding="utf-8") as obo:
        for line in obo:
            line = line.strip()
            if line == "[Term]":
                term = None
            elif line.startswith("id: HP:"):
                term = line[4:]
            elif line.startswith("is_a: ") and term:
                parents[term].add(line[6:].split(" ")[0])
    return parents


def hpo_ancestors(term: str, parents: dict[str, set[str]]) -> set[str]:
    seen, stack = set(), [term]
    while stack:
        current = stack.pop()
        for parent in parents.get(current, ()):
            if parent not in seen:
                seen.add(parent)
                stack.append(parent)
    return seen


def symptom_term_expansion(crosswalk_path: Path, obo_path: Path) -> dict[str, tuple[set[str], set[str]]]:
    """symptom -> (root term ids, excluded subtree root ids), the crosswalk's reading of HPO."""
    expansion = {}
    with open(crosswalk_path, encoding="utf-8") as crosswalk:
        for row in csv.DictReader(crosswalk):
            roots = {term for term in row.get("hpo_ids", "").split(";") if term}
            excluded = {term for term in row.get("excluded_hpo_ids", "").split(";") if term}
            expansion[row["target_symptom"]] = (roots, excluded)
    return expansion


def term_read_as(term: str, symptom: str, expansion: dict[str, tuple[set[str], set[str]]], parents: dict[str, set[str]]) -> bool | None:
    """Whether the crosswalk reads the HPO term as the symptom: the term or an ancestor is a root and none is an excluded root."""
    if symptom not in expansion or not term:
        return None
    roots, excluded = expansion[symptom]
    lineage = hpo_ancestors(term, parents) | {term}
    return bool(lineage & roots) and not bool(lineage & excluded)


def validation_set_checks(rows: list[dict], nodes: pd.DataFrame, evidence: pd.DataFrame, expansion: dict | None = None, parents: dict | None = None) -> dict:
    gene_symbols = set(nodes.loc[nodes.node_type.isin(["gene", "protein"]), "gene_symbol"].dropna().astype(str))
    labels_by_name: dict[str, pd.DataFrame] = {name: group for name, group in evidence.assign(_name=evidence.perturbation_label.astype(str).str.lower()).groupby("_name")}
    results = []
    for row in rows:
        if row.get("kind") == "withdrawn":
            results.append({**row, "verdict": "withdrawn"})
            continue
        name = row["perturbation"].strip().lower() if row["perturbation_type"] == "drug" else row["perturbation"].strip()
        key = name.lower()
        group = labels_by_name.get(key)
        result = {**row, "target_in_graph": row["target_gene"] in gene_symbols, "perturbation_in_table": group is not None}
        if row.get("hpo_term") and expansion is not None and parents is not None:
            result["term_read_by_crosswalk"] = term_read_as(row["hpo_term"], row["symptom"], expansion, parents)
        if group is not None:
            targets: dict[str, float] = {}
            for value in group.perturbation_nodes.dropna().astype(str).unique():
                targets.update(perturbation_targets(value))
            result["table_targets"] = {symbol: sign for symbol, sign in targets.items()}
            result["target_in_table"] = row["target_gene"] in targets
            expected = EXPECTED_SIGN.get(row["direction"])
            actual = targets.get(row["target_gene"])
            result["sign_agrees"] = None if expected is None or actual is None else (expected * actual > 0)
            pair = group[(group.symptom == row["symptom"]) & (group.relation == (row["relation"] if row["relation"] != "absent" else "induces"))]
            positive_label = pair[(pair.positive_report_count > 0) & pair.grade.isin(LABEL_GRADES)] if len(pair) and "positive_report_count" in pair.columns else pair
            result["pair_rows"] = int(len(pair))
            result["pair_is_positive_label"] = bool(len(positive_label))
            result["pair_negative_reports"] = int(pair.negative_report_count.sum()) if len(pair) and "negative_report_count" in pair.columns else 0
            result["pair_grades"] = sorted(set(pair.grade.astype(str))) if len(pair) else []
            if row["kind"] == "chain":
                if result["pair_is_positive_label"]:
                    result["verdict"] = "supported"
                elif result.get("term_read_by_crosswalk") is False:
                    result["verdict"] = "term_outside_crosswalk"  # the construct decision excludes the term (a developmental intellectual disability is not cognitive impairment)
                else:
                    result["verdict"] = "soft_only" if len(pair) else "missing"
            else:
                result["verdict"] = "contradicted" if result["pair_is_positive_label"] else ("respected_with_negative" if result["pair_negative_reports"] else "respected")
        else:
            result["verdict"] = "perturbation_absent"
        results.append(result)
    summary: dict[str, Counter] = {"chain": Counter(), "negative": Counter(), "withdrawn": Counter()}
    by_tier: dict[str, Counter] = defaultdict(Counter)
    for result in results:
        summary[result["kind"]][result["verdict"]] += 1
        by_tier[f"{result['kind']}:{result['tier']}"][result["verdict"]] += 1
    sign_disagreements = [r for r in results if r.get("sign_agrees") is False]
    targets_missing = [r for r in results if r.get("perturbation_in_table") and not r.get("target_in_table")]
    return {"rows": len(rows), "summary": {kind: dict(counter) for kind, counter in summary.items()}, "by_tier": {key: dict(counter) for key, counter in sorted(by_tier.items())},
            "sign_disagreements": [{"perturbation": r["perturbation"], "target_gene": r["target_gene"], "direction": r["direction"], "table_targets": r.get("table_targets")} for r in sign_disagreements],
            "targets_missing_from_table": [{"perturbation": r["perturbation"], "target_gene": r["target_gene"], "table_targets": r.get("table_targets")} for r in targets_missing],
            "results": results}


# ----------------------------------------------------------------------------------------------------------------------
# 5. coverage


def coverage(evidence: pd.DataFrame) -> dict:
    rows = []
    for symptom, group in evidence.groupby("symptom"):
        labels = group[group.grade.isin(LABEL_GRADES) & (group.positive_report_count > 0)] if "positive_report_count" in group.columns else group
        rows.append({"symptom": symptom, "label_pairs": int(len(labels)), "label_genes": int((labels.perturbation_type == "gene").sum()),
                     "label_drugs": int((labels.perturbation_type == "drug").sum()),
                     "soft_pairs": int((~group.grade.isin(LABEL_GRADES)).sum()),
                     "pairs_with_a_negative_report": int((group.negative_report_count > 0).sum()) if "negative_report_count" in group.columns else 0,
                     "pairs_all_negative": int((group.positive_report_count == 0).sum()) if "positive_report_count" in group.columns else 0,
                     "classes": sorted(set(group.evidence_class.astype(str)))})
    rows.sort(key=lambda row: -row["label_pairs"])
    return {"symptoms": len(rows), "by_symptom": rows, "symptoms_without_drug_labels": [row["symptom"] for row in rows if row["label_drugs"] == 0],
            "symptoms_without_gene_labels": [row["symptom"] for row in rows if row["label_genes"] == 0],
            "symptoms_with_any_negative": sum(1 for row in rows if row["pairs_with_a_negative_report"] or row["pairs_all_negative"])}


def therapeutic_direction(evidence: pd.DataFrame, indications_path: Path, crosswalk_path: Path) -> dict:
    """How much of the therapeutic direction (relation relieves) the symptom crosswalk can label.

    The study predicts symptoms; a drug label states its indications as disorders ("Attention deficit/hyperactivity
    disorder", "Dementia Alzheimer's type"), and the crosswalk maps symptom-level MedDRA preferred terms. So an
    indication reaches a label only where the label happens to name the symptom. This counts the loss at each step, so
    a missing relieves chain in the validation set can be read as a construct decision rather than as a data gap.
    """
    result = {"relieves_pairs_in_the_evidence": int((evidence.relation == "relieves").sum()) if "relation" in evidence.columns else 0}
    if not indications_path.exists() or not crosswalk_path.exists():
        return {**result, "indications_read": False}
    terms: dict[str, str] = {}
    with open(crosswalk_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            for term in str(row.get("meddra_preferred_terms", "")).split(";"):
                if term.strip():
                    terms[term.strip().lower()] = row["target_symptom"]
    rows = []
    with gzip.open(indications_path, "rt", encoding="utf-8") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if len(row) > 6 and row[4] == "PT" and row[2] == "NLP_indication"]
    study_drugs = set(evidence[evidence.perturbation_type == "drug"].perturbation_id) if "perturbation_type" in evidence.columns else set()
    on_a_study_drug = [row for row in rows if row[0] in study_drugs]
    mapped = [row for row in on_a_study_drug if row[6].strip().lower() in terms]
    unmapped_terms = Counter(row[6].strip() for row in on_a_study_drug if row[6].strip().lower() not in terms)
    drugs_without_an_indication = len(study_drugs - {row[0] for row in rows})
    return {**result, "indications_read": True, "crosswalk_preferred_terms": len(terms),
            "indication_rows": len(rows), "indication_rows_on_a_study_drug": len(on_a_study_drug),
            "indication_rows_naming_a_study_symptom": len(mapped),
            "distinct_pairs_naming_a_study_symptom": len({(row[0], terms[row[6].strip().lower()]) for row in mapped}),
            "study_drugs_with_no_indication_row": drugs_without_an_indication, "study_drugs": len(study_drugs),
            "most_common_indication_terms_outside_the_crosswalk": unmapped_terms.most_common(10)}


# ----------------------------------------------------------------------------------------------------------------------
# rendering


RULES_APPLIED = [
    ("Imaging-derived findings are not admitted as perturbation-to-symptom evidence.", "T1-weighted grey-matter volume and cortical thickness move with iron, myelin, hydration, "
     "lithium's T1 shortening and antipsychotic exposure rather than with tissue alone; BOLD depends on neurovascular coupling that drugs such as caffeine change without a "
     "neural change; PET occupancy rests on a reference region free of specific binding. Only PET target occupancy with a valid reference region would qualify, as tier-B "
     "physiology, and none is in the sources."),
    ("A mouse readout enters only with a stated validity.", "configs/impc_symptom_crosswalk.csv names, per symptom, the proxy and how far it stands for the symptom; weak proxies "
     "are excluded unless asked for, prepulse inhibition is never read as psychosis, and the symptoms with no defensible proxy have no rows."),
    ("Spontaneous reports enter at the lowest tier.", "FAERS counts have no exposure denominator and one report names several drugs; a signal is disproportion, an "
     "absence is under-reporting with a large exposure, and neither is a rate."),
    ("A combination arm is attributed to no drug.", "A trial arm naming two study drugs would carry one arm's events into two records across the held-out split, "
     "whose groups are built on targets; such arms are left out of both sides."),
    ("A curated negative interaction is a contested pair, not a deleted edge.", "IntAct's non-interacting pairs include textbook interactions reported absent in one "
     "assay; the audit counts them per source and flags the edges, and nothing is removed on their account."),
]


def render(stamp: str, sections: dict, graph_dir: Path, evidence_dir: Path) -> str:
    integrity, compartments, balance, validation, cover, negatives = (sections[name] for name in ("integrity", "compartments", "balance", "validation", "coverage", "edge_negatives"))
    therapeutic = sections["therapeutic_direction"]
    lines = ["# Physiology audit (generated by experiments/audit_physiology.py — do not edit by hand)", "", stamp, "",
             f"Graph `{graph_dir}`, evidence `{evidence_dir}`. The question asked of every section is whether the graph and the evidence behind it represent the "
             "underlying physiology as truly as the sources allow, which the owner set as the first priority on 15 October 2026.", "",
             "## 1. Graph integrity", "",
             f"{integrity['nodes']:,} nodes ({', '.join(f'{k} {v:,}' for k, v in integrity['nodes_by_type'].items())}) and {integrity['edges']:,} edges "
             f"({', '.join(f'{k} {v:,}' for k, v in integrity['edges_by_relation'].items())}); by source: {', '.join(f'{k} {v:,}' for k, v in integrity['edges_by_source'].items())}.", "",
             f"- exact duplicate edge rows: **{integrity['exact_duplicate_edge_rows']:,}** over {integrity['exact_duplicate_edge_pairs']:,} pairs (a duplicate doubles a node's in-degree in the adjacency)",
             f"- pairs carrying both signs under one relation: **{integrity['contradictory_sign_pairs']:,}**" + (f" (first: {', '.join(e['source_id'] + '->' + e['target_id'] for e in integrity['contradictory_examples'][:5])})" if integrity["contradictory_examples"] else ""),
             f"- unsigned (sign 0) fraction by relation: {', '.join(f'{k} {v:.0%}' for k, v in integrity['unsigned_fraction_by_relation'].items())}",
             f"- self-loops {integrity['self_loops']:,}; edges with an endpoint not in the node table {integrity['dangling_endpoint_edges']:,}; isolated nodes {integrity['isolated_nodes']:,}",
             f"- stoichiometric coefficients in the graph the model reads: **{'present' if integrity['stoichiometry']['column_present'] else 'ABSENT'}**"
             + (f", {integrity['stoichiometry']['edges_with_coefficient_other_than_one']:,} of {integrity['stoichiometry']['metabolic_edges']:,} metabolic edges with a coefficient other than 1"
                if integrity["stoichiometry"]["column_present"] else " (every reaction edge counts as one copy; item 20's conserved-currency arithmetic cannot hold)"),
             "", "Highest-degree nodes (a currency metabolite should be tagged; an untagged hub is a short cut through the graph):", "",
             "| node | name | degree | currency |", "| --- | --- | --- | --- |"]
    for node in integrity["top_degree_nodes"]:
        lines.append(f"| {node['node']} | {node['name']} | {node['degree']:,} | {'yes' if node['is_currency'] else 'no'} |")
    lines += ["", "## 2. Compartments", "",
              f"{compartments['reactions']:,} reactions; **{compartments['reactions_spanning_compartments_without_transport_flag']:,}** have substrates and products in different "
              f"compartments without being flagged as transport" + (f" (first: {', '.join(compartments['examples'][:6])})" if compartments["examples"] else "") + ". "
              f"Metabolites by compartment: {', '.join(f'{k} {v:,}' for k, v in compartments['metabolites_by_compartment'].items())}; without a compartment: {compartments['metabolites_without_compartment']}.", ""]
    lines += ["## 3. Mass and charge balance (Human-GEM SBML through cobra)", ""]
    if balance.get("skipped"):
        lines.append(f"Skipped: {balance['skipped']}.")
    else:
        lines += [f"{balance['reactions']:,} reactions, {balance['boundary_reactions']:,} boundary (not balanced by construction), {balance['internal_reactions']:,} internal. "
                  f"Element imbalances: **{balance['element_imbalances']:,}** ({balance['element_imbalance_fraction_of_internal']:.2%} of internal); proton-and-charge-only "
                  f"imbalances {balance['proton_and_charge_only_imbalances']:,} (a protonation-state convention, tolerable); charge-only {balance['charge_only_imbalances']:,}.", ""]
        if balance["examples"]:
            lines += ["| reaction | imbalance |", "| --- | --- |"]
            for example in balance["examples"]:
                lines.append(f"| {example['reaction']} | {example['imbalance']} |")
    lines += ["", "## 4. The cited validation set", "",
              f"{validation['rows']} rows of `configs/physiology_validation_set.csv` (chains the evidence should hold, negatives it should not contradict; the citations are "
              "the research pass's of 15 October 2026 and the ones marked as pharmacology texts are for the owner to replace with a primary source).", "",
              f"Chains: {', '.join(f'{k} {v}' for k, v in validation['summary']['chain'].items())}. Negatives: {', '.join(f'{k} {v}' for k, v in validation['summary']['negative'].items())}. "
              f"Withdrawn after the first run (the label records the event the research pass called absent): {sum(validation['summary']['withdrawn'].values())}; "
              "a chain whose HPO term the crosswalk excludes is reported as term_outside_crosswalk, a construct decision and not a missing label.", "",
              "By tier: " + "; ".join(f"{key}: {', '.join(f'{k} {v}' for k, v in counts.items())}" for key, counts in validation["by_tier"].items()) + ".", ""]
    if validation["sign_disagreements"]:
        lines += ["**Sign disagreements** (the table's perturbation sign at the target is opposite to the established direction):", ""]
        for item in validation["sign_disagreements"]:
            lines.append(f"- {item['perturbation']} at {item['target_gene']} ({item['direction']}): table targets {item['table_targets']}")
        lines.append("")
    if validation["targets_missing_from_table"]:
        lines += ["**Targets the table does not carry for a perturbation it holds** (the drug is in the table with other targets):", ""]
        for item in validation["targets_missing_from_table"]:
            lines.append(f"- {item['perturbation']}: expected {item['target_gene']}, table targets {item['table_targets']}")
        lines.append("")
    lines += ["| kind | perturbation | target | symptom | relation | tier | verdict | pair grades | negatives |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for result in validation["results"]:
        lines.append(f"| {result['kind']} | {result['perturbation']} | {result['target_gene']} | {result['symptom']} | {result['relation']} | {result['tier']} | "
                     f"**{result['verdict']}** | {', '.join(result.get('pair_grades', []))} | {result.get('pair_negative_reports', 0)} |"
                     + (f" {result['status']}" if result.get("status") else ""))
    lines += ["", "## 5. Evidence coverage per symptom", "",
              f"{cover['symptoms']} symptoms; without a drug label: {', '.join(cover['symptoms_without_drug_labels']) or 'none'}; without a gene label: "
              f"{', '.join(cover['symptoms_without_gene_labels']) or 'none'}; with any negative report: {cover['symptoms_with_any_negative']}.", "",
              "| symptom | label pairs | genes | drugs | soft pairs | pairs with a negative | all-negative pairs | classes |", "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in cover["by_symptom"]:
        lines.append(f"| {row['symptom']} | {row['label_pairs']} | {row['label_genes']} | {row['label_drugs']} | {row['soft_pairs']} | {row['pairs_with_a_negative_report']} | "
                     f"{row['pairs_all_negative']} | {', '.join(row['classes'])} |")
    lines += ["", "## 6. The therapeutic direction: how much of it the crosswalk can label", ""]
    if not therapeutic.get("indications_read"):
        lines.append(f"{therapeutic['relieves_pairs_in_the_evidence']:,} relieves pairs in the evidence; the indication table or the crosswalk is not on disk, so the loss is not measured.")
    else:
        lines += [f"{therapeutic['relieves_pairs_in_the_evidence']:,} pairs of the evidence table carry relation `relieves`. SIDER states an indication as a MedDRA preferred term, and the "
                  f"crosswalk maps {therapeutic['crosswalk_preferred_terms']} symptom-level terms, so an indication reaches a label only where the label happens to name the symptom rather "
                  "than the disorder:", "",
                  f"- NLP_indication rows at preferred-term level: **{therapeutic['indication_rows']:,}**, of which **{therapeutic['indication_rows_on_a_study_drug']:,}** are on one of the "
                  f"{therapeutic['study_drugs']} study drugs",
                  f"- of those, **{therapeutic['indication_rows_naming_a_study_symptom']:,}** name a term the crosswalk maps, giving "
                  f"**{therapeutic['distinct_pairs_naming_a_study_symptom']:,}** distinct (drug, symptom) pairs",
                  f"- **{therapeutic['study_drugs_with_no_indication_row']}** study drugs have no indication row in SIDER at all (the release predates several of them)", "",
                  "The commonest indication terms the crosswalk does not map are disorders, not symptoms: "
                  + ", ".join(f"{term} ({count})" for term, count in therapeutic["most_common_indication_terms_outside_the_crosswalk"][:8]) + ". "
                  "So a `relieves` chain the validation set states and the table lacks (methylphenidate for hyperactivity, donepezil for cognitive impairment) is a construct decision of the "
                  "same kind as `term_outside_crosswalk` above, not a gap in the source; a chain whose drug has no indication row at all (ketamine, esketamine, phenelzine, suvorexant) is a "
                  "gap in SIDER 4.1, which the trial-arm table can speak to and the label sources cannot."]
    lines += ["", "## 7. Curated non-interacting pairs against the protein edges", ""]
    if negatives.get("skipped"):
        lines.append(f"Skipped: {negatives['skipped']}.")
    else:
        lines += [f"{negatives['human_negative_pairs']:,} human pairs from IntAct name {negatives['proteins_named']:,} proteins; **{negatives['edges_on_a_negative_pair']:,}** graph edges "
                  "fall on a pair curators recorded as non-interacting. Per source, against the edges that join two named proteins:", "",
                  "| source | edges on named proteins | on a negative pair | fraction |", "| --- | --- | --- | --- |"]
        for source, counts in negatives["by_source"].items():
            fraction = counts["edges_on_a_negative_pair"] / counts["edges_on_named_proteins"] if counts["edges_on_named_proteins"] else 0.0
            lines.append(f"| {source} | {counts['edges_on_named_proteins']:,} | {counts['edges_on_a_negative_pair']:,} | {fraction:.1%} |")
        lines += ["", "A curated negative is one assay's null in one context, and several of the flagged pairs are textbook interactions (MYD88-IRAK4, BAD-BCL2L1 on the "
                  "release graph), so the figure bounds a source's contradiction rate on the pairs the set can speak to and removes nothing; the flagged edges are "
                  "listed in the JSON beside this document as contested."]
    lines += ["", "## 8. Rules applied to what did not enter", ""]
    for rule, why in RULES_APPLIED:
        lines.append(f"- **{rule}** {why}")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders_drugs_v4"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--sbml", type=Path, default=Path("data/raw/Human-GEM/model/Human-GEM.xml"))
    parser.add_argument("--intact", type=Path, default=Path("data/raw/intact/intact_negative.txt"))
    parser.add_argument("--validation-set", type=Path, default=DEFAULT_VALIDATION_SET)
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--sider-indications", type=Path, default=Path("data/raw/sider_4.1/meddra_all_indications.tsv.gz"))
    parser.add_argument("--hpo-obo", type=Path, default=Path("data/raw/hpo/hp.obo"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/audit"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/physiology_audit.md"))
    parser.add_argument("--skip-balance", action="store_true", help="skip the cobra mass-balance pass (minutes)")
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet")
    evidence = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet")
    sections = {"integrity": graph_integrity(nodes, edges), "compartments": compartment_consistency(nodes, edges)}
    if arguments.skip_balance or not arguments.sbml.exists():
        sections["balance"] = {"skipped": "--skip-balance" if arguments.skip_balance else f"{arguments.sbml} is not on disk"}
    else:
        sections["balance"] = mass_charge_balance(arguments.sbml)
    expansion = symptom_term_expansion(arguments.crosswalk, arguments.hpo_obo) if arguments.crosswalk.exists() and arguments.hpo_obo.exists() else None
    parents = hpo_parents(arguments.hpo_obo) if expansion is not None else None
    sections["validation"] = validation_set_checks(read_validation_set(arguments.validation_set), nodes, evidence, expansion, parents)
    sections["coverage"] = coverage(evidence)
    sections["therapeutic_direction"] = therapeutic_direction(evidence, arguments.sider_indications, arguments.crosswalk)
    if arguments.intact.exists():
        from mechanistic_pathway_learning.graph.edge_negatives import audit as negatives_audit, read_human_negative_pairs
        sections["edge_negatives"] = negatives_audit(read_human_negative_pairs(arguments.intact), edges, nodes)
    else:
        sections["edge_negatives"] = {"skipped": f"{arguments.intact} is not on disk"}
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    (arguments.output_dir / "physiology_audit.json").write_text(json.dumps(sections, indent=1, default=str) + "\n", encoding="utf-8")
    stamp = provenance_stamp("experiments/audit_physiology.py", [str(arguments.graph_dir), str(arguments.evidence_dir), str(arguments.validation_set), str(arguments.intact)])
    arguments.markdown_output.write_text(render(stamp, sections, arguments.graph_dir, arguments.evidence_dir), encoding="utf-8")
    print(f"wrote {arguments.markdown_output} and {arguments.output_dir / 'physiology_audit.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
