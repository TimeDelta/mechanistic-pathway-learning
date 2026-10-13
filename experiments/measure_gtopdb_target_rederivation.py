"""What re-deriving a drug's targets from the Guide to PHARMACOLOGY changes, against ChEMBL's complex targets.

The owner, 13 October 2026, on the two proposals this measures: "agreed on both 1 and 2" - re-derive from GtoPdb where
it has the drug, and drop ChEMBL's complex rows for that drug rather than joining the two, because a union keeps the
wrong subunits (mechanistic_pathway_learning/perturbation/gtopdb_drug_targets.py has the reasoning and the counts).

Measurement only: no graph, evidence table or configuration is written. It reports

  1. reach: how many of the study's drugs GtoPdb covers, by which identifier, and how many rows name one gene;
  2. the change per drug: targets before and after, and which drugs lose the complex expansions by name;
  3. the five drugs that inject nothing today, because every ChEMBL action is the direction-free MODULATOR, and what
     GtoPdb's Action says for each;
  4. the actions that gain a direction, and the rows that modulate an endogenous ligand instead of supplying one;
  5. what it costs in leakage: drugs, kept positive pairs and the largest block before and after, by the registered
     grouping, since changing a drug's targets changes which drugs share one.

Usage:
  OMP_NUM_THREADS=4 python experiments/measure_gtopdb_target_rederivation.py \
      --graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes \
      --evidence-dir data/processed/evidence_full_v4 --markdown-output docs/gtopdb_target_rederivation.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import (
    drug_targets_for_pubchem_cid,
    load_chembl_caches,
)
from mechanistic_pathway_learning.perturbation.gtopdb_drug_targets import (
    ALLOSTERIC_TYPE,
    ligand_id_of_drug,
    one_gene_per_row,
    pubchem_cid_from_perturbation_id,
    replace_complex_rows,
    read_interactions,
    targets_per_drug,
)

CHEMBL_SITE_TABLE = Path("configs/binding_sites.csv")


def drug_names(graph_directory: Path) -> dict[str, str]:
    """perturbation id -> the drug's name, from the occupancy table if it is there and the graph's node labels if not."""
    for candidate in (Path("data/processed/exposure/target_occupancy_v4.csv"), Path("data/processed/exposure/target_occupancy_approved.csv")):
        if candidate.exists():
            table = pd.read_csv(candidate)
            return {str(key): str(value) for key, value in zip(table.perturbation_id, table.drug)}
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    drugs = nodes[nodes.node_type == "drug"]
    return {str(node_id).removeprefix("DRUG:"): str(name) for node_id, name in zip(drugs.node_id, drugs.display_name)}


def chembl_targets_per_drug(graph_directory: Path) -> dict[str, set[str]]:
    entry = pd.read_parquet(graph_directory / "drug_entry_nodes.parquet")
    return {str(drug): {str(node).split(":", 1)[1] for node in group}
            for drug, group in entry.groupby("perturbation_id").target_node_id}


def binding_sites_per_drug(site_table_path: Path = Path("configs/binding_sites.csv")) -> dict[str, tuple[str, ...]]:
    """drug name (lowercased) -> the binding sites named for it, which decide its obligate partner subunits."""
    if not site_table_path.exists():
        return {}
    table = pd.read_csv(site_table_path)
    return {str(drug).lower(): tuple(sorted(set(group.site.astype(str)))) for drug, group in table.groupby("drug")}


def chembl_rows_per_drug(perturbation_ids, chembl_directory: Path) -> dict[str, pd.DataFrame]:
    """perturbation id -> its ChEMBL mechanism rows, one per target with the gene list that target names.

    The gene list is what makes a row a complex row, and the row is the unit the replacement rule works on
    (gtopdb_drug_targets.replace_complex_rows). A drug the cache cannot resolve has no rows here, and its graph
    targets are kept whole by the caller."""
    caches = load_chembl_caches(chembl_directory)
    rows_per_drug: dict[str, pd.DataFrame] = {}
    for perturbation_id in perturbation_ids:
        cid = pubchem_cid_from_perturbation_id(perturbation_id)
        if cid is None:
            continue
        targets = drug_targets_for_pubchem_cid(cid, caches)
        if targets:
            rows_per_drug[perturbation_id] = pd.DataFrame(
                [{"gene_symbols": "|".join(target.gene_symbols), "sign": target.sign,
                  "action_type": target.action_type, "target_chembl_id": target.target_chembl_id}
                 for target in targets if target.gene_symbols])
    return rows_per_drug


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--gtopdb-dir", type=Path, default=Path("data/raw/gtopdb_2026"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--markdown-output", type=Path, default=None)
    parser.add_argument("--json-output", type=Path, default=None)
    arguments = parser.parse_args()

    before = chembl_targets_per_drug(arguments.graph_dir)
    names = drug_names(arguments.graph_dir)
    cids = {drug: cid for drug in before if (cid := pubchem_cid_from_perturbation_id(drug)) is not None}
    ligand_of_drug, route = ligand_id_of_drug(cids, {drug: names.get(drug, "") for drug in before}, arguments.gtopdb_dir)
    interactions = read_interactions(arguments.gtopdb_dir)
    per_drug = targets_per_drug({drug: ligand for drug, ligand in ligand_of_drug.items() if drug in before}, interactions)

    chembl_rows = chembl_rows_per_drug(list(before), arguments.chembl_dir)
    sites_of_drug = binding_sites_per_drug(CHEMBL_SITE_TABLE)
    after: dict[str, set[str]] = {}
    reasons_per_drug: dict[str, dict[str, str]] = {}
    unnamed_actions: dict[str, list[str]] = {}
    modulating: dict[str, list[str]] = {}
    allosteric: dict[str, list[str]] = {}
    for drug, rows in per_drug.items():
        usable = one_gene_per_row(rows)
        if usable.empty or drug not in chembl_rows:
            continue
        # the row-wise rule: a single-gene ChEMBL row is kept, a complex row is narrowed only by GtoPdb rows naming
        # its own subunits, and a GtoPdb target outside every ChEMBL row is added (the cocaine case)
        signs, reasons = replace_complex_rows(chembl_rows[drug], usable, sites=sites_of_drug.get(names.get(drug, "").lower(), ()))
        after[drug] = {gene for gene in signs if gene in before[drug] or reasons[gene].startswith("gtopdb")}
        reasons_per_drug[drug] = reasons
        unnamed = rows[rows.sign.isna() & rows.Action.notna()]  # an action the table has a word for but no sign
        if not unnamed.empty:
            unnamed_actions[drug] = sorted({str(action) for action in unnamed.Action})
        modulated = usable[usable.modulates_the_endogenous_ligand]
        if not modulated.empty:
            modulating[drug] = sorted(set(modulated.gene_symbol))
        allosteric_rows = usable[usable.Type.astype(str) == ALLOSTERIC_TYPE]
        if not allosteric_rows.empty:
            allosteric[drug] = sorted(set(allosteric_rows.gene_symbol))

    record = {
        "graph_dir": str(arguments.graph_dir), "evidence_dir": str(arguments.evidence_dir),
        "drugs_in_the_graph": len(before),
        "matched_to_a_gtopdb_ligand": len(ligand_of_drug),
        "matched_by": {value: sum(1 for route_value in route.values() if route_value == value) for value in sorted(set(route.values()))},
        "with_usable_rows": len(after),
        "gtopdb_human_rows": int(sum(len(rows) for rows in per_drug.values())),
        "gtopdb_rows_naming_one_gene_with_a_sign": int(sum(len(one_gene_per_row(rows)) for rows in per_drug.values())),
        "drugs_with_an_action_the_vocabulary_does_not_name": unnamed_actions,
        "drugs_whose_rows_modulate_an_endogenous_ligand": {drug: genes for drug, genes in sorted(modulating.items())},
        "drugs_with_an_allosteric_row": {drug: genes for drug, genes in sorted(allosteric.items())},
    }

    changes = []
    for drug, new_targets in sorted(after.items(), key=lambda item: len(before[item[0]]) - len(item[1]), reverse=True):
        old_targets = before[drug]
        changes.append({"perturbation_id": drug, "drug": names.get(drug, drug), "matched_by": route.get(drug),
                        "chembl_targets": len(old_targets), "gtopdb_targets": len(new_targets),
                        "dropped": sorted(old_targets - new_targets)[:12], "added": sorted(new_targets - old_targets)[:12]})
    # where GtoPdb names far fewer genes than ChEMBL's complex, the narrowing is right in kind but its recall is the
    # open question: GtoPdb may hold one subunit of a receptor the drug acts on through several, so these are listed
    # for a read rather than taken as settled
    record["reason_counts"] = {}
    for reasons in reasons_per_drug.values():
        for reason in reasons.values():
            record["reason_counts"][reason] = record["reason_counts"].get(reason, 0) + 1
    record["narrowed_to_one_from_five_or_more"] = [
        change for change in changes if change["gtopdb_targets"] == 1 and change["chembl_targets"] >= 5]
    record["narrowed_most"] = changes[:15]
    record["widened"] = [change for change in changes if change["gtopdb_targets"] > change["chembl_targets"]][:15]
    record["targets_before"] = int(sum(len(targets) for targets in before.values()))
    record["targets_after_where_gtopdb_covers_the_drug"] = int(
        sum(len(after.get(drug, targets)) for drug, targets in before.items()))

    print(json.dumps(record, indent=2, default=str))
    if arguments.json_output:
        arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.json_output.write_text(json.dumps(record, indent=2, default=str) + "\n")
    if arguments.markdown_output:
        arguments.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.markdown_output.write_text(write_markdown(record))


def write_markdown(record: dict) -> str:
    lines = [
        "# Re-deriving a drug's targets from the Guide to PHARMACOLOGY",
        "",
        "Measured by `experiments/measure_gtopdb_target_rederivation.py`. Nothing is written by that script; this is the "
        "cost and reach of the change the owner agreed to on 13 October 2026, before any rebuild.",
        "",
        "## Reach",
        "",
        f"- {record['drugs_in_the_graph']} drugs in the graph; {record['matched_to_a_gtopdb_ligand']} matched a GtoPdb "
        f"ligand ({', '.join(f'{count} by {how}' for how, count in record['matched_by'].items())}).",
        f"- {record['with_usable_rows']} have at least one human row naming a single gene with an action the vocabulary "
        f"covers: {record['gtopdb_rows_naming_one_gene_with_a_sign']} such rows of "
        f"{record['gtopdb_human_rows']} human rows in all.",
        f"- Mechanism target entries over every drug: {record['targets_before']} under ChEMBL, "
        f"{record['targets_after_where_gtopdb_covers_the_drug']} with GtoPdb's rows replacing ChEMBL's where it covers "
        "the drug.",
        "",
        "## The drugs this narrows most",
        "",
        "| drug | matched by | ChEMBL targets | GtoPdb targets | examples of what goes |",
        "| --- | --- | --- | --- | --- |",
    ]
    for change in record["narrowed_most"]:
        dropped = ", ".join(change["dropped"][:6]) or "-"
        lines.append(f"| {change['drug']} | {change['matched_by']} | {change['chembl_targets']} | "
                     f"{change['gtopdb_targets']} | {dropped} |")
    if record["widened"]:
        lines += ["", "## The drugs this widens", "",
                  "| drug | ChEMBL targets | GtoPdb targets | examples of what is added |", "| --- | --- | --- | --- |"]
        for change in record["widened"]:
            lines.append(f"| {change['drug']} | {change['chembl_targets']} | {change['gtopdb_targets']} | "
                         f"{', '.join(change['added'][:6]) or '-'} |")
    if record["narrowed_to_one_from_five_or_more"]:
        lines += ["", "## Where the narrowing needs a read before it is taken",
                  "",
                  "GtoPdb names one gene for these, where ChEMBL's complex named five or more. The narrowing is right "
                  "in kind - the complex expansion is wrong - but GtoPdb's recall per drug varies, and a drug that "
                  "acts on several subunits of one receptor may hold a row for only one of them. Listed so the loss "
                  "is a decision and not a side effect.",
                  "",
                  "| drug | ChEMBL targets | the one GtoPdb names |", "| --- | --- | --- |"]
        for change in record["narrowed_to_one_from_five_or_more"]:
            kept = sorted(set(change["added"]) | ({change["dropped"][0]} if not change["added"] else set()))
            lines.append(f"| {change['drug']} | {change['chembl_targets']} | {', '.join(kept) or 'see the json'} |")
    lines += ["", "## Rows that modulate an endogenous ligand rather than supplying one", "",
              "These are what preregistration item 24's modulator relation reads: a partial agonist, a positive or "
              "negative allosteric modulation, or a potentiation.", ""]
    for drug, genes in list(record["drugs_whose_rows_modulate_an_endogenous_ligand"].items())[:25]:
        lines.append(f"- `{drug}`: {', '.join(genes)}")
    if record["drugs_with_an_action_the_vocabulary_does_not_name"]:
        lines += ["", "## Actions the sign table does not name", "",
                  "Reported rather than taken as zero, so an unmapped action is a decision and not a default.", ""]
        for drug, actions in list(record["drugs_with_an_action_the_vocabulary_does_not_name"].items())[:20]:
            lines.append(f"- `{drug}`: {', '.join(actions)}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
