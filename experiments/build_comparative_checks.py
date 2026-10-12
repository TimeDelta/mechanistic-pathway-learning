"""Write the table of comparative checks: a drug with a second, varied factor against the drug alone.

Four sources, one row per statement (mechanistic_pathway_learning/evaluation/comparative_checks.py says how a row is
put to a model, and that nothing is trained on one):

  carrier_level     configs/carrier_and_displacement_statements.csv, written by hand from drug labels and two hospital
  displacement      cohorts, each row with its quotation and where it was read (docs/label_source_review.md, section 2)
  enzyme_activity   configs/drug_clearing_enzymes.csv, the rows of the FDA Table of Pharmacogenetic Associations that
                    give a clearance edge (experiments/build_drug_clearing_enzymes.py)
  second_drug       the CRESCENDDI reference set of drug-drug-event triples (Kontsioti et al. 2022,
                    doi:10.1038/s41597-022-01159-y; the two files are released under CC0), where both drugs are study
                    drugs and the event maps to a study symptom. A positive control gives two rows, one with each drug
                    as the first, so the pair is compared with either drug alone; a negative control is reported and
                    not judged. The event's MedDRA name is not copied, only its code and the study symptom, because
                    MedDRA terms are not redistributed from this repository.

The drugs are those of --evidence-dir; a statement that names a drug outside it is listed and left out.

Usage:
  PYTHONPATH=. python experiments/build_comparative_checks.py --evidence-dir data/processed/evidence_full_v4
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evaluation.comparative_checks import CHECK_COLUMNS, carrier_genes

STATEMENTS = Path("configs/carrier_and_displacement_statements.csv")
CLEARING_ENZYMES = Path("configs/drug_clearing_enzymes.csv")
CRESCENDDI_DIRECTORY = Path("data/raw/label_review/crescenddi")
CROSSWALK = Path("docs/symptom_crosswalk.csv")
OUTPUT = Path("configs/comparative_checks.csv")
CRESCENDDI_SOURCE = "CRESCENDDI reference set of drug-drug interactions (Kontsioti et al. 2022)"
CRESCENDDI_IDS = {"positive": "doi:10.6084/m9.figshare.14847093.v1", "negative": "doi:10.6084/m9.figshare.14847099.v1"}


def study_drugs(evidence_directory: Path) -> tuple[dict[str, str], dict[str, set[tuple[str, float]]]]:
    """(drug name in lower case -> perturbation id, perturbation id -> its (target node, sign) pairs)."""
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet", columns=["perturbation_id", "perturbation_type", "perturbation_label", "perturbation_nodes"])
    drugs = records[records.perturbation_type == "drug"].drop_duplicates("perturbation_id")
    perturbation_id_of = {label.lower(): perturbation_id for label, perturbation_id in zip(drugs.perturbation_label, drugs.perturbation_id)}
    targets_of = {perturbation_id: {(triple[0], float(triple[1])) for triple in json.loads(nodes)} for perturbation_id, nodes in zip(drugs.perturbation_id, drugs.perturbation_nodes)}
    return perturbation_id_of, targets_of


def carrier_and_displacement_rows(statements: pd.DataFrame, perturbation_id_of: dict[str, str]) -> tuple[list[dict], list[str]]:
    rows, left_out = [], []
    for statement in statements.itertuples(index=False):
        names = [statement.drug] + ([statement.displacing_drug] if statement.displacing_drug else [])
        missing = [name for name in names if name.lower() not in perturbation_id_of]
        if missing:
            left_out.append(f"{statement.statement_id}: {', '.join(missing)} not a study drug")
            continue
        rows.append({"statement_id": statement.statement_id, "kind": statement.kind, "drug_perturbation_id": perturbation_id_of[statement.drug.lower()], "drug": statement.drug,
                     "lowered_genes": ";".join(carrier_genes(statement.carrier)),
                     "second_drug_perturbation_id": perturbation_id_of[statement.displacing_drug.lower()] if statement.displacing_drug else "", "second_drug": statement.displacing_drug,
                     "symptom": statement.symptom, "expectation": "at_least_as_high", "source": statement.source, "source_id": statement.source_id,
                     "source_version": statement.source_version, "quote": statement.quote, "second_mechanism_or_limit": statement.second_mechanism_or_limit})
    return rows, left_out


def enzyme_rows(clearing_enzymes: pd.DataFrame, perturbation_id_of: dict[str, str]) -> list[dict]:
    rows = []
    with_edge = clearing_enzymes[clearing_enzymes.gives_an_edge.astype(str).str.lower() == "true"]
    for enzyme in with_edge.itertuples(index=False):
        if enzyme.drug.lower() not in perturbation_id_of:
            continue
        rows.append({"statement_id": f"enzyme_{enzyme.drug.lower().replace(' ', '_')}_{enzyme.enzyme_gene}", "kind": "enzyme_activity",
                     "drug_perturbation_id": perturbation_id_of[enzyme.drug.lower()], "drug": enzyme.drug, "lowered_genes": enzyme.enzyme_gene,
                     "second_drug_perturbation_id": "", "second_drug": "", "symptom": "", "expectation": "at_least_as_high",
                     "source": f"FDA Table of Pharmacogenetic Associations, section {enzyme.fda_table_section}", "source_id": enzyme.source,
                     "source_version": f"page saved {enzyme.page_saved_on}", "quote": f"{enzyme.fda_table_drug}, {enzyme.fda_table_gene}, {enzyme.fda_table_subgroup}: {enzyme.fda_table_description}",
                     "second_mechanism_or_limit": "The row gives a direction for concentrations or adverse reactions as a whole and names no symptom, so the check reads the drug's own kept positive symptoms."})
    return rows


def second_drug_rows(crescenddi_directory: Path, crosswalk: pd.DataFrame, perturbation_id_of: dict[str, str], targets_of: dict[str, set[tuple[str, float]]]) -> list[dict]:
    symptom_of_term = {term.strip().lower(): row.target_symptom for row in crosswalk.itertuples() for term in str(row.meddra_preferred_terms).split(";") if term.strip() and term.strip() != "nan"}
    rows = []
    for polarity, file_name in (("positive", "positive_controls.xlsx"), ("negative", "negative_controls.xlsx")):
        controls = pd.read_excel(crescenddi_directory / file_name)
        controls = controls.assign(first=controls.DRUG_1_CONCEPT_NAME.str.lower(), second=controls.DRUG_2_CONCEPT_NAME.str.lower(),
                                   symptom=controls.EVENT_CONCEPT_NAME.str.lower().map(symptom_of_term))
        controls = controls[controls["first"].isin(perturbation_id_of) & controls["second"].isin(perturbation_id_of) & controls.symptom.notna()]
        seen = set()
        for control in controls.itertuples(index=False):
            pair = tuple(sorted((control.first, control.second)))
            if (pair, control.symptom) in seen or pair[0] == pair[1]:  # several event terms of one symptom, or one drug under two names
                continue
            seen.add((pair, control.symptom))
            first_targets, second_targets = targets_of[perturbation_id_of[pair[0]]], targets_of[perturbation_id_of[pair[1]]]
            shared_nodes = {node for node, _ in first_targets} & {node for node, _ in second_targets}
            same_sign = bool(first_targets & second_targets)
            relation = "share no target" if not shared_nodes else ("share a target with the same sign" if same_sign and len(first_targets & second_targets) == len(shared_nodes) else "share a target with opposite or mixed signs")
            evidence = "" if polarity == "negative" else "; ".join(
                f"{name} {value}" for name, value in (("event source", control.EVENT_SOURCE), ("Micromedex evidence level", control.MICROMEDEX_EVID_LEVEL), ("BNF evidence level", control.BNF_EVID_LEVEL))
                if isinstance(value, str) and value)
            for order, (first, second) in enumerate((pair, pair[::-1])):
                if polarity == "negative" and order == 1:
                    continue  # a negative control is reported once
                rows.append({"statement_id": f"pair_{polarity}_{first.replace(' ', '_')}_with_{second.replace(' ', '_')}_{control.symptom}", "kind": "second_drug",
                             "drug_perturbation_id": perturbation_id_of[first], "drug": first, "lowered_genes": "",
                             "second_drug_perturbation_id": perturbation_id_of[second], "second_drug": second, "symptom": control.symptom,
                             "expectation": "at_least_as_high" if polarity == "positive" else "none", "source": f"{CRESCENDDI_SOURCE}, {polarity} controls",
                             "source_id": CRESCENDDI_IDS[polarity], "source_version": "files downloaded 2026-10-10",
                             "quote": f"{polarity} control: {pair[0]} with {pair[1]}, MedDRA code {int(control.MDR_CODE)} ({control.symptom})" + (f"; {evidence}" if evidence else ""),
                             "second_mechanism_or_limit": f"The two drugs {relation}. " + ("A pair that shares a target with the same sign adds inputs at one node, so the comparison holds for any encoder that is monotone in its input."
                                                                                       if relation.endswith("same sign") else "")})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--statements", type=Path, default=STATEMENTS)
    parser.add_argument("--clearing-enzymes", type=Path, default=CLEARING_ENZYMES)
    parser.add_argument("--crescenddi-dir", type=Path, default=CRESCENDDI_DIRECTORY)
    parser.add_argument("--crosswalk", type=Path, default=CROSSWALK)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    perturbation_id_of, targets_of = study_drugs(arguments.evidence_dir)
    rows, left_out = carrier_and_displacement_rows(pd.read_csv(arguments.statements, dtype=str).fillna(""), perturbation_id_of)
    rows += enzyme_rows(pd.read_csv(arguments.clearing_enzymes, dtype=str).fillna(""), perturbation_id_of)
    rows += second_drug_rows(arguments.crescenddi_dir, pd.read_csv(arguments.crosswalk), perturbation_id_of, targets_of)
    table = pd.DataFrame(rows, columns=list(CHECK_COLUMNS))
    if table.statement_id.duplicated().any():
        raise SystemExit(f"duplicate statement ids: {sorted(table.statement_id[table.statement_id.duplicated()])[:5]}")
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(arguments.output, index=False)
    print(f"wrote {arguments.output}: {len(table)} checks, by kind and expectation {table.groupby(['kind', 'expectation']).size().to_dict()}")
    for line in left_out:
        print("left out:", line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
