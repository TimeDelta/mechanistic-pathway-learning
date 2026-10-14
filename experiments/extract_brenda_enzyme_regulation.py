"""Pull BRENDA's enzyme-regulation records for human proteins into one table.

Decision 9.14 and docs/missing_physiology_review.md section 2: the graph can receive a metabolite response but
almost cannot act on one, 600 of 266,966 edges running metabolite to protein and 4 `inhibits` in the whole graph.
BRENDA's inhibitor, activating-compound, Ki, IC50, cofactor and metals/ions fields are that missing layer.

Sign convention matches the graph: an inhibitor arrives at -1, an activating compound at +1, and a cofactor or a
metal ion carries no sign because `cofactor_of` is unsigned.
"""
import argparse
import collections
import csv
import json
import re
from pathlib import Path

import ijson

REGULATION_FIELDS_AND_RELATIONS = {
    "inhibitor": ("inhibits", -1.0),
    "ki_value": ("inhibits", -1.0),
    "ic50_value": ("inhibits", -1.0),
    "activating_compound": ("activates", 1.0),
    "cofactor": ("cofactor_of", 0.0),
    "metals_ions": ("cofactor_of", 0.0),
}
UNIPROT_ACCESSION_PATTERN = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\b")


def human_accessions_of_protein_entry(protein_entry: dict) -> list[str]:
    """BRENDA keys `protein` by index; each entry carries `organism` and, for 40 percent of them, `accessions`.

    Only Homo sapiens entries are kept, because the graph is human. An entry with no accession still identifies the
    enzyme through the EC number, which the caller joins to the local UniProt dump's EC column, so returning an
    empty list here does not lose the record - it moves it to the EC route.
    """
    if "homo sapiens" not in str(protein_entry.get("organism", "")).lower():
        return []
    accession_field = protein_entry.get("accessions") or ""
    text_to_search = accession_field if isinstance(accession_field, str) else " ".join(map(str, accession_field))
    return sorted({match.group(1) for match in UNIPROT_ACCESSION_PATTERN.finditer(text_to_search)})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("brenda_json", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    rows_written = 0
    counts_by_relation = collections.Counter()
    distinct_ligand_names: set[str] = set()
    distinct_human_accessions: set[str] = set()
    enzyme_entries_seen = 0

    with open(arguments.brenda_json, "rb") as handle, open(arguments.output, "w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=["ec_number", "field", "relation_type", "sign", "ligand_name",
                                                         "uniprot_accession", "organism_entry", "brenda_comment",
                                                         "brenda_reference_count"])
        writer.writeheader()
        for ec_number, enzyme_record in ijson.kvitems(handle, "data"):
            enzyme_entries_seen += 1
            protein_entries = enzyme_record.get("protein") or {}
            accessions_by_protein_index = {}
            organism_by_protein_index = {}
            human_protein_indices_without_an_accession = set()
            for protein_index, protein_entry in protein_entries.items():
                if not isinstance(protein_entry, dict):
                    continue
                protein_index = str(protein_index)
                is_human = "homo sapiens" in str(protein_entry.get("organism", "")).lower()
                accessions = human_accessions_of_protein_entry(protein_entry)
                if accessions:
                    accessions_by_protein_index[protein_index] = accessions
                    organism_by_protein_index[protein_index] = str(protein_entry.get("organism", ""))[:60]
                elif is_human:
                    # no accession on this entry; the EC number still identifies the enzyme, so the row is written
                    # with an empty accession and joined through the local UniProt EC column instead
                    human_protein_indices_without_an_accession.add(protein_index)
                    accessions_by_protein_index[protein_index] = [""]
                    organism_by_protein_index[protein_index] = str(protein_entry.get("organism", ""))[:60]
            if not accessions_by_protein_index:
                continue
            for field_name, (relation_type, sign) in REGULATION_FIELDS_AND_RELATIONS.items():
                for record in enzyme_record.get(field_name) or []:
                    ligand_name = str(record.get("value", "")).strip()
                    if not ligand_name or ligand_name.lower() in ("more", "-"):
                        continue
                    for protein_index in record.get("proteins") or []:
                        for accession in accessions_by_protein_index.get(str(protein_index), []):
                            writer.writerow({"ec_number": ec_number, "field": field_name, "relation_type": relation_type,
                                             "sign": sign, "ligand_name": ligand_name, "uniprot_accession": accession,
                                             "organism_entry": organism_by_protein_index.get(str(protein_index), ""),
                                             "brenda_comment": str(record.get("comment", ""))[:200],
                                             "brenda_reference_count": len(record.get("references") or [])})
                            rows_written += 1
                            counts_by_relation[relation_type] += 1
                            distinct_ligand_names.add(ligand_name)
                            distinct_human_accessions.add(accession)

    summary = {"enzyme_entries_seen": enzyme_entries_seen, "rows_written": rows_written,
               "counts_by_relation": dict(counts_by_relation),
               "distinct_ligand_names": len(distinct_ligand_names),
               "distinct_human_accessions": len(distinct_human_accessions)}
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
