"""Rows of the FDA Table of Pharmacogenetic Associations that name a study drug: by gene, by what the row says happens,
and whether the gene has a node in the confirmatory graph.

Part of the label source review of 10 October 2026 (docs/label_source_review.md). Counting only: no label, table or
configuration of the study is changed.

Usage:
  python experiments/label_review/fda_pgx_overlap.py data/raw/label_review/fda_pharmacogenetic_associations/table_pharmacogenetic_associations.html data/processed/label_review/fda_pharmacogenetic_rows.csv
"""
import re, sys
from html.parser import HTMLParser
import pandas as pd

class TableRows(HTMLParser):
    def __init__(self):
        super().__init__(); self.tables, self.row, self.cell, self.in_cell = [], None, [], False
    def handle_starttag(self, tag, attrs):
        if tag == "table": self.tables.append([])
        elif tag == "tr": self.row = []
        elif tag in ("td", "th"): self.in_cell, self.cell = True, []
    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.in_cell:
            self.row.append(" ".join("".join(self.cell).split())); self.in_cell = False
        elif tag == "tr" and self.row is not None and self.tables:
            self.tables[-1].append(self.row); self.row = None
    def handle_data(self, data):
        if self.in_cell: self.cell.append(data)

parser = TableRows()
parser.feed(open(sys.argv[1], errors="replace").read())
rows = []
for section, table in enumerate(parser.tables, start=1):
    for row in table[1:]:
        if len(row) >= 4:
            rows.append({"section": section, "drug": row[0], "gene": row[1], "subgroup": row[2], "description": row[3]})
table = pd.DataFrame(rows)
print("rows:", len(table), "by section:", table.section.value_counts().sort_index().to_dict())
graph_nodes = set(pd.read_parquet("data/processed/graph_full_neuronal_split_binders/nodes.parquet", columns=["node_id"]).node_id)
def study_names(evidence_directory):
    records = pd.read_parquet(f"{evidence_directory}/evidence_records.parquet", columns=["perturbation_type", "perturbation_label"])
    return {name.lower() for name in records[records.perturbation_type == "drug"].perturbation_label.unique()}
name_sets = {"current rule (154 drugs)": study_names("data/processed/evidence_full_v3_parkinsonism"),
             "up to three targets (222)": study_names("data/processed/label_review/evidence_n_only_cap3"),
             "any ATC, up to three targets (924)": study_names("data/processed/label_review/evidence_any_atc_cap3")}
def first_drug_names(cell):  # "Amitriptyline", "Codeine", "Valproic Acid"; combination rows name two
    return [part.strip().lower() for part in re.split(r" and |/|,", cell) if part.strip()]
table["drug_names"] = table.drug.map(first_drug_names)
table["gene_in_graph"] = table.gene.map(lambda gene: any(f"GENE:{part.strip()}" in graph_nodes for part in re.split(r"[ ,/]|and", gene) if part.strip()))
table["says_adverse"] = table.description.str.contains("adverse|toxicit|side effect|QT|sedation|respiratory depression|seizure", case=False)
table["says_higher_concentration"] = table.description.str.contains("higher systemic|higher plasma|higher concentration|increased exposure|higher exposure", case=False)
for label, names in name_sets.items():
    matched = table[table.drug_names.map(lambda parts: any(part in names for part in parts))]
    print(f"\n{label}: rows {len(matched)}, drugs {matched.drug.nunique()}, by gene {matched.gene.value_counts().to_dict()}")
    print(f"   gene has a graph node: {int(matched.gene_in_graph.sum())} | says adverse reaction or toxicity: {int(matched.says_adverse.sum())} | says higher concentration: {int(matched.says_higher_concentration.sum())}")
    if label.startswith("up to three"):
        print(matched[["drug", "gene", "subgroup", "description"]].to_string(index=False, max_colwidth=110))
table.to_csv(sys.argv[2], index=False)
