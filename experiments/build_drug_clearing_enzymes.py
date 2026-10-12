"""Which enzyme clears which study drug, from the FDA Table of Pharmacogenetic Associations.

Why. The user, 11 October 2026: "Add the enzyme labels as the carrier twin instead then". A plasma carrier and a clearing
enzyme stand in the same relation to a drug: more of either leaves less free drug. The carrier's statement holds for few
drugs (docs/drug_exposure_table.md) and has four to six sources (docs/label_source_review.md, section 2). The
enzyme's holds whatever the route, and a regulator has reviewed each row. This script reads the regulator's table and
writes one row per study drug and enzyme, with the table's own sentence, for two uses:

  - the graph builder adds a clearance edge, relation clears and sign -1, from the enzyme to the drug's node
    (mechanistic_pathway_learning/graph/drug_entry_nodes.py);
  - experiments/build_comparative_checks.py turns each row into a comparative statement: the drug with its enzyme
    lowered should score at least as high as the drug alone.

Which rows give an edge. The sign -1 says "less enzyme, more drug". A row supports that when it is about a metabolising
enzyme and says that poor (or intermediate or poor) metabolisers have higher concentrations of the drug or of its active
metabolite, or a higher risk of adverse reactions. Rows left out, each with its reason in the output:

  - an immune gene (HLA-A, HLA-B): the row is about hypersensitivity and not about clearance;
  - a row that states no direction ("May alter systemic concentrations.");
  - a row in which the enzyme forms the active compound, so less enzyme gives less effect (tramadol with CYP2D6:
    "Results in lower systemic active metabolite concentrations and may result in reduced efficacy.");
  - a row about another ingredient of a combination (dextromethorphan under bupropion).

Source. https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations, the page saved on
10 October 2026 under data/raw/label_review/fda_pharmacogenetic_associations. FDA's website policies (read 11 October
2026; content current as of 01/19/2024): "Unless otherwise noted, the contents of the FDA website (www.fda.gov) — both
text and graphics — are not copyrighted. They are in the public domain and may be republished, reprinted and otherwise
used freely by anyone without the need to obtain permission from FDA. Credit to the U.S. Food and Drug Administration as
the source is appreciated but not required."

Usage:
  PYTHONPATH=. python experiments/build_drug_clearing_enzymes.py --evidence-dir data/processed/evidence_full_v4
"""
from __future__ import annotations

import argparse
import re
from html.parser import HTMLParser
from pathlib import Path

import pandas as pd

FDA_TABLE = Path("data/raw/label_review/fda_pharmacogenetic_associations/table_pharmacogenetic_associations.html")
OUTPUT = Path("configs/drug_clearing_enzymes.csv")
FDA_TABLE_ADDRESS = "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations"
# the genes of the table that code for an enzyme which metabolises drugs; a gene outside this set gives no clearance edge
METABOLISING_ENZYMES = frozenset({"CYP2D6", "CYP2C19", "CYP2C9", "CYP2B6", "CYP3A5", "NAT2", "UGT1A1", "DPYD", "TPMT", "NUDT15", "BCHE"})
LOWERED_ENZYME_SUBGROUP = re.compile(r"poor metaboli[sz]ers", re.IGNORECASE)
HIGHER_WITH_LESS_ENZYME = re.compile(r"(?:results? in|have) higher systemic(?: active metabolite)? concentrations|higher adverse reaction risk", re.IGNORECASE)
ENZYME_FORMS_THE_ACTIVE_COMPOUND = re.compile(r"lower systemic active metabolite", re.IGNORECASE)
NAMES_ANOTHER_COMPOUND = re.compile(r"concentrations of ([a-z][a-z-]+)", re.IGNORECASE)
SECTION_NAMES = {1: "data support therapeutic management recommendations", 2: "data indicate a potential impact on safety or response",
                 3: "data demonstrate a potential impact on pharmacokinetic properties only"}


class TableRows(HTMLParser):
    """The rows of every table of a page, as lists of cell texts."""

    def __init__(self) -> None:
        super().__init__()
        self.tables: list[list[list[str]]] = []
        self.row: list[str] | None = None
        self.cell: list[str] = []
        self.in_cell = False

    def handle_starttag(self, tag, attrs) -> None:
        if tag == "table":
            self.tables.append([])
        elif tag == "tr":
            self.row = []
        elif tag in ("td", "th"):
            self.in_cell, self.cell = True, []

    def handle_endtag(self, tag) -> None:
        if tag in ("td", "th") and self.in_cell and self.row is not None:
            self.row.append(" ".join("".join(self.cell).split()))
            self.in_cell = False
        elif tag == "tr" and self.row is not None and self.tables:
            self.tables[-1].append(self.row)
            self.row = None

    def handle_data(self, data) -> None:
        if self.in_cell:
            self.cell.append(data)


def table_rows(page_path: Path = FDA_TABLE) -> pd.DataFrame:
    """One row per association: section (1 to 3), drug, gene, subgroup and description, as the page prints them."""
    parser = TableRows()
    parser.feed(page_path.read_text(errors="replace"))
    rows = [{"section": section, "drug": row[0], "gene": row[1], "subgroup": row[2], "description": row[3]}
            for section, table in enumerate(parser.tables, start=1) for row in table[1:] if len(row) >= 4]
    return pd.DataFrame(rows)


def drugs_named(cell: str) -> list[str]:
    """The drug names of a table cell, lower case: "Bupropion and Dextromethorphan" names two."""
    return [part.strip().lower() for part in re.split(r" and |/|,", cell) if part.strip()]


def edge_decision(study_drug: str, gene: str, subgroup: str, description: str) -> tuple[bool, str]:
    """(gives a clearance edge, the reason) for one row that names a study drug."""
    genes = [part for part in re.split(r"[ ,/]+|and", gene) if part]
    enzymes = [part for part in genes if part in METABOLISING_ENZYMES]
    if not enzymes or len(genes) != 1:
        return False, "not a single metabolising enzyme" + (" (an immune gene: the row is about hypersensitivity)" if gene.startswith("HLA") else "")
    if ENZYME_FORMS_THE_ACTIVE_COMPOUND.search(description) or re.search(r"ultrarapid metaboli[sz]ers\d*$", subgroup.strip(), re.IGNORECASE):
        return False, "the enzyme forms the active compound, or the row is about more enzyme: the sign of a clearance edge would be wrong"
    other_compound = NAMES_ANOTHER_COMPOUND.search(description)
    if other_compound and other_compound.group(1).lower() not in (study_drug, "the", "drug", "parent"):
        return False, f"the row is about {other_compound.group(1).lower()}, another ingredient"
    if not LOWERED_ENZYME_SUBGROUP.search(subgroup):
        return False, "the subgroup is not poor metabolisers"
    if not HIGHER_WITH_LESS_ENZYME.search(description):
        return False, "no direction stated"
    return True, "poor metabolisers have higher concentrations or a higher risk of adverse reactions"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"), help="names the study drugs")
    parser.add_argument("--fda-table", type=Path, default=FDA_TABLE)
    parser.add_argument("--page-saved-on", default="2026-10-10", help="the date the page was saved, recorded with every row")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet", columns=["perturbation_id", "perturbation_type", "perturbation_label"])
    drugs = records[records.perturbation_type == "drug"].drop_duplicates("perturbation_id")
    perturbation_id_of = {label.lower(): perturbation_id for label, perturbation_id in zip(drugs.perturbation_label, drugs.perturbation_id)}
    rows = []
    for row in table_rows(arguments.fda_table).itertuples(index=False):
        for name in drugs_named(row.drug):
            if name not in perturbation_id_of:
                continue
            gives_an_edge, reason = edge_decision(name, row.gene, row.subgroup, row.description)
            rows.append({"perturbation_id": perturbation_id_of[name], "drug": name, "enzyme_gene": row.gene if gives_an_edge else "", "gives_an_edge": gives_an_edge, "reason": reason,
                         "fda_table_section": row.section, "fda_table_section_title": SECTION_NAMES.get(row.section, ""), "fda_table_drug": row.drug, "fda_table_gene": row.gene,
                         "fda_table_subgroup": row.subgroup, "fda_table_description": row.description, "source": FDA_TABLE_ADDRESS, "page_saved_on": arguments.page_saved_on})
    table = pd.DataFrame(rows).sort_values(["gives_an_edge", "drug", "fda_table_gene"], ascending=[False, True, True])
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(arguments.output, index=False)
    with_edge = table[table.gives_an_edge]
    print(f"rows naming a study drug: {len(table)}; giving a clearance edge: {len(with_edge)} over {with_edge.drug.nunique()} drugs, by enzyme {with_edge.enzyme_gene.value_counts().to_dict()}")
    print("left out:", table[~table.gives_an_edge].reason.value_counts().to_dict())
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
