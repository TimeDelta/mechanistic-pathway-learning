"""The notice that travels with every table derived from the Comparative Toxicogenomics Database (CTD).

CTD's data files open with the terms of their use. Four are listed there, and this module is how the repository meets
the two that a file can meet: the citation and the links. The other two, the notice to CTD and CTD's access to the
publication of its data, are acts of the repository's owner; the notice text says so and docs/data_sources.md records
whether they were done.

The terms are quoted from the header of CTD_chemicals_diseases.tsv.gz (report created 29 September 2026), read on
11 October 2026, and from the data availability statement of Davis et al. 2025. CTD's own pages for the legal notice,
the citation guidelines and the linking instructions decline automated readers and were not read
(docs/label_source_review.md, section 11.2), so the wording of the citation follows the paper and the form of the
links follows the Bioregistry's records for CTD (bioregistry.io/ctd.chemical and bioregistry.io/ctd.disease, read on
11 October 2026).
"""
from __future__ import annotations

from pathlib import Path

CTD_TERMS_FILE_NAME = "CTD_TERMS.md"
CTD_HOME = "https://ctdbase.org/"
CTD_LEGAL_NOTICE = "http://ctdbase.org/about/legal.jsp"
CTD_CITATION_GUIDELINES = "http://ctdbase.org/about/publications/#citing"
CTD_LINKING_INSTRUCTIONS = "http://ctdbase.org/help/linking.jsp"
CTD_CONTACT = "http://ctdbase.org/help/contact.go"
CTD_CITATION = ("Davis AP, Wiegers TC, Sciaky D, Barkalow F, Strong M, Wyatt B, Wiegers J, McMorran R, Abrar S, Mattingly CJ. "
                "Comparative Toxicogenomics Database's 20th anniversary: update 2025. Nucleic Acids Res. 2025;53(D1):D1328-D1334. "
                "https://doi.org/10.1093/nar/gkae883")
# the four terms, as the header of the data file prints them
CTD_TERMS_AS_PRINTED = (
    "All forms of publication (e.g., web sites, research papers, databases, software applications, etc.) that use or rely on CTD data must cite CTD.",
    "All electronic or online applications must include hyperlinks from contexts that use CTD data to the applicable CTD data pages.",
    "You must notify CTD, and describe your use of our data",
    "For quality control purposes, you must provide CTD with periodic access to your publication of our data.",
)
CTD_COMMERCIAL_USE_AS_PRINTED = ("All aspects of CTD are freely available to non-commercial users; commercial users are asked to purchase a license "
                                 "in order to download CTD data files to their systems.")


def ctd_chemical_page(chemical_mesh_id: str) -> str:
    """CTD's page for a chemical, from its MeSH identifier (the form the Bioregistry records for ctd.chemical)."""
    return f"http://ctdbase.org/detail.go?type=chem&acc={chemical_mesh_id.removeprefix('MESH:')}"


def ctd_disease_page(disease_mesh_id: str) -> str:
    """CTD's page for a disease, from its MeSH identifier (the form the Bioregistry records for ctd.disease)."""
    return f"http://ctdbase.org/detail.go?type=disease&db=MESH&acc={disease_mesh_id.removeprefix('MESH:')}"


def ctd_terms_notice(report_created: str | None, file_sha256: str | None, rows_derived: int | None, retrieved: str | None = None,
                     derived_table: str = "ctd_relations_filtered") -> str:
    """The text of CTD_TERMS.md for a directory that holds a table derived from CTD."""
    lines = [
        "# Data from the Comparative Toxicogenomics Database (CTD) in this directory",
        "",
        f"The table `{derived_table}` is derived from CTD's file `CTD_chemicals_diseases.tsv.gz` ({CTD_HOME}),"
        + (f" report created {report_created}" if report_created else "")
        + (f", SHA-256 `{file_sha256}`" if file_sha256 else "") + "."
        + (f" It holds {rows_derived:,} rows: the curated chemical-disease statements (the DirectEvidence column) for the study's drugs and symptoms, one row per cited paper." if rows_derived is not None else ""),
        "Rows that came through the literature reports into other tables of this repository carry the source name `CTD-curated`.",
        "",
        f"Curated chemical-disease data were retrieved from the Comparative Toxicogenomics Database (CTD), {CTD_HOME}"
        + (f", on {retrieved}" if retrieved else "") + ". The file's header reads: \"Copyright 2002-2012 MDI Biological Laboratory. All rights reserved. Copyright 2012-2026 NC State University. All rights reserved.\"",
        "",
        "## Citation",
        "",
        CTD_CITATION,
        "",
        "## CTD's terms, as the data file prints them",
        "",
        f"\"Use is subject to the terms set forth at {CTD_LEGAL_NOTICE}. These terms include:\"",
        "",
    ]
    lines += [f"{number}. \"{term}\"" for number, term in enumerate(CTD_TERMS_AS_PRINTED, start=1)]
    lines += [
        "",
        f"Davis et al. 2025, under data availability: \"{CTD_COMMERCIAL_USE_AS_PRINTED}\"",
        "",
        "So this table is for non-commercial use. A commercial user needs CTD's licence before using it.",
        "",
        "## How this repository meets them",
        "",
        "1. **Citation.** Above, in `README.md` and as reference `davis2025ctd` in `docs/references.bib`. Any publication that uses this table has to cite CTD as well.",
        "2. **Links.** Every row names its CTD chemical in `chemical_mesh_id` and its CTD disease in `mesh_descriptor`. Their CTD pages are",
        f"   `{ctd_chemical_page('<chemical_mesh_id>')}` and `{ctd_disease_page('<mesh_descriptor>')}`,",
        f"   for example {ctd_chemical_page('D000077443')} and {ctd_disease_page('D001008')}.",
        "3. **Notice to CTD.** Sent by the repository's owner through CTD's contact form; `docs/data_sources.md` records the date, or that it has not been sent yet.",
        "4. **Access.** The repository is public, and the notice gives CTD its address.",
        "",
        "## What was not read",
        "",
        f"CTD's legal notice ({CTD_LEGAL_NOTICE}), citation guidelines ({CTD_CITATION_GUIDELINES}) and linking instructions ({CTD_LINKING_INSTRUCTIONS}) decline automated readers.",
        "The terms above are the ones the data file and the paper print; the header says the legal notice's terms \"include\" them, so the list may not be complete.",
        "The form of the page addresses is the one the Bioregistry records for CTD (https://bioregistry.io/ctd.chemical, https://bioregistry.io/ctd.disease).",
        "",
    ]
    return "\n".join(lines)


def write_ctd_terms(directory: Path, report_created: str | None, file_sha256: str | None, rows_derived: int | None, retrieved: str | None = None) -> Path:
    """Write CTD_TERMS.md into a directory that holds a table derived from CTD; returns its path."""
    path = directory / CTD_TERMS_FILE_NAME
    path.write_text(ctd_terms_notice(report_created, file_sha256, rows_derived, retrieved), encoding="utf-8")
    return path
