# Data from the Comparative Toxicogenomics Database (CTD) in this directory

The table `ctd_relations_filtered` is derived from CTD's file `CTD_chemicals_diseases.tsv.gz` (https://ctdbase.org/), report created Tue Sep 29 13:17:30 EDT 2026, SHA-256 `e11e6dc36a27d36ed88576e529a92211186e3e0b3cea4271cc78dcf3d1b821ee`. It holds 5,191 rows: the curated chemical-disease statements (the DirectEvidence column) for the study's drugs and symptoms, one row per cited paper.
Rows that came through the literature reports into other tables of this repository carry the source name `CTD-curated`.

Curated chemical-disease data were retrieved from the Comparative Toxicogenomics Database (CTD), https://ctdbase.org/, on 2026-10-06. The file's header reads: "Copyright 2002-2012 MDI Biological Laboratory. All rights reserved. Copyright 2012-2026 NC State University. All rights reserved."

## Citation

Davis AP, Wiegers TC, Sciaky D, Barkalow F, Strong M, Wyatt B, Wiegers J, McMorran R, Abrar S, Mattingly CJ. Comparative Toxicogenomics Database's 20th anniversary: update 2025. Nucleic Acids Res. 2025;53(D1):D1328-D1334. https://doi.org/10.1093/nar/gkae883

## CTD's terms, as the data file prints them

"Use is subject to the terms set forth at http://ctdbase.org/about/legal.jsp. These terms include:"

1. "All forms of publication (e.g., web sites, research papers, databases, software applications, etc.) that use or rely on CTD data must cite CTD."
2. "All electronic or online applications must include hyperlinks from contexts that use CTD data to the applicable CTD data pages."
3. "You must notify CTD, and describe your use of our data"
4. "For quality control purposes, you must provide CTD with periodic access to your publication of our data."

Davis et al. 2025, under data availability: "All aspects of CTD are freely available to non-commercial users; commercial users are asked to purchase a license in order to download CTD data files to their systems."

So this table is for non-commercial use. A commercial user needs CTD's licence before using it.

## How this repository meets them

1. **Citation.** Above, in `README.md` and as reference `davis2025ctd` in `docs/references.bib`. Any publication that uses this table has to cite CTD as well.
2. **Links.** Every row names its CTD chemical in `chemical_mesh_id` and its CTD disease in `mesh_descriptor`. Their CTD pages are
   `http://ctdbase.org/detail.go?type=chem&acc=<chemical_mesh_id>` and `http://ctdbase.org/detail.go?type=disease&db=MESH&acc=<mesh_descriptor>`,
   for example http://ctdbase.org/detail.go?type=chem&acc=D000077443 and http://ctdbase.org/detail.go?type=disease&db=MESH&acc=D001008.
3. **Notice to CTD.** Sent by the repository's owner through CTD's contact form; `docs/data_sources.md` records the date, or that it has not been sent yet.
4. **Access.** The repository is public, and the notice gives CTD its address.

## What was not read

CTD's legal notice (http://ctdbase.org/about/legal.jsp), citation guidelines (http://ctdbase.org/about/publications/#citing) and linking instructions (http://ctdbase.org/help/linking.jsp) decline automated readers.
The terms above are the ones the data file and the paper print; the header says the legal notice's terms "include" them, so the list may not be complete.
The form of the page addresses is the one the Bioregistry records for CTD (https://bioregistry.io/ctd.chemical, https://bioregistry.io/ctd.disease).
