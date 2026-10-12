# Benet, Broccatelli and Oprea 2011: BDDCS class and percent excreted unchanged in urine for 927 drugs

> Copied on 11 October 2026 from `data/raw/exposure/benet2011_bddcs/sources.md`, the note written beside the raw files (which are not in git). Paths such as `raw/...` are relative to that directory. The transcription scripts are in `experiments/exposure_sources/`.

Retrieved on 12 October 2026 (UTC; the evening of 11 October 2026 local time). Transcription only: nothing here is
matched to the study's drugs.

## Citation

Benet LZ, Broccatelli F, Oprea TI. BDDCS applied to over 900 drugs. AAPS J 2011;13(4):519-547.
doi:10.1208/s12248-011-9290-9. PMID 21818695. PMCID PMC3231854.

The tables transcribed are Tables I to V of the article, one per BDDCS class: "Measured and In Silico Data for 351
BDDCS Class 1 Drugs", "... 265 BDDCS Class 2 Drugs", "... 247 BDDCS Class 3 Drugs", "... 53 BDDCS Class 4 Drugs",
"... 11 BDDCS Class 0 Drugs".

## Route

The article page in PubMed Central, read as HTML with `experiments/read_public_page.py` (the PMC robots file allows
`/articles/`). The page holds the five tables in full. The programmatic routes and the supplementary spreadsheet were
closed (see "Refused").

## Addresses read

| Address | How | Result |
| --- | --- | --- |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:"10.1208/s12248-011-9290-9"&resultType=core&format=json` | curl | record: PMC3231854, `isOpenAccess` N, no licence field, `hasSuppl` Y (`raw/metadata/epmc_10.1208_s12248-011-9290-9.json`) |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3231854/fullTextXML` | curl, twice | HTTP 500 both times: `{"status":500,"error":"Internal Server Error", ...}` |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3231854/supplementaryFiles` | curl | HTTP 200 with an error document: "Article with id PMC3231854 is not open access one" |
| `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=3231854&retmode=xml` | curl | front matter only, with the comment "The publisher of this article does not allow downloading of the full text in XML form." (`raw/metadata/efetch_PMC3231854.xml`) |
| PubMed connector, full text of PMC3231854 | connector | abstract only, empty full text |
| `https://api.crossref.org/works/10.1208/s12248-011-9290-9` | curl | licence metadata (`raw/metadata/crossref_benet2011.json`) |
| `https://pmc.ncbi.nlm.nih.gov/articles/PMC3231854/` | `experiments/read_public_page.py` | **the source page**, HTTP 200, 994369 bytes (`raw/2f21927cccf25526.html`) |
| `https://pmc.ncbi.nlm.nih.gov/about/copyright/` | `experiments/read_public_page.py` | PMC copyright notice (`raw/metadata/d630fe4ee6d42383.html`) |
| `https://pmc.ncbi.nlm.nih.gov/articles/instance/3231854/bin/12248_2011_9290_MOESM2_ESM.xlsx` and `..._MOESM1_ESM.doc` | `experiments/read_public_page.py` | **refused**: a page "Preparing to download ..." that sets a proof-of-work challenge for the browser; not solved, files not saved |
| `https://static-content.springer.com/esm/art%3A10.1208%2Fs12248-011-9290-9/MediaObjects/12248_2011_9290_MOESM2_ESM.xlsx` and `..._MOESM1_ESM.doc` | `experiments/read_public_page.py` | **not read**: "the site's robots file disallows this address for automated readers" |
| `https://escholarship.org/uc/item/74q7552x` (a search result of unknown content) | `experiments/read_public_page.py` | **not read**: "the site did not serve a robots file to an automated reader, which is taken as a refusal" |
| GitHub code search for a redistributed BDDCS table | GitHub connector | nothing found |

No address waited on a permission request.

## Raw files and SHA-256

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `raw/2f21927cccf25526.html` (PMC article page) | 994369 | `2f21927cccf2552647cba9779fe8c5f56ae5449ccb822a891d6e62db4be65c1e` |
| `raw/2f21927cccf25526.html.json` (address, time, SHA-256) | 465 | `e6ba7afdb243a15c5c1f74cc1645d4452031cfda14385df1cd067b013d0e8006` |
| `raw/metadata/efetch_PMC3231854.xml` | 8288 | `2fb41ec7f022ee9131ef011f3df6df569e4907016a178734d221f9ab8a0f9f40` |
| `raw/metadata/epmc_10.1208_s12248-011-9290-9.json` | 6461 | `531c45573a938061bdb380f50d92431e31ba9cdd4aa333ec7ff4da4091fd5e19` |
| `raw/metadata/crossref_benet2011.json` | 21891 | `58abcf416936a69e987cca7cccf553ff172e4c59fb74e74ff62d444d3d035553` |
| `raw/metadata/d630fe4ee6d42383.html` (PMC copyright notice) | 85464 | `d630fe4ee6d42383d868faede01ffa0e9e37d1630dc3e19d1727a235976c2f4e` |

`benet2011_bddcs.csv`: 94804 bytes, SHA-256 `2954651b18411e46531f1549bd9b337745a53adf4998412b95340858a82e83ab`.

## Licence

**No open licence. The article is free to read in PubMed Central under ordinary copyright.**

- Copyright statement of the article, read in the NCBI record (`raw/metadata/efetch_PMC3231854.xml`,
  `<copyright-statement>`) and on the PMC page: "© American Association of Pharmaceutical Scientists 2011".
- Europe PMC record: `isOpenAccess` "N", no licence. The NCBI record: "The publisher of this article does not allow
  downloading of the full text in XML form."
- PMC copyright notice (`https://pmc.ncbi.nlm.nih.gov/about/copyright/`): "In the absence of a copyright statement,
  users should assume that standard copyright protection applies, unless the article contains an explicit statement
  to the contrary." and "Users of PMC are directly and solely responsible for compliance with copyright restrictions
  and are expected to adhere to the terms and conditions defined by the copyright holder." and "Systematic
  downloading of batches of articles from the main PMC web site, in any way, is prohibited because of copyright
  restrictions." (Two article pages were read for this work, one for this source and one for the 2016 update.)
- Crossref licence field: `http://www.springer.com/tdm` (the publisher's text-and-data-mining terms; that page was
  not read).
- The article states of its supplement: "The online version of this article (doi:10.1208/s12248-011-9290-9) contains
  supplementary material, which is available to authorized users."

The CSV therefore holds facts transcribed from a copyrighted article; it should not be redistributed as a file
without the owner's decision.

## How the CSV was made

Script: `experiments/exposure_sources/make_benet2011_csv.py` (run with `python3 -I`, BeautifulSoup 4.15.0 with lxml).
The page has exactly five `<table>` elements, each with one header row, 17 columns and no merged cells.

- `bddcs_class` is taken from the caption of the table in which the drug is listed (Table I = 1, II = 2, III = 3,
  IV = 4, V = 0); `source_table` and `source_row` (position within the table) say where the row was read.
- Source column to CSV column: `Generic name` → `name`; `Maximum strength dose value` →
  `maximum_strength_dose_value`; `Maximum strength dose unit` → `maximum_strength_dose_unit`; `Formulation` →
  `formulation`; `Route` → `route`; `Measured solubility (mg/mL)` → `measured_solubility_mg_per_ml`; `Dose number`
  → `dose_number`; `% Excreted unchanged in urine` → `percent_excreted_unchanged_urine`; `MW drug` → `mw_drug`;
  `Measured LogS molar` → `measured_log_s_molar`; `Measured LogP` → `measured_log_p`; `Measured LogD74` →
  `measured_log_d_7_4`; `minVSLgS 3–7.5` → `min_vs_lgs_3_7_5`; `CLogP` → `clogp`; `HBA` → `hba`; `HBD` → `hbd`;
  `PSA` → `psa`.
- Cell text is kept as printed (for example `1.20` stays `1.20`). Runs of whitespace, which are line wraps of the
  page source, are collapsed to one space. In numeric columns the typographic minus sign (U+2212) is written as the
  ASCII hyphen-minus. Nothing is computed.
- The tables have no column for excretion in bile and none for the extent of metabolism as a number; the class is
  the only statement about metabolism.
- No cell carries a footnote marker, and the tables have no footnotes. No numeric cell carries a qualifier.
- Two rows have a `note`: `Penicillin G; benzylpenicillin`, whose dose cells read `39 (60,000` and
  `mg/mL units/mL)` and are kept as printed; and `SN-38; 7-ethyl-10- hydroxycamptothecin`, where the page source
  breaks the line after a hyphen, so the space after `10-` comes from the line break.
- 30 rows are active metabolites: their `maximum_strength_dose_unit` reads `Active metabolite` and their dose,
  formulation and route cells are empty, as printed.

## Rows and non-missing values

927 rows: class 1, 351; class 2, 265; class 3, 247; class 4, 53; class 0, 11.

| Column | Non-missing | Distinct |
| --- | --- | --- |
| `name` | 927 | 927 |
| `bddcs_class` | 927 | 5 |
| `maximum_strength_dose_value` | 897 | 110 |
| `maximum_strength_dose_unit` | 927 | 12 |
| `formulation` | 897 | 31 |
| `route` | 897 | 12 |
| `measured_solubility_mg_per_ml` | 635 | 338 |
| `dose_number` | 536 | 157 |
| `percent_excreted_unchanged_urine` | 834 | 147 |
| `mw_drug` | 927 | 839 |
| `measured_log_s_molar` | 634 | 399 |
| `measured_log_p` | 430 | 322 |
| `measured_log_d_7_4` | 478 | 341 |
| `min_vs_lgs_3_7_5` | 922 | 573 |
| `clogp` | 924 | 571 |
| `hba` | 927 | 28 |
| `hbd` | 927 | 26 |
| `psa` | 927 | 683 |
| `source_table` | 927 | 5 |
| `source_row` | 927 | 351 |
| `note` | 2 | 2 |

`percent_excreted_unchanged_urine` by class: class 1, 302 of 351; class 2, 228 of 265; class 3, 243 of 247;
class 4, 50 of 53; class 0, 11 of 11. Values run from 0 to 100.

## Checks

- Each table has as many rows as its caption says (351, 265, 247, 53, 11); the script stops if it does not.
- The CSV reproduces the counts given in the abstract: "927 drugs, which include 30 active metabolites" (927 rows,
  30 with `Active metabolite`); "Of the 897 parent drugs, 78.8% (707) are administered orally" (897 rows with a
  route, 707 with route `Oral`); "this value is reported for 72.7% (513) of these orally administered drugs" (513
  oral rows with a measured solubility); class 4 "(N = 53 compared with over 200 each in classes 1–3)" (53 rows).
- Twelve spot checks, cell text of the page beside the CSV (`experiments/exposure_sources/spot_checks.py`):

| Name | Table and row | Class | % excreted unchanged in urine (page / CSV) | Dose (page) | CLogP (page / CSV) |
| --- | --- | --- | --- | --- | --- |
| Amiodarone hydrochloride | Table II row 10 | 2 | `0` / `0` | `400` `mg` | `8.95` / `8.95` |
| Diazepam | Table I row 106 | 1 | `0.5` / `0.5` | `10` `mg` | `2.96` / `2.96` |
| Warfarin | Table II row 261 | 2 | `1` / `1` | `10` `mg` | `2.90` / `2.90` |
| Haloperidol | Table II row 108 | 2 | `1` / `1` | `20` `mg` | `3.85` / `3.85` |
| Fluoxetine | Table I row 146 | 1 | `1.25` / `1.25` | `20` `mg` | `4.57` / `4.57` |
| Lithium carbonate | Table III row 140 | 3 | `95` / `95` | `600` `mg` | `−3.65` / `-3.65` |
| Metformin | Table III row 145 | 3 | `99` / `99` | `1000` `mg` | `−1.63` / `-1.63` |
| Gabapentin | Table III row 110 | 3 | `100` / `100` | `800` `mg` | `−0.66` / `-0.66` |
| Acetazolamide | Table IV row 1 | 4 | `90` / `90` | `250` `mg` | `−0.98` / `-0.98` |
| Amphetamine sulfate | Table V row 1 | 0 | `40` / `40` | `10` `mg` | `1.74` / `1.74` |
| Atenolol | Table III row 15 | 3 | `94` / `94` | `100` `mg` | `−0.11` / `-0.11` |
| Penicillin G; benzylpenicillin | Table III row 180 | 3 | empty / empty | `39 (60,000` `mg/mL units/mL)` | `1.75` / `1.75` |

## Refused, and other things to know

- The supplementary files (ESM 1, a 5.6 MB `.doc`; ESM 2, a 274.9 KB `.xlsx`, as listed on the PMC page) were not
  obtained: PubMed Central answers the download address with a proof-of-work page for a browser, the publisher's
  file host disallows the address in its robots file, and the Europe PMC service answers "Article with id
  PMC3231854 is not open access one". None of these was worked around. Addresses for the owner to open:
  `https://pmc.ncbi.nlm.nih.gov/articles/PMC3231854/` (section "Supplementary Materials") or
  `https://doi.org/10.1208/s12248-011-9290-9`. What the spreadsheet adds to the five printed tables is not known.
- `escholarship.org` did not serve a robots file to the reader, which the reader takes as a refusal.
- The 2016 update (`../hosey2016_bddcs/`) changes the class of 13 of these drugs and the percent excreted unchanged
  of four; those changes are **not** applied in this file.
- No text addressed to an AI assistant was seen in the page.
