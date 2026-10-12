# Hosey, Chan and Benet 2016: BDDCS class for 180 further drugs, and corrections to the 2011 table

> Copied on 11 October 2026 from `data/raw/exposure/hosey2016_bddcs/sources.md`, the note written beside the raw files (which are not in git). Paths such as `raw/...` are relative to that directory. The transcription scripts are in `experiments/exposure_sources/`.

Retrieved on 12 October 2026 (UTC; the evening of 11 October 2026 local time). Transcription only: nothing here is
matched to the study's drugs, and the corrections are not applied to the 2011 table.

## Citation

Hosey CM, Chan R, Benet LZ. BDDCS predictions, self-correcting aspects of BDDCS assignments, BDDCS assignment
corrections, and classification for more than 175 additional drugs. AAPS J 2016;18(1):251-260.
doi:10.1208/s12248-015-9845-2. PMID 26589308. PMCID PMC4706282. (Published online 20 November 2015.)

## Route

The article page in PubMed Central, read as HTML with `experiments/read_public_page.py`. The page holds three
tables: Table I "BDDCS Class Changes from Initial Publication", Table II "Parameter Changes from Initial
Publication", Table III "Newly BDDCS Classified Drugs". The page lists no supplementary file.

## Addresses read

| Address | How | Result |
| --- | --- | --- |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:"10.1208/s12248-015-9845-2"&resultType=core&format=json` | curl | record: PMC4706282, `isOpenAccess` N, no licence field, `hasSuppl` N (`raw/metadata/epmc_10.1208_s12248-015-9845-2.json`) |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4706282/fullTextXML` | curl | HTTP 500 |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4706282/supplementaryFiles` | curl | HTTP 200 with an error document: "Article with id PMC4706282 is not open access one" |
| `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=4706282&retmode=xml` | curl | front matter only, with the comment "The publisher of this article does not allow downloading of the full text in XML form." (`raw/metadata/efetch_PMC4706282.xml`) |
| `https://api.crossref.org/works/10.1208/s12248-015-9845-2` | curl | licence metadata (`raw/metadata/crossref_hosey2016.json`) |
| `https://pmc.ncbi.nlm.nih.gov/articles/PMC4706282/` | `experiments/read_public_page.py` | **the source page**, HTTP 200, 198829 bytes (`raw/6e9eae3c00ba1e65.html`) |
| `https://pmc.ncbi.nlm.nih.gov/about/copyright/` | `experiments/read_public_page.py` | PMC copyright notice (`raw/metadata/d630fe4ee6d42383.html`) |
| `https://escholarship.org/uc/item/82n593zd` (the open-access copy named in the Europe PMC record) | `experiments/read_public_page.py` | **not read**: "the site did not serve a robots file to an automated reader, which is taken as a refusal" |

No address waited on a permission request.

## Raw files and SHA-256

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `raw/6e9eae3c00ba1e65.html` (PMC article page) | 198829 | `6e9eae3c00ba1e65bed753ede3dc4d95ac27742849437366c48011e138e09466` |
| `raw/6e9eae3c00ba1e65.html.json` (address, time, SHA-256) | 465 | `8d7d6baed76f0e56c3c525794d66c432845aad61c925ad13c46b7d129219479f` |
| `raw/metadata/efetch_PMC4706282.xml` | 6734 | `fe31497fbdba9af3be6b501e3a7a09062df817524e06b21db869066acdcd056c` |
| `raw/metadata/epmc_10.1208_s12248-015-9845-2.json` | 6316 | `8fa16a5740823fcc46107d12549da5a3781b655ea2d24e9eab7d6c8e26a22173` |
| `raw/metadata/crossref_hosey2016.json` | 20561 | `6db514b514e219d1e243f6228a8b7cae51519d9a79b164b26a47b77fac5f6ab9` |
| `raw/metadata/d630fe4ee6d42383.html` (PMC copyright notice) | 85464 | `d630fe4ee6d42383d868faede01ffa0e9e37d1630dc3e19d1727a235976c2f4e` |

CSV files: `hosey2016_bddcs.csv` 5166 bytes, SHA-256
`4a3794bbf7d716f882db1f333c6dab33ffede5dcfe61b45ebba9fc12dd846e81`; `hosey2016_table1_class_changes.csv` 2011
bytes, `5613b8a35955c164f368b353e4b898972be0afc3a30d27228c3940d6e6d1e6f2`;
`hosey2016_table2_parameter_changes.csv` 689 bytes,
`1b8f2d5b07f95184f0b9bb4bc88e375e6d7a3d854b052743c5dcdd486b45acd6`.

## Licence

**No open licence was verified. The article is free to read in PubMed Central under ordinary copyright.**

- Copyright statement, read in the NCBI record (`raw/metadata/efetch_PMC4706282.xml`) and on the PMC page:
  "© American Association of Pharmaceutical Scientists 2015".
- Europe PMC record: `isOpenAccess` "N", no licence; it names an open-access copy at
  `https://escholarship.org/uc/item/82n593zd` (source "Unpaywall", "Open access"). That page could not be read, so
  the terms of that copy are **not verified**; it is the address for the owner to open.
- PMC copyright notice (`https://pmc.ncbi.nlm.nih.gov/about/copyright/`): "Users of PMC are directly and solely
  responsible for compliance with copyright restrictions and are expected to adhere to the terms and conditions
  defined by the copyright holder."
- Crossref licence field: `http://www.springer.com/tdm` (not read).

## How the CSV files were made

Script: `experiments/exposure_sources/make_hosey2016_csv.py` (run with `python3 -I`, BeautifulSoup 4.15.0 with lxml).
Cell text is kept as printed, with runs of whitespace collapsed to one space. Nothing is computed.

- `hosey2016_bddcs.csv` is Table III: `Drug` → `name`, `Class` → `bddcs_class`, plus `source_table`, `source_row`
  and `note`. The table gives the class only: no percent excreted unchanged, no route, no solubility.
- `hosey2016_table1_class_changes.csv` is Table I: `Compound` → `name`, `Listed class` → `listed_bddcs_class`,
  `Updated class` → `updated_bddcs_class`, `Major route of elimination` → `major_route_of_elimination`. The header
  `Listed class` carries footnote marker a; the footnote, "As listed in Benet et al. (30)", is quoted in `note`.
- `hosey2016_table2_parameter_changes.csv` is Table II: `Compound` → `name`, `Parameter value listed` →
  `parameter_value_listed`, `Updated parameter value` → `updated_parameter_value`. The table has three sub-headings
  that span the value columns ("% dose excreted as unchanged drug in urine", "Solubility", "Administration route");
  each is written into the `parameter` column of the rows beneath it. Values are kept as printed, units included
  (`0.0000204 mg/mL`).

## Rows and non-missing values

| File | Rows | Columns, all without missing values unless stated |
| --- | --- | --- |
| `hosey2016_bddcs.csv` | 180 | `name` 180 (180 distinct), `bddcs_class` 180, `source_table` 180, `source_row` 180, `note` 0 |
| `hosey2016_table1_class_changes.csv` | 13 | `name`, `listed_bddcs_class`, `updated_bddcs_class`, `major_route_of_elimination`, `source_table`, `source_row`, `note`: 13 each |
| `hosey2016_table2_parameter_changes.csv` | 9 | `name`, `parameter`, `parameter_value_listed`, `updated_parameter_value`, `source_table`, `source_row`: 9 each; `note` 0 |

Table III by class: class 1, 56; class 2, 69; class 3, 41; class 4, 11; class 0, 3. Table II holds four changes of
"% dose excreted as unchanged drug in urine" (memantine 71 → 48, pravastatin 20 → 47, ranitidine 30 → 69,
rosuvastatin 5 → 30), two of solubility and three of administration route.

## Checks

- The abstract: "we incorporate more than 175 drugs into the system and amend the classification of 13 drugs".
  Table III has 180 rows and Table I has 13.
- Ten spot checks of Table III, cell text of the page beside the CSV (every eighteenth row,
  `experiments/exposure_sources/spot_checks.py`):

| Name | Class on page | Class in CSV | Row |
| --- | --- | --- | --- |
| Afatinib | 3 | 3 | 1 |
| Bidisomide | 3 | 3 | 19 |
| Chlorpropamide | 0 | 0 | 37 |
| Dihydralazine | 1 | 1 | 55 |
| Flavoxate | 2 | 2 | 73 |
| Isoxepac | 2 | 2 | 91 |
| Methysergide | 1 | 1 | 109 |
| Paritapravir | 2 | 2 | 127 |
| Practolol | 3 | 3 | 145 |
| Tesaglitazar | 2 | 2 | 163 |

- Tables I and II are small enough to have been read in full against the page: 13 and 9 rows, as printed.

## Refused, and other things to know

- `escholarship.org` did not serve a robots file to the reader, which the reader takes as a refusal; not read by
  other means.
- Names are as printed, including the article's own spellings (for example `Paritapravir`).
- `Chlorpropamide` is listed here in class 0 and is also in the 2011 Table V (class 0); other overlaps with the 2011
  table were not examined, because matching is outside this task.
- No text addressed to an AI assistant was seen in the page.
