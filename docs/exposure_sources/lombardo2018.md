# Lombardo, Berellini and Obach 2018: human intravenous pharmacokinetic parameters for 1352 compounds

> Copied on 11 October 2026 from `data/raw/exposure/lombardo2018/sources.md`, the note written beside the raw files (which are not in git). Paths such as `raw/...` are relative to that directory. The transcription scripts are in `experiments/exposure_sources/`.

Retrieved on 12 October 2026 (UTC; the evening of 11 October 2026 local time). Transcription only: nothing here is
matched to the study's drugs.

## Citation

Lombardo F, Berellini G, Obach RS. Trend analysis of a database of intravenous pharmacokinetic parameters in humans
for 1352 drug compounds. Drug Metab Dispos 2018;46(11):1466-1477. doi:10.1124/dmd.118.082966. PMID 30115648. Not in
PubMed Central.

The table transcribed is the article's supplementary Table S1. A text box in the workbook reads: "SUPPORTING
INFORMATION for: "TREND ANALYSIS OF A DATABASE OF INTRAVENOUS PHARMACOKINETIC PARAMETERS IN HUMANS FOR 1352 DRUG
COMPOUNDS" Franco Lombardo, Giuliano Berellini, R. Scott Obach Table S1. Summary of human iv PK parameters for 1352
compounds. REVISED and CORRECTED October 4 2018". The workbook properties give creation on 2018-07-15 and last
modification on 2018-10-09.

## Route

The publisher's copy was not read: the journal site's robots file disallows the supplementary-file address (see
"Refused"). The file used is the copy that the United States Environmental Protection Agency keeps in the public
development repository of its `httk` R package, where the package's build script reads it
(`datatables/load_package_data_tables.R`, lines 545 to 547: `read_excel(` `"Lombardo2018-Supplemental_82966_revised_corrected.xlsx",` `skip=8))`).

Repository: https://github.com/USEPA/CompTox-ExpoCast-httk, branch `main`, commit
`f6605eebc6eb5a483a8824dce6cbddb36c8fb091`.

## Addresses read

| Address | How | Result |
| --- | --- | --- |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:"10.1124/dmd.118.082966"&resultType=core&format=json` | curl | record; `isOpenAccess` N, `inEPMC` N, no PMC id, `hasSuppl` N (`raw/metadata/epmc_10.1124_dmd.118.082966.json`) |
| `https://api.crossref.org/works/10.1124/dmd.118.082966` | curl | licence and copyright metadata (`raw/metadata/crossref_lombardo2018.json`) |
| `https://dmd.aspetjournals.org/content/suppl/2018/08/16/dmd.118.082966.DC1` | `experiments/read_public_page.py` | **not read**: "the site's robots file disallows this address for automated readers" |
| `https://doi.org/10.1124/dmd.118.082966` | `experiments/read_public_page.py` | a redirect page of the publisher's linking service, 2670 bytes (`raw/metadata/d96da739fd918515.html`); not followed further |
| `https://www.elsevier.com/tdm/userlicense/1.0/` | `experiments/read_public_page.py` | read (`raw/metadata/98d7bac9e7114f07.html`) |
| `https://www.elsevier.com/legal/tdmrep-license` | `experiments/read_public_page.py` | read (`raw/metadata/4b1de4fc446a4a68.html`) |
| GitHub code and repository search, directory listings and README of `USEPA/CompTox-ExpoCast-httk`, `srijitseal/PKSmart`, `bartongroup/DDUpred`, `dreadlesss/pharmacokinetic_prediction`, `CDDLeiden/PK-in-generative-drug-design` | GitHub connector | found the copies listed below |
| `https://raw.githubusercontent.com/USEPA/CompTox-ExpoCast-httk/f6605eebc6eb5a483a8824dce6cbddb36c8fb091/datatables/Lombardo2018-Supplemental_82966_revised_corrected.xlsx` | curl | **the source file** |
| same commit: `datatables/Lombardo2018-Supplemental_82966_revised_corrected.csv`, `datatables/README.md`, `datatables/load_package_data_tables.R`, `README.md`, `httk/DESCRIPTION`, `httk/LICENSE` | curl | saved in `raw/usepa_httk/` (`/` in the path written as `_` in the file name) |
| `https://raw.githubusercontent.com/dreadlesss/pharmacokinetic_prediction/master/{dataset.xlsx,README.md,LICENSE}` (commit `5839c3033157143a81a16e56ff7083a2d5b5916c`) | curl | cross-check copy, `raw/crosscheck_dreadlesss/` |
| `https://raw.githubusercontent.com/bartongroup/DDUpred/2ccf4535227d7d41f54a75a5fb8f35219aba66cd/training_data/{lombardo.csv.gz,lombardo.info.csv}` and `README.md` | curl | cross-check copy (a 1141-row subset), `raw/crosscheck_ddupred/` |
| `https://raw.githubusercontent.com/srijitseal/PKSmart/661aa31617a3fd4c05ccaaff1ff5ed168163cf68/data/Human_PK_data.csv` | curl | inspected, not kept: a derived table of 1283 standardised structures without names or CAS numbers |
| `https://zenodo.org/api/records/10611606` | curl | metadata only (PKSmart archive, CC BY 4.0, one 45 MB zip); not downloaded |

Web searches (the search tool) found no further copy. No address waited on a permission request.

## Raw files and SHA-256

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `raw/usepa_httk/datatables_Lombardo2018-Supplemental_82966_revised_corrected.xlsx` | 561132 | `4c462bf75759b657e0c90542434a83610fe65d8c00ca164decab3955da75a7fb` |
| `raw/usepa_httk/datatables_Lombardo2018-Supplemental_82966_revised_corrected.csv` | 841946 | `b2ce65f836cbea105bf131aeeb620d2bb9146077677835ce91faeb10f09a27c0` |
| `raw/usepa_httk/datatables_load_package_data_tables.R` | 131060 | `af400fc0c6ed17ce35f37a4cc3b85d950088badbe15f602007221e2609cb4bce` |
| `raw/usepa_httk/datatables_README.md` | 5612 | `7e44012351f952e092b5c9479ecd6ff493381685b3fa679b27463b41d524ee36` |
| `raw/usepa_httk/README.md` | 8755 | `6834a24164d8b580b6f193a4398f91c592b310182bc687d5c4e6ae11b67cd335` |
| `raw/usepa_httk/httk_DESCRIPTION` | 9129 | `d897914beb5b5c23a2e12963d563485e8de490557ff41a51844f1aa1730d8bbe` |
| `raw/usepa_httk/httk_LICENSE` | 43 | `f5fbaf0a52c80fd0b440420236b8fa95412d2cc7b87bdc40113046470900e9ec` |
| `raw/crosscheck_dreadlesss/dataset.xlsx` | 516908 | `8b61fd458e710ca25ea99c50906f457e4c940b120dfe85ce940ceea9bb83be78` |
| `raw/crosscheck_dreadlesss/README.md` | 1921 | `2aa79dba6bce892a1e56cc4e596aab58753ce303b5a516a46575f2f2a88f4583` |
| `raw/crosscheck_dreadlesss/LICENSE` | 35149 | `3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986` |
| `raw/crosscheck_ddupred/lombardo.csv.gz` | 1077641 | `8d85cd61d202d491367d784ac556a26f078516a56ecc8ff6ed9a2fbb8106ae60` |
| `raw/crosscheck_ddupred/lombardo.info.csv` | 820 | `f7e166de46109228cefe1a3353b753bb96683c1687bb3ca542f90cf5f55b02a0` |
| `raw/crosscheck_ddupred/README.md` | 1520 | `4a4d3eb592e77e5a92f70a48a1758576ba6e772e2a2261d3cef6c18072ade277` |
| `raw/metadata/crossref_lombardo2018.json` | 18535 | `cc2975bbd62cb99f14e4170dfc85a488a6934ddb6f36a78e4b1b05c2e5fd3f9a` |
| `raw/metadata/epmc_10.1124_dmd.118.082966.json` | 6698 | `1615fb722166dc231da760f0f47a9db9b6a9de3f1b24504159ae55574cfa2cc4` |
| `raw/metadata/d96da739fd918515.html` (doi.org redirect page) | 2670 | `d96da739fd9185154794429c49e74124a1a92c7e9e54e4ae74e17fc5311e4fc3` |
| `raw/metadata/98d7bac9e7114f07.html` (Elsevier text and data mining licence page) | 209290 | `98d7bac9e7114f07140c4942e0089a8c83b97d4d1c3bd93868d1528863e1481c` |
| `raw/metadata/4b1de4fc446a4a68.html` (Elsevier TDM reservation page) | 198744 | `4b1de4fc446a4a689d71fcdd57bfdda19f5f6a986a25abe0412035f1335e9153` |

Each `.html` file saved by `read_public_page.py` has a `.json` record beside it with address, time and SHA-256.
The doi.org redirect page embeds the network address of the proxy that made the request (not the owner's).
`lombardo2018.csv`: 860771 bytes, SHA-256 `3d9e0507b284a61d6f7536fd81682e1f727fcfacc33a95904e2829cd33c79a13`.

## Licence: checked again on 12 October 2026, and permission is required

The sections below were written on 11 October 2026 and concluded "not verified". A second pass on 12 October 2026
tried the routes they had not, and the position is now negative rather than unknown. The additions, each quoted with
where it was read; `docs/label_source_review.md` section 13.9 is the summary.

- **The only licence the article has forbids this use.** Crossref records two licence URLs for the DOI, both with
  `"content-version": "tdm"`, and no article user licence. From the first: "The TDM click-through agreement does not
  allow you to: ... Substantially or systematically reproduce, retain or redistribute the Dataset." and "Text mining
  access for subscription content is provided to subscribers for non-commercial research purposes."
  (`https://www.elsevier.com/tdm/userlicense/1.0/`, read with `experiments/read_public_page.py`, 12 October 2026).
  The publisher's copyright policy: a user may mine subscription articles "without sharing any adaptation of the
  original content with others", and "additional uses/sharing beyond those listed would require permissions to be
  requested from Elsevier" (`https://www.elsevier.com/about/policies-and-standards/copyright`, same date).
- **Free to read is not free to redistribute.** OpenAlex reports `is_oa: true`, `oa_status: "bronze"` and
  `license: null`, with the publisher's own PDF as the open location (`api.openalex.org`, curl, 12 October 2026). The
  note below that records Europe PMC's `isOpenAccess: N` should be read with this beside it: the full text is free to
  read at the publisher and was unread here only because of the robots refusal.
- **No post-embargo licence will attach.** The journal is not among the titles in the publisher's open-archive
  programme, whose page states that articles in it "are free for everyone to read and download, subject to Elsevier's
  user license" after an embargo (`https://www.elsevier.com/open-access/open-archive`, 12 October 2026). The society's
  journals "will remain hybrid" under the partnership announcement (same date).
- **The data-as-facts route is narrower than it looks.** "Creating an original figure or table from data or factual
  information that was not previously in figure or table format typically does not require permission, subject to
  proper acknowledgement" — and Table S1 was already a table
  (`https://www.elsevier.com/about/policies-and-standards/copyright/permissions`, 12 October 2026). The same page:
  "As a general rule, permission should be sought from the rights holder to reproduce any substantial part of a
  copyrighted work. This includes any text, illustrations, charts, tables, photographs, or other material from
  previously published sources."
- **The authors keep the data.** "Elsevier does not require, nor ever request, rights to the underlying research
  itself. Authors (and institutions or employers) will continue to retain patent, trademark and other intellectual
  property rights relating to the research, as well as the right to use their research data freely and without
  restriction." (copyright page, as above). That is the shorter of the two routes: a deposit of Table S1 by its
  authors under CC BY or CC0. No deposit exists today; Zenodo, figshare and Dryad were all searched on 12 October 2026.
- **The httk route does not carry.** A search of the whole `USEPA/CompTox-ExpoCast-httk` repository finds the table
  named in eight files, every one a citation or build code, with no permission, provenance or rights statement; the
  per-dataset R documentation lists the citation under `\references` and nothing else. The EPA's data licence covers
  data "produced by the U.S EPA" (`pasteur.epa.gov/license/sciencehub-license.html`, 12 October 2026), which this is
  not, and the package's own `LICENSE` names a package author as copyright holder, which speaks to the package.
- **A warning for a later pass.** The httk team's figshare deposit, doi:10.23645/epacomptox.6062791, is licensed
  "Creative Commons Zero v1.0 Universal" (DataCite and figshare APIs, curl, 12 October 2026) but its file list holds
  one 0-byte `index.html`. It is a metadata stub and does not contain this table. Do not read that CC0 as covering it.
- **The permissions route, if it is wanted.** "ASPET handles all permissions requests through the Copyright Clearance
  Center's Republication Service." (`https://www.aspet.org/aspet/journals/permissions-request`, 12 October 2026); the
  journal's link from that page is `https://marketplace.copyright.com/rs-ui-web/mp/search/journal/1521009X`. A request
  has to name the revised supplement explicitly, the whole of it rather than an excerpt, and ask for non-exclusive
  worldwide electronic redistribution for non-commercial research with the society's credit line. Worth asking
  outright whether the society considers the factual values themselves, as distinct from the selection and
  arrangement of the table, to fall within its claim.
- **No precedent found.** A full-text search of the open-access literature for a paper reusing this table with a
  stated permission returns nothing, though 29 articles carry the society's conventional permission credit line for
  other material, so the route is used.

**Correction to the "Route" and "Refused" sections below, and it matters for every future retrieval of anything on
this host.** Those sections say the robots file disallows the supplementary-file address. It is broader than that: the
host's robots file carries `Disallow: /` for the user agents an AI assistant's reader announces itself under, so
**every address on `dmd.aspetjournals.org` is refused**, the article landing page and the free full-text PDF included.
Trying another path there is futile; the user has to open them personally.

## Licence (the first pass, 11 October 2026)

**The article and its supplement carry no open licence that could be verified. The copy used comes from a repository
whose stated licence is MIT; whether that licence can cover a publisher's supplementary file is not established.**
The owner should decide whether this is sufficient before the table is redistributed.

What was read, and where:

- Article copyright, Crossref record (`raw/metadata/crossref_lombardo2018.json`, field `assertion`): "© 2018 by The
  American Society for Pharmacology and Experimental Therapeutics". PubMed (PubMed connector, PMID 30115648):
  "Copyright © 2018 by The American Society for Pharmacology and Experimental Therapeutics." Europe PMC:
  `isOpenAccess` "N", no licence field.
- Licence addresses in the Crossref record: `https://www.elsevier.com/tdm/userlicense/1.0/` and
  `https://www.elsevier.com/legal/tdmrep-license` (the journal is now hosted by Elsevier). From the second page:
  "For all other content, TDM rights are reserved and consent for TDM needs to be requested via the email contact".
  From the first: "Text mining access for subscription content is provided to subscribers for non-commercial
  research purposes." The doi.org redirect page carries `<meta name="tdm-reservation" content="1" />`.
- The workbook itself: no licence or terms statement (only the text box quoted above).
- The redistributing repository, `README.md` of `USEPA/CompTox-ExpoCast-httk` at the commit above, section
  "License": "License: MIT <https://github.com/ULRI-CIRI/httk-dev/blob/dev/httk/LICENSE>"; `httk/DESCRIPTION`:
  "License: MIT + file LICENSE"; `httk/LICENSE`: "YEAR: 2025" / "COPYRIGHT HOLDER: John Wambaugh". The same README,
  section "Disclaimer": "This software/application was initially developed by the U.S. Environmental Protection
  Agency (USEPA). No warranty expressed or implied is made regarding the accuracy or utility of the system, nor shall
  the act of distribution constitute any such warranty." The repository says nothing specific about the licence of
  the third-party data files in `datatables/`.
- Cross-check copy `dreadlesss/pharmacokinetic_prediction`: `LICENSE` is the GNU General Public License version 3;
  its README says "dataset.xlsx (download from supporting information)".
- Cross-check copy `bartongroup/DDUpred`: no licence file or statement found.

Routes that were closed: the publisher's article and supplement pages (robots rule), the article's data
availability statement (subscription full text, not read). Address for the owner to open:
https://doi.org/10.1124/dmd.118.082966 (supplement: "Supplemental Data" link on that page).

## How the CSV was made

Script: `experiments/exposure_sources/make_lombardo2018_csv.py` (run with `python3 -I`, openpyxl 3.1.5), on sheet
`Data_sheet` of the EPA workbook. Row 9 is the header, rows 10 to 1361 are the 1352 compounds; rows 1 to 8 hold only
the text box.

- Source column to CSV column: `Name` → `name`; `CAS #` → `cas_number`; `SMILES` → `smiles`; `human VDss (L/kg)` →
  `vdss_l_per_kg`; `human CL (mL/min/kg)` → `clearance_ml_per_min_per_kg`; `fraction unbound in plasma (fu)` →
  `fraction_unbound_plasma`; `MRT (h)` → `mrt_h`; `terminal  t1/2 (h)` → `half_life_h`; `Reference` → `reference`;
  `Comments ` → `source_comments`; `Notes` → `source_notes`; `Year of first disclosure` →
  `year_of_first_disclosure`; `MW` → `mw`; `HBA` → `hba`; `HBD` → `hbd`; `TPSA_NO` → `tpsa_no`; `RotBondCount` →
  `rot_bond_count`; `moka_ionState7.4` → `moka_ion_state_7_4`; `MoKa.LogP` → `moka_logp`; `MoKa.LogD7.4` →
  `moka_logd_7_4`. Added: `source_row` (row number in the sheet) and `note` (transcription remarks; the source's own
  remarks are in `source_comments` and `source_notes`).
- Numbers are written with the shortest decimal text that reads back to the number stored in the workbook (the
  workbook stores, for example, `8.6999999999999993` for 8.7). Nothing is rounded, computed or filled in.
- 45 `fraction_unbound_plasma` cells hold a long decimal in the workbook (for example alfentanil
  `8.5999999999999896E-2`, written here as `0.0859999999999999`; candoxatrilat `0.8723436361917`). They are written
  unrounded, and `note` quotes how the EPA's comma-separated export of the same workbook shows the cell (`0.086`,
  `0.872`).
- 13 compounds have the text `n/a` in `MoKa.LogP` and `MoKa.LogD7.4`; those cells are empty in the CSV and quoted in
  `note`.
- 9 names end in a space in the workbook (for example `'Amiodarone '`); the space is stripped and the fact is in
  `note`.
- No pharmacokinetic cell carries a qualifier such as `<` or `~`; an empty cell is an empty CSV field.
- `reference`, `source_comments` and `source_notes` keep the line breaks and the zero-width characters of the
  workbook, so the CSV has quoted fields that span lines (1352 records on 1645 lines). Read it with a CSV parser.

## Rows and non-missing values

1352 rows (the abstract: "a data set of 1352 drugs").

| Column | Non-missing | Distinct |
| --- | --- | --- |
| `name` | 1352 | 1352 |
| `cas_number` | 1352 | 1352 |
| `smiles` | 1352 | 1351 |
| `vdss_l_per_kg` | 1315 | 366 |
| `clearance_ml_per_min_per_kg` | 1350 | 569 |
| `fraction_unbound_plasma` | 920 | 175 |
| `mrt_h` | 1309 | 525 |
| `half_life_h` | 1335 | 502 |
| `reference` | 1352 | 1329 |
| `source_comments` | 1021 | 949 |
| `source_notes` | 112 | 110 |
| `year_of_first_disclosure` | 1352 | 114 |
| `mw` | 1352 | 1208 |
| `hba` | 1352 | 57 |
| `hbd` | 1352 | 40 |
| `tpsa_no` | 1352 | 799 |
| `rot_bond_count` | 1352 | 57 |
| `moka_ion_state_7_4` | 1352 | 5 |
| `moka_logp` | 1339 | 269 |
| `moka_logd_7_4` | 1339 | 314 |
| `source_row` | 1352 | 1352 |
| `note` | 66 | 37 |

Ranges in the file: VDss 0.03 (Secretin (porcine)) to 700 L/kg (Hydroxychloroquine); CL 0.0037
(Hydroxystaurosporine, 7-) to 1070 mL/min/kg (Artesunate); fu 0.0002 (Amiodarone) to 1; MRT 0.01 to 1400 h;
t1/2 0.02 to 1344 h; year of first disclosure 1878 to 2014.

## Checks

`experiments/exposure_sources/verify_lombardo2018.py`:

- All 6229 numeric cells of the five pharmacokinetic columns equal the numbers in the sheet XML of the workbook
  (read without openpyxl); 0 differences.
- Against `raw/crosscheck_dreadlesss/dataset.xlsx` (an independently uploaded copy, 1352 rows, header in row 1,
  last saved 2020-04-15): all 1352 compounds match on name and CAS number, and all five pharmacokinetic columns
  agree; 0 differences.
- Against the EPA's own CSV export: 1352 rows; the only differences are the 45 long-decimal fu cells, which the
  export shows as displayed.
- Against `raw/crosscheck_ddupred/lombardo.csv.gz` (1141 rows): 1130 match on CAS number and agree except for the
  same long-decimal fu cells; the other 11 have CAS numbers that a spreadsheet turned into dates in that copy (for
  example Ajmaline `07/12/4360`).
- The article's full text is behind the subscription, so its summary statistics were not compared; the count of
  1352 is confirmed by the abstract and by the workbook's own title.

Twelve spot checks (text of the sheet XML, then the CSV value where it is written differently):

| Name in workbook | CAS # | VDss (L/kg) | CL (mL/min/kg) | fu | MRT (h) | t1/2 (h) | Sheet row |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Amiodarone ` | 1951-25-3 | 60 | 1.9 | 2.0000000000000001E-4 (0.0002) | 620 | 820 | 75 |
| `Diazepam` | 439-14-5 | 1 | 0.38 | 2.3E-2 (0.023) | 44 | 42 | 391 |
| `Warfarin` | 81-81-2 | 0.13 | 5.5E-2 (0.055) | 1.4999999999999901E-2 (0.0149999999999999) | 39 | 29 | 1332 |
| `Haloperidol` | 52-86-8 | 17 | 7.8 | 0.08 | 36 | 35 | 608 |
| `Morphine` | 57-27-2 | 2.2999999999999998 (2.3) | 26 | 0.65 | 1.5 | 2 | 843 |
| `Caffeine` | 58-08-2 | 0.63 | 1.4 | 0.64 | 7.5 | 4.9000000000000004 (4.9) | 193 |
| `Midazolam` | 59467-70-8 | 1.1000000000000001 (1.1) | 5.3 | 1.7000000000000001E-2 (0.017) | 3.5 | 3.1 | 820 |
| `Propranolol` | 525-66-6 | 3.1 | 12 | 0.13 | 4.3 | 3.4 | 1038 |
| `Hydroxychloroquine` | 118-42-3 | 700 | 11 | 0.56999999999999995 (0.57) | 1300 | 850 | 618 |
| `Citalopram` | 59729-33-8 | 12 | 4.3 | 0.2 | 47 | 33 | 302 |
| `Venlafaxine` | 93413-69-5 | 4.4000000000000004 (4.4) | 14 | 0.73 | 5.2 | 5 | 1308 |
| `Artesunate` | 88495-63-0 | 15 | 1070 | 0.25 | 0.32 | 0.22 | 104 |

## Refused, and other things to know

- `dmd.aspetjournals.org`: the robots file disallows the supplementary-file address for automated readers. Not read
  by any other means.
- The publisher marks the article with a text-and-data-mining reservation (see "Licence").
- The EPA build script corrects eight CAS numbers and one name's CAS number after reading this workbook (lines
  548 to 567 of `raw/usepa_httk/datatables_load_package_data_tables.R`, for example "85650-52-8" to "61337-67-5").
  Those corrections are the EPA's, are **not** applied here, and may be useful when CAS numbers are matched.
- Obach, Lombardo and Waters 2008 (670 compounds) was not retrieved, because the 2018 table was had in full. The same
  EPA directory holds `datatables/Obach2008.xlsx` (50543 bytes) should it be wanted.
- No text addressed to an AI assistant was seen in any file or page.
