# Human blood-to-plasma concentration ratios: tables held in the US EPA `httk` repository

> Copied on 11 October 2026 from `data/raw/exposure/blood_to_plasma/sources.md`, the note written beside the raw files (which are not in git). Paths such as `raw/...` are relative to that directory. The transcription scripts are in `experiments/exposure_sources/`.

Retrieved on 12 October 2026 (UTC; the evening of 11 October 2026 local time), in a short look for an openly
distributed table. Transcription only: nothing here is matched to the study's drugs.

**These are the EPA's typed copies of published tables, not the publishers' files.** They were not compared with the
articles, which are behind subscriptions.

## What is here

| CSV | Rows | What it is |
| --- | --- | --- |
| `uchimura2010.csv` | 96 | `datatables/Uchimura 2010 cas.xlsx`: human and rat blood-to-plasma ratio and fraction unbound, per compound with CAS number |
| `paixao2012_table2.csv` | 164 | sheet `Table 2` of `datatables/Paixao-2012.xlsx`: per drug, among others, plasma clearance and `Rb` |
| `httk_human_rblood2plasma.csv` | 280 | the rows of `datatables/HTTK-Chem-Props.txt` (the package's compiled chemical table, 18476 records) that have a value in `Human.Rblood2plasma`, with the reference the table gives for it |

The third file overlaps the first two: of its 280 values, 124 are referenced "Paixao 2012" and 42 "Uchimura 2010";
the rest are referenced "Sternbeck 2012" (28), "ECVAM" (20), "Obach 1999" (17), "Poothong 2017" (14) and 20 other
references with one to five values each (35 values). It also holds environmental chemicals, not only drugs.

## Citations

- Uchimura T, Kato M, Saito T, Kinoshita H. Prediction of human blood-to-plasma drug concentration ratio. Biopharm
  Drug Dispos 2010;31(5-6):286-297. doi:10.1002/bdd.711.
- Paixão P, Gouveia LF, Morais JAG. Prediction of the human oral bioavailability by using in vitro and in silico
  drug related parameters in a physiologically based absorption model. Int J Pharm 2012;429(1-2):84-98.
  doi:10.1016/j.ijpharm.2012.03.019. The workbook's `References` sheet says: "Rb values are from Paixao et al.
  (2009)", that is Paixão P, Gouveia LF, Morais JAG. Prediction of drug distribution within blood. Eur J Pharm Sci
  2009;36(4-5):544-554. doi:10.1016/j.ejps.2008.12.011.
- The compiled table: US EPA, `httk` (high-throughput toxicokinetics) R package, development repository
  https://github.com/USEPA/CompTox-ExpoCast-httk, branch `main`, commit
  `f6605eebc6eb5a483a8824dce6cbddb36c8fb091`.

The article citations were confirmed with the Crossref service (`https://api.crossref.org/works?query.bibliographic=...`,
curl). That the workbooks hold the tables of those articles rests on the file names and on the EPA build script
(`raw/usepa_httk/datatables_load_package_data_tables.R`, lines 672 to 709 and 792 to 808).

## Addresses read

All with curl, from `https://raw.githubusercontent.com/USEPA/CompTox-ExpoCast-httk/f6605eebc6eb5a483a8824dce6cbddb36c8fb091/`:
`datatables/Uchimura%202010%20cas.xlsx`, `datatables/Paixao-2012.xlsx`, `datatables/HTTK-Chem-Props.txt`,
`datatables/load_package_data_tables.R`, `README.md`, `httk/DESCRIPTION`, `httk/LICENSE`. The directory listing was
read with the GitHub connector. One Europe PMC search for open-access articles with "blood-to-plasma" in the title
(curl) returned nothing useful. No address waited on a permission request, and no site refused.

## Raw files and SHA-256

| File | Bytes | SHA-256 |
| --- | --- | --- |
| `raw/usepa_httk/datatables_Uchimura 2010 cas.xlsx` | 17241 | `5e16abfea0d5025db55536765ebb3310443e9f9b40b51d3e22cb0a2d995945bc` |
| `raw/usepa_httk/datatables_Paixao-2012.xlsx` | 33421 | `9210c9d52aecfa0abafe98f83bc2dc0066f46a04fdc0ef53aa1ec24abe0f8fcc` |
| `raw/usepa_httk/datatables_HTTK-Chem-Props.txt` | 9218791 | `09acad773cc0bac380a989751ab8844c73f01de14ebbdbcde35e4ef159c7844a` |
| `raw/usepa_httk/datatables_load_package_data_tables.R` | 131060 | `af400fc0c6ed17ce35f37a4cc3b85d950088badbe15f602007221e2609cb4bce` |
| `raw/usepa_httk/README.md` | 8755 | `6834a24164d8b580b6f193a4398f91c592b310182bc687d5c4e6ae11b67cd335` |
| `raw/usepa_httk/httk_DESCRIPTION` | 9129 | `d897914beb5b5c23a2e12963d563485e8de490557ff41a51844f1aa1730d8bbe` |
| `raw/usepa_httk/httk_LICENSE` | 43 | `f5fbaf0a52c80fd0b440420236b8fa95412d2cc7b87bdc40113046470900e9ec` |

CSV files: `uchimura2010.csv` 9160 bytes, `400b368cd3cfaa197a078b4d46c4359af8dd9886f6ed10e207ff48833b57abea`;
`paixao2012_table2.csv` 11034 bytes, `65f346225ce8301ff03d5d8574c7aa6f5502c498c726ca94273d9862f15805a8`;
`httk_human_rblood2plasma.csv` 16911 bytes, `c76374f08e3a74e64a156623bdec8eb9b26b033dc4f06ead7fda8af8a2fbef90`.

## Licence

- The repository, `README.md` at the commit above, section "License": "License: MIT
  <https://github.com/ULRI-CIRI/httk-dev/blob/dev/httk/LICENSE>". `httk/DESCRIPTION`: "License: MIT + file LICENSE".
  `httk/LICENSE`: "YEAR: 2025" / "COPYRIGHT HOLDER: John Wambaugh". Section "Disclaimer" of the README: "This
  software/application was initially developed by the U.S. Environmental Protection Agency (USEPA). No warranty
  expressed or implied is made regarding the accuracy or utility of the system, nor shall the act of distribution
  constitute any such warranty."
- The compiled table is the package's own data set and falls under that statement. The two workbooks are typed
  copies of tables from copyrighted journal articles; the repository says nothing specific about them, and the
  licences of the articles themselves were **not verified** (not looked up in this short pass). Addresses for the
  owner: https://doi.org/10.1002/bdd.711 and https://doi.org/10.1016/j.ijpharm.2012.03.019.

## How the CSV files were made

Script: `experiments/exposure_sources/make_blood_to_plasma_csv.py` (run with `python3 -I`, openpyxl 3.1.5). Numbers are
written with the shortest decimal text that reads back to the stored number; text cells are written as they stand;
nothing is computed.

- `uchimura2010.csv`, first sheet, header in row 1: `Name` → `name`; `cas` → `cas_number`; `Polarity` → `polarity`;
  `logD` → `log_d`; `Human.fup` → `fraction_unbound_plasma_human`; `Human.Rb2p` → `blood_to_plasma_ratio_human`;
  `Rat.fup` → `fraction_unbound_plasma_rat`; `Rat.Rb2p` → `blood_to_plasma_ratio_rat`. 31 rows hold the spreadsheet
  marker `#N/A` in the two rat columns; those cells are empty in the CSV and quoted in `note`. Five further rows have
  empty rat cells.
- `paixao2012_table2.csv`, sheet `Table 2`, header in row 1: `Drug` → `name`; `CAS` → `cas_number` (empty for five
  drugs); `Foral` → `f_oral`; `frenal` → `f_renal`; `fp` → `fp`; `Clplasma (ml.min−1.kg−1)` →
  `cl_plasma_ml_per_min_per_kg`; `Rb` → `blood_to_plasma_ratio_rb`; `Drug Class` → `drug_class`; `ClH` → `cl_h`;
  `Fabs` → `f_abs`; `Fmet` → `f_met`; `Reference` → `reference` (letter codes; the key is the workbook's
  `References` sheet, for example "(A) Goodman et al. (1996)", "(B) Goodman et al. (2006)").
  **`fp` is left under its source name on purpose: it appears to be the fraction bound, not the fraction unbound.**
  Ziprasidone has `fp` 1 and zolpidem 0.92, and the EPA script computes `Paixao2012.table2$fup <- 1-Paixao2012.table2$fp`
  (line 675). `cl_h` is likewise left under its source name; its definition was not checked against the article.
- `httk_human_rblood2plasma.csv`: `Compound` → `name`; `CAS` → `cas_number`; `DTXSID` → `dtxsid`;
  `Human.Rblood2plasma` → `blood_to_plasma_ratio_human`; `Human.Rblood2plasma.Reference` → `reference`;
  `source_row` is the line of the file. The file is tab-separated without quoting. Three records (Pcb136, Pcb153,
  Pcb155) are broken over two lines inside a name field; the script puts them together, and none of them has a
  blood-to-plasma value. Compound names are in the package's own capitalisation (`Cyclosporin a`, `2,4-d`).

## Rows and non-missing values

| File | Column: non-missing |
| --- | --- |
| `uchimura2010.csv` (96 rows) | `name` 96; `cas_number` 96; `polarity` 96; `log_d` 96; `fraction_unbound_plasma_human` 96; `blood_to_plasma_ratio_human` 96; `fraction_unbound_plasma_rat` 60; `blood_to_plasma_ratio_rat` 60; `note` 31 |
| `paixao2012_table2.csv` (164 rows) | `name` 164; `cas_number` 159; every other source column 164; `note` 0 |
| `httk_human_rblood2plasma.csv` (280 rows) | `name` 280; `cas_number` 280; `dtxsid` 280; `blood_to_plasma_ratio_human` 280 (0.256 to 10); `reference` 280; `note` 0 |

## Spot checks (cell in the raw file, then the CSV)

| File | Name | Raw | CSV |
| --- | --- | --- | --- |
| Uchimura | Alfentanil (71195-58-9) | Human.Rb2p 0.63, Human.fup 0.08 | 0.63, 0.08 |
| Uchimura | Amiodarone (1951-25-3) | Human.Rb2p 0.73, Human.fup 0.0002 | 0.73, 0.0002 |
| Uchimura | Amitriptyline (50-48-6) | Human.Rb2p 0.86; Rat.Rb2p `#N/A` | 0.86; empty, noted |
| Uchimura | Verapamil (52-53-9) | Human.Rb2p 0.77 | 0.77 |
| Uchimura | Zolpidem (82626-48-0) | Human.Rb2p 0.66 | 0.66 |
| Paixão | Acebutolol (37517-30-9) | Rb 1, fp 0.26, Clplasma 6.8 | 1, 0.26, 6.8 |
| Paixão | Acetaminophen (103-90-2) | Rb 1.04, fp 0.2, Clplasma 5 | 1.04, 0.2, 5 |
| Paixão | Alendronate (no CAS) | Rb 1.7, fp 0.78, Clplasma 1.1 | 1.7, 0.78, 1.1 |
| Paixão | Ziprasidone (146939-27-7) | Rb 0.81, fp 1, Clplasma 11.7 | 0.81, 1, 11.7 |
| Paixão | Zolpidem (82626-48-0) | Rb 0.76, fp 0.92, Clplasma 4.5 | 0.76, 0.92, 4.5 |
| httk | Caffeine (58-08-2) | 1.08, Sternbeck 2012 | 1.08, Sternbeck 2012 |
| httk | Warfarin (81-81-2) | 0.55, Obach 1999 | 0.55, Obach 1999 |
| httk | Diazepam (439-14-5) | 0.78, Sternbeck 2012 | 0.78, Sternbeck 2012 |
| httk | Amiodarone (1951-25-3) | 0.73, Paixao 2012 | 0.73, Paixao 2012 |

## Other things to know

- The two tables disagree for some drugs (zolpidem 0.66 in Uchimura, 0.76 in Paixão); both are kept as given.
- In the compiled table one value per chemical survives the package's loading order. The EPA build script loads
  the Uchimura values with `overwrite=FALSE` (lines 794 to 799) and later the "Sternbeck 2012" values with
  `overwrite=TRUE` for every chemical whose value is not referenced "Paixao 2012" (lines 867 to 878, comment
  "overwrite all Rb2p but paixao"). Consistent with that, only 42 of the 96 Uchimura values appear there; the two
  workbooks are the fuller record of their sources.
- No text addressed to an AI assistant was seen in any file.
