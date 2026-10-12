# Schulz et al. 2020: therapeutic and toxic blood concentrations (Additional file 1)

> Copied on 11 October 2026 from `data/raw/exposure/schulz2020/sources.md`, the note written beside the raw files (which are not in git). Paths such as `raw/...` are relative to that directory. The transcription scripts are in `experiments/exposure_sources/`.

Retrieved and transcribed on 12 October 2026 (UTC; the evening of 11 October 2026 local time). Nothing here is
matched to the study's drugs, and no value was converted, estimated or filled in.

## Citation

Schulz M, Schmoldt A, Andresen-Streichert H, Iwersen-Bergmann S. Revisited: Therapeutic and toxic blood
concentrations of more than 1100 drugs and other xenobiotics. Crit Care 2020;24(1):195.
doi:10.1186/s13054-020-02915-5. PMID 32375836. PMCID PMC7201985.

The table is not in the article body (the full-text XML holds no table). It is "Additional file 1", a 208-page PDF:
pages 1-82 the table, pages 82-84 units and the meaning of the brackets, pages 85-118 the 633 numbered annotations
and the abbreviations, pages 119-208 the reference list.

## Addresses read

| Address | How | Saved as |
|---|---|---|
| `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:"10.1186/s13054-020-02915-5"&resultType=core&format=json` | curl (Europe PMC service) | `raw/europepmc_search_schulz2020.json` |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7201985/fullTextXML` | curl (Europe PMC service) | `raw/PMC7201985_fullTextXML.xml` |
| `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7201985/supplementaryFiles` | curl (Europe PMC service; answered 200 on the first try) | `raw/PMC7201985_supplementaryFiles.zip` |
| `https://static-content.springer.com/esm/art%3A10.1186%2Fs13054-020-02915-5/MediaObjects/13054_2020_2915_MOESM1_ESM.pdf` | curl (publisher's supplementary-file host) | `raw/springer_static_content/13054_2020_2915_MOESM1_ESM.pdf` |
| PubMed connector, `get_copyright_status` for PMID 32375836 | connector | not saved; it answered licence "CC BY 4.0", source "pmc" |

No address waited on a permission request. No site refused. No text addressed to an AI assistant was seen in the
XML or in the PDF (searched for "ignore previous", "AI assistant", "language model", "system prompt" and similar).

## SHA-256 of the raw files

| File | SHA-256 |
|---|---|
| `raw/europepmc_search_schulz2020.json` | `93adff3e73f3774025757b41c88eeb2319601e8320093ab46fadb9bbfe90afeb` |
| `raw/PMC7201985_fullTextXML.xml` | `f3229f74770f20d53236d96f5fb35eae3c964236dc2ec7df584b2a31dfb5e275` |
| `raw/PMC7201985_supplementaryFiles.zip` | `6917dec68ed6d1e25d1de67f4bfea04fda55bac4191ce8f3ecc4a4fb6750af4a` |
| `raw/supplementaryFiles_extracted/13054_2020_2915_MOESM1_ESM.pdf` (unzipped from the archive) | `758e3b6fc8dade1e4c47adba02c2375bf91f63da16e7daa7d8ff3d4125a0a1e6` |
| `raw/supplementaryFiles_extracted/13054_2020_2915_Fig1_HTML.jpg` | `62b6be20820279b7ace6debb6db7ed1d771fba8d53736df666122565e3dab322` |
| `raw/supplementaryFiles_extracted/13054_2020_2915_Fig1_HTML.gif` | `54cefe1e06bf351c9e8384252e4ba8ecb39ca56d8403315a39739551f1e73c2d` |
| `raw/springer_static_content/13054_2020_2915_MOESM1_ESM.pdf` | `758e3b6fc8dade1e4c47adba02c2375bf91f63da16e7daa7d8ff3d4125a0a1e6` |

The PDF from Europe PMC and the PDF from the publisher's host are byte-identical.

## Licence (verified)

Read in `raw/PMC7201985_fullTextXML.xml`, element `<permissions><license>` (the PubMed Central full text served by
Europe PMC). Copyright statement: "© The Author(s) 2020". Licence paragraph, verbatim:

> Open Access This article is licensed under a Creative Commons Attribution 4.0 International License, which
> permits use, sharing, adaptation, distribution and reproduction in any medium or format, as long as you give
> appropriate credit to the original author(s) and the source, provide a link to the Creative Commons licence, and
> indicate if changes were made. The images or other third party material in this article are included in the
> article's Creative Commons licence, unless indicated otherwise in a credit line to the material. If material is
> not included in the article's Creative Commons licence and your intended use is not permitted by statutory
> regulation or exceeds the permitted use, you will need to obtain permission directly from the copyright holder.
> To view a copy of this licence, visit http://creativecommons.org/licenses/by/4.0/. The Creative Commons Public
> Domain Dedication waiver (http://creativecommons.org/publicdomain/zero/1.0/) applies to the data made available
> in this article, unless otherwise stated in a credit line to the data.

The same XML lists the additional file with the caption: "Additional file 1. Therapeutic (“normal”), toxic, and
comatose-fatal blood-plasma concentrations (mg/L, if not stated otherwise) in man. A table containing data for more
than 1100 drugs and other xenobiotics in man including annotations and references." The Europe PMC record gives
`license: cc by` and `isOpenAccess: Y`. The PDF itself carries no licence line; every page carries this note:
"Note on liability: The data and information have been extracted from scientific sources (primary and secondary
literature) and have been carefully reviewed. Nevertheless, transcription errors cannot be excluded. No
responsibility is taken for the correctness of the information or data included. As a matter of principle, this
compilation should be used with expert knowledge."

## Units, as the source states them

- Table header: "Blood-plasma concentration (mg/L)" over "therapeutic (“normal”)", "toxic (from)" and
  "comatose-fatal (from)"; "t½ (h)"; "References".
- Page 82: "Drug concentrations are mentioned as mg/L (= μg/mL), if not otherwise stated."
- Page 82-83: "To convert to SI units: mg/L (= μg/mL) x CF = μmol/L. To convert from SI units: μmol/L ÷ CF = mg/L
  (= μg/mL)" (CF from the molecular weight; nothing was converted here).
- Page 83: "Therapeutic: blood-plasma/serum concentrations (in general, trough at steady state) observed following
  therapeutically effective doses; no or only minimal side effects (drugs); “normal”: concentrations associated
  with no or only minimal toxic effects (other xenobiotics, including ‘recreational’ use). Toxic:
  blood-plasma/serum concentrations which produce toxicity/clinically relevant side effects/symptoms.
  Comatose-fatal: blood-plasma/serum (comatose) concentrations and whole blood (fatal) concentrations reported to
  have caused coma and death, respectively."
- Page 84, brackets: "Meaning of brackets e.g., therapeutic plasma concentration = (0.02-) 0.08-0.25 or t½ = (4-)
  8-22 (-35): These brackets may have slightly different meanings and depend on the substance. [...] So, if not
  specifically mentioned or annotated, the data in brackets widen the range but do not cover most of the cases."
- Page 118: "Abbreviations: appr. approximately, AUC area under the (blood plasma concentration-time) curve, Cmax
  maximum (peak) plasma/serum concentration, Cmin minimum (trough) plasma/serum concentration (usually at steady
  state), h hour/hours, [...] t½ terminal elimination half-life (if not stated otherwise), tmax time to peak
  concentration (Cmax)".
- 64 concentration cells state another unit in the cell (ng/mL, ng/L, ng/g, µmol/L, pg/mL, mg/kg, %); 80 half-life
  cells state days, weeks, months, years or minutes. These keep their text, have no numeric columns, carry the
  form `not_parsed_other_unit` and a note. 68 further cells carry an annotation whose text mentions another unit
  (for example lithium, annotation 79, in mmol/L); they are flagged in `note` and the annotation is in `remarks`.

## Files written

| File | Content | SHA-256 |
|---|---|---|
| `schulz2020.csv` | the table, 1,246 rows | `894eb5fa7864e3502a3ca60461a7abcec8e2b50bc9124af92d8b161ef57266b5` |
| `schulz2020_annotations.csv` | the 633 numbered annotations (number, PDF page, text) | `979dc9f88b78292b0abc97bdb79fe1c5cd81db7b2c3fac94975416821a2e17e3` |
| `derived/schulz2020_cells.json` | every table cell as printed lines with character runs (body, raised, lowered) | `cb20f9e0505ec8da3dbc64ccbfe2c0d0a6f9db53803c889e6f991455822e2812` |
| `derived/schulz2020_annotations.json` | annotations and the abbreviations paragraph | `aea35d5abadb07877d6825b610fb570d204a304c112b11828631d808a8df9374` |
| `derived/MOESM1_pdftotext_layout.txt` | `pdftotext -layout` of the whole PDF (used only for checking; it also holds the reference list, pages 119-208) | `350dd3ed98091a645174e7524c9317517b555f28cdd394d2f8780dc3729bf005` |
| `derived/schulz2020_summary.json` | counts printed by the build script | `eaaa107784a7395c781620c4d8e387aa4ba7ea02bbe922b75755caf2db40202c` |

## How the CSV was made

Scripts (copies in `experiments/exposure_sources/`; they were run with `python3 -I` from a directory that holds no downloaded
file; pdfplumber 0.11.10, pdfminer.six 20260107, poppler 24.02.0, Python 3.13):

1. `schulz_cells.py <pdf> 1 82 derived/schulz2020_cells.json`. A table row is the band between two horizontal rules
   that the PDF draws under each row; the six columns have their vertical rules at the same x positions on every
   page. Characters at 6.5 pt that sit above the baseline are the raised annotation numbers; 6.5 pt characters on
   the baseline are subscripts (only in "Vitamin B6" and "Vitamin B12", four cells). 1,323 bands = 82 header bands + 1,241 printed
   rows; all 66,903 non-blank characters between the first and last rule of the 82 pages were assigned to a cell.
2. `schulz_annotations.py <pdf> 85 118 derived/schulz2020_annotations.json`. 633 annotations, numbered 1 to 633
   without a gap.
3. `schulz_build_csv.py derived/schulz2020_cells.json derived/schulz2020_annotations.json schulz2020.csv
   schulz2020_annotations.csv derived/schulz2020_summary.json`.

Columns of `schulz2020.csv`:

- `source_row` (1 to 1,241, the printed row in page order), `source_line`, `page` (page of the PDF).
- `name`: the substance name as printed, without the raised annotation numbers. `name_printed`: the same with them,
  written `^{...}`. `name_annotations`: the numbers.
- For `therapeutic`, `toxic`, `comatose_fatal`: `<column>_printed` (the cell as printed, raised annotation numbers
  written `^{...}`), `<column>_text` (the same without them), `<column>_annotations`, `<column>_low_mg_per_l`,
  `<column>_high_mg_per_l`, `<column>_parse`.
- `half_life_printed`, `half_life_text`, `half_life_annotations`, `half_life_low_h`, `half_life_high_h`,
  `half_life_parse`.
- `references_printed`: the reference numbers as printed (the list itself is in the PDF, pages 119-208).
- `see`: the target when the therapeutic cell reads "see ..." (113 rows are such cross-references).
- `remarks`: the text of every annotation cited in the row, as `[number] text`, joined by ` || `.
- `note`: what the transcription did or found in this row.

Rules, all visible in the `_parse` columns:

- A number goes into a numeric column only when the cell text has one of these forms (N a number as printed):
  `range` "N-N"; `single` "N" (written to both low and high; in the toxic and comatose-fatal columns the header
  reads "(from)"); `up_to` "-N" and `less_than` "< N" (high only); `greater_than` "> N" and `at_least` "≥ N" (low
  only); `approx_range`, `approx_single`, `approx_up_to` (the same, preceded by "appr."); and
  `range_with_bracketed_extension` / `single_with_bracketed_extension` for "(N-) N-N", "N-N (-N)", "(N-) N-N (-N)",
  "(N-) N", "N (-N)", where **only the unbracketed core** is written to low and high and the bracketed outer
  bounds stay in the text. A consumer that wants strictly plain cells keeps `range` and `single`.
- Everything else keeps its text only: `not_parsed` (for example carbamazepine "2-8 (4-12)", phenytoin "(5-15)
  10-20"), `not_parsed_several_values` (lists such as "0.32; 0.5"), `not_parsed_other_unit`, `dash` (a printed
  "-"), `cross_reference`, `annotation_only` (mitragynine, toxic: the cell prints only annotation 339).
- Printed line breaks inside a cell are joined with one space; a line that ends in a hyphen with no space after it
  is joined without a space and the hyphen is kept (so "See Bromobenzo-" / "difuranylisopropylamine" reads
  "See Bromobenzo-difuranylisopropylamine"). `note` says when a cell was printed on several lines.
- Five printed rows hold two entries, one per printed line, each with its own values on that line (seen on the
  rendered pages): amitriptyline / "Amitriptyline + Nortriptyline" (row 65), aripiprazole / "Aripiprazole +
  Dehydroaripiprazole" (88), clomipramine / "Clomipramin + N-Desmethylclomipramin" (297), doxepin / "Doxepin plus
  N-desmethyldoxepin" (427), imipramine / "plus Desipramine" (626). Each line is its own CSV row (`source_line` 1
  and 2), which is why the CSV has 1,246 rows for 1,241 printed rows. Cells printed once for the row
  (comatose-fatal, t½, references) are in the row of line 1.
- Two "≥" signs (enoximone, mebendazole) are set in the Symbol font and come out of the PDF as a private-use code
  point; they are written "≥" (checked on the rendered page). The same holds for ±, ≤, →, µ and α in the
  annotations.
- Oddities kept as printed: "0.0015-0-005" (nickel), "363, 364, [1084]]" (fluoxetine references), "[1989, 1090]" (lenalidomide references; the list ends at 1232), "9.1^{8, 396;}"
  (ocfentanil: the semicolon is printed raised), names printed twice in the source ("5-ABP" on pages 1 and 6,
  "5F-AMB (5F-AMB-PINACA)" on pages 4 and 36, "5F-Cumyl-PEGACLONE" on pages 22 and 38).

## Rows and non-missing values

1,241 printed rows; 1,246 CSV rows; 1,243 distinct names; 113 cross-reference rows; 1,133 rows with a value.

| Column | Non-missing | Column | Non-missing |
|---|---|---|---|
| `name` | 1,246 | `name_annotations` | 175 |
| `therapeutic_printed` | 1,184 (1,017 without the 113 "see" and 54 "-") | `therapeutic_annotations` | 288 |
| `therapeutic_low_mg_per_l` | 808 | `therapeutic_high_mg_per_l` | 871 |
| `toxic_printed` | 581 | `toxic_annotations` | 284 |
| `toxic_low_mg_per_l` | 460 | `toxic_high_mg_per_l` | 460 |
| `comatose_fatal_printed` | 433 | `comatose_fatal_annotations` | 286 |
| `comatose_fatal_low_mg_per_l` | 269 | `comatose_fatal_high_mg_per_l` | 269 |
| `half_life_printed` | 1,049 | `half_life_annotations` | 143 |
| `half_life_low_h` | 895 | `half_life_high_h` | 897 |
| `references_printed` | 1,145 | `see` | 113 |
| `remarks` | 716 | `note` | 272 |

Forms of the therapeutic cell: range 504, range_with_bracketed_extension 201, approx_range 68, up_to 48,
approx_single 19, less_than 14, approx_up_to 12, greater_than 10, single 5, at_least 1; cross_reference 113,
not_parsed 99, dash 54, not_parsed_other_unit 28, not_parsed_several_values 8; empty 62.
Forms of the half-life cell: range 652, range_with_bracketed_extension 127, single 57, approx_single 47,
approx_range 9, single_with_bracketed_extension 3, less_than 1, up_to 1; not_parsed_other_unit 80, not_parsed 44,
not_parsed_several_values 23, dash 5; empty 197.
Toxic: single 295, range 131, approx_single 14, range_with_bracketed_extension 10,
single_with_bracketed_extension 8, approx_range 2; not_parsed_several_values 89, not_parsed_other_unit 12, dash 11,
not_parsed 7, annotation_only 1, cross_reference 1; empty 665.
Comatose-fatal: single 173, range 81, range_with_bracketed_extension 7, single_with_bracketed_extension 6,
approx_single 2; not_parsed_several_values 130, not_parsed_other_unit 24, not_parsed 10; empty 813.

143 rows cite reference [4], which the PDF's list gives as "Hiemke C, Bergemann N, Clement HW, [...] Consensus
guidelines for therapeutic drug monitoring in neuropsychopharmacology: update 2017. Pharmacopsychiatry.
2018;51(1/02):9-62."; 169 toxic cells carry annotation 311 ("“laboratory alert level”; mostly according to AGNP
Consensus Guidelines for therapeutic drug monitoring in psychiatry: update 2017 [4] [...]").

## Checks

- `schulz_cell_check.py`: every printed line of every cell (5,809 lines, pdfplumber) occurs character for
  character on the same page of the `pdftotext -layout` text (poppler, an independent reader). Not found: 0.
- `schulz_count_check.py`: the name lines read by the two tools agree on all 82 pages (1,283 lines each, in the
  same order).
- Pages 5, 31 and 54 were rendered and read by eye (the two-entry rows, the "≥" sign, the mitragynine row).
- Spot checks with `schulz_spot_check.py`. "PDF" is the `pdftotext -layout` line of the raw file, in which a raised
  annotation number runs on from the value (so "0.3311" is 0.3 with annotation 311).

| Substance (page) | PDF line | CSV: therapeutic; toxic; comatose-fatal; t½ |
|---|---|---|
| Amitriptyline (5) | `Amitriptyline7, 48  0.05-0.3  0.5-0.6  1.5-2  (10-) 30-50` and `Amitriptyline + Nortriptyline  0.08-0.2  0.3311` | `0.05-0.3` [0.05, 0.3]; `0.5-0.6`; `1.5-2`; `(10-) 30-50` [30, 50]. Line 2: `0.08-0.2`; `0.3^{311}` |
| Clozapine (21) | `Clozapine136  (0.1-) 0.35-0.6  0.6-1 (9.58)  1.28; 28; 5.28  (6-) 12-16` | `(0.1-) 0.35-0.6` [0.35, 0.6]; `0.6-1 (9.5^{8})` text only; `1.2^{8}; 2^{8}; 5.2^{8}` text only; `(6-) 12-16` [12, 16] |
| Diazepam (25) | `Diazepam19  0.1-2 (-2.5)459  3311-5459  24-48` | `0.1-2 (-2.5)^{459}` [0.1, 2]; `3^{311}-5^{459}` [3, 5]; empty; `24-48` [24, 48] |
| Haloperidol (40) | `Haloperidol  (0.001-) 0.005-0.017  0.006/0.023521;  0.188, 74; 0.5  10-35153` and `0.015311; 0.05` | `(0.001-) 0.005-0.017` [0.005, 0.017]; `0.006/0.023^{521}; 0.015^{311}; 0.05` text only; `0.18^{8, 74}; 0.5` text only; `10-35^{153}` [10, 35] |
| Lithium (46) | `Lithium  4-879  8311; 13  16  8-5028` | `4-8^{79}` [4, 8] (annotation 79 gives mmol/L); `8^{311}; 13` text only; `16` [16, 16]; `8-50^{28}` [8, 50] |
| Valproic acid (79) | `Valproic acid  (40-) 50-100 (-150)  120311; 150-200  5568; 7208  7-17 (-20)` | `(40-) 50-100 (-150)` [50, 100]; `120^{311}; 150-200` text only; `556^{8}; 720^{8}` text only; `7-17 (-20)` [7, 17] |
| Phenytoin (63) | `Phenytoin  (5-15) 10-2081  20-25311  388; 438  10-6037` | `(5-15) 10-20^{81}` text only; `20-25^{311}` [20, 25]; `38^{8}; 43^{8}` text only; `10-60^{37}` [10, 60] |
| Fluoxetine (36) | `Fluoxetine  0.12-0.5130  1311  2.2; 3.8507; 68  2-6 days130` | `0.12-0.5^{130}` [0.12, 0.5]; `1^{311}` [1, 1]; `2.2; 3.8^{507}; 6^{8}` text only; `2-6 days^{130}` text only (days) |
| Olanzapine (58) | `Olanzapine548  (0.001-) 0.02-0.08548  0.1311; 0.128  0.258; 18; 2.58; 4.98  (20-) 30-60` | `(0.001-) 0.02-0.08^{548}` [0.02, 0.08]; `0.1^{311}; 0.12^{8}` text only; `0.25^{8}; 1^{8}; 2.5^{8}; 4.9^{8}` text only; `(20-) 30-60` [30, 60] |
| Methadone (50) | `Methadone  (0.05-) 0.1-0.6135  0.3339; 0.6311  0.4 (0.05-1)362  24-48 (-55)` | `(0.05-) 0.1-0.6^{135}` [0.1, 0.6]; `0.3^{339}; 0.6^{311}` text only; `0.4 (0.05-1)^{362}` text only; `24-48 (-55)` [24, 48] |
| Carbamazepine (14) | `Carbamazepine12  2-8 (4-12)  10  20  12-60 (7-35)140` | `2-8 (4-12)` text only; `10` [10, 10]; `20` [20, 20]; `12-60 (7-35)^{140}` text only |
| Quetiapine (67) | `Quetiapine549  0.1-0.5239, 600  1311; 1.8  0.958; 1.48; 1.98;  appr. 5-7240` and `98; 12.78` | `0.1-0.5^{239, 600}` [0.1, 0.5]; `1^{311}; 1.8` text only; `0.95^{8}; 1.4^{8}; 1.9^{8}; 9^{8}; 12.7^{8}` text only; `appr. 5-7^{240}` [5, 7] |
| Morphine (54) | `Morphine288  0.01-0.1  0.1339  0.1339  1-4` | `0.01-0.1` [0.01, 0.1]; `0.1^{339}` [0.1, 0.1]; `0.1^{339}` [0.1, 0.1]; `1-4` [1, 4] |
| Digoxin (26) | `Digoxin  0.0005-0.0009  (0.0012-)  0.004-0.005  40-70` and `0.002383` | `0.0005-0.0009` [0.0005, 0.0009]; `(0.0012-) 0.002^{383}` [0.002, 0.002]; `0.004-0.005` [0.004, 0.005]; `40-70` [40, 70] |
| Paracetamol (60) | `Paracetamol  (5-) 10-25  100-150595  200-300  2-4` | `(5-) 10-25` [10, 25]; `100-150^{595}` [100, 150]; `200-300` [200, 300]; `2-4` [2, 4] |
| Risperidone (69) | `Risperidone48, 550  0.02-0.06272  0.12311  1.88  2-4159` | `0.02-0.06^{272}` [0.02, 0.06]; `0.12^{311}` [0.12, 0.12]; `1.8^{8}` [1.8, 1.8]; `2-4^{159}` [2, 4] |
| Sertraline (71) | `Sertraline423  0.01-0.15 (-0.5)  0.298; 0.3311  1.1-1.5; 1.68; 38  (22-) 24-28 (36)` | `0.01-0.15 (-0.5)` [0.01, 0.15]; `0.29^{8}; 0.3^{311}` text only; `1.1-1.5; 1.6^{8}; 3^{8}` text only; `(22-) 24-28 (36)` text only |
| Zolpidem (82) | `Zolpidem  0.08-0.16629 (-0.2)  0.32311; 0.5  1.5-4  1-5` | `0.08-0.16^{629} (-0.2)` [0.08, 0.16]; `0.32^{311}; 0.5` text only; `1.5-4` [1.5, 4]; `1-5` [1, 5] |
| Lamotrigine (45) | `Lamotrigine  1-15  20  368, 347; 508  (15-) 23-37109` | `1-15` [1, 15]; `20` [20, 20]; `36^{8, 347}; 50^{8}` text only; `(15-) 23-37^{109}` [23, 37] |
| Theophylline (74) | `Theophylline  (5-) 8-15 (-20)82  20  50  6-941` | `(5-) 8-15 (-20)^{82}` [8, 15]; `20` [20, 20]; `50` [50, 50]; `6-9^{41}` [6, 9] |

## What went wrong or was left out

- Nothing was refused and nothing failed for this source.
- The reference list (1,232 entries, pages 119-208) was not turned into a CSV; it is in the PDF and in
  `derived/MOESM1_pdftotext_layout.txt`.
- The bracketed outer bounds, lists of several values and cells with their own unit have no numeric columns by
  design; their text is complete.
- Schulz et al. 2012 (doi:10.1186/cc11441) was not fetched, because the 2020 table was had in full.
