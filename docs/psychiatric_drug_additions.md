# Targets of the eighteen psychiatric drugs absent from the perturbation set

14 October 2026, at the owner's request ("source the missing psychiatric drugs"). Section 5 of
`docs/missing_physiology_review.md` named eighteen psychiatric drugs that the 222-drug perturbation set does not
hold. This is where their molecular targets come from, what each one would enter as and which of them the graph can
hold today.

**Status.** The table is `configs/psychiatric_drug_additions_v4_draft.csv`, 130 rows over 18 drugs, one row per
drug and target gene. It is a draft and enters at the v4 rebuild, which the owner has not ordered. **No table in
use was changed**, the graph was not rebuilt and nothing is staged or committed. A shorter account of the same work
is section 8 of `docs/missing_physiology_review.md`; this document is the full one.

## The schema, and where it was taken from

The columns follow **`data/processed/graph_full_neuronal_split_binders_drugs_v4/drug_entry_nodes.parquet`**, which
is the drug-target table the v4 graph is built from. Its five columns are `perturbation_id`, `drug_node_id`,
`target_node_id`, `sign` and `magnitude`, and the draft carries all five under those names with those conventions:

| convention | in use | in the draft |
|---|---|---|
| perturbation id | `CID1` + an 8-digit PubChem CID, or `RXCUI:<id>` | same |
| drug node id | `DRUG:<perturbation id>` | same |
| target node id | `PROTEIN:<UniProt entry name>` | same |
| sign | `map_drug_targets_to_graph_nodes.ACTION_TYPE_SIGN` and `gtopdb_drug_targets.GTOPDB_ACTION_SIGN` | same |
| magnitude | 1 / number of graph nodes under the target entry (`perturbation_nodes`) | same |

The draft adds eight columns the in-use table does not need, because the in-use table is the output of a build and
the draft is its input: `drug`, `gene`, `role`, `action_type`, `p_affinity`, `p_affinity_type`,
`p_affinity_source`, `target_chembl_id`, `target_is_graph_node`, `source` and `note`.

`p_affinity` is the affinity column name, and it was taken from the two places that already carry it:
`data/processed/exposure/target_occupancy_v4.csv` (written by `experiments/measure_target_occupancy.py`, which also
fixes `role` as `mechanism` or `other`) and `occupancy_calibration.py`, whose
`fit_concentration_calibration(best_mechanism_p_affinity, ...)` reads it. The unit is the one both use: **-log10 of
a molar constant**, so a pKi, a pKd, a pIC50, a pEC50 or a pChEMBL value goes in directly, and
`p_affinity_type` records which it was.

## Three of the eighteen are already in the set, under other names

Checked by PubChem CID and by the OnSIDES ingredient bridge rather than by display name. Three of the eighteen are
perturbations of the v4 graph already, and the names hid it:

| drug | perturbation id | label on the drug node | targets it carries today |
|---|---|---|---|
| vortioxetine | `CID109966051` | `Lu` (from the development code "Lu AA21004") | HTR1A +1, SLC6A4 -1, HTR3A -1 |
| sodium oxybate | `CID100010413` | `hydroxybutyrate` | GABBR1 +1, GABBR2 +1, each at magnitude 0.5 |
| escitalopram | `CID100002771` | `citalopram` | SLC6A4 -1 |

`data/processed/label_review/drug_names.csv` settles each: `CID109966051` has the PubChem title "Vortioxetine",
`CID100010413` has the title "4-Hydroxybutyric Acid" with the RxNorm name "oxybate", and `CID100002771` lists the
RxNorm names "citalopram;escitalopram". The rows above were read back out of
`drug_entry_nodes.parquet` to confirm they are seeded and not merely present.

**Escitalopram is the hardest of the three.** The project's identifier scheme gives the S-enantiomer and the
racemate one perturbation id, so escitalopram cannot enter as a perturbation of its own without splitting that id,
and its one target (SLC6A4) is the target citalopram already carries. So the genuinely absent count is **fifteen**,
not eighteen. `configs/centrally_acting_drugs.csv` separately holds `CID100003519` guanfacine and `CID100003008`
dextromethorphan as candidates, both outside the 222.

## Sources, and the rule each row follows

Everything came from files the repository had already downloaded, plus the session's ChEMBL web service for what
those files do not hold:

| source | what it supplied | licence |
|---|---|---|
| `data/raw/chembl` (ChEMBL_37, 7,561 mechanism rows) | the mechanism target and its action type | already in `docs/data_sources.md` |
| `data/raw/gtopdb_2026` (Guide to PHARMACOLOGY 2026.3) | per-subunit rows, the `Primary Target` flag, affinities | already in `docs/data_sources.md` |
| `data/raw/drugcentral_dti_2021_09_01` (release 2021_09_01) | affinities and the `MOA = 1` mechanism-of-action flag | already in `docs/data_sources.md` |
| `data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json` | the perturbation id and ChEMBL parent of each drug the bridge holds | already in `docs/data_sources.md` |
| the session's ChEMBL MCP service | the gene symbols of `CHEMBL2096662`, which the local `targets.json` does not cache, and the affinities no local file holds | already in `docs/data_sources.md` |

Each drug was resolved to a **ChEMBL parent id confirmed by two independent local sources**: the `ChEMBl ID` column
of the Guide to PHARMACOLOGY `ligand_id_mapping.csv` and the `chembl_parent` field of the OnSIDES bridge. The two
agree wherever both hold the drug. Resolving by guessed salt-form ids was tried and rejected: it pulled in other
molecules' mechanisms, because a salt id shared with another parent carries that parent's rows.

The rules are the ones already registered, not new ones:

- a ChEMBL mechanism row naming one gene is kept, always;
- a ChEMBL mechanism complex row is narrowed by the Guide to PHARMACOLOGY where that source names at least two of
  its subunits, and stays whole otherwise (`gtopdb_drug_targets.replace_complex_rows`, the owner's decision of 13
  October 2026);
- a Guide to PHARMACOLOGY gene belonging to no ChEMBL row of the drug is added;
- `p_affinity` is the maximum over every DrugCentral and Guide to PHARMACOLOGY measurement of that drug and gene,
  under the owner's rule of 9 October 2026, with one documented extension below;
- `magnitude` is 1 / number of graph nodes under the target entry;
- `role` is `mechanism` where ChEMBL files a mechanism row naming the gene, where DrugCentral sets `MOA = 1`, or
  where the Guide to PHARMACOLOGY marks the row `Primary Target`, and `other` otherwise. The draft keeps both so
  the owner can filter. **Only the `mechanism` rows are candidates for a seed.**

Two choices are not measurements and are marked as such in the table:

- **a GABA-A pentamer's EC50 is carried by every subunit of that pentamer.** For zuranolone and brexanolone no
  source holds a per-subunit affinity, and the concentration calibration needs a best-bound target per drug, so the
  pentamer value is the only number there is. `p_affinity_source` names the pentamer, for example
  `ChEMBL:CHEMBL3430899 alpha4/beta3/delta`.
- **where the Guide to PHARMACOLOGY files several rows for one gene with different actions, the row it marks
  `Primary Target` sets the sign.** Guanfacine is why: that source gives ADRA2A both "Agonist" (`Primary Target`
  false) and "Partial agonist" (`Primary Target` true), and a dictionary built row by row would have been decided
  by file order. The draft takes the partial agonism, so ADRA2A enters at +0.5 while ADRA2B and ADRA2C enter at +1.

**The one extension to the affinity rule.** Four values come from the ChEMBL service rather than from DrugCentral
or the Guide to PHARMACOLOGY, because those two hold nothing better for the pair, and each carries its assay id in
`p_affinity_source` so the owner can check it:

| drug and gene | value | assay | why |
|---|---|---|---|
| viloxazine SLC6A2 | pKd 6.81 (Kd 155 nM) | `CHEMBL751379` | DrugCentral holds pIC50 6.70 for the pair; this is a Kd rather than an IC50, so it is the value closer to a dissociation constant, which is the ground the maximum rule rests on |
| zuclopenthixol DRD2 | pAC50 7.35 (AC50 45 nM) | `CHEMBL5291781` | DrugCentral's DRD1 and DRD2 rows for this drug carry no value at all, and the Guide to PHARMACOLOGY holds no human row for it |
| zuranolone alpha1/beta2/gamma2 | pEC50 6.82 (EC50 152 nM) | `CHEMBL5609126` | neither local source holds any zuranolone affinity |
| zuranolone alpha4/beta3/delta | pEC50 6.87 (EC50 135 nM) | `CHEMBL5609171` | the same |

Without this extension viloxazine reads 6.70 rather than 6.81, and zuclopenthixol and zuranolone have no
`p_affinity` at all, so six of the fifteen genuinely absent drugs would be imputed-concentration drugs with nothing
to impute from instead of four. **The extension is visible per row and is the owner's to ratify or refuse.**

## Targets per drug

`m` is rows with role `mechanism`, `o` rows with role `other`. "ChEMBL targets" counts ChEMBL **mechanism target
entries**, which is what the preregistered cap counts: a family or a complex counts as one however many genes it
expands to. "best pK" is the highest `p_affinity` over the mechanism rows, which is the number
`occupancy_calibration.py` would read.

| drug | perturbation id | m | o | ChEMBL targets | mechanism genes | best pK | where the rows came from |
|---|---|---|---|---|---|---|---|
| lithium | `CID100011125` | 3 | 0 | 2 | GSK3A, GSK3B, IMPA1 | none | ChEMBL only |
| escitalopram | `CID100002771` | 1 | 0 | 1 | SLC6A4 | 9.00 | ChEMBL single-gene row |
| vortioxetine | `CID109966051` | 6 | 2 | 3 | HTR1A, HTR1B, HTR1D, HTR3A, HTR7, SLC6A4 | 8.80 | 3 ChEMBL, 3 DrugCentral `MOA = 1`, 2 GtoPdb as `other` |
| topiramate | `CID100005514` | 37 | 3 | 6 | CA2, CA4, 16 GABA-A subunits, GRIA1-4, GRIK1-5, 10 SCN genes | 8.30 | 37 ChEMBL (4 complexes kept whole), 3 GtoPdb as `other` |
| guanfacine | `CID100003519` | 3 | 0 | 1 | ADRA2A, ADRA2B, ADRA2C | 8.70 | GtoPdb narrowing one ChEMBL family |
| viloxazine | `CID100005666` | 1 | 0 | 1 | SLC6A2 | 6.81 | ChEMBL single-gene row; GtoPdb holds no human row |
| levomilnacipran | `CID106917779` | 2 | 1 | 2 | SLC6A2, SLC6A4 | 7.95 | 2 ChEMBL, 1 GtoPdb as `other` |
| zuranolone | `CID186294073` | 16 | 0 | 1 | the 16 GABA-A subunits of `CHEMBL2093872` | 6.87 | ChEMBL complex kept whole |
| brexanolone | `CID100092786` | 6 | 0 | 1 | GABRA1 to GABRA6 | 7.22 | GtoPdb narrowing the same ChEMBL complex |
| tianeptine | `CID100068870` | 0 | 1 | 0 | none; OPRM1 measured only | none | DrugCentral measurement, no action named |
| flupentixol | `CID105281881` | 2 | 2 | 0 | DRD1, DRD2 | 8.82 | DrugCentral `MOA = 1` and GtoPdb `Primary Target`; DRD3 and DRD5 fall to `other` |
| zuclopenthixol | `CID100012453` | 2 | 0 | 0 | DRD1, DRD2 | 7.35 | DrugCentral `MOA = 1` only |
| sertindole | `CID100060149` | 8 | 6 | 3 | ADRA1A, ADRA1B, ADRA1D, DRD2, DRD3, DRD4, HTR2A, HTR2C | 9.70 | 4 ChEMBL, 2 GtoPdb narrowing `CHEMBL2096662`, 2 GtoPdb `Primary Target` |
| solriamfetol | `RXCUI:2121751` | 0 | 3 | 0 | none; SLC6A2, SLC6A3, SLC6A4 measured only | none | GtoPdb and DrugCentral measurements |
| trimipramine | `CID100005584` | 10 | 0 | 6 | ADRA1A, ADRA1B, ADRA1D, DRD2, HRH1, HTR2A, HTR2C, SLC6A2, SLC6A3, SLC6A4 | 8.85 | 8 ChEMBL, 2 DrugCentral `MOA = 1` |
| mianserin | `CID100004184` | 3 | 6 | 0 | HTR2A, HTR2B, HTR2C | 8.75 | DrugCentral `MOA = 1` and GtoPdb `Primary Target`; 6 GtoPdb rows fall to `other` |
| dextromethorphan | `CID100003008` | 2 | 1 | 2 | GRIN2A, SIGMAR1 | 6.28 | 2 ChEMBL, 1 GtoPdb as `other` |
| sodium oxybate | `CID100010413` | 2 | 1 | 1 | GABBR1, GABBR2 | none | ChEMBL complex kept whole, 1 GtoPdb as `other` |

Drugs in the OnSIDES bridge take the bridge's perturbation id, which is how the 222 were named; the seven the
bridge does not hold (viloxazine, levomilnacipran, zuranolone, brexanolone, tianeptine, flupentixol, sertindole)
take `CID1` plus the PubChem CID of their Guide to PHARMACOLOGY ligand record. Where the two sources disagree about
a CID the bridge wins, which is why zuclopenthixol is `CID100012453` and not `CID105311507`, topiramate
`CID100005514` and not `CID105284627`, and dextromethorphan `CID100003008` and not `CID105360696`.

**Five drugs have no ChEMBL mechanism row at all.** `docs/drug_targets_any_type.md` had already recorded two of
them: mianserin and zuclopenthixol are among its 17 drugs with no human mechanism target. Rechecked against the
ChEMBL service on 14 October 2026: mianserin `CHEMBL6437`, zuclopenthixol `CHEMBL53904`, tianeptine
`CHEMBL1289110` and flupentixol `CHEMBL54661` return no mechanism row, and solriamfetol `CHEMBL4297620` returns one
row whose `target_chembl_id` is null and whose `mechanism_of_action` reads "Unknown". Their targets therefore rest
on DrugCentral's mechanism-of-action flag and the Guide to PHARMACOLOGY's primary-target flag, and the draft's
`source` column says so per row.

## Which targets are graph nodes

Checked row by row against the `node_id` column of
`data/processed/graph_full_neuronal_split_binders_drugs_v4/nodes.parquet` (49,796 nodes, of which 12,709 `protein`
and 12,810 `gene`).

**Only two gene symbols in the whole draft are not nodes: GABRG1 and GABRP.** Both are subunits of ChEMBL's
16-subunit GABA-A complex `CHEMBL2093872`, and both are absent as a `protein` node, as a `gene` node and under a
`GENE:` identifier. They affect the two drugs that carry that complex whole:

| drug | mechanism genes | of them nodes | magnitude per node |
|---|---|---|---|
| topiramate | 37 | 35 | 1/14 over the GABA-A complex, 1/4 over AMPA, 1/5 over kainate, 1/10 over the sodium channels, 1 for CA2 and CA4 |
| zuranolone | 16 | 14 | 1/14 |

Neither drug is blocked by it, and the two rows carry an empty `target_node_id` and an empty `magnitude` rather
than a guessed one. **No drug in the draft needs a new node to be added.** Every other target is a `protein` node,
including lithium's GSK3A, GSK3B and IMPA1, tianeptine's OPRM1, sodium oxybate's HCAR1, dextromethorphan's TAS2R1
and mianserin's LPAR1.

**IMPA2 and INPP1 are nodes but are not in the draft.** Section 5 named them among lithium's targets and they are
indeed nodes, but neither carries a row in ChEMBL, in DrugCentral or in the Guide to PHARMACOLOGY. Adding them
would be hand curation from the literature, not a read of these three sources, so they were left out.

## The licence rows for `docs/data_sources.md`

**No new row is needed.** Every source the draft rests on is already in that table with a verified licence. The
quotes, and where each was read:

- **ChEMBL**: "The ChEMBL data is made available on a Creative Commons Attribution-Share Alike 3.0 Unported
  License." Read from `chembl.gitbook.io/chembl-interface-documentation/about`, section "Data Licensing", on 14
  October 2026 with `experiments/read_public_page.py`; the saved page is
  `data/raw/public_pages/697ab689ae23925f.html`, sha256
  `697ab689ae23925f02dc08d645ed6527a08519a3a0f79cccb523645fe92b7f90`. Checked again here against that saved copy,
  whose digest still matches what the table records.
- **Guide to PHARMACOLOGY**: "The Guide to PHARMACOLOGY database is licensed under the Open Data Commons Open
  Database License (ODbL)."; "Its contents are licensed under a Creative Commons Attribution-ShareAlike 4.0
  International License" (`guidetopharmacology.org/about.jsp`, read 11 October 2026).
- **DrugCentral**: "DrugCentral is available under Creative Commons license, download and use of this resource
  evidences your agreement to all the terms and conditions of license." (`drugcentral.org/privacy`, read 11
  October 2026).
- **OnSIDES v3.1.1**: CC BY 4.0, with the note that the MedDRA term names inside are licensed and not
  redistributed. Only the ingredient bridge was read here, which holds identifiers and no MedDRA term.

Nothing is redistributed: the draft holds gene symbols, signs, affinities and the identifiers that point back at
each source.

**One amendment worth making, which is a pin and not a licence.** The ChEMBL-web-service row names
`CHEMBL5291781` and `CHEMBL5609126` as the assays that showed the service's release is not v34. A third assay used
by the draft, **`CHEMBL5609171`** (Eur J Med Chem 2024, zuranolone at the alpha4/beta3/delta pentamer), belongs in
that list for the same reason. The sentence to add to the row's version cell:

> Also `CHEMBL5609171` (Eur J Med Chem 2024), the source of zuranolone's alpha4/beta3/delta pEC50.

## The four hard cases

**Lithium is not a receptor occupancy, and a target row for it does not mean what the other rows mean.**

- *What the row would represent*: ChEMBL's two curated mechanisms, inositol monophosphatase 1 inhibition
  (`CHEMBL1786`, one gene) and glycogen synthase kinase-3 inhibition (`CHEMBL2095188`, a `PROTEIN FAMILY` that
  expands to GSK3A and GSK3B at magnitude 0.5 each), both at sign -1. All three genes are nodes.
- *What it would not represent*: **any affinity**, because no source holds one. Lithium acts at millimolar plasma
  concentrations by competing with magnesium at the enzymes' metal sites. That is an uncompetitive,
  magnesium-dependent inhibition rather than a one-site equilibrium binding, so `occupancy_calibration.py` can
  neither fit nor read it, and lithium would be an imputed-concentration drug with no `p_affinity` to impute from.
  Nor does the row represent the phosphoinositide depletion the mechanism is actually about, since the graph has
  one node per enzyme and no inositol cycle to deplete. Nor IMPA2 and INPP1, which the sources do not file.
- Lithium is also absent for a reason the targets do not fix: the OnSIDES bridge gives lithium carbonate
  `nervous_system_atc: false` and an empty ATC list, so the admission rule has no code to read.

**Sodium oxybate is already in the set, and its high-affinity site is the one that is missing.** It enters as
`CID100010413` on GABBR1 and GABBR2 at sign +1 from ChEMBL's GABA-B complex row, with no affinity. **The GHB
receptor has no node and no row**: no gene symbol for it appears in any of the three sources. The Guide to
PHARMACOLOGY's only affinity for the drug is HCAR1 at pEC50 1.82, a 15-millimolar full agonism of the lactate
receptor, which is not the psychiatric mechanism and is marked `other`; and that source's one GABA-B row names the
auxiliary subunits as well (`GABBR1|GABBR2|KCTD8|KCTD12|KCTD16`) and is dropped as a multi-gene row, which is the
behaviour `one_gene_per_row` is for. So the drug's representation is GABA-B agonism alone, with no number attached.

**Brexanolone and zuranolone are positive allosteric modulators of GABA-A, which the graph represents poorly.**
Both carry sign +0.5 in the draft, from ChEMBL `POSITIVE ALLOSTERIC MODULATOR` and the Guide to PHARMACOLOGY
`Potentiation`, which asserts that the receptor's activity rises. But a modulator with no agonist present does
nothing, and **the graph has no edge from GABA to its receptor to modulate**: section 2 of
`docs/missing_physiology_review.md` counts 4 `inhibits` edges metabolite-to-protein in the whole graph. A +0.5 on a
subunit therefore stands in for a potentiation the graph cannot express. The two also diverge for a reason that is
an artefact rather than pharmacology: **zuranolone keeps all 16 subunits because the Guide to PHARMACOLOGY holds no
zuranolone row at all, while brexanolone narrows to GABRA1 to GABRA6 because that source names exactly those six**
(all with `Primary Target` false, so they enter through ChEMBL's complex row). The narrowing drops every beta,
gamma and delta subunit, and the neurosteroid site sits at the alpha/beta interface, so it removes half a binding
site. That is the error `OBLIGATE_SUBUNITS_OF_A_SITE` exists to prevent, except that `configs/binding_sites.csv`
names no neurosteroid site, so no partner is added back. Two drugs with the same mechanism enter with 6 and 16
targets. Either both should be narrowed, with a neurosteroid-site entry in `configs/binding_sites.csv` giving
GABRB1, GABRB2 and GABRB3 as obligate partners, or neither should be. **This is a decision for the owner, and the
draft records both as they fall out of the registered rule.**

**Dextromethorphan's targets are sourced here, but whether it is admitted is a separate open question.** It
carries GRIN2A at -1 and SIGMAR1 at +1 from ChEMBL, with SIGMAR1's pKi 6.28 the only affinity. The Guide to
PHARMACOLOGY adds GRIN2C at pIC50 6.33 with no `Action`, so it takes no sign, and TAS2R1 as an agonism, a bitter
taste receptor, marked `other`. It is already on the item 31 shortlist for an unrelated reason: **its US label is
an over-the-counter Drug Facts with no pharmacology section**, so the central-action rule of item 14(a) has nothing
to read and item 31(b) sends it to a second medicines agency. The two questions are independent.

## What this does not settle

- **The target counts are far from uniform**, and admission turns on the cap rather than on the gene counts:
  trimipramine 10 genes, sertindole 8, zuranolone 16, topiramate 37 against viloxazine 1. Counting ChEMBL mechanism
  target entries, which is what the cap counts, **4 of the 15 genuinely absent drugs pass a cap of one**
  (guanfacine, viloxazine, zuranolone, brexanolone), 6 fail it with 2 to 6 entries (levomilnacipran,
  dextromethorphan, lithium, sertindole, topiramate, trimipramine) and **5 fail for a different reason, having no
  ChEMBL mechanism target at all** (tianeptine, flupentixol, zuclopenthixol, solriamfetol, mianserin). Admitting
  these drugs is therefore not separable from the target-cap decision `docs/drug_targets_any_type.md` already puts
  to the owner, and for those 5 it is a question about which table a mechanism may be read from.
- **Four drugs have no mechanism affinity at all** (lithium, sodium oxybate, solriamfetol, tianeptine), so each
  would be an imputed-concentration drug with nothing to impute from. The 14 drugs that do have one all fall inside
  the calibration's fitted support of pK 4.27 to 10.40, so none of them would be clamped.
- **Two drugs have no `mechanism` row at all** (tianeptine, solriamfetol), so under the registered rule neither has
  a seed, only measurements.
- **Two places where the draft is narrower than the pharmacology**, both visible in the `source` column:
  viloxazine's modern label describes serotonin modulation at 5-HT2C and 5-HT2B that none of the three sources
  files as a mechanism row, and topiramate's GABA-A row has a DrugCentral WOMBAT-PK entry naming three subunits
  (`GABRA1|GABRG2|GABRB2`) that the registered rule ignores, because narrowing reads the Guide to PHARMACOLOGY
  only.
- **Solriamfetol is the clearest case where the rule and the prose disagree.** Its one ChEMBL row has a null
  target, but the same row's comments name the mechanism: "Known to function as a dopamine and norepinephrine
  reuptake inhibitor and has low affinity for the dopamine transporter (Ki 14.2 uM) and norepinephrine transporter
  (Ki 3.7 uM)." and "Does not inhibit serotonin reuptake (IC50 > 100 uM). No binding to the dopamine, serotonin,
  norepinephrine, GABA, adenosine, histamine, orexin, benzodiazepine, muscarininc acetylcholine or nicotinic
  acetylcholine receptors." The two quoted Ki values are the draft's SLC6A2 pKi 5.43 and SLC6A3 pKi 4.85, so the
  rows are right and only their `role` is wrong: ChEMBL names SLC6A2 and SLC6A3 as the mechanism in prose while
  filing no mechanism target, and the registered rule reads the target field. SLC6A4 is kept as `other` for the
  opposite reason: DrugCentral's drug-label row carries a Ki of 4.088 for it, which the same comment says is not a
  mechanism. **Reading a mechanism out of a comment is a new rule and the draft does not do it**, so solriamfetol
  stands at 0 mechanism rows and 3 measurements.

## What was checked, and how

The table was rebuilt from the sources independently and compared row by row, and the following all hold:

| check | result |
|---|---|
| `target_is_graph_node` recomputed against the v4 `nodes.parquet` | 0 mismatches of 130 |
| `role` recomputed from the three stated criteria | 0 rows differ of 130 |
| `drug_node_id` equals `DRUG:` + `perturbation_id` | holds for all 130 |
| duplicate drug-and-gene pairs | none |
| `magnitude` equals 1 / graph nodes under the target entry | holds; the only empty values are GABRG1 and GABRP, which have no node |
| `sign` against `ACTION_TYPE_SIGN` and `GTOPDB_ACTION_SIGN` | every action type maps to the registered sign; the one empty sign is tianeptine's OPRM1, whose DrugCentral row names no action |
| the four ChEMBL-service affinities | each reproduced from the service with its assay id: `CHEMBL751379` Kd 155 nM, `CHEMBL5291781` AC50 45 nM, `CHEMBL5609126` EC50 152 nM, `CHEMBL5609171` EC50 135 nM |
| the ChEMBL mechanism rows of all 18 drugs | reproduced from the local ChEMBL_37 cache keyed on the confirmed parent, and spot-checked against the service for viloxazine, levomilnacipran and sertindole, which agree |
| the ChEMBL licence sentence | re-read from the saved page, whose sha256 matches what `docs/data_sources.md` records |
