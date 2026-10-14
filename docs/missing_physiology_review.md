# Critical review: physiology that should be represented and is not

14 October 2026, at the owner's request ("a thorough critical review on missing physiology that should be
represented"). It extends `docs/physiological_adherence_review.md`, which asks whether what **is** represented is
faithful; this asks what is **absent**.

A separate pass read the repository without seeing this session's work, and its report is the basis of what follows.
**Every count and quoted line below that is marked (verified) was re-measured in this session against
`data/processed/graph_full_neuronal`, the confirmatory graph.** Findings marked (relayed) are that pass's and were
not independently re-measured; they are stated as its claims, not as facts of the repository, and each should be
checked before it is acted on. Nothing here is a decision.

---

## 1. The chloride charge sum is broken, and it is a defect rather than a scope limit (verified)

**Fixed and queued: `tests/test_chloride_charge_sum.py`, decision 1.2.**

The graph carries chloride under two identifiers. `MAM01442*` is Human-GEM's chloride, CHEBI:17996, SBML charge
-1. `MAM02144*` is display-named **"Cl-"** in three of its four compartments and is in **neither** the Human-GEM
SBML nor its `metabolites.tsv`. `INORGANIC_IONS` keys chloride on `MAM01442` alone
(`experiments/build_neuronal_graph_variant.py` line 43) and `ion_charges` filters SBML charges to those bases, so
`reactome_import.py`'s `inward_charge = sum(ion_charges.get(base, 0) * amount for ...)` ignores chloride carried as
`MAM02144`, and `ions_moved_inward` does not even count it as moved.

Three consequences, each measured:

| reaction | what it does | what the graph says |
| --- | --- | --- |
| `RCT_975340` | GABA_A heteropentamer, Cl- inward | **no electrical edge at all** |
| `RCT_975389` | glycine receptor, Cl- inward | **no electrical edge at all** |
| `RCT_426086` | NKCC1, Na+ + K+ + 2Cl- inward, **electroneutral** | `changes_membrane_potential` **+1.0** |
| `RCT_426155` | KCC2, K+ + Cl- outward, **electroneutral** | `changes_membrane_potential` **-1.0** |

The signs are exactly what is left when chloride is dropped: NKCC1's Na+ and K+ leave +1, KCC2's K+ leaves -1. So
two electroneutral cotransporters are read as **opposing currents**, and the electrical layer can be depolarised by
a ligand-gated receptor and not hyperpolarised by one. No reaction touching `MAM02144*` reaches
`NEURON_MAM01442c`, so neither end of chloride homeostasis is connected.

**Why it matters.** Every benzodiazepine, barbiturate, propofol, etomidate and neurosteroid in the drug set acts
through GABA_A, and none of them has an electrical route. It bears on somnolence_or_hypersomnia,
cognitive_impairment, disinhibition_or_impulsivity, psychomotor_retardation, catatonia and anxiety.

**The fix is a mapping change and needs the graph rebuilt**, so it is queued for decision 1.2 rather than applied.
The tests assert the defect as it stands and each names what it should say afterwards, so whoever rebuilds can
invert them and watch them pass. **This is the highest-value item in this review**: it is cheap, checkable and
affects a large, named group of drugs.

---

## 2. The metabolic layer can receive but almost cannot act (verified)

Of **266,966** edges, **600** run from a metabolite to a gene or protein entity: 436 `binds`, 83 `cofactor_of`, 41
`oxidatively_modifies`, 36 `activates` (the curated presynaptic ligands) and **4 `inhibits` in the whole graph**.
**93 of 8,580 metabolites** carry any such edge. The reverse direction is 23,901 `catalyzed_by` edges plus the
complex route, so the two layers are joined at a ratio of about one to forty.

Specific absences, each checked: cortisol -> `NR3C1` **0**, triiodothyronine -> `THRA` **0**, any metabolite ->
`PPARG` **0**. `NR3C1` is in the graph with 291 `regulates_transcription_of` edges and is reachable only from other
proteins.

**Why it matters more than the ratio suggests.** The study's motivating class is the monogenic metabolic disorders,
and the mechanism in every one of them is a metabolite acting on something: ammonia in the urea-cycle defects,
homocysteine at NMDA receptors, quinolinate at NMDA, copper in Wilson disease, phenylalanine competing for
transport. The model can move those metabolite levels and then very little happens. The relayed pass reports that
quinolinate, L-kynurenine, homocysteine and ammonia have **no** protein edge at all; I verified the pattern on the
nuclear receptors and did not check those four individually.

**Addable, partly, with sources the project already uses**: Guide to PHARMACOLOGY and ChEMBL carry
endogenous-ligand-to-nuclear-receptor pairs, and the curated presynaptic table is the template for the relation. A
table of order 100 pairs would cover the hormones and the neuroactive metabolites in the 24-symptom set.

---

## 3. One membrane potential, and no mitochondrial one (verified)

`node_type == "membrane_potential"` has exactly **1** node, `VM_c`, the plasma membrane. There is no inner
mitochondrial membrane potential, so the proton-motive force is not a state variable. The project has built pool
nodes for precisely this kind of problem elsewhere - `NEURON_*` for four ions, `REDOX_NADPH_*` and `REDOX_GSH_*` -
so the pattern exists and has been applied to the antioxidant arm and not the energy arm. The relayed pass adds
that complex I's only propagating output is the ubiquinone pool, because H+, NADH, NAD+, ATP, ADP and Pi are all
currency metabolites that receive and do not send; I did not re-measure that.

**Why it matters.** fatigue is the second-largest symptom in the label set, and complex-I inhibition is the
canonical molecular route to parkinsonism.

---

## 4. The eleven relayed findings, now verified (14 October 2026)

A second delegated pass checked all eleven against the repository; the most consequential corrections were
re-measured again in this session. **Three were refuted or materially corrected, and one of the refutations would
have made a disclosure sentence false.**

**A structural fact that conditions three claims.** `graph_full_neuronal`, the confirmatory graph, has **no drug
nodes and no `targets`, `clears` or `sequesters` relation**; the drug layer lives in
`graph_full_neuronal_split_binders_drugs_v4_enzymes`. So the subunit-composition and pharmacokinetic claims, and
the drug-target half of the G-protein claim, are statements about the v4 drug graph only and must be re-measured if
the drug layer is rebuilt. The non-drug subgraphs are identical between the two, which was checked rather than
assumed.

### Corrected, and the correction matters

- **No closed adaptation loop — the defect is about 20 times larger than reported, and the loop is closable.** The
  claim said 7 of 22 GRK/arrestin-to-receptor edges carry the wrong sign. Re-measured here: **94 of 116**
  GRK/arrestin edges into GPCR-family genes read `activates` at **+1** (13 `inhibits`, 9 `binds`); the delegated
  pass counts 147 of 196 on a wider receptor set. Either way a desensitising kinase is recorded as **activating**
  its receptor across most of the population, which is the single largest data defect any of the eleven points at.
  And five receptor-to-arrestin `activates` edges do exist (CCR7, PTH1R, CHRNA7, GPR182, LHCGR), so the loop is
  closable in the graph - just not for any psychiatric receptor. The honest disclosure is "no closed adaptation
  loop **for any psychiatric receptor**", which is a different statement and a different fix.
- **Conductance: the premises hold and the conclusion is false, so the drafted disclosure would have been wrong.**
  `LEAK_CONDUCTANCE_*` really is a module constant and there really are only 4 `driving_force` edges. But "the only
  route from a channel to V_m is the slow concentration route" is false: **405 `changes_membrane_potential` edges
  run reaction to `VM_c` directly**, 353 at +1 and 52 at -1 (verified here). The real problem is already recorded
  in `membrane_potential.py`: the in-degree divisor means each of the 405 arrives divided by 405, with a measured
  median field of 2.9e-15. That is a normalisation defect, not a scope limit, and disclosing it as a scope limit
  would have concealed something fixable.
- **GABA_A subunit composition: about four times smaller than reported.** `RCTE_975268` has **13 members, not 4**,
  and there are at least six GABA_A-related entities rather than one, so **2 of the 13 drug-target subunits reach
  no channel (GABRD, GABRE), not 9**. The parts worth keeping are exact: 164 of 222 drugs in 36 identical-input
  groups, and cocaine's only graph input being ten sodium-channel subtypes with no monoamine transporter at all.
- **Astrocytic K+ buffering does have a route.** KCNJ10 catalyses no reaction directly, but it is `member_of`
  `RCTE_1299206` (Kir 4.1/5.1), which catalyses `RCT_1296045`, which moves `NEURON_MAM02200c` to
  `NEURON_MAM02200e` and carries `changes_membrane_potential` **-1.0**. The true limitation is that pool nodes
  carry cell-class weight 1.0 in every class, so the route is not astrocyte-specific.
- **Second messengers**: every number exact, but cAMP reaches the PKA holoenzyme in 2 hops through a Reactome
  binding reaction, and the neuronal calcium pool has 14 signed edges into the transmitter-release reactions. So
  "bypassed" should read "bypassed on the direct metabolite-to-protein relation".
- **G-protein coupling**: every graph number right (the Gq count is 24 and 13, not 23 and 12). But the claim that a
  per-receptor polarity is impossible in the encoder is **wrong**: sign and weight are both per-edge and decision
  4.11 now propagates the weight, so the obstacle is that the graph writes +1 on all 224 receptor-to-G-alpha edges.
  **The defect is in the data, not the model**, which makes it cheaper to fix than the review first said.

### Verified clean, no correction worth making

No temporal field (28 columns, and `chronic`, `weeks`, `onset`, `latency`, `tolerance`, `desensit`, `downregulat`
all zero in the preregistration); pharmacokinetics (1,182 drug edges = 1,058 `targets` + 101 `sequesters` + 23
`clears`, four enzymes only, both prodrugs edge-free, ABCB1/ABCG2/CYP3A4 all present with zero drug edges); the
endocrine genes (24 of 24 present, no organ column, DRD2 inhibits PRL at -1); Complex I (its only propagating
output is the ubiquinone pool, because an edge leaving a currency node is dropped at encoder lines 771-774); and no
regions in the propagation, which the repository already says in its own words at `node_descriptors.py` lines 55-57.

## 5. Not physiology but found while checking: psychiatric drugs absent from the set

The 222-drug set has good modern coverage - cariprazine, brexpiprazole, lumateperone, pimavanserin, iloperidone,
asenapine, lurasidone, esketamine, the three orexin antagonists - and these are absent:

**lithium**, escitalopram, vortioxetine, topiramate, guanfacine, viloxazine, levomilnacipran, zuranolone,
brexanolone, tianeptine, flupentixol, zuclopenthixol, sertindole, solriamfetol, trimipramine, mianserin,
dextromethorphan, sodium oxybate.

**Every one of the five checked has its canonical targets already in the graph**: lithium's GSK3B, IMPA1, IMPA2 and
INPP1; topiramate's GRIK1, SCN1A, CA2 and GABRA1; vortioxetine's SLC6A4, HTR1A, HTR3A and HTR7; escitalopram's
SLC6A4; guanfacine's ADRA2A. So they are absent by **selection**, not because the graph cannot hold them.
**Lithium is the conspicuous one**: it is the archetypal mood stabiliser, its symptom profile is as well documented
as any drug in the set, and a study predicting psychiatric symptoms that contains no lithium will be asked about it.

## 6. Superseded: the original relayed list

Stated as the delegated pass's claims. Each is specific enough to check, and should be checked before it is acted
on.

- **G-protein coupling is promiscuous.** 12 of the 42 drug-target GPCRs are said to `activate` both a Gs-family
  and a Gi-family subunit at +1, with HTR1A activating GNAS; and because the encoder holds one gain per relation,
  the Gs/Gi axis is claimed to rest on a single global ratio with a recoverable window of (1.447, 4.32). If that is
  right it voids any cAMP-mediated mechanism card, and it is fixable from Guide to PHARMACOLOGY's transduction
  field, which is already a project source.
- **Conductance is not a state variable**, only ion concentration is, so a synaptic conductance change - the fast
  route - has no representation and only the slow concentration route exists. This one I think is right and is a
  genuine scope limit rather than a defect: a product of two state variables is outside the registered linearity.
- **Second messengers are bypassed.** cAMP, cGMP, IP3 and the calcium pools exist but are claimed to send almost
  nothing (cAMP 1 outgoing protein edge, cGMP 0), while OmniPath gives a 2-hop protein route where the
  physiological one is 10 hops.
- **No closed adaptation loop**, although GRKs and arrestins are nodes with 22 edges back to receptors and no
  edges forward to them - and 7 of those 22 are said to carry sign +1, the wrong direction for a desensitising
  kinase. Tolerance, antidepressant latency and D2 upregulation all live here.
- **Receptor subunit composition is absent**: 164 of 222 drugs are claimed to fall into 36 groups sharing an
  identical input, so three quarters of the drug perturbations are 36 distinct inputs rather than 222. The GABA_A
  group is 22 drugs hitting 13 subunits uniformly, where the only pentamer entity has 4 members.
- **No circuits and no regions in the propagation**, so one D2 input must explain both psychosis and parkinsonism.
  This is an in-principle limit, not a quantitative one, and should be stated next to the equifinality endpoint.
- **No temporal field anywhere.** The evidence table has no dose, duration, time-to-onset or acute/chronic column,
  so a steady-state model is read against chronic-dose labels and the mismatch cannot even be stratified.
- **Pharmacokinetics**: no active metabolites as entities, four clearing enzymes only, no efflux transporter edges;
  and tramadol and dihydrocodeine are prodrugs whose `clears` edge, if added at the registered -1, would invert.
- **Glia are expression channels, not a compartment**, so astrocytic K+ buffering has no reaction-level route.
- **Endocrine axes cannot close** without organ labels, which the graph does not have.

---

## What I would put to the owner

1. **The chloride fix, at the rebuild.** Cheap, checkable, and it affects every GABAergic drug in the set.
2. **A metabolite-to-protein table**, order 100 curated pairs, because without it the study's founding examples
   have no mechanism inside the model. This is the finding I would most want checked and then acted on.
3. **Verify the G-protein coupling claim**, because if it holds it is the difference between a mechanism card that
   means something and one that does not, and the fix is a file the project already downloads.
4. **Four disclosure sentences**, which cost nothing and are currently missing: conductance is not a state
   variable; there are no circuits, so a one-target perturbation with two opposite symptoms cannot have their
   relative strength predicted even in principle; the labels are not steady-state observations; no endocrine loop
   is representable.


## 7. Correction: BRENDA is CC BY 4.0, so the allosteric gap is a download rather than a curation

This session told the owner that BRENDA "restricts bulk download". **That was asserted from memory and it is wrong.**
`brenda-enzymes.org/download.php`, read 14 October 2026 with `experiments/read_public_page.py`: "All copyrightable
parts of BRENDA are licensed under **Creative Commons Attribution License 4.0 (CC BY 4.0)**", and the whole database
is published as a 71.59 MB text tarball and a 79.35 MB JSON tarball. The only gate is "Before downloading any files
from BRENDA, you have to actively accept the license", a checkbox the owner clicks.

So there is no prohibition to work around, and the enzyme-regulation gap - 173 regulatory edges against 60,838 into
reactions - is addressable by downloading a released file and mapping its inhibitor and activator records onto
Human-GEM, not by hand curation. CC BY 4.0 also permits redistribution with attribution, so a derived table may
ship in this repository, which the licence of a curated-by-hand table would not have settled.

**This is what the licence-verification ruling of 11 October 2026 exists to catch** ("No. You have to verify the
license"), and it caught a claim that would otherwise have sent the owner into unnecessary curation. SABIO-RK, which
this review proposed as the free alternative, is moot; its terms page served no robots file and was recorded as a
refusal rather than read.

## 8. Targets sourced for the eighteen absent psychiatric drugs

14 October 2026, at the owner's request ("source the missing psychiatric drugs"). The draft table is
`configs/psychiatric_drug_additions_v4_draft.csv`, 130 rows over 18 drugs, one row per drug and target gene. It is a
draft and enters at the v4 rebuild, which the owner has not ordered. **No table in use was changed**; the only other
file touched is `docs/data_sources.md`, where the ChEMBL licence now carries its quoted sentence.

### Three of the eighteen are already in the 222-drug set

Section 5 named eighteen drugs as absent. Checked by PubChem CID and by the OnSIDES ingredient bridge rather than by
display name, **three of them are in the set already** and the names hid it:

| drug | perturbation id in the set | label on the drug node | targets it carries today |
| --- | --- | --- | --- |
| vortioxetine | `CID109966051` | `Lu` (from "Lu AA21004") | HTR1A, SLC6A4, HTR3A |
| sodium oxybate | `CID100010413` | `hydroxybutyrate` | GABBR1, GABBR2 |
| escitalopram | `CID100002771` | `citalopram` | SLC6A4 |

`data/processed/label_review/drug_names.csv` carries the RxNorm names that settle each: `CID109966051` is
"Vortioxetine", `CID100010413` is "4-Hydroxybutyric Acid" with the RxNorm name "oxybate", and `CID100002771` lists
"citalopram;escitalopram". **Escitalopram is the hardest of the three**: the project's identifier scheme gives the
enantiomer and the racemate one perturbation id, so escitalopram cannot enter as a perturbation of its own without
splitting that id, and its one target (SLC6A4) is the target citalopram already carries. The genuinely absent count
is fifteen, not eighteen. `configs/centrally_acting_drugs.csv` also already holds `CID100003519` guanfacine and
`CID100003008` dextromethorphan as candidates, both outside the 222.

### Sources, and the rule each row follows

Everything came from files the repository had already downloaded, plus the session's ChEMBL web service for what
those files do not hold:

- `data/raw/chembl` (ChEMBL_37, 7,561 mechanism rows) for the mechanism target and its action;
- `data/raw/gtopdb_2026` (Guide to PHARMACOLOGY 2026.3) for per-subunit rows and affinities;
- `data/raw/drugcentral_dti_2021_09_01` for affinities and the mechanism-of-action flag;
- the ChEMBL MCP service for five drugs the local dump leaves without a mechanism row, for the gene symbols of
  `CHEMBL2096662` (a SELECTIVITY GROUP the local `targets.json` does not cache: HTR2A and HTR2C), and for the
  pentamer EC50 values of zuranolone and brexanolone.

The rules are the ones already registered, not new ones:

- a ChEMBL mechanism row naming one gene is kept;
- a ChEMBL mechanism complex row is narrowed by GtoPdb where GtoPdb names at least two of its subunits, and stays
  whole otherwise (`gtopdb_drug_targets.replace_complex_rows`, the owner's decision of 13 October 2026);
- `p_affinity` is the maximum over every DrugCentral and GtoPdb measurement of that drug and gene, as -log10 of a
  molar constant, under the owner's rule of 9 October 2026;
- `magnitude` is 1 / number of graph nodes under the target entry (`perturbation_nodes`);
- `role` is `mechanism` where ChEMBL files a mechanism row, DrugCentral sets `MOA = 1` or GtoPdb marks `Primary
  Target`, and `other` otherwise. The draft keeps both so the owner can filter; only the `mechanism` rows are
  candidates for a seed.

One choice is not a measurement and is marked as such in the table: **a GABA-A pentamer's EC50 is carried by every
subunit of that pentamer**. For zuranolone and brexanolone no source holds a per-subunit affinity, and the
concentration calibration (`occupancy_calibration.py`) needs a best-bound target per drug, so the pentamer value is
the only number there is. A second choice: where GtoPdb files several rows for one gene with different actions, the
row GtoPdb marks `Primary Target` sets the sign. Guanfacine is why - GtoPdb gives ADRA2A both "Agonist" and "Partial
agonist", and a dictionary built row by row would have been decided by file order.

### Targets per drug

`m` is mechanism rows, `o` rows with role `other`. "best pK" is the highest `p_affinity` over the mechanism rows,
which is the number the concentration calibration would read.

| drug | perturbation id | m | o | mechanism targets | best pK |
| --- | --- | --- | --- | --- | --- |
| lithium | `CID100011125` | 3 | 0 | GSK3A, GSK3B, IMPA1 | none |
| escitalopram | `CID100002771` | 1 | 0 | SLC6A4 | 9.00 |
| vortioxetine | `CID109966051` | 6 | 2 | HTR1A, HTR1B, HTR1D, HTR3A, HTR7, SLC6A4 | 8.80 |
| topiramate | `CID100005514` | 37 | 3 | CA2, CA4, 16 GABA-A subunits, GRIA1-4, GRIK1-5, 10 SCN genes | 8.30 |
| guanfacine | `CID100003519` | 3 | 0 | ADRA2A, ADRA2B, ADRA2C | 8.70 |
| viloxazine | `CID100005666` | 1 | 0 | SLC6A2 | 6.81 |
| levomilnacipran | `CID106917779` | 2 | 1 | SLC6A2, SLC6A4 | 7.95 |
| zuranolone | `CID186294073` | 16 | 0 | the 16 GABA-A subunits of `CHEMBL2093872` | 6.87 |
| brexanolone | `CID100092786` | 6 | 0 | GABRA1 to GABRA6 | 7.22 |
| tianeptine | `CID100068870` | 0 | 1 | none; OPRM1 measured only | none |
| flupentixol | `CID105281881` | 2 | 2 | DRD1, DRD2 | 8.82 |
| zuclopenthixol | `CID100012453` | 2 | 0 | DRD1, DRD2 | 7.35 |
| sertindole | `CID100060149` | 8 | 6 | ADRA1A, ADRA1B, ADRA1D, DRD2, DRD3, DRD4, HTR2A, HTR2C | 9.70 |
| solriamfetol | `RXCUI:2121751` | 0 | 3 | none; SLC6A2, SLC6A3, SLC6A4 measured only | none |
| trimipramine | `CID100005584` | 10 | 0 | ADRA1A, ADRA1B, ADRA1D, DRD2, HRH1, HTR2A, HTR2C, SLC6A2, SLC6A3, SLC6A4 | 8.85 |
| mianserin | `CID100004184` | 3 | 6 | HTR2A, HTR2B, HTR2C | 8.75 |
| dextromethorphan | `CID100003008` | 2 | 1 | GRIN2A, SIGMAR1 | 6.28 |
| sodium oxybate | `CID100010413` | 2 | 1 | GABBR1, GABBR2 | none |

Five drugs have **no ChEMBL mechanism row at all**, which `docs/drug_targets_any_type.md` had already recorded for
two of them (mianserin and zuclopenthixol are among its 17 with no human mechanism target). Checked again through
the ChEMBL service on 14 October 2026: mianserin `CHEMBL6437`, zuclopenthixol `CHEMBL53904`, tianeptine
`CHEMBL1289110` and flupentixol `CHEMBL54661` return `{"count":0}`; solriamfetol `CHEMBL4297620` returns one row
whose `target_chembl_id` is null and whose `mechanism_of_action` reads "Unknown". Their targets therefore come from
DrugCentral's mechanism-of-action flag and GtoPdb's primary-target flag, and the draft says so per row.

### Which targets are graph nodes

**Only two gene symbols in the whole draft are not nodes of
`data/processed/graph_full_neuronal_split_binders_drugs_v4`: GABRG1 and GABRP.** Both are subunits of ChEMBL's
16-subunit GABA-A complex `CHEMBL2093872`, and both are absent as a `protein` node and as a `gene` node. They affect
topiramate and zuranolone, which carry that complex whole. Neither drug is blocked by it: topiramate keeps 35 of its
37 mechanism genes and zuranolone 14 of its 16, so each still has graph nodes and a `magnitude` of 1/14 over the
complex rather than 1/16. **No drug in the draft needs a new node to be added.** Every other target - including
lithium's GSK3A, GSK3B and IMPA1, tianeptine's OPRM1, sodium oxybate's HCAR1, dextromethorphan's TAS2R1 and
mianserin's LPAR1 - is a `protein` node.

IMPA2 and INPP1, which section 5 named among lithium's targets and which are indeed nodes, carry **no row in
ChEMBL, DrugCentral or the Guide to PHARMACOLOGY**, so they are not in the draft. Adding them would be hand
curation from the literature, not a read of these three sources.

### The four hard cases

- **Lithium** is not a receptor occupancy and a target row for it does not mean what the other rows mean. What the
  row would represent: ChEMBL's two curated mechanisms, inositol monophosphatase 1 inhibition (`CHEMBL1786`) and
  glycogen synthase kinase-3 inhibition (`CHEMBL2095188`, which expands to GSK3A and GSK3B), both at sign -1. What
  it would not represent: any affinity, because **no source holds one** - lithium acts at millimolar plasma
  concentrations by competing with magnesium at the enzymes' metal sites, which is an uncompetitive,
  magnesium-dependent inhibition and not a one-site equilibrium Kd, so `occupancy_calibration.py` can neither fit
  nor read it and lithium would be an imputed-concentration drug with no `p_affinity` to impute from; the
  phosphoinositide depletion and the downstream effects that the mechanism is actually about, since the graph has
  one node per enzyme and no inositol cycle to deplete; and IMPA2 and INPP1, which the sources do not file.
- **Sodium oxybate** is already in the set on GABBR1 and GABBR2 at sign +1 from ChEMBL's GABA-B complex row, with no
  affinity; **the GHB receptor has no node and no row** - no gene symbol for it appears in any of the three sources,
  GtoPdb's only affinity for the drug is HCAR1 at pEC50 1.82 (a 15-millimolar full agonism, which is the lactate
  receptor and not the psychiatric mechanism), and GtoPdb's one GABA-B row names the KCTD auxiliary subunits
  (`GABBR1|GABBR2|KCTD8|KCTD12|KCTD16`) and is dropped as a multi-gene row. So the drug's representation is GABA-B
  agonism alone, and the high-affinity GHB site is simply absent.
- **Brexanolone and zuranolone** are positive allosteric modulators of GABA-A, which the graph represents poorly. In
  the draft both carry sign +0.5 (`POSITIVE ALLOSTERIC MODULATOR` and GtoPdb `Potentiation`), which asserts that the
  receptor's activity rises - but a modulator with no agonist present does nothing, and the graph has no edge from
  GABA to its receptor to modulate: section 2 of this review records only 4 `inhibits` edges metabolite-to-protein
  in the whole graph. The two also diverge in the draft for a reason that is an artefact, not pharmacology:
  **zuranolone keeps all 16 subunits because GtoPdb holds no zuranolone row at all, while brexanolone narrows to
  GABRA1 to GABRA6 because GtoPdb names exactly those six**. That narrowing drops every beta, gamma and delta
  subunit, and the neurosteroid site sits at the alpha/beta interface, so it removes half a binding site - the error
  `OBLIGATE_SUBUNITS_OF_A_SITE` exists to prevent, except that `configs/binding_sites.csv` names no neurosteroid
  site and so no partner is added back. Two drugs with the same mechanism therefore enter with 6 and 16 targets.
  Either both should be narrowed with a neurosteroid-site entry in `configs/binding_sites.csv` giving GABRB1,
  GABRB2 and GABRB3 as obligate partners, or neither should be narrowed. **This is a decision for the owner and the
  draft records both as they fall out of the registered rule.**
- **Dextromethorphan** carries GRIN2A at -1 and SIGMAR1 at +1 from ChEMBL, with SIGMAR1's pKi 6.28 the only
  affinity; GtoPdb adds GRIN2C at pIC50 6.33 with no Action, so it takes no sign, and TAS2R1 as an agonism, which is
  a bitter taste receptor and is marked `other` rather than mechanism. It is already on the item 31 shortlist for an
  unrelated reason - **its US label is an over-the-counter Drug Facts with no pharmacology section**, so the
  central-action rule of item 14(a) has nothing to read and item 31(b) sends it to a second medicines agency. Its
  targets are therefore sourced here, but whether it is admitted at all is still open under item 31, and the two
  questions are independent.

### What this does not settle

- The targets per drug are far from uniform: trimipramine 10 genes, sertindole 8, zuranolone 16, topiramate 37
  against viloxazine 1. The preregistered cap counts **ChEMBL mechanism targets**, not genes, and a family or a
  complex counts as one, so the figures that decide admission are these: **4 of the 15 genuinely absent drugs pass
  the cap of one** (guanfacine 1 target, viloxazine 1, zuranolone 1, brexanolone 1), 6 fail it with 2 to 6 targets
  (levomilnacipran 2, dextromethorphan 2, lithium 2, sertindole 3, topiramate 6, trimipramine 6), and **5 fail for a
  different reason - they have no ChEMBL mechanism target at all** (tianeptine, flupentixol, zuclopenthixol,
  solriamfetol, mianserin), so no cap admits them and their rows rest on DrugCentral and GtoPdb instead. Admitting
  these drugs is therefore not separable from the target-cap decision that `docs/drug_targets_any_type.md` already
  puts to the owner, and for those 5 it is a question about which table a mechanism may be read from.
- Four drugs have **no mechanism affinity at all** (lithium, sodium oxybate, solriamfetol and tianeptine), so each
  would be an imputed-concentration drug with nothing to impute from.
- Viloxazine's modern label describes serotonin modulation (5-HT2C and 5-HT2B) that none of the three sources files
  as a mechanism row, and topiramate's GABA-A row has a DrugCentral WOMBAT-PK entry naming three subunits
  (`GABRA1|GABRG2|GABRB2`) that the registered rule ignores because narrowing reads GtoPdb only. Both are places
  where the draft is narrower than the pharmacology, and both are visible in the table's `source` column.
