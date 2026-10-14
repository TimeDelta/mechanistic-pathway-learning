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
