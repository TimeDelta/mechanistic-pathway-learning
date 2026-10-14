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

## 4. Findings relayed and not re-measured

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
