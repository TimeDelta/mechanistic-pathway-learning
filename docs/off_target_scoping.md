# Off-target binding: scoping (8 October 2026)

Scoping only. Nothing is trained, and no evidence table, graph or configuration is changed. The numbers come from
`experiments/scope_off_target_binding.py` (summary in `data/processed/off_target_scoping/summary.json`) and
`experiments/scope_open_targets_symptom_labels.py` (`open_targets_symptom_labels.json` beside it). Both read the
cap-3 evidence (`evidence_full_v3_max3_targets`, 222 drugs), the selection `better_v2_full_v3_max3_targets` and
`graph_full_neuronal`. The downloads sit in their own directories under `data/raw`, which git ignores.

## Why

ChEMBL lists a drug's therapeutic mechanisms, not the other proteins it binds (docs/drug_targets_any_type.md).
Amitriptyline's seeds are its two transporters, so its antihistamine and antimuscarinic effects have no input. The
safety-pharmacology literature treats binding at a fixed panel of off-targets as the main predictor of adverse
effects (Bowes et al. 2012; Brennan et al. 2024; Lynch et al. 2017).

## Sources

| source | release | what it gives |
|---|---|---|
| DrugCentral `drug.target.interaction.tsv.gz` (Avram et al. 2023) | 2021_09_01, the latest flat file; newer releases ship only as a PostgreSQL dump | one row per drug, target and measurement; `ACT_VALUE` as −log10 molar; `ACTION_TYPE` mostly on mechanism rows; `MOA` = 1 for the mechanism of action; rows from ChEMBL, DrugMatrix, PDSP, IUPHAR, WOMBAT-PK and drug labels |
| IUPHAR/BPS Guide to PHARMACOLOGY `interactions.csv`, `ligand_id_mapping.csv` (Harding et al. 2024) | 2026.3 (16 September 2026) | affinities with the ligand's type and action: full, partial and inverse agonist, antagonist, channel blocker, allosteric modulator and others |

Drugs are matched by the PubChem CID in the perturbation identifier, the normalised name and the INN: 216 of 222
match. Benztropine, deprenyl, dothiepin, esketamine, eslicarbazepine and propoxyphene do not.

## Rules

**One value per drug-target pair.** The affinity of drug d at gene g is the maximum of every measurement of d at g
in both sources (the user's decision of 9 October 2026, "Oh. In that case, I guess Max is fine."; the median until
then) (types Ki, Kd, IC50, EC50, Kb, A2 and AC50, all as −log10 molar). A value is never borrowed from
another ligand: noradrenaline at the noradrenaline transporter and amphetamine at the same transporter are two
different pairs. The value belongs to the pair, so it enters as the magnitude of that drug's seed at that gene, not
as a node descriptor (a node descriptor would give every drug the same value at the gene).

Why the maximum rather than the median. A pair with more than one value does not hold one assay run twice: DrugCentral
keeps one value per pair 96.5 percent of the time and GtoPdb already medians its own papers, so the repeats are two
curators' aggregates of different literature. 14.7 percent of the 26,597 pairs with an affinity have more than one
value (median 2, maximum 20); of the 2,367 pairs in both sources the median absolute difference is 0.10 log units, the
90th percentile 1.13, and 12.1 percent differ by more than 1 log unit, with GtoPdb higher in 35.7 percent. With two
values the median is the midpoint, so the maximum sits half a gap above it (median half-gap 0.05, 90th percentile 0.57,
which moves an off-target's occupancy at the 0.25 threshold from 0.25 to between about 0.35 and 0.70).

Pooling Ki with IC50 and EC50 is a simplification. IC50 depends on the substrate or radioligand concentration in the
assay, and EC50 measures potency, not binding (Hulme and Trevethick 2010). By Cheng-Prusoff an IC50 is at or above the
Ki, so an IC50-derived value is biased downward and the maximum is the one closer to a dissociation constant, which is
the pharmacological argument for the decision; 34.3 percent of DrugCentral pairs carry no Ki or Kd. The untried
alternative, recorded as a note: the type-aware rule, the median over Ki and Kd where they exist and the maximum over
IC50 and EC50 otherwise.

**Occupancy.** At equilibrium and one binding site, the fraction of target bound at free drug concentration C is
θ = C / (C + K), where K is the equilibrium dissociation constant (Hulme and Trevethick 2010). Affinity is a property
of the pair; occupancy depends on the pair and on the concentration the drug reaches. What matters in vivo is the
free concentration: "in vivo efficacy is determined by the free (unbound) drug concentration surrounding the
therapeutic target, not by the free drug fraction" (Smith et al. 2010). The occupancy formula takes C as an input
and does not compute it.

The scoping uses a placeholder for C: at the therapeutic dose the drug's best-bound mechanism target is 90 percent
occupied, so C = 9 K_primary and θ_target = 1 / (1 + 10^(pK_primary − pK_target) / 9). The primary target is the
mechanism gene (ChEMBL seed) with the highest median affinity; 183 drugs have one, and the other 33 (no measured
affinity at any mechanism gene) get no off-targets under this rule. A threshold on θ then reads as a maximum
affinity gap to the primary target:

| occupancy threshold | off-target may bind this much more weakly than the primary target |
|---|---|
| 0.5 | 0.95 log units (9-fold) |
| 0.25 | 1.43 log units (27-fold) |
| 0.1 | 1.91 log units (81-fold) |

A measured C would replace the placeholder: free C = therapeutic total plasma concentration (Schulz et al. 2020,
more than 1,100 drugs) × fraction unbound in plasma (Lombardo et al. 2018, 1,352 drugs). That route needs a molar
conversion and a match for each drug; it is not built.

**Sign.** Each pair's action strings from both sources map to −1 (inhibitor, antagonist, blocker, inverse agonist,
negative modulator, gating inhibitor), +1 (agonist, activator, opener, positive modulator, potentiator, releasing
agent) or none (no action recorded, or actions of both signs). Partial agonism is recorded by GtoPdb and flagged
apart: its net sign depends on endogenous tone (an antagonist where the endogenous agonist is high, an agonist where
it is low).

## Results

| | occupancy ≥ 0.1 | ≥ 0.25 | ≥ 0.5 |
|---|---|---|---|
| off-target pairs kept | 743 | 540 | 357 |
| of which with an action (a sign) | 224 | 165 | 119 |
| of which partial agonist | 25 | 17 | 12 |
| of which the gene is a graph node | 731 | 532 | 349 |
| drugs with at least one | 127 | 117 | 105 |
| off-targets per drug, median (max) | 1 (31) | 1 (24) | 0 (17) |
| drugs with an off-target that is a labelled gene perturbation | 46 | 36 | 24 |

Most kept pairs (375 of 540 at 0.25) carry no action. Affinity panels record binding, not function.

These counts are lower than the median rule gave (540 pairs over 117 drugs at 0.25 against 622 over 119, and 357 against 408 at 0.5): the maximum raises the primary target's affinity along with the off-targets', and the primary sets the concentration an off-target has to compete with, so fewer pairs clear a threshold.

| drug | mechanism genes (primary pK) | highest-occupancy off-targets: gene, pK, occupancy, sign |
|---|---|---|
| amitriptyline | SLC6A2, SLC6A4 (8.37) | HRH1 9.30 0.99 −1; ADRA1D 8.25 0.87; ADRA1A 8.20 0.86 −1; CHRM4 8.14 0.84 −1; HRH4 8.14 0.84; HTR2C 8.10 0.83 |
| sertraline | SLC6A4 (9.54) | SLC6A3 7.60 0.09; SIGMAR1 7.24 0.04; ADRA2B 7.08 0.03; ADRA2A 6.96 0.02; ADRA1A 6.86 0.02; SLC6A2 6.38 0.01 (none kept at 0.25) |
| haloperidol | DRD2, DRD3, DRD4, HTR2A (9.42) | SIGMAR1 8.92 0.74; ADRA1B 8.10 0.30 |
| methylphenidate | SLC6A2, SLC6A3 (7.21) | ADRA2C 6.07 0.39 |

**Measurement coverage is uneven.** Human genes with an affinity per matched drug: median 6, quartiles 3 and 14,
maximum 68. Only 44 drugs have rows from DrugMatrix (a fixed panel run on every compound it holds) and 58 from PDSP.
A drug studied more has more measured off-targets, so off-target count partly measures research attention. Equal
loss shares per drug (`--equal-drug-shares`) stop a drug with many labels from taking a larger share of the gradient,
and input normalisation (`--normalise-drug-input`) stops a drug with many seeds from carrying a larger input. Neither
fills in the binding that was never measured. The evaluation stratified by measurement coverage, which the user asked
for, is built: see **Stratified evaluation by measurement coverage** below.

## Stratified evaluation by measurement coverage

The overstudied-drug reading, approved by the user on 8 October 2026. The point is not to correct for measurement
coverage but to report the metric inside strata of it: if a model's advantage over the baselines sits only among the
drugs that were screened against a panel, the advantage is a reading of how well the drug was studied.

`experiments/scope_off_target_binding.py --coverage-table` writes one row per perturbation
(`data/processed/off_target_scoping/measurement_coverage.parquet`): the drug's measured human genes with and without an
affinity, its primary affinity, which sources it came from, whether it carries DrugMatrix or PDSP rows, and its
off-target count at each occupancy threshold.
`mechanistic_pathway_learning/evaluation/measurement_coverage_strata.py` turns that into one stratum per scored row,
and `macro_auprc_by_stratum` in the metrics module reads macro AUPRC inside each. Both
`experiments/run_baselines.py` and `experiments/aggregate_main_model_runs.py` take `--measurement-coverage PATH` and
then carry `macro_auprc_by_measurement_coverage` in their results and a table in their markdown. Without the flag
nothing changes, and `experiments/score_confirmatory.py` is untouched: this is a development reading, not part of the
registered confirmatory decision.

Strata on `evidence_full_v3_max3_targets` with the better-label selection (1,619 perturbations, kept positives of the
label selection):

| stratum | perturbations | kept positive pairs | symptoms with at least 5 positives |
|---|---|---|---|
| `gene_perturbation` | 1,397 | 1,793 | 16 |
| `drug_without_binding_measurements` | 6 | 16 | 0 |
| `measured_genes_bin_0_[0,3]` | 72 | 273 | 14 |
| `measured_genes_bin_1_[4,11]` | 72 | 374 | 17 |
| `measured_genes_bin_2_[11,68]` | 72 | 417 | 16 |

Three readings of the drug terciles are available, so the comparison is not starved of labels. Points to keep in mind:

- A gene perturbation has no drug measurement and a drug with no row in either source has none either, so both are
  their own strata rather than being pooled with the measured drugs.
- The terciles are equal-count, ties broken by position, so a count that straddles a boundary (11 here) appears in two
  bins. The bin label carries its range and size, so this is visible in every table.
- Kept positives rise with coverage (273, 374, 417). That is the confounding this reading is for: the better-measured
  drugs also carry more labels, so a model and a popularity baseline both have more to work with there.
- The strata are drawn from the whole evidence table, not per fold, so the same drug is in the same stratum in every
  fold and the reading is comparable across runs.

## Leakage groups

`merge_drugs_with_their_targets` joins two perturbations when they share a target node. Off-targets join drugs
through receptors that most CNS drugs bind (adrenergic α1 and α2, H1, muscarinic, σ1, 5-HT2), so the groups collapse.
Three joining rules, at each threshold (largest group by drugs; share of all 222 drugs; share of kept positive pairs
in the group holding most):

| rule | ≥ 0.1 | ≥ 0.25 | ≥ 0.5 |
|---|---|---|---|
| every off-target joins, as a mechanism target does | 209 (0.94); 0.51 | 209 (0.94); 0.51 | 208 (0.94); 0.50 |
| an off-target joins its drug only to that gene's own perturbation | 204 (0.92); 0.50 | 204 (0.92); 0.50 | 203 (0.91); 0.50 |
| drugs joined only to the labelled genes they bind (mechanism or off-target); drugs sharing a target stay apart | 144 (0.65); 0.40 | 142 (0.64); 0.40 | 129 (0.58); 0.37 |

For comparison, without off-targets: the current rule gives 17 drug-holding groups, the largest holding 114 drugs
(0.51) and 0.27 of kept positive pairs; drugs joined only to labelled genes gives 67 drug-holding groups, the largest
holding 70 drugs (0.32) and 0.20 of kept positive pairs.

Reading: if off-targets join groups, about 94 percent of drugs fall into one group. A grouped five-fold split then
puts nearly every drug in one test fold, and a lockbox of 20 percent of drugs cannot be drawn as whole groups. Even
the narrowest rule leaves 58 to 65 percent of drugs in one group. The options are:

1. group by mechanism targets as now, add off-targets as input only, and report as a diagnostic the test drugs that
   share an off-target (occupancy ≥ 0.25) with a training drug or a labelled training gene;
2. join drugs only to the labelled genes they bind through mechanism targets (67 drug-holding groups, largest 32
   percent), with the same diagnostic for off-targets;
3. join on off-targets and give up grouped evaluation of drugs.

The current rule already concentrates drugs at cap 3: the two largest drug-holding groups hold 114 (51 percent) and 71
(32 percent) of 222 drugs (docs/drug_targets_any_type.md has the largest by perturbations). A lockbox of 20 percent
of drugs drawn as whole groups can then come only from the 71-drug group or from the small groups.

## Cap and input

The cap of 3 counts mechanism targets and decides which drugs enter. Off-targets would not count against it; they
would enter as extra seeds with magnitude equal to their occupancy and their sign (or the unsigned channel below).
Total input per drug = the mechanism magnitudes (1 each) + the sum of off-target occupancies. Normalising that total
to 1 (`--normalise-drug-input`) would shrink a drug's mechanism seeds as more of its off-targets are measured, so
input strength would follow measurement coverage. Normalising the mechanism part alone, or not normalising once
occupancy sets magnitudes, avoids that.

Pairs without an action could enter through a separate unsigned input channel: the graph learns that the drug binds
the gene without being told the direction. Partial agonists would go there too. This changes the encoders' input
layer, not the graph.

## Albumin and the other plasma substrate binders

Decided on 9 October 2026 (the user: "All of the plasma substrate binders should be modeled. Once AAG and albumin are
added, there will be state to attach that separate extension to."). The scope, the binders, the data audit and the two
build stages are in docs/plasma_protein_binding.md, which also records why orosomucoid rather than albumin carries
this study's drugs, most of them bases. Nothing is built yet: the carriage edges need one UniProt fetch and the
fraction unbound per drug is in no pinned source here.

`GENE:ALB` has 49 edges in `graph_full_neuronal`: 32 regulates_transcription_of, 10 binds (protein partners such as
LRP2, FCGRT, VCAM1 and B2M), 3 activates, 2 member_of and 2 catalyzed_by (the Human-GEM reactions that make albumin
from aminoacyl-tRNAs and break it down). No edge binds albumin to a small molecule, so the graph does not model
albumin carrying tryptophan, fatty acids, bilirubin, calcium or thyroxine, or the free fraction of any metabolite.
For drugs, albumin binding enters through the fraction unbound in the free concentration above. Albumin-metabolite
binding edges would be a graph change.

## Open Targets target safety curation as symptom labels

The Open Targets target safety data hold two curations. `adverse_effects.tsv`: 4,939 rows, 122 targets and 580
effect terms, each with a direction (activation or inhibition) and a dosing (acute, general, chronic or
developmental), from Lynch et al. 2017 (3,823 rows), Bowes et al. 2012 (647) and Urban et al. 2012 (469).
`secondary_pharmacology.json`: 210 records, 33 targets and 167 events (Brennan et al. 2024). These give
target-level labels: "inhibiting gene X causes symptom Y". The study has none of that kind today; its gene rows are
loss of function.

Mapping the effect terms to the 24 target symptoms by keyword (`experiments/scope_open_targets_symptom_labels.py`,
then each caught term checked by hand). The counts below are after the two pattern changes of 9 October 2026, the
`psychomotor_retardation` narrowing and the stereotypy flag, and after parkinsonism joined the symptom list:

| | count |
|---|---|
| caught (symptom, target, direction) triples | 200 |
| set aside: opposite direction ("decreased anxiety", "anxiolysis", "decreased aggression", "increased memory") | 44 |
| set aside: both directions in one term ("increased/decreased sleep") or no direction ("cognitive effects") | 8 |
| set aside: another thing ("respiratory depression", "Parkinson's disease" as a diagnosis) | 2 |
| positive triples | 146 (61 activation, 85 inhibition) |
| of which rodent readouts (locomotor activity, catalepsy) | 17 |
| of which wording the rodent assay and the clinic share (stereotypy) | 11 |
| symptom-target pairs | 141 (5 positive under both directions) |
| targets (graph gene nodes) | 52 (50) |
| target symptoms reached | 18 of 24 |

Per symptom: cognitive impairment 31, insomnia 18, psychosis 14, hyperactivity 13, compulsive behaviour 11,
anxiety 9, somnolence 9, depressed mood 8, psychomotor agitation 7, irritability or aggression 6, increased
appetite 5, catatonia 4, fatigue 3, decreased libido 2, elevated mood 2, parkinsonism 2, suicidality 1, abnormal
dreams 1. Not reached: anhedonia, apathy, disinhibition, emotional lability, psychomotor retardation, self-injury.
`psychomotor_retardation` left the list when its pattern stopped matching the rodent open-field count, which is
the whole of what the curation held for it; 4,344 of the 4,939 curation rows match no symptom pattern at all.

What a build would need:
- a perturbation kind "pharmacological inhibition (or activation) of gene X", seeded at X with sign −1 or +1, kept
  apart from the loss-of-function gene rows (a lifelong gene loss and an acute inhibition need not give the same
  symptoms);
- a crosswalk from effect terms to symptoms reviewed term by term, not keywords;
- its own evidence class and weight, below grade B, with the rodent readouts lower still;
- care with leakage: Lynch et al. and Bowes et al. draw partly on what known drugs do in patients, possibly including
  lockbox drugs. These rows would be training-only, kept out of the confirmatory metric, with a sensitivity analysis
  with and without them.

## Decisions for the user

1. Leakage rule with off-targets (options 1 to 3 above). Joining on off-targets puts about 94 percent of drugs in one
   group.
2. The cap-3 lockbox: with 51 and 32 percent of drugs in two groups, how to draw 20 percent of drugs.
3. The occupancy anchor: the 90-percent placeholder, or measured free concentration (Schulz et al. 2020 × Lombardo
   et al. 2018) where available.
4. The occupancy threshold: 0.25 as default, or a prespecified development grid {0.1, 0.25, 0.5}.
5. Input normalisation once off-targets enter: mechanism part only, or none.
6. Open Targets target-level labels: build them with the requirements above, or not.
7. Albumin-metabolite binding edges in the graph, or not.

Decided since this document was written (the user, 8 and 9 October 2026; docs/preregistration.md carries each decision
with its wording):

- the stratified evaluation by measurement coverage is built (section above), and the human-evidence review of the
  rodent-only curation rows is in docs/rodent_readout_human_evidence.md;
- one value per pair is the maximum, not the median (section "Rules"), which settles nothing about item 3 or 4;
- item 5, input normalisation once off-targets enter: the summary of the off-target occupancies enters the drug input
  normalisation, and no drug-level scalar is added beside the graph ("You can drop the off target summary into the
  normalization but mark it down as a note in case a situation comes up where it is needed to try"). The note holds the
  three feature forms that were offered instead: the mean occupancy over the kept off-targets, log(number of
  off-targets) and log(sum of the occupancies), the last two count-driven and closer to a research-attention proxy;
- item 7, the plasma binders: yes, in the two stages of docs/plasma_protein_binding.md;
- item 1 is measured rather than decided. docs/leakage_group_rules.md and
  docs/leakage_group_rules_evidence_full_v2.md test the target-family rules with both tie-breaks and report what
  leakage each one leaves; the registered grouping stays for the confirmatory reading.

## Parkinsonism, split from psychomotor slowing and then scored in its own right

The user, 9 October 2026: "it's more important to split the bradykinesia one since it has strong clinical treatment
implications between psychomotor slowing and drug-induced Parkinsonism." The split is made in
`experiments/scope_open_targets_symptom_labels.py`, where `psychomotor_retardation` matched `\bbradykine` until
then. It cost nothing to make, because the live labels never held the merger: `docs/symptom_crosswalk.csv` defines
`psychomotor_retardation` by two MedDRA terms, `Psychomotor retardation` and `Bradyphrenia`, gives it no HPO term
(HPO retired HP:0025356 into Global developmental delay) and no MeSH descriptor. The parkinsonian terms were in the
scoping pattern only, so the pattern was wider than the symptom it mapped onto.

Two things the narrowing showed:

- Not one row of the Open Targets curation names `psychomotor_retardation`. All 77 rows the old pattern caught read
  "decreased locomotor activity" (74) or "increased/decreased locomotor activity" (3) - the rodent open-field count,
  over 30 targets. The curation could never have supplied an honest label for that symptom, and after the narrowing
  it supplies none.
- The bradykinesia that appeared in docs/rodent_readout_human_evidence.md came from the human-evidence sentence,
  not from the curation: the DRD2-inhibition claim, curated as "decreased locomotor activity", has human evidence
  naming drug-induced parkinsonism, and `\bbradykine` sent that relabel straight back to the claimed symptom.

**Then the symptom was added** (the user, the same day: "it's still a potential MAJOR side effect of some drugs and
it would be very bad clinically if that side effect was not predicted so add it as a symptom"). The reach below is
what the decision was made on; `docs/parkinsonism_symptom.md` is what the rebuild delivered.

| source | reach |
|---|---|
| SIDER preferred terms, wide union including muscle rigidity, hypertonia and Parkinson's disease | 67 of this study's 142 drugs |
| SIDER preferred terms, the nine the crosswalk ended up giving the symptom | 53 of this study's 142 drugs |
| HPO, the union of the crosswalk's six roots, direct annotations | 181 of this study's gene perturbations |
| HPO HP:0001300 Parkinsonism alone | 119 |
| HPO HP:0002067 Bradykinesia alone | 103 |
| for comparison, `psychomotor_retardation` in SIDER (Psychomotor retardation, Bradyphrenia) | 9 drugs, of which 8 also carry a parkinsonian term |

Those are raw annotation counts before the better_v2 selection masks low-frequency positives and grade-C-only
pairs. After it, the symptom holds 138 kept positives, 106 gene and 32 drug. The comparison with
`psychomotor_retardation` is sharper than an earlier draft of this document said: that symptom has 5 evidence rows
but only 1 kept positive in the v2 selection and 2 in v3, against a macro floor of 5 kept positives, so it does not
enter the macro average at all, while 8 of the 9 drugs behind it in SIDER also carry a parkinsonian term. The
construct the study did score here was below its own reporting floor and nearly nested in one it did not score. The
narrower SIDER count is the one the crosswalk earns: Hypertonia, Muscle rigidity and Tremor are left out as
nonspecific and Parkinson's disease as diagnosis level under A7, which costs 14 drugs and is the price of the terms
meaning what they say.

The curation now supplies the symptom two target-direction pairs, DRD1 inhibition and SLC6A3 inhibition, from the
curated terms "parkinsonism" and "parkinsonian symptoms (tremors)"; neither is a rodent readout. All 7 parkinsonian
rows in the curation fall inside the symptom's pattern, so the crosswalk's exclusion of bare rigidity costs nothing
here. The DRD2 relabel lands on `parkinsonism` rather than on nothing.

**What still has to follow.** Adding a symptom changes the evidence table, the label selection and the outcome
matrix, and `configs/lockbox_v2.json` pins `evidence_records_sha256` and `label_selection_sha256`. The rebuild went
into new directories, so nothing the lockbox describes was touched. Moving to it means redrawing the lockbox from
the same rule and repeating the development runs, and it means adopting the salt-form and any-target-type change of
`docs/drug_targets_any_type.md` at the same time, because the current code no longer reproduces the v2 tables the
lockbox pins. That is the user's call and it has to be made before the lockbox is finalised, not after.

## References

See docs/references.bib: avram2023drugcentral, harding2024iuphar, bowes2012reducing, brennan2024state,
lynch2017potential, urban2012screening, smith2010effect, schulz2020revisited, lombardo2018trend,
hulme2010ligand.
