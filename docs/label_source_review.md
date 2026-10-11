# Label sources: what the dataset uses, what it leaves out and what could be added

Written 10 October 2026, at the user's request before the held-out set is fixed: "Before finalizing the lockbox, I
want to be sure we have added all types of training examples that make sense and ones that vary the carrier with the
drug have yet to be added. Can you please also review the project for more potentially useful labels to add to the
dataset?"

Sections 1 to 8 change no label, table or configuration. Section 9, added later the same day, records the user's
decisions on them, what was built and what the builds showed; it changes no study table either. Every count marked
"measured" was made on the tables of this container by a script under `experiments/label_review/` (last section). Counts from outside sources are quoted
from the source or marked as keyword counts, which are upper bounds. Article metadata and abstracts were retrieved
through PubMed and each article is given with its DOI.

## The short answer

- **No dataset varies a plasma carrier with a psychotropic or an opioid and records a symptom.** What exists is label
  text for valproate and phenytoin, two small hospital cohorts and a handful of displacement pairs, each of which has
  a second mechanism. That is four to six comparative statements, enough for a named check and not for training.
- **The same question has an answer for another factor that changes a drug's exposure: the enzyme that clears it.**
  The FDA's table of pharmacogenetic associations has 26 rows on the study's 154 drugs and 44 on the 222 drugs of the
  three-target rule, every gene has a node in the graph, and the direction is the one pharmacokinetics supports for
  every route. An enzyme-to-drug clearance edge would be the carrier edge's twin, with examples behind it.
- **The largest gain in labels is in drugs the study already excludes by rule.** Allowing up to three mechanism
  targets adds 68 drugs and raises the kept drug pairs from 754 to 1,137. The drugs it adds are the ones psychiatry
  uses most: aripiprazole, risperidone, olanzapine, quetiapine, clozapine, haloperidol, venlafaxine, duloxetine,
  bupropion, mirtazapine, amitriptyline and methylphenidate among them.
- **Two corrections to labels already in use cost nothing.** 17 genes whose only causal association is a gain of
  function are seeded as a loss of function, and 40 kept drug pairs are never more frequent on the drug than on
  placebo in any label that reports both arms.

## 1. What the dataset uses now

Measured on `evidence_full_v3_parkinsonism` with the `better_v2` selection: 1,568 perturbations (1,414 genes and 154
drugs) and 24 symptoms.

| class | source | pairs | kept positives | note |
| --- | --- | --- | --- | --- |
| gene, loss of function (grade A) | HPO annotations through OMIM and Orphanet | 3,470 | 1,899 | every gene one seed of sign -1 |
| drug (grade B) | SIDER 4.1 and OnSIDES v3.1.1 label events | 1,503 | 754 over 106 drugs | one mechanism target, a nervous-system ATC code |
| grade C | non-causal gene associations | 804 | masked | neither positive nor negative |

Built into the table and read by no run: 85 `relieves` pairs over 55 drugs. Absent from the table: the literature
class (214,375 reports, a soft prior in the design). The full inventory, with the sentence and line that records
each item, is in the review notes (`repository_label_inventory.md`; see the last section).

## 2. A carrier varied with the drug: what exists

**Label text.** openFDA's 262,895 label documents were searched by phrase (keyword counts, upper bounds):

- "lower baseline albumin" occurs in 428 documents and every one is a valproate product. The Depakote label ties it
  to a named symptom, under "Somnolence in the Elderly": "A significantly higher proportion of valproate patients had
  somnolence compared to placebo... There was a trend for the patients who experienced these events to have a lower
  baseline albumin concentration, lower valproate clearance, and a higher BUN."
- Phenytoin labels name hypoalbuminaemia and attach an instruction, not a symptom: "Because the fraction of unbound
  phenytoin is increased in patients with renal or hepatic disease, or in those with hypoalbuminemia, the monitoring
  of phenytoin serum levels should be based on the unbound fraction".
- Alpha-1-acid glycoprotein appears in 2,461 documents as a pharmacokinetic description and in none beside an
  adverse event.

**Cohorts.** The nearest thing to a dataset is two hospital studies: 144 phenytoin and 383 valproate patients in
which the albumin-adjusted concentration went with neurotoxicity (Madan et al. 2021,
[doi:10.1111/bcp.14853](https://doi.org/10.1111/bcp.14853)), and 41 patients in which unbound valproate correlated
with neurological symptoms (Doré et al. 2017, [doi:10.1002/phar.1965](https://doi.org/10.1002/phar.1965)). A
systematic review of free valproate monitoring concludes: "This review demonstrates a paucity of data informing the
clinical utility of free VPA serum levels" (Lin et al. 2022,
[doi:10.1007/s40262-022-01171-w](https://doi.org/10.1007/s40262-022-01171-w)).

**Evidence against a free-drug effect for the other carrier.** A study of 1,106 clozapine samples states the caveat
itself: "if the level of the drug-binding plasma protein increases, the total amount of drug in plasma will increase,
whereas the free concentration of the drug will remain relatively unchanged" (Espnes et al. 2026,
[doi:10.1111/bcpt.70208](https://doi.org/10.1111/bcpt.70208)). This agrees with the pharmacokinetic reading in
`docs/drug_entry_nodes.md`.

**Displacement pairs.** No curated list exists. The classic examples were reattributed: "recent evaluation of these
DDIs revealed that inhibition of metabolism of the victim drugs by the corresponding perpetrator drugs, rather than
protein binding displacement, to be the culprit in these cases" (Bohnert and Gan 2013,
[doi:10.1002/jps.23614](https://doi.org/10.1002/jps.23614)). The pairs with a central nervous
system consequence that survive are valproate with phenytoin, valproate with diazepam or carbamazepine and aspirin
with valproate or phenytoin, and each of them has a second, metabolic mechanism.

**The one carrier-mediated neurological outcome written into labels is for another class and another direction.**
"Ceftriaxone can displace bilirubin from its binding to serum albumin, leading to a risk of bilirubin encephalopathy
in these patients." The drug occupies the carrier and the carrier's own cargo does the harm. That is the direction
`docs/drug_entry_nodes.md` found the literature to support, and neither the drugs nor the outcome are in this study.

**What follows.** Four to six comparative statements can be written down with a citation each (valproate with low
albumin and somnolence; phenytoin or valproate with low albumin and neurotoxicity; the displacement pairs). They can
serve as a named check of the carrier edge's direction, reported apart from every score. They cannot train it.

## 3. Other factors that vary with a drug

A training example of this kind needs a second factor the graph can represent and an outcome that differs from the
drug's own. Four were looked at.

| factor | source | what it holds | on this study's drugs (measured) | outcome form |
| --- | --- | --- | --- | --- |
| metaboliser status | FDA Table of Pharmacogenetic Associations (147 rows as parsed here, content of 10 September 2026) | drug, gene, affected subgroup, one sentence | 26 rows over 22 drugs under the current rule; 44 rows over 39 drugs with up to three targets (CYP2D6 30, CYP2C19 5, HLA-B 5, CYP2C9 2, NAT2 1, HLA-A 1) | free text with a direction: "Results in higher systemic concentrations and higher adverse reaction risk"; a symptom in few rows |
| metaboliser status and other variants | ClinPGx (PharmGKB) variant annotations, 14,560 rows; CC BY-SA 4.0, and CC BY-NC-SA 4.0 for the 2026 variant release | direction of effect, "risk of" a side effect, significance | 224 annotations touch one of the 24 symptoms by keyword; about 13 are significant, positive and anchored to a metaboliser phenotype | structured sentence; the variant's effect on the gene is known for star alleles only |
| a second drug | CRESCENDDI, CC0 (10,286 positive and 4,544 negative drug-drug-event controls) | a MedDRA event per pair | both drugs in the study and a study symptom: 2 pairs under the current rule, 43 pairs with up to three targets (somnolence and cognitive impairment), and 19 negative-control pairs | one event term per pair |
| a second drug | TWOSIDES (5,729,992 triplets from 2014 reports) | an event and a reporting ratio per pair | not counted: the authors call it out of date and the licence is not stated | disproportionality signal |
| a patient population | label sections for children, the elderly and pregnancy | free text | not structured; the paediatric file in the OnSIDES archive on disk is a manual annotation of 200 labels for testing an extractor, not a table of events | text |
| dose | none found across drugs | | | |

Three things follow from the table.

- **The enzyme is the better-supported twin of the carrier.** Lower clearance raises exposure whatever the route,
  which is not true of binding, and a regulator has reviewed the rows. All the genes named have nodes in the graph.
  An edge from a clearing enzyme to its drug, of sign -1 and carrying the enzyme's change, is the edge built for
  sequestration with another relation name.
- **Most of these rows say "more of everything the drug does", not which symptom.** As a yes or no label the
  combined example (the drug with its enzyme lowered) has the drug's own positives. What it adds is an ordering: the
  combination should score at least as high as the drug alone on the drug's positive symptoms. Using that needs a
  comparative loss, or the examples serve as a check. In the linear-response encoder the ordering holds by
  construction for any edge of that sign, so there it can fit one gain and test nothing.
- **Two-drug examples become available only with the multi-target drugs.** Under the current rule two pairs remain.

## 4. More examples of the kinds the study already has

`experiments/label_review/build_label_variants.sh` rebuilt the evidence table and the selection under six drug rules.

| rule | drugs | kept drug pairs | drugs with a kept positive | leakage groups holding a drug | largest group, perturbations (drugs in it) |
| --- | --- | --- | --- | --- | --- |
| nervous-system code, one target (current) | 154 | 754 | 106 | 30 | 236 (66) |
| nervous-system code, up to three targets | 222 | 1,137 | 156 | 17 | 241 (71) |
| nervous-system code, any number of targets | 230 | 1,173 | 163 | 15 | 423 (194) |
| any ATC code, one target | 720 | 1,733 | 413 | 196 | 300 (129) |
| any ATC code, up to three targets | 924 | 2,390 | 538 | 174 | 747 (447) |
| any ATC code, any number of targets | 955 | 2,464 | 561 | 168 | 804 (491) |

The current rule reproduces the registered counts, which checks the method.

**Up to three targets.** 68 drugs enter with 373 kept pairs. By kept positives: aripiprazole 18, risperidone 16,
venlafaxine 16, olanzapine 16, duloxetine 15, quetiapine 15, bupropion 14, trazodone 14, mirtazapine 13,
buprenorphine 12, levetiracetam 11, methylphenidate 10, then pimozide, ziprasidone, nortriptyline, clozapine,
imipramine, amitriptyline and haloperidol at 9 each. Catatonia gets its first kept drug pair and parkinsonism goes
from 32 to 57. `docs/drug_targets_any_type.md` gave two reasons against: a drug with k targets "carries total
magnitude k", and any cap above one "merges lockbox_v2 groups with v2 development perturbations". The first still
holds (corrected 10 October: an earlier version of this sentence said drug nodes removed it; the confirmatory graph has
no drug nodes, and with them each mechanism edge still delivers magnitude 1 per target; section 9). The second is a
reason to settle this before the held-out set is fixed, which is now; section 9 measures it. The cost that remains is in the last columns:
the drugs fall into 17 leakage groups instead of 30, and one group holds 71 of them, so a score for drugs rests on
few independent groups under either rule.

**Any number of targets** adds 8 drugs over the three-target rule and joins 194 of the 230 drugs into one group of
423 perturbations. That is a poor trade.

**Without the nervous-system rule.** 566 drugs enter at one target with 954 kept pairs, and the pairs are mostly
fatigue (309 of the 1,733 kept pairs), insomnia, anxiety and somnolence. By ATC group of the entering drugs:
cardiovascular 117 drugs and 237 kept pairs, antineoplastic and immunomodulating 151 and 226, alimentary 86 and 140,
respiratory 50 and 109, musculo-skeletal 40 and 102, sex hormones and genito-urinary 35 and 81, systemic hormones
30 and 47. Among them are drugs with a known central action that the ATC letter hides (guanfacine, baclofen,
tizanidine, promethazine, nabilone, rimonabant, and metoclopramide at three targets) and hormones whose receptors
are in the graph (testosterone 11 kept pairs, progesterone 8, cortisol 6, the gonadotropin-releasing hormone
agonists). Beside them are proton-pump inhibitors, anti-inflammatory drugs and cytotoxic drugs whose fatigue and
insomnia are more likely the illness or the treatment's general burden. The ATC code was a stand-in for central
action, and the record already says so ("a measured blood-brain barrier flag is still to replace the ATC N proxy",
`docs/onsides_label_slice_spec.md`). Lifting it whole would more than double the drug labels with pairs of that
second kind.

## 5. Corrections to labels already in use

**Gain of function seeded as loss of function.** Measured: 117 reports over 38 genes carry Orphanet's association
type "Disease-causing germline mutation(s) (gain of function) in", and 17 genes have causal reports of that type
only (C1R, C1S, CASR, DNAJC13, FGF23, GRIN2D, KCNA2, KCNE2, KCNQ1, LRRK2, NFKB2, PIEZO1, SCNN1A, SCNN1B, SCNN1G,
TBX20, TRPC6). Every gene is seeded with sign -1. The preregistration discloses this and leaves it ("so the signed
encoder reads the wrong sign for them"), because the table could not change then. It can now. The problem is larger
than Orphanet's field shows: "dominant-negative and gain-of-function mechanisms account for 48% of phenotypes in
dominant genes" (Badonyi and Marsh 2025,
[doi:10.1038/s41467-025-63234-3](https://doi.org/10.1038/s41467-025-63234-3)), and of the 1,663 genes with a row,
408 have a disease annotated autosomal dominant. Gene2Phenotype records a mechanism (loss of function, dominant
negative, gain of function) on the same record as the phenotype terms; it was not fetched and its licence was not
verified.

**Events no more frequent on the drug than on placebo.** Measured on SIDER's frequency table: 373 of the 754 kept
drug pairs have a drug arm and a placebo arm in the same label. In 314 every such label has the drug arm above its
placebo arm, in 19 the labels disagree, and in **40 no label does**. The 40 include fluoxetine with depressed mood,
alprazolam with insomnia, buspirone with insomnia, cabergoline and tolcapone with parkinsonism, and memantine with
psychosis and agitation: the symptom of the condition being treated, listed as an event. The evidence
specification records the gap (`docs/evidence_reports_spec.md`: "A treatment frequency at or below the placebo
frequency is not counted as evidence against (left for the appraisal work)"). The other 381 kept pairs have no
placebo arm to check.

## 6. Other single-perturbation label types

| type | source | size | reading |
| --- | --- | --- | --- |
| target inhibited or activated | Open Targets target-safety curation (scoped in `docs/off_target_scoping.md`) | 146 positive triples over 52 targets, 61 of them activations | a drug-like perturbation without a drug; an open decision of that document. Not the only source of positive-sign examples, as an earlier version of this row said: 388 of the 791 mechanism edges of the study's drugs are positive. Section 9 counts what it adds |
| relief of a symptom | SIDER indications | 85 pairs over 55 drugs in the table | recorded at the level of a diagnosis; a second relation with its own link matrix |
| copy-number syndromes | HPO entries of DECIPHER | 8 entries with a target symptom, 11 pairs | too few; each spans many genes |
| dosage gain | ClinGen triplosensitivity | 1.5 percent of 1,461 genes at the top score | scarce and at the level of a disease |
| mouse knockouts | IMPC release 23.0, CC BY 4.0 | 9,277 genes, 113,803 significant calls | "Some of these symptoms (e.g., depressed mood, suicidality) cannot be assessed in mice" (Nestler and Hyman 2010, [doi:10.1038/nn.2647](https://doi.org/10.1038/nn.2647)); usable for hyperactivity, anxiety-like behaviour, cognition, sleep, appetite and motor signs, and better as a test |
| chemicals that are not drugs, nutrients | CTD; no structured source links a nutrient level to a symptom | manganese with parkinsonism and vitamin B12 verified at the level of a syndrome | would be hand curation |
| spontaneous reports | OffSIDES | 329 off-label events a drug on average | it "recovers 38.8% (18,842 associations) of SIDER's associations", so reports and labels disagree on most pairs; a test, as `docs/pathway_test_options.md` already lists it |
| symptom-level genetics | item-level genome-wide studies, loss-of-function burden | found for about half of the symptoms; none found for catatonia, compulsive behaviour, libido, impulsivity, mania, adult hyperactivity or parkinsonian signs | effects are small and gene assignment is ambiguous; a test |

## 7. What the code assumes about a perturbation

A combined example (a drug with a gene, a carrier or a second drug) meets these assumptions, each found in the
code:

- A perturbation's type is "gene" or "drug" and nothing else; the selection, the loss shares, the validation strata
  and both lockbox scripts branch on "drug".
- A gene has one seed of sign -1.
- `merge_drugs_with_their_targets` joins drug rows through their seed nodes, so two drugs that shared a carrier seed
  would be chained into one leakage group.
- With drug nodes the loader replaces all of a drug's seeds by its one node.
- `rebalance_lockbox.py` stops when two perturbations share a seed node.

The seeds themselves are a general list of node, sign and magnitude, and batching handles any number. A third
perturbation type with its own group rule (a combination belongs to the union of its components' groups, so it
never sits on both sides of a split) is the change this needs.

## 8. Recommendations, for the user to decide

| | change | effect on the dataset | needs | recommendation |
| --- | --- | --- | --- | --- |
| 1 | drugs with up to three mechanism targets | 154 to 222 drugs, 754 to 1,137 kept drug pairs | the drug-node graph rebuilt on the new table; a held-out set drawn on it | **add** |
| 2 | sign +1 for the 17 genes with gain-of-function associations only | 17 perturbations change sign | one rule in the assembler; the 21 genes with both types stay as they are unless split | **fix** |
| 3 | mask kept drug pairs never above placebo | 40 pairs masked, neither positive nor negative | one rule in the selection | **fix** |
| 4 | enzyme-to-drug clearance edges and the pharmacogenetic rows as comparative examples | 26 to 44 rows; a new edge type | a clearance relation, a combined perturbation type, a comparative loss or a reported check | add as a check first; training use is a design decision |
| 5 | carrier statements for valproate and phenytoin | 4 to 6 statements | the combined perturbation type | a named check, reported apart |
| 6 | two-drug examples from CRESCENDDI | 43 positive and 19 negative pairs after change 1 | the combined perturbation type and group rule | a check first; small |
| 7 | target-level labels from Open Targets | 146 positive triples, 61 activations | its own evidence class, kept apart from gene loss | the user's open decision in `docs/off_target_scoping.md`; training-only as that document proposes |
| 8 | drugs outside the nervous-system code | up to 566 more drugs | a rule for central action that the ATC letter does not give | not as a whole; hormones and the centrally acting drugs named above are the defensible subset |
| 9 | relieves, mouse knockouts, spontaneous reports, genetics | | | not as training labels; tests |

Changes 1 to 3 alter the table that a held-out set is drawn on, so they come before it. Changes 4 to 6 add a kind of
perturbation and can follow, because a combination inherits the groups of its parts.

## 9. Decisions of 10 October, what was built and what the builds showed

The user's answers, through the question prompt: "Up to three targets", "Fix gain-of-function sign", "Mask
never-above-placebo"; combined perturbations "Build, report as checks"; "Add hormones and named central drugs". Then,
on the Open Targets statements: "can you measure mechanistic overlap with the held out ones just to be sure they don't
overlap and then include only the ones from that subset that pass the overlap test in training and put the others in
the lockbox?"

### 9.1 The three label changes, as flags that are off by default

| change | flag | where |
| --- | --- | --- |
| up to three mechanism targets | `--max-drug-targets 3` (existed) | `build_onsides_reports.py` and the assembler |
| sign +1 for a gene whose every causal association is a gain of function | `--seed-gain-of-function-genes-positive` | the assembler (`seed_monogenic_reports`) |
| a kept drug pair never more frequent on the drug than on placebo is set aside | `--mask-never-above-placebo` | `build_label_selection.py` (`placebo_arm_comparison.py`) |

Without the flags the builders reproduce the tables of the current rule byte for byte (evidence records, reports,
OnSIDES reports and selection, compared after the change). Tests: `tests/test_label_review_changes.py`.

The placebo comparison is made inside one label and one preferred term at a time. That is stricter than the count of
section 5, which pooled the terms of a symptom within a label, and it gives the same 40 pairs on the current table.

Measured on a counting build with all three (`data/processed/label_review/evidence_approved`, not a study table):

| | current rule | with the three changes |
| --- | --- | --- |
| perturbations | 1,568 | 1,636 |
| drugs | 154 | 222 |
| kept pairs | 2,653 | 2,975 |
| kept drug pairs | 754 | 1,076 |
| drugs with a kept pair | 106 | 156 |
| kept drug pairs set aside for placebo | 0 | 61 |
| genes seeded +1 | 0 | 17 (18 kept pairs over 13 of them) |

Apart from the seeds of the 17 genes, the evidence table equals the three-target table of section 4.

**The input of a multi-target drug.** On the three-target table 46 drugs carry a total seed magnitude of 2 and 23 a
total of 3: each target enters at magnitude 1. That is the pharmacological reading (occupancy at one target does not
fall because the drug has others), and it makes the size of a drug's field follow its number of targets. Drug nodes
do not change it, because each mechanism edge carries the magnitude of the seed it replaces. The user's decision
(`docs/confirmatory_architecture_review.md`, question 6).

### 9.2 The three-target rule and the held-out set

`experiments/label_review/held_out_overlap.py` reads the membership of `configs/lockbox_v2.json` (never a prediction)
and applies the user's rule: an added example that overlaps a held-out perturbation is held out with it. Overlap is
tested on nodes, one step: an added drug is held out when it acts on a node that a held-out perturbation acts on or is.

The held-out set carried to the current table by its leakage groups is 320 perturbations, 31 of them drugs. It holds
the dopamine D2 group (bromocriptine, cabergoline, fluphenazine, perphenazine, pramipexole, prochlorperazine,
ropinirole, rotigotine, sulpiride, benperidol, levodopa), the alpha-2 agonists (clonidine, dexmedetomidine,
lofexidine) and two H1 antagonists (diphenhydramine, hydroxyzine).

| of the 68 added drugs | drugs | kept pairs | which |
| --- | --- | --- | --- |
| act on a held-out node | 29 | 160 | the antipsychotics and other drugs with a D2, alpha-2 or H1 action: aripiprazole 18, olanzapine 15, quetiapine 14, mirtazapine 13, risperidone 13, asenapine 10, clozapine 9, haloperidol 9 and 21 others |
| do not | 39 | 192 | venlafaxine 16, duloxetine 15, bupropion 11, buprenorphine 11, levetiracetam 11, methylphenidate 10, trazodone 10 and 32 others |

Under the rule 349 of 1,636 perturbations are held out, 60 of 222 drugs, and 307 of 1,076 kept drug pairs.

Two things the rule leaves open.

- **The registered group rule chains, and chaining is not usable here.** It joins every perturbation that shares a
  node with another. The 29 drugs act on a held-out node and on development nodes, so on the three-target table the
  registered rule would pull 149 more perturbations into the held-out set (96 drugs, 715 kept pairs). `lockbox_v2`
  cannot be carried over to this table under the registered rule without that.
- **The 29 share their other targets with development drugs.** Their targets outside the held-out set are 23 nodes
  (serotonin receptors, the serotonin and noradrenaline transporters, alpha-1 receptors, H3). 43 development drugs act
  on one of those nodes (254 kept pairs), the reuptake inhibitors, the tricyclics and the triptans among them. A model
  trained on those drugs has seen part of the mechanism of each of the 29.

Three ways to place the 29, for the user:

| | the 29 | cost |
| --- | --- | --- |
| a | held out and reported as their own stratum, outside the primary score | the primary score keeps the registered separation; the antipsychotics do not train and do not decide H1 |
| b | held out and in the primary score | the primary score then includes drugs that share a target with 43 training drugs |
| c | left out of the table | 160 kept pairs unused, among them most of the added parkinsonism pairs |

Recommended: a. It keeps `lockbox_v2` and its separation as registered, uses the 29 as held-out evidence, and needs
no redraw. None of the 68 was ever in a pilot.

Catatonia's one kept pair is in the held-out set, so catatonia has no training positive under any of these.

### 9.3 The Open Targets statements against the held-out set

The 52 statements that no study drug label makes, in human wording, are 49 once the two whose target is not a graph
gene (GABRG1, SCN4B) are left out. The same script places them, with the held-out set of 9.2 (the 29 included). A
statement is held out when its target is a node of a held-out perturbation. Statements with the same direction,
symptom and references whose targets are subunits or members of one receptor (one ChEMBL target record lists both
genes, or the gene symbols share their stem) are one sentence of the source, so they go together.

| | statements | targets | independent sentences |
| --- | --- | --- | --- |
| held out | 16 (14 by their own target, 2 with their receptor) | 10 | 13 |
| train | 33 (17 activations, 16 inhibitions; 20 from acute dosing only) | 19 | 25 |

Held out: AR inhibition (four symptoms; AR is a held-out gene), HTR2A activation (three), HTR2C activation (two),
ADRA2A, DRD1, HTR3A, HTR7 and the three sodium-channel beta subunits. Training: the GABA-A subunits with insomnia
(eight statements, one sentence), CNR1 activation (six symptoms), NR3C1, GRIN1, CCKAR, the muscarinic receptors and
others.

What the measurement says about the plan.

- **The overlap test is the node rule.** A graph distance does not work as a test: on the confirmatory graph the
  development perturbations already sit a median of 2 edges from a held-out seed (515 at one edge, 606 at two). The
  statements that would train sit at 1 edge (17), 2 (12) or 3 (4), no closer than ordinary training examples.
- **It is small.** 33 statements are 1.4 percent of the training kept pairs, and 25 independent sentences. What they
  add that the table lacks is the activation direction on 17 of them.
- **Five statements are contradicted by a study drug** that acts on the target in the statement's direction and does
  not list the symptom; one of the five would train (CNR1 inhibition with depressed mood). The 52 are by construction
  the statements no study drug confirms.
- **Six of the inhibition statements that would train name a target that is also a gene perturbation**, the same node
  lowered; the gene's own labels list the symptom in two.
- **Eleven would train on a target that no other training example touches.**
- **In the lockbox they would be scored against the weakest labels of the study**: review statements about a target
  class, with no frequency and no recorded absence. Recommended: held-out statements are reported as their own
  stratum, as in 9.2 a, and do not enter the primary score.

The weight in training is a separate question from the split. The study weights every example by its evidence grade;
these statements need a grade. The user's decision.

### 9.4 Hormones and named centrally acting drugs: the draft list

`configs/admitted_atc_prefixes_draft.csv`, for the user's sign-off. It changes no table until it is passed to the
builders with `--admitted-atc-prefixes`. Each entry has its reason and source; a named drug has the sentence of its
label that states the central action. The rule for an entry is pharmacology, not the adverse events the drug is known
for, since admitting a drug because it causes the symptoms would select on the outcome.

| entry | what it is | drugs added | kept pairs | largest |
| --- | --- | --- | --- | --- |
| H | systemic hormonal preparations | 31 | 45 | desmopressin 7, cortisol 6, dexamethasone 5, octreotide 5, paricalcitol 5 |
| G03 | sex hormones and modulators of the genital system | 18 | 42 | testosterone 11, progesterone 8, norethisterone 4, levonorgestrel 4 |
| C02A | antiadrenergic agents, centrally acting | 5 | 14 | guanfacine 8, moxonidine 3, methyldopa 2 |
| M03B | muscle relaxants, centrally acting agents | 5 | 21 | baclofen 10, tizanidine 5, cyclobenzaprine 3 |
| A08AA | centrally acting antiobesity products | 5 | 1 | phentermine 1 |
| R05DA | opium alkaloids and derivatives | 3 | 8 | dextromethorphan 4, hydrocodone 3, codeine 1 |
| named | metoclopramide 8, rimonabant 8, promethazine 7, nabilone 6, aprepitant 4, dronabinol 3 | 6 | 36 | |
| all | | 73 | 167 | 43 of the 73 have a kept pair |

With the list and the three changes the table has 295 drugs and 1,243 kept drug pairs. Against the held-out set, 13
of the 73 act on a held-out node and would be held out (46 kept pairs: testosterone, guanfacine, metoclopramide,
promethazine, tizanidine, danazol, methyldopa, cinacalcet, chlorzoxazone and four with no kept pair); 60 would train
(121 kept pairs).

Considered and left off the draft, each for the user to overrule (counts from the any-ATC three-target table of
section 4, before the placebo mask):

| group or drug | drugs, kept pairs | why not |
| --- | --- | --- |
| L02, endocrine therapy (the gonadotropin-releasing hormone agonists, anti-oestrogens, anti-androgens, aromatase inhibitors) | 24, 90 | hormone antagonists, so within "hormones" on one reading; given to cancer patients, whose fatigue, insomnia and low mood have another cause. The largest single addition available, and the one most open to that objection |
| A10A, insulins | 3, 1 | hormones by any definition; one kept pair |
| C07A, beta blockers | 18, 64 | not hormones and not named as centrally acting; outside the user's wording |
| first-generation antihistamines other than promethazine (cyproheptadine, cyclizine, clemastine, carbinoxamine, chlorpheniramine, meclizine) | 6, 12 | their labels state sedation as an effect, which is close to selecting on the outcome |
| 5-HT3 antagonists (granisetron, ondansetron, palonosetron) | 3, 13 | act at the area postrema and on vagal afferents, outside the blood-brain barrier |
| domperidone | 1, 2 | the peripheral counterpart of metoclopramide |
| montelukast, isotretinoin, ribavirin | 3, 31 | known for psychiatric adverse events, with no established central target: selecting on the outcome |

One inconsistency found on the way and not changed: the SIDER route tests a drug's SIDER ATC codes only, the OnSIDES
route the union of the RxNav and SIDER codes. Fampridine has a nervous-system code in RxNav and none in SIDER, so its
OnSIDES statements qualify and its SIDER events do not (2 kept pairs).

Sources of the quotations in the list: US labels through openFDA (set identifiers in the file; responses cached under
`data/raw/label_review/openfda_pharmacology`); the 2006 label of Cesamet and the European public summary of Acomplia,
both read through a fetch tool, which returns a model's reading of the page and not the page. No page contained text
addressed to an AI assistant.

## How the counts were made

Scripts, all under `experiments/label_review/`, none of which writes to a study table:

| script | what it counts |
| --- | --- |
| `build_label_variants.sh` with `run_with_any_atc.py` | the six evidence tables of section 4, built into `data/processed/label_review` by the study's own builders; the runner replaces the nervous-system test before the builders import it and edits no file |
| `summarise_label_variants.py` | drugs, kept pairs and leakage groups of each variant |
| `atc_breakdown.py` | the ATC groups of the drugs that enter without the nervous-system rule |
| `placebo_check.py` | drug arm against placebo arm within each SIDER label |
| `hpo_unused.py` | the DECIPHER entries of the HPO annotation file |
| `fda_pgx_overlap.py` | rows of the FDA pharmacogenetic table that name a study drug |
| `ddi_overlap.py` | CRESCENDDI controls with both drugs and a symptom in the study |
| `open_targets_overlap.py` | the Open Targets statements a study drug label already makes |
| `held_out_overlap.py` | which side of the held-out set each added drug and each Open Targets statement falls on (section 9) |

Downloads, under `data/raw/label_review` (not in git):

| file | source | SHA-256, first 16 |
| --- | --- | --- |
| `fda_pharmacogenetic_associations/table_pharmacogenetic_associations.html` | https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations, fetched 10 October 2026 | `94b93d6dc36cfdd2` |
| `crescenddi/positive_controls.xlsx` | https://ndownloader.figshare.com/files/28584873 (doi:10.6084/m9.figshare.14847093.v1, CC0) | `9a6f2ab3db42e2e6` |
| `crescenddi/negative_controls.xlsx` | https://ndownloader.figshare.com/files/28584876 (doi:10.6084/m9.figshare.14847099.v1, CC0) | `f66a06d89960f0e7` |

The outside sources were read by three research passes whose notes hold every quote with its source: the
repository inventory (333 quotes checked against the files), the combined-perturbation sources and the
single-perturbation sources. The notes are outside the repository. Limits they record: several licences were not
retrieved (Gene2Phenotype, OffSIDES and TWOSIDES, CTD, MGI); the ClinPGx counts are keyword matches; the ClinGen and
IMPC statistics pages could not be opened; no page read contained text addressed to an AI assistant.

Sources named above and not linked there: CRESCENDDI, Kontsioti et al. 2022
([doi:10.1038/s41597-022-01159-y](https://doi.org/10.1038/s41597-022-01159-y)); TWOSIDES, Tatonetti et al. 2012
([doi:10.1126/scitranslmed.3003377](https://doi.org/10.1126/scitranslmed.3003377)); ClinGen dosage curation,
[doi:10.1002/humu.24291](https://doi.org/10.1002/humu.24291); IMPC,
[doi:10.1093/nar/gkaf1148](https://doi.org/10.1093/nar/gkaf1148).
