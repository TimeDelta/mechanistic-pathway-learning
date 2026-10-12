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
- **Built on the user's decisions of 11 October (section 13):** the evidence table v4 with the three label changes; the
  carrier, displacement, enzyme and two-drug statements as 214 checks that train nothing; a clearance edge from an
  enzyme to its drug; a table per drug of clearance, extraction and unbound concentration in plasma and brain; and
  CTD's terms beside the data derived from it. Section 13 also records three findings that differ from what was
  expected: the labels do not reward a drug for having many targets, so that confound is small; displacement is not
  shown to raise unbound concentrations; and no source supports competitive against non-competitive binding as a label.

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

**The other carrier: what is expected and what has been measured.** Rewritten on 11 October 2026, after the user
questioned the weight put on one paper. An earlier version of this paragraph was headed "Evidence against a free-drug
effect for the other carrier" and rested on one sentence of Espnes et al. 2026
([doi:10.1111/bcpt.70208](https://doi.org/10.1111/bcpt.70208)): "if the level of the drug-binding plasma protein
increases, the total amount of drug in plasma will increase, whereas the free concentration of the drug will remain
relatively unchanged". That sentence is from the paper's introduction and states an expectation. The study measured
total clozapine, orosomucoid and C-reactive protein in 1,106 samples of 369 patients and no unbound concentration
("analysing free (unbound) plasma clozapine concentrations would clearly have added valuable information"). Its own
data speak against binding as the whole account: "the unproportional high increase in clozapine concentrations compared
to the increase in AGP levels clearly indicates that increased protein binding is not the sole explanation". The two
markers rose together (correlation 0.591), so the study cannot separate more binding from slower clearance. It was not
evidence for the heading it was given.

What has been measured, in two small studies:

- **The unbound fraction falls as orosomucoid rises.** 1.22 percent against 1.65 percent in 6 samples with raised
  and 20 with normal orosomucoid (Man et al. 2019,
  [doi:10.1007/s40262-019-00744-6](https://doi.org/10.1007/s40262-019-00744-6)), and in 319 samples of 44 patients "A
  statistically significant (negative) correlation between AGP concentrations and clozapine (P = 0.000) and
  norclozapine (P = 0.028) unbound fractions was found" (Beex-Oosterhuis et al. 2025,
  [doi:10.1097/JCP.0000000000001979](https://doi.org/10.1097/JCP.0000000000001979)). Espnes et al. call the first "The
  only study that has measured free clozapine concentrations"; the second appeared in March 2025.
- **The unbound concentration has one comparison behind it**, in the same 26 samples: "7.2 ng/mL vs. 6.3 ng/mL (mean
  difference = 0.9 ng/mL, 95% confidence interval − 4.45 to 2.54" (Man et al. 2019). An interval about as wide as the
  mean does not show that the unbound concentration is unchanged. The authors had no doses, no smoking status and
  several samples from one patient.
- **Symptoms have not been studied against the carrier.** Case reports go both ways: a trough concentration "almost
  4-fold the upper limit of the reference range" in a patient who "did not report any adverse effects" (Espnes et al.
  2012, [doi:10.1097/FTD.0b013e3182666c62](https://doi.org/10.1097/FTD.0b013e3182666c62)) and a raised concentration
  "Without Signs of Intoxication" (Spangemacher et al. 2024,
  [doi:10.1097/jcp.0000000000001791](https://doi.org/10.1097/jcp.0000000000001791)), against a systematic review of 40
  cases in which "Adverse events, typically sedation, were associated with raised clozapine levels during infection"
  (Clark et al. 2018, [doi:10.1016/j.schres.2017.03.045](https://doi.org/10.1016/j.schres.2017.03.045)).

So the statement that more orosomucoid leaves the free drug and its effects unchanged is the pharmacokinetic
expectation for an oral drug cleared by the liver (`docs/drug_entry_nodes.md`), repeated by these authors and tested
by none of them. The conclusion of this section does not rest on it: no dataset varies this carrier with a drug and
records a symptom.

No review of the 2026 paper exists yet. PubMed and Crossref index no comment, letter or correction (11 October 2026).
Three works cite it; the one read in full text cites it in passing for "Inflammatory states have also been associated
with altered clozapine levels" (Özdin and Böke 2026,
[doi:10.1186/s12888-026-08386-2](https://doi.org/10.1186/s12888-026-08386-2)), and the other two, a systematic review
of clozapine and B-lymphocytes and a case report of myocarditis, were not read. PubPeer declines automated readers and
was not checked.

**Displacement pairs.** No curated list exists. The classic examples were reattributed: "recent evaluation of these
DDIs revealed that inhibition of metabolism of the victim drugs by the corresponding perpetrator drugs, rather than
protein binding displacement, to be the culprit in these cases" (Bohnert and Gan 2013,
[doi:10.1002/jps.23614](https://doi.org/10.1002/jps.23614)). The pairs with a central nervous
system consequence that survive are valproate with phenytoin, valproate with diazepam or carbamazepine and aspirin
with valproate or phenytoin, and each of them has a second, metabolic mechanism. The regulators' guideline on
interaction studies leaves the topic out of its scope: "Other types of pharmacokinetic interactions, such as impact
on absorption (e.g., gastric pH change, gastric motility change, formation of chelation or complexation), food
effects, or protein binding displacement, are not part of this document and may be covered by regional guidelines"
(ICH M12, section 1.3, EMA/CHMP/ICH/652460/2022). Its section 7.2, "Protein Binding", is about measuring the unbound
fraction of the drug that causes an interaction, not about displacement.

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
| a second drug | TWOSIDES (5,729,992 triplets from 2014 reports) | an event and a reporting ratio per pair | not counted: the authors' site says of OffSIDES and TWOSIDES "Note that these resources are quite a bit out of date.", and the data files carry no licence (section 11) | disproportionality signal |
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
the drugs fall into 17 leakage groups instead of 30. Corrected on 11 October: an earlier version of this sentence
said "one group holds 71 of them". That is the group with the most perturbations (241). The group with the most
drugs is another one, of 167 perturbations and 114 drugs, so two groups hold 185 of the 222 drugs and 89 percent of
the kept drug pairs (section 12). A score for drugs rests on few independent groups under either rule.

**Any number of targets** adds 8 drugs over the three-target rule and joins 194 of the 230 drugs into one group of
423 perturbations. An earlier version called that a poor trade. Section 12 measures it again: the join is the work of
one drug, and the three-target table was already concentrated in two groups.

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
408 have a disease annotated autosomal dominant. Gene2Phenotype records a mechanism on the same record as the
phenotype terms: "These include ‘Loss of function’, ‘Dominant negative’ and ‘Gain of function’. There are also
categories for ‘Undetermined’ and ‘Undetermined non-loss-of-function’" (Yates et al. 2024,
[doi:10.1186/s13073-024-01398-1](https://doi.org/10.1186/s13073-024-01398-1)). Its data were not fetched. Its
licence is CC BY 4.0, read on 11 October 2026 (section 11).

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
| dosage gain | ClinGen triplosensitivity | 3 of 1,532 genes (LMNB1, PLP1, TPSAB1) and 20 of 513 regions at the top score, counted on 11 October 2026 in ClinGen's own files (section 11). An earlier version of this row said "1.5 percent of 1,461 genes"; the paper's figure covers genes and regions together: "few genes/regions (1.5%) have sufficient evidence for TS (TS score 3)" (Thaxton et al. 2022, [doi:10.1002/humu.24291](https://doi.org/10.1002/humu.24291)) | scarcer for single genes than that row said, and at the level of a disease |
| mouse knockouts | IMPC, CC BY 4.0 | release 23.0: 9,277 genes, 113,803 significant calls; the current release, 24.0 of 16 March 2026, states "Number of phenotyped genes: 9,605" and "Number of phenotype calls: 111,664" (section 11.1) | "Some of these symptoms (e.g., depressed mood, suicidality) cannot be assessed in mice" (Nestler and Hyman 2010, [doi:10.1038/nn.2647](https://doi.org/10.1038/nn.2647)); usable for hyperactivity, anxiety-like behaviour, cognition, sleep, appetite and motor signs, and better as a test |
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

## 10. A prediction at the peak of a drug's concentration: the question, and what data exist

The user, 10 October 2026: "What do you think about adding a prediction at a calculated point of peak concentration?",
then "Or peak ispersion", "I meant pharmacokinetically" and "I'd like the count please". The standard terms for what
is meant are distribution and the effect-site (biophase) concentration: the unbound concentration of the drug where
its targets are.

### 10.1 The assessment given

- **A peak effect exists apart from exposure.** With immediate-release against extended-release quetiapine, "One hour
  after dosing on day 1, sedation was significantly greater with quetiapine IR than with quetiapine XR (mean VAS score,
  33.2 vs 11.3", at peak concentrations of 689 and 382 ng/mL and exposures of 2,836 and 2,515 ng h/mL (Datto et al.
  2009, abstract, [doi:10.1016/j.clinthera.2009.03.002](https://doi.org/10.1016/j.clinthera.2009.03.002)). A
  steady-state reading cannot express that.
- **The quantity that matters is at the site of action.** "the unbound (free) concentration of a drug at the site of
  action, such as the brain, is the driving force for pharmacological responses" (Loryan et al. 2022, abstract,
  [doi:10.1007/s11095-022-03246-6](https://doi.org/10.1007/s11095-022-03246-6)). It has two parts. The extent is the
  unbound brain-to-plasma partition coefficient, Kp,uu, which scales every reading, the steady state included, and is
  the measured quantity the record still owes in place of the ATC proxy for central action. The rate decides how much
  of a plasma peak the brain sees: "intrinsic brain equilibrium half-life [t(1/2eq,in) = V(b)ln2/(PS . f(u,brain))]"
  (Liu et al. 2005, abstract, [doi:10.1124/jpet.104.079319](https://doi.org/10.1124/jpet.104.079319)). A drug that
  equilibrates fast sees about the plasma peak times Kp,uu, and a slow one sees the average.
- **The model has no concentration, and the placeholder hides the cases that matter.** A drug enters at magnitude 1 per
  target, and the off-target scoping assumes the best-bound target 90 percent occupied. Measured occupancy can swing far
  more: quetiapine gives "transiently high (58%-64%) D2 occupancy 2 to 3 hours after a single dose that then decreases
  to minimal levels in 12 hours" (Kapur et al. 2000, abstract,
  [doi:10.1001/archpsyc.57.6.553](https://doi.org/10.1001/archpsyc.57.6.553)). A peak reading needs concentrations.
- **The labels carry no timing**, so a peak reading cannot be trained or scored on its own. It would be a reported
  check. The one source of the project with a time dimension is the Open Targets statements: 83 of the 146 record
  acute dosing only and 11 chronic dosing only.
- **It adds nothing for the carrier edges beyond the plasma peak.** Kp,uu "is independent of protein binding in blood
  and in brain" (Bohnert and Gan 2013, quoted in `docs/drug_entry_nodes.md`), and the half-life above has no
  plasma-binding term.
- **A peak cannot be read off the propagation.** The graph has no rate constants, so a step is not a unit of time. A
  scratch probe of the untrained confirmatory linear response (8 knockouts and 8 drugs, not in the repository) looked
  for the step of widest spread. With the input held on, the spread grows to its settled value and no node passes its
  settled value by 5 percent in more than 0.5 percent of entries. With the input given once, the spread counted against
  a fixed threshold peaks at step 4 and counted against the current largest value never peaks, so the step is set by
  the threshold.

### 10.2 The count

`docs/effect_site_data_coverage.md`, generated by `experiments/label_review/effect_site_coverage.py` from tables that
three research passes copied out of retrieved sources. For the 222 drugs of the three approved label changes (1,076
kept drug pairs):

| data | drugs | kept drug pairs on them |
| --- | --- | --- |
| extent: a measured Kp,uu in any species | 95 | 686 |
| extent: a human value (cerebrospinal fluid or PET) | 24 | 230 |
| rate: permeability and brain binding both, an equilibration half-life or a human effect-compartment rate | 60 | 427 |
| rate: a human effect-compartment rate from a journal | 13 | |
| extent and rate | 45 | 366 |
| occupancy: a human PET or SPECT study, any | 62 | |
| occupancy: at more than one time within hours to days of a dose | 24 | 215 |
| extent, rate and an occupancy time course | 10 | 111 |
| none of the three in any retrieved source | 97 | |

The ten with all three are bupropion, citalopram, clozapine, haloperidol, methylphenidate, olanzapine, paroxetine,
risperidone, sertraline and sulpiride; five of them (the antipsychotics) fall on the held-out side under the node
rule of section 9.2, as do 11 of the 24 with an occupancy time course.

How to read the count.

- **It is a lower bound on what is published.** The classic tables could not be downloaded: Fridén 2009 was read only
  as ChEMBL's transcription, 26 sources for the extent yielded no row (Kalvass 2007, Maurer 2005 and Doran 2005 among
  them) and the six primary sources named for the rate (Summerfield 2007 among them) were read only as abstracts. The
  occupancy studies were read as abstracts, one study per drug. The PubMed keyword counts in the generated
  document are the upper bounds: 68 drugs for an unbound brain concentration, 95 for an effect compartment, 62 for an
  occupancy with words for time. Two page requests of this count were withdrawn because a permission request went
  unanswered; both were repeated on 11 October (section 11) and neither adds a row: Spreafico and Jacobson 2013 had
  been read by another route and tabulates no value per drug, and Mandema and Danhof 1992, a review of
  electroencephalogram effect measures wanted for the rate, is a closed article.
- **Most extent and rate values are from rats or mice.** 38 of the 95 have a microdialysis value, the reference
  method. Homogenate, slice and microdialysis values can differ several-fold for one drug.
- **The rate inputs are not on one scale.** Permeability is reported in different units and under different assumed
  surface areas, and older in situ values are limited by flow for lipophilic drugs, so half-lives calculated from them
  come out too long.
- **The human rate constants are all for anaesthetics, opioids and sedatives** (alfentanil, dexmedetomidine,
  esketamine, etomidate, fentanyl, lorazepam, methadone, morphine, oxycodone, propofol, remifentanil, remimazolam,
  sufentanil). None was found in open text for an antidepressant, an antipsychotic or an anticonvulsant.
- **No occupancy time course was found for a benzodiazepine, a hypnotic or an antihistamine**, and no human occupancy
  study at all for a cannabinoid, an orexin antagonist or a melatonin agonist.
- **Three time courses are not courses within a dose** and are left out of the 24: cariprazine (one time per dose
  across days of dosing) and the irreversible enzyme inhibitors rasagiline and selegiline (recovery over weeks). Depot
  injections are not counted for the oral drug.

What follows from it, for the user's decision. The extent could be supplied for under half of the drugs, mostly from
rodents, so using it as a model input needs a rule for the drugs without one. An effect-site peak could be calculated
for the 45 with extent and rate and checked against measurement for 10. That supports a reported check on those drugs
and does not support a model input.

Text addressed to an AI assistant: the repository of the dosing-simulation program whose defaults one pass read
(stanpumpR) holds `AGENTS.md` and `CLAUDE.md`, files that orient AI coding agents working on that code; they were read
as text and not followed. The PubMed connector attaches a notice requiring attribution and DOI links, which is met
here. No journal article contained such text.

## 11. Pages and licences read again on 11 October 2026

Twelve page requests made by the research passes behind this document ended because a permission request waited
five minutes for the user and was then withdrawn, and the passes wrote the page down as not retrieved. The user's
ruling of 11 October: that is not an acceptable reason to leave a page unread, and a licence has to be verified. All
twelve were requested again that morning, and each licence that had been left unverified was then read from the
provider's own page, data file, source repository or paper (ClinGen's from its terms page as the user read it). No count in sections 1 to 10 changed except the ClinGen
row of section 6. None of the pages contained text addressed to an AI assistant.

### 11.1 The twelve pages

| page | what the second request returned | what changed here |
| --- | --- | --- |
| ClinGen gene-dosage statistics (`search.clinicalgenome.org/kb/gene-dosage`) | the site declines automated readers | the counts were taken from the files ClinGen publishes for download (`ftp.clinicalgenome.org`, gene list dated 11 October 2026, region list dated 28 April 2026): 1,532 genes, of which 3 have triplosensitivity score 3 (LMNB1, PLP1, TPSAB1), 1 score 2, 10 score 1, 1,306 score 0, 2 are scored unlikely to be dosage sensitive, 1 is marked as not to be evaluated and 209 are not yet evaluated; 513 regions, of which 20 have score 3. Section 6 corrected |
| ClinGen terms of use (`www.clinicalgenome.org/docs/terms-of-use/`) | the site declines automated readers; the user opened the page the same morning and pasted its text | licence verified (11.2) |
| DECIPHER downloads (`www.deciphergenomics.org/about/downloads/data`) | HTTP 403 to the fetch tool and to the shell | terms read from the DECIPHER paper instead (11.2) |
| IMPC, "Is IMPC data freely available?" | read | licence confirmed with a quote (11.2) |
| MGI statistics (`informatics.jax.org/mgihome/homepages/stats/all_stats.shtml`) | read; the page as returned shows no counts, only "last database update 09/29/2026" | nothing; no MGI count is used in this document |
| MGI copyright notice | read | licence verified (11.2) |
| `nsides.io` | read | the out-of-date statement is now quoted in section 3; no licence terms on the site (11.2) |
| openFDA label fields (`open.fda.gov/apis/drug/label/searchable-fields/` and the field file behind it, `open.fda.gov/fields/druglabel.yaml`) | read | the label sections named in section 3 exist as fields, each with a `_table` twin: `pediatric_use`, `geriatric_use`, `pregnancy`, `nursing_mothers`, `labor_and_delivery`, `teratogenic_effects`, `nonteratogenic_effects`, `pregnancy_or_breast_feeding`, `use_in_specific_populations`. There is no `lactation` field and no field for hepatic or renal impairment. The fields hold the section's text, so the reading "not structured" stands |
| ICH M12 guideline (EMA, 63 pages) | read | the scope sentence is quoted in section 2 |
| Smith and Waters 2019, Table I (`pmc.ncbi.nlm.nih.gov/articles/PMC7089466/table/Tab1/`) | the seven tables were read through NCBI's E-utilities | none here. Table II gives the concentrations in healthy humans: albumin 35 to 50 mg/ml (500 to 750 µM), orosomucoid 0.4 to 1.0 mg/ml (9 to 24 µM). Table I gives orosomucoid in health and disease by species; for humans, 0.5 to 1.0 mg/ml healthy and up to 3.0 mg/ml in inflammation, infection and cancer, a 3 to 6-fold change. The note of the carrier review that described Table I as the albumin and orosomucoid table was wrong about which table it is |
| Spreafico and Jacobson 2013 (`PMC4005893`) | a browser challenge at the PMC site | none. The rate pass had already read the article through NCBI's E-utilities; it has one table, of correlation coefficients, and no value per drug. The extent pass's "not retrieved" was wrong, and there is no row to add to section 10 |
| Mandema and Danhof 1992 (`link.springer.com/article/10.2165/00003088-199223030-00003`) | read: the publisher shows the summary and the reference list only | none. The article is closed and its tables stay unread; section 10 now says so |

A thirteenth page was read later the same morning, after the user added a rule to `.claude/settings.json` that
lets the fetch tool read any domain: the IMPC release summary (`mousephenotype.org/data/release`), which a pass had
reported it could not open. It gives release 24.0, dated 16 March 2026, with 9,605 phenotyped genes, 10,341 mutant
lines and 111,664 phenotype calls. The row of section 6 now carries both releases; its reading is unchanged.

### 11.2 Licences

| source | used by the study | terms, in the provider's words | where it was read |
| --- | --- | --- | --- |
| Gene2Phenotype | no (candidate of section 5) | "The G2P dataset is released under the Creative Commons Attribution 4.0 (CC BY 4.0) licence. If publishing work using G2P data, we request that you cite the latest G2P publication." | the page `www.ebi.ac.uk/gene2phenotype/responsible-handover`, read in the script file the site served on 11 October 2026 (the site draws its pages with a script, so the fetch tool sees an empty page), and in the site's source, `EBI-G2P/gene2phenotype_web`, commit 30dae72 of 9 October 2026 |
| CTD | yes (literature reports, evidence class E3) | "All aspects of CTD are freely available to non-commercial users; commercial users are asked to purchase a license in order to download CTD data files to their systems." (Davis et al. 2025, [doi:10.1093/nar/gkae883](https://doi.org/10.1093/nar/gkae883)). The header of the pinned file adds: "All forms of publication (e.g., web sites, research papers, databases, software applications, etc.) that use or rely on CTD data must cite CTD."; "You must notify CTD, and describe your use of our data"; "For quality control purposes, you must provide CTD with periodic access to your publication of our data." | the paper's data availability statement, and the first lines of `CTD_chemicals_diseases.tsv.gz` (report created 29 September 2026, the file `docs/data_sources.md` pins). The legal notice itself (`ctdbase.org/about/legal.jsp`) sits behind a check that the visitor is human and was not read; the header says its terms "include" the ones quoted, so the list may not be complete |
| MGI | no | "MGI data and annotations are licensed under a Creative Commons Attribution 4.0 International License" | `informatics.jax.org/mgihome/other/copyright.shtml` |
| IMPC | no (candidate of section 6) | "The IMPC makes this dataset publicly available under the Creative Commons Attribution 4.0 International license"; "You must give appropriate credit, provide a link to the license, and indicate if changes were made." | `mousephenotype.org/help/faqs/is-impc-data-freely-available/` |
| OffSIDES and TWOSIDES | no | none stated for the data. The authors' site, the release notes and the release repository (`tatonetti-lab/nsides-release`) carry no licence. The project's code repository (`tatonetti-lab/nsides`) is under "Attribution-NonCommercial-ShareAlike 4.0 International", and the 2022 rebuild of the OffSIDES method (`tatonetti-lab/offsides`) is under the MIT licence, which covers its software | `nsides.io`; the three repositories, read on 11 October 2026 |
| DECIPHER | no | "The openly consented patient data is available for bulk download for research purposes, subject to a data access agreement."; "DECIPHER also shares the data in bulk for display, subject to a Data Display Agreement." (Foreman et al. 2022, [doi:10.1002/humu.24340](https://doi.org/10.1002/humu.24340)) | the paper; the downloads page refused the request. The DECIPHER entries counted in section 6 come from the HPO annotation file, not from DECIPHER |
| ClinGen dosage curation | no (candidate of section 6) | "All curated content published by ClinGen is available free of restriction under the CC0 1.0 Universal (CC0 1.0) Public Domain Dedication. However, ClinGen requests that you give attribution to ClinGen and provide the date accessed whenever possible and appropriate." | the terms page (`www.clinicalgenome.org/docs/terms-of-use/`, dated January 26, 2018), opened by the user on 11 October 2026 and pasted into the session; the page declines automated readers, and the readme beside the download files states no licence |

Attribution, in the form ClinGen asks for: the following curated content was obtained from the Clinical Genome
Resource (www.clinicalgenome.org): Gene Dosage Sensitivity classification [11 October 2026] (Rehm et al. 2015,
[doi:10.1056/NEJMsr1406261](https://doi.org/10.1056/NEJMsr1406261)).

Two consequences for the study, both about CTD, the one source in this table that the study uses.

- `docs/data_sources.md` described the terms as "free for non-commercial research with citation". That is correct
  and incomplete; the row now carries the notification and access terms as well.
- The repository is public and its release directory holds a table derived from CTD
  (`data/releases/v0.4/literature/ctd_relations_filtered.csv`, 5,191 rows). CTD is cited (reference
  `davis2025ctd`). Whether CTD has been told of the use is not recorded anywhere in the repository; only the user
  can send that notice (`ctdbase.org/help/contact.go`). The licence sentence of the release manifest names SIDER and
  MedDRA and does not name CTD.

## 12. The limit on a drug's targets, measured again (11 October 2026)

The user, on 11 October: "I thought I had previously said to remove any limit on the number of targets", and a few
minutes later: "Oh. Was there something about the ones having more than three targets being not so relevant and concern
about confounds?"

**The record.** No instruction to remove the limit is in this session's record, in the documents or in the user's
other conversations that could be searched; the earliest part of the session survives only as a summary. What is
recorded: on 8 October, "fix that bug and incorporate drugs for any type of target" (any type of target, not any
number; `docs/drug_targets_any_type.md`, which measured the limits and left the choice to the user); on 10 October
the answer "Up to three targets" through the question prompt, and then "please continue implementing the three target
rule, gain of function sign, and placebo mask". The reasons given at the time match the user's second message:
"Most of that gain comes at a cap of 3", and "Total input magnitude grows with the number of targets (up to 7). A model
can then learn 'large input, many symptoms' from target count alone" (`docs/drug_targets_any_type.md`); and in
section 4 of this document, that any number of targets joins the drugs into one leakage group.

**Measured again**, on counting builds that carry the approved sign and placebo changes
(`build_approved_variants.sh`, `target_limit_comparison.py`, `bridging_drugs.py`). The groups are those of the
registered rule, which joins perturbations that share a node.

| drug rule | drugs | kept drug pairs | largest total input of a drug | drugs in the group with most drugs (its share of kept drug pairs) | second group |
| --- | --- | --- | --- | --- | --- |
| three targets | 222 | 1,076 | 3 | 114 (58.6%) | 71 (30.6%) |
| four targets | 224 | 1,085 | 4 | 187 (89.2%) | 6 |
| six targets | 226 | 1,102 | 6 | 190 (89.4%) | 6 |
| no limit | 230 | 1,112 | 7 | 194 (89.5%) | 6 |
| no limit, orphenadrine left out | 229 | 1,111 | 7 | 116 (57.5%) | 77 (32.0%) |
| draft list, three targets | 295 | 1,243 | 3 | 218 (84.8%) | 6 |
| draft list, three targets, dextromethorphan left out | 294 | 1,239 | 3 | 128 (54.0%) | 89 (30.8%) |
| draft list, no limit | 305 | 1,279 | 7 | 227 (85.2%) | 6 |
| draft list, no limit, both left out | 303 | 1,274 | 7 | 130 (53.2%) | 95 (31.9%) |

The eight drugs that enter only without a limit, with their mechanism targets and kept pairs: topiramate (6 targets,
16 kept pairs), nefazodone (4, 8), sevoflurane (7, 6), desflurane (7, 2), isoflurane (7, 2), trimipramine (6, 1),
orphenadrine (4, 1) and enflurane (7, 0). That is 36 kept pairs, 3.3 percent more, and 24 of them are topiramate's
and nefazodone's.

Three readings.

- **The main reason I gave against lifting the limit was wrong.** The drugs do not join into one group because drugs
  with many targets enter. They join because of one drug: orphenadrine acts on the histamine H1 receptor and the
  noradrenaline transporter, which lie in the block of 114 drugs, and on the NMDA receptor and sodium channels, which
  lie in the block of 71. It has one kept pair. With it left out, the table without a limit keeps the two blocks.
  Dextromethorphan (sigma-1 receptor and NMDA receptor, 4 kept pairs) does the same to the draft list at three targets,
  so the draft list and the limit have to be judged together. And the three-target table was never in many groups: its
  two blocks hold 89 percent of the kept drug pairs, which section 4 misreported.
- **What stands of the earlier reasons.** The gain is small, and four of the eight drugs are volatile anaesthetics,
  whose label events follow anaesthesia and surgery (10 kept pairs; this is my reading of what such labels describe,
  not a measurement). The input of a drug rises with its targets, to 7 for the anaesthetics and 6 for topiramate and
  trimipramine, so the question of how a multi-target drug's input is scaled
  (`docs/confirmatory_architecture_review.md`, question 6) weighs more without a limit.
- **Against the held-out set**, under the user's overlap rule of 10 October (`held_out_overlap.py`): six of the eight
  fall on the held-out side (the four anaesthetics, orphenadrine and trimipramine, 12 kept pairs) and topiramate and
  nefazodone train (24 kept pairs). In all 355 of 1,644 perturbations are held out, 66 of 230 drugs and 319 of 1,112
  kept drug pairs, against 349, 60 and 307 under three targets. The development drugs that share a target with a
  held-out added drug rise from 43 (254 kept pairs) to 103 (555 kept pairs), because the held-out anaesthetics act on
  the GABA-A receptor, as the benzodiazepines and many other development drugs do.

**Standing.** Three targets is the recorded decision and stays the rule of the approved builds until the user says
otherwise. The choice is now between three options, and none of them is registered: keep three; lift the limit and
set aside, by a written rule, a drug whose targets lie in both blocks; or lift it and accept one group.

## 13. Decisions of 11 October, 20:00, and what was built on them

The user's message of 11 October, 20:00, in the order of its lines, with what each led to. Nothing here draws or
registers a held-out set and nothing was trained.

| the user | what was done | where |
| --- | --- | --- |
| "'large input, many symptoms' from target count alone is counterable … What options did you come up with?" | the options and a measurement of the confound | 13.1 |
| "Sounds like phenytoin and valproate can be added for the carrier variance labels but not much else?" | confirmed: five statements over the two drugs, albumin only | 13.3 |
| "Add the valproate, etc displacements as checks then since those are important and survived" | five displacement statements, as checks | 13.3 |
| "Add the enzyme labels as the carrier twin instead then" | a clearance edge from an enzyme to its drug, and 23 statements as checks; read as a check, as row 4 of section 8 recommends, with training use left to the user | 13.3 |
| "Sounds like one of the label options you found can be used for competitive vs non competitive? That would be huge!" | no source supports it as a label; what exists is counted | 13.4 |
| "I accept all the suggestions under section 8 of the label review." "Please do the recommended label fixes" | the evidence table `evidence_full_v4` and its selection | 13.2 |
| "I accept put the three mechanism target pairs being mostly in the held out set" | the placement of 9.2, option a, as a draft table | 13.2 |
| "Why would working no better than the placebo not be counted as evidence against a label?" | an option that counts those pairs as negatives; the mask stays the default until the user chooses | 13.2 |
| "Has selectivity been modeled (different binding affinities for receptor subtypes, etc)?" | it has not; a first measurement of occupancy per target | 13.5 |
| "Can you please add the CTD reference to be compliant?" | citation, terms and links beside the derived table; a draft of the notice only the user can send | 13.6 |
| "Peak tissue unbound concentration can be relegated to confirmation only" | recorded; nothing built | section 10 |
| "Were you able to find clearance per drug?" "Seems worth adding data on a drug's concentration at its site of action." | a table per drug of clearance, extraction ratio and unbound concentration in plasma and brain | 13.5 |

### 13.1 The input of a multi-target drug

The question is number 6 of section 8 of `docs/confirmatory_architecture_review.md`; this review only pointed to it.

| option | what it does | cost |
| --- | --- | --- |
| 1 per target (as built) | a drug with three targets carries a total input of 3 | the size of the input follows the number of targets curated for the drug |
| a total of 1 (`--normalise-drug-input`) | each of k targets gets 1/k | says that a drug engages a target less because it has others, which pharmacology does not support |
| occupancy per target | the input at a target is C / (C + K) at the therapeutic unbound concentration | both numbers exist for 115 of the 222 drugs in plasma and 66 in brain (13.5), so the other drugs would sit on another scale |
| 1 per target, countered in the evaluation | a baseline that sees the number of targets and the total input only; results by number of targets; a permutation that exchanges target sets among drugs with the same number | adds no model; it shows the size of the confound and does not remove it |
| equal loss share per drug (`--equal-drug-shares`) | acts on the loss | does not change the input |

Measured on `evidence_full_v4` and the confirmatory graph: 155 drugs carry a total input of 1 (one of them 1.3), 45
carry 2 and 22 carry 3. Their kept positive pairs average 4.7, 5.1 and 5.6, and the rank correlation between a drug's
total input and its number of kept positive pairs is 0.06 (p = 0.41). Without a limit the eight further drugs carry 4,
6 or 7 and the correlation is 0.05. So the labels do not reward "large input, many symptoms": a model that read the
size of the input alone would gain almost nothing on them. The user's reading holds, and on stronger grounds than
the argument of section 12 gave. Occupancy does not remove the count either: the summed occupancy of a drug's
mechanism genes correlates 0.84 with their number.

Recommended, for the user: 1 per target with the three evaluation controls, and occupancy reported as a sensitivity
analysis on the drugs that have it.

### 13.2 The evidence table v4

`scripts/rebuild_data_layer.sh` builds it with the study's builders: `--max-drug-targets 3` and
`--seed-gain-of-function-genes-positive` for `data/processed/evidence_full_v4`, and `--mask-grades C
--mask-never-above-placebo` for `data/processed/label_selection/better_v2_full_v4.parquet`. The records, the reports and
the selection equal the counting build of 9.1 byte for byte (SHA-256, first 16: records `e0206dbe181747e8`, reports
`50d847913fa8ef9f`, selection `1d0cde3b1292217c`): 1,636 perturbations, 222 drugs, 2,975 kept pairs, 1,076 of them
drug pairs.

**Still open on this table**, each for the user: the limit on a drug's targets (section 12; three stays until the user
decides), the rows of the draft list of admitted drug groups (9.4; not in v4), and the placebo pairs below.

**Placebo pairs as negatives.** The objection is sound. An unlabelled pair is a weak negative in training and a
negative in scoring, on no evidence at all; a pair set aside is neither. So the 61 pairs for which a label shows no
excess over placebo count for less than pairs that no source mentions. `build_label_selection.py
--never-above-placebo-as-negative` writes the same 61 rows with `counts_as_negative`, and the loader turns each into a
labelled negative with the weight of any other negative (`better_v2_full_v4_placebo_negative.parquet`; drug negatives
go from 3,066 to 3,127 and nothing else changes). Two limits on reading them as evidence of absence: 48 of the 61
rest on one label entry with both arms, and SIDER gives percentages without the sizes of the arms; 7 of the 61 are
the drug's own recorded indication (fluoxetine, bupropion and buspirone with depressed mood, cabergoline and tolcapone
with parkinsonism, risperidone with irritability, hydroxybutyrate with somnolence). Recommended: count them as ordinary
negatives. The mask stays the default until the user says so.

**Held-out placement.** `held_out_overlap.py` on v4 gives the counts of 9.2: 29 of the 68 added drugs act on a held-out
node and are held out (160 of their 352 kept pairs), 39 train (192). Under the rule 349 of 1,636 perturbations are held
out, 60 of 222 drugs and 307 of 1,076 kept drug pairs. `configs/held_out_additions_v4_draft.csv` lists the 68 with the
nodes that place each one. It is a draft: no held-out set is drawn or registered on v4. The user's wording was
"mostly in the held out set"; by kept pairs it is 45 percent, the antipsychotics among them.

**Drug nodes.** `configs/drug_plasma_carriers_v4.csv` names a carrier for 62 of the 222 drugs (13 more than before:
five by the database alone, eight with a sentence of their label), and `graph_full_neuronal_split_binders_drugs_v4` holds 222 drug nodes
and 101 sequestration edges. Reading the labels of the added drugs showed two kinds of sentence the carrier reading
took wrongly, both now refused with a test: another ingredient's binding in a combination label (naloxone under
buprenorphine) and a sentence that names albumin without a word for binding (a dose rule of a radioligand). The table
of the earlier 154 drugs is reproduced unchanged.

### 13.3 Carrier, displacement, enzyme and two-drug statements, as checks

`configs/comparative_checks.csv` (214 rows; `experiments/build_comparative_checks.py`) holds every statement with its
quotation and where it was read. `mechanistic_pathway_learning/evaluation/comparative_checks.py` defines what each kind
compares, and `experiments/run_comparative_checks.py` puts them to a graph or to a finished run. Nothing is trained on
a check, and a check is reported apart from every score.

| kind | statements | source | what is compared | development, held out |
| --- | --- | --- | --- | --- |
| carrier level | 5 (valproate 3, phenytoin 2; albumin) | two drug labels and two hospital cohorts (`configs/carrier_and_displacement_statements.csv`) | the drug with its carrier lowered, against the drug | 5, 0 |
| displacement | 5 (valproate on phenytoin, diazepam and carbamazepine; aspirin on valproate and phenytoin) | the valproate label and the epilepsy guideline of Patsalos et al. 2008 | the two drugs with the carrier lowered, against the two drugs | 5, 0 |
| enzyme activity | 23 over 22 drugs (CYP2D6 16, CYP2C19 4, CYP2C9 2, NAT2 1) | FDA Table of Pharmacogenetic Associations (`configs/drug_clearing_enzymes.csv`) | the drug with the enzyme lowered, against the drug | 12, 11 |
| second drug | 80 positive triples, each asked with either drug as the reference, and 21 negative controls | CRESCENDDI | the two drugs, against one of them | 22 of 160 rows, 138; negative controls 15, 6 |

All 214 can be expressed on the graph with drug nodes and clearance edges (`docs/comparative_checks.md`). A finished
run is asked the development checks only, because its data have the lockbox removed.

**What these statements can and cannot carry.**

- **The carrier and displacement statements show a higher unbound fraction, and the sources say that the unbound
  concentration need not follow.** The guideline of the International League Against Epilepsy: "Although unbound drug
  concentrations will not be affected by changes in serum proteins, changes in protein binding need to be taken into
  account when interpreting total serum drug concentrations", and of displacement, "the unbound concentration remains
  unchanged or may even increase when the displacing drug also inhibits phenytoin metabolism" (Patsalos et al. 2008,
  [doi:10.1111/j.1528-1167.2008.01561.x](https://doi.org/10.1111/j.1528-1167.2008.01561.x)). Every displacement
  statement names inhibition of metabolism in the same sentence, and for valproate with phenytoin the rise in unbound
  drug is attributed to it (Lai and Huang 1993, [doi:10.1002/bdd.2510140409](https://doi.org/10.1002/bdd.2510140409)).
  So "survived" means that displacement is documented for these pairs. It does not mean that displacement is what
  raises the effect. Valproate is the exception with support of its own: its binding and metabolism are saturable,
  and in 41 patients "high unbound VPA concentrations are associated with clinically significant neurologic adverse
  symptoms" (Doré et al. 2017, [doi:10.1002/phar.1965](https://doi.org/10.1002/phar.1965)).
- **Displacement is written as a lowered carrier.** The graph has no edge from a drug to the carrier it occupies. Such
  an edge would make every drug with a known carrier move the carrier when taken alone, which is the extra input the
  change-carrying edges were built to avoid. A displacement check therefore tests the same edge as a carrier check,
  with two drugs present.
- **The enzyme rows are the better twin.** 23 of the 44 rows that name a study drug give an edge. Left out: 12 that
  state no direction ("May alter systemic concentrations."), 6 on an immune gene, 2 in which the enzyme forms the
  active compound (tramadol, where a clearance edge would have the wrong sign) and 1 about another ingredient. A row
  names concentrations or adverse reactions as a whole and no symptom, so a check reads the drug's own kept positive
  symptoms.
- **In the linear-response encoder a carrier or an enzyme reaches its drug through one signed edge with one gain per
  relation.** The comparison then holds or fails for every drug at once. There the checks test one learned sign per
  relation and not each statement.
- **The two-drug statements are mostly weak.** 74 of the 80 positive triples carry the evidence level "Theoretical"
  in the reference set (5 "Probable", 1 "Established"). In 31 the two drugs share a target with the same sign, which
  any encoder that is monotone in its input passes; the 47 that share no target are the informative ones, and most
  of them include a held-out drug.

**The clearance edge.** Relation `clears`, sign -1, from the enzyme's protein node to the drug's node
(`graph/drug_entry_nodes.py`; `experiments/build_drug_entry_node_variant.py --clearing-enzymes`). In the
message-passing encoder it carries the enzyme's change from its unperturbed state, as the carrier's edge does, each
relation through its own weights; a drug taken alone hears neither. The variant
`graph_full_neuronal_split_binders_drugs_v4_enzymes` is the v4 drug graph with 23 edge rows and one relation
appended; its node table is the same. No registered configuration reads it.

### 13.4 Competitive against non-competitive binding

No source found supports this as a label. What was counted:

- Guide to PHARMACOLOGY, release 2026.3 (24,934 interaction rows), marks 98 rows "Competitive" and 4
  "Non-competitive". None of the 647 human rows of the study's drugs carries either. What those rows do hold is the
  action (antagonist 226, inhibition 99, full agonist 92, partial agonist 60, agonist 56, positive 28, inverse agonist
  24 and others) and, for 40 rows over 19 drugs, a mode of binding (allosteric modulator 31, pore blocker 7,
  irreversible inhibition 3).
- Of the 80 positive two-drug triples, 2 have opposing actions at one target (aripiprazole with loxapine at the D2
  receptor, a partial agonist with an antagonist).
- The displacement pairs are competition at a carrier, not at a target, and 13.3 says what they can carry.

Competition at a target depends on the two concentrations and affinities, θ_A = (A/K_A) / (1 + A/K_A + B/K_B). The
table of 13.5 holds those numbers for about half the drugs, so the quantity could be computed for a pair of drugs
at one target. That would be a model input built from pharmacology, not a label, and no outcome data exist to test
it on beyond the two triples.

### 13.5 Clearance, concentration at the site of action and occupancy

`docs/drug_exposure_table.md` (`experiments/build_drug_exposure_table.py`), for the 222 drugs of v4:

| quantity | drugs |
| --- | --- |
| intravenous plasma clearance (Lombardo et al. 2018) | 137 |
| hepatic extraction ratio, computed | 101 (low 45, intermediate 29, high 27; above 1 for 15, where an assumption fails) |
| fraction unbound in plasma | 163 (125 from Lombardo et al.; 38 from the drug's FDA label, read by hand and checked against the cached sentence) |
| therapeutic total plasma concentration (Schulz et al. 2020) | 161 |
| unbound plasma concentration at therapeutic levels | 128 |
| unbound brain concentration (× Kp,uu,brain) | 69; the ratio is from humans for 12, from a non-human primate for 1, from other animals for 53 and of unstated species for 3 |

Two results.

- **The carrier edge's direction holds on pharmacokinetic grounds for few drugs.** More carrier lowers unbound
  exposure for a drug of high hepatic extraction given past the first pass (`docs/drug_entry_nodes.md`). Of the 62
  drugs with a sequestration edge, 45 have an extraction ratio, and 6 of those are of that kind (chlorpromazine,
  cocaine, diphenhydramine, lidocaine, morphine, sufentanil).
- **Occupancy** (`docs/target_occupancy_measured.md`). At the upper end of the therapeutic range the best-bound
  mechanism gene of a drug has a median occupancy of 0.92 in plasma (115 drugs; below 0.5 for 19) and 0.85 in brain
  (66 drugs; below 0.5 for 11). The placeholder of `docs/off_target_scoping.md`, 0.9 for every drug, is close to the
  median; in plasma 63 of the 115 drugs are at or above it and 19 are below 0.5. Within a drug the occupancies of its mechanism genes differ by a median
  of 0.17 in plasma and 0.33 in brain. Selectivity among a drug's targets is therefore measurable for about half the
  drugs and is not in the model: every mechanism target enters at 1.

Limits: active metabolites are not added to the parent's concentration; an EC50 or IC50 is pooled with Ki and Kd;
the concentrations are at steady state, not at a peak; local anaesthetics act at the site of injection, so their
plasma occupancy says little. The joined table stays out of git because the tables of Lombardo et al. and Benet et
al. are not licensed for redistribution; the licence of the copy of Lombardo's table that was read (the US EPA's
`httk` repository, MIT for the package) does not settle whether the table may be redistributed, which is the user's
call. The consensus guideline of Hiemke et al. 2018 would be the better source of therapeutic ranges for
psychotropic drugs and could not be read (`docs/exposure_sources/hiemke2018_agnp.md`).

### 13.6 CTD

- `data/releases/v0.4/literature/CTD_TERMS.md` quotes the four terms of the data file and the non-commercial term of
  Davis et al. 2025, cites the paper and gives the form of the links to CTD's chemical and disease pages. The fetch
  script writes the same file beside every table it derives, so a later release carries it
  (`mechanistic_pathway_learning/evidence/ctd_terms.py`). The release manifest lists the file and its licence sentence
  now names CTD; no table of the release changed and the release verifies.
- `README.md` carries the citation and the terms in short. `docs/experiment_design.md` called CTD "public"; corrected.
- `docs/ctd_notice_draft.md` is the notice for the user to send through CTD's contact form. It asks CTD whether a
  derived table may be redistributed in a public repository, because CTD's legal notice, citation guidelines and
  linking instructions decline automated readers and could not be read. The form of the links is the one the
  Bioregistry records for CTD.

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
| `effect_site_keyword_counts.py` | PubMed keyword counts per drug for occupancy, effect-compartment rates and unbound brain concentrations (section 10) |
| `effect_site_coverage.py` | the study drugs named in the retrieved tables of brain exposure, equilibration rate and occupancy over time (section 10) |
| `build_approved_variants.sh` | the four counting builds that carry the approved sign and placebo changes: three targets or no limit, with or without the draft list (sections 9 and 12) |
| `target_limit_comparison.py` | drugs, kept drug pairs, the group that holds the most drugs and each drug's total input, by variant (section 12) |
| `bridging_drugs.py` | the added drugs whose targets lie in both large drug blocks, and the groups with those drugs left out (section 12) |

Scripts of section 13, under `experiments/`:

| script | what it does |
| --- | --- |
| `read_label_protein_binding.py` | the percentage bound to plasma proteins of 123 drugs, read by hand from the cached label sentences and checked against them (`configs/drug_label_protein_binding.csv`) |
| `build_drug_exposure_table.py`, `measure_target_occupancy.py` | the exposure table and the occupancies of 13.5; the transcription scripts of the source tables are in `experiments/exposure_sources/` and the notes on each source in `docs/exposure_sources/` |
| `build_drug_clearing_enzymes.py` | the rows of the FDA pharmacogenetic table that give a clearance edge (`configs/drug_clearing_enzymes.csv`) |
| `build_comparative_checks.py`, `run_comparative_checks.py` | the table of checks and the report of what a graph can express (`docs/comparative_checks.md`) |
| `scope_label_plasma_binder_coverage.py --evidence-dir`, `tabulate_drug_routes.py` | the carrier and route tables of the v4 drugs |

Downloads, under `data/raw/label_review` (not in git):

| file | source | SHA-256, first 16 |
| --- | --- | --- |
| `fda_pharmacogenetic_associations/table_pharmacogenetic_associations.html` | https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations, fetched 10 October 2026 | `94b93d6dc36cfdd2` |
| `crescenddi/positive_controls.xlsx` | https://ndownloader.figshare.com/files/28584873 (doi:10.6084/m9.figshare.14847093.v1, CC0) | `9a6f2ab3db42e2e6` |
| `crescenddi/negative_controls.xlsx` | https://ndownloader.figshare.com/files/28584876 (doi:10.6084/m9.figshare.14847099.v1, CC0) | `f66a06d89960f0e7` |
| `clingen_dosage/ClinGen_gene_curation_list_GRCh38.tsv` | https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv, file dated 11 October 2026 | `e61cbc55bac615d8` |
| `clingen_dosage/ClinGen_region_curation_list_GRCh38.tsv` | https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv, file dated 28 April 2026 | `06bd0333e2c94585` |
| `clingen_dosage/README` | https://ftp.clinicalgenome.org/README, dated 18 March 2026 | `db2b953aa846ce30` |
| `ctd_terms/CTD_chemicals_diseases.header.txt` | the comment lines at the top of https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz (first 64 kB requested, 11 October 2026) | `1dfd31a9d620d6e4` |
| `ctd_terms/PMC11701581_fullTextXML.xml` | Europe PMC full text of Davis et al. 2025 | `2f1a859ba261356d` |
| `gene2phenotype_licence/ResponsibleHandoverView_commit_30dae72.vue` | https://github.com/EBI-G2P/gene2phenotype_web, `src/views/ResponsibleHandoverView.vue` at commit 30dae72; the script file the live site served (`index-DQD2xSP4.js`, SHA-256 beginning `f86c6741b989f04b`, 3.1 MB) holds the same sentence and was not kept | `ee554cfdb6c9c573` |
| `ich_m12/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf` | https://www.ema.europa.eu/en/documents/scientific-guideline/ich-m12-guideline-drug-interaction-studies-step-5_en.pdf | `126762154af51347` |
| `openfda_fields/druglabel.yaml` | https://open.fda.gov/fields/druglabel.yaml | `21a6aa9ae2a03a46` |
| `smith_waters_2019/PMC7089466_efetch.xml` | NCBI E-utilities, `efetch.fcgi?db=pmc&id=7089466` | `2aee3c5b8a952a8c` |
| `openfda_carrier_statements/{depakote,phenytoin,diazepam,carbamazepine}.json` | openFDA label queries of 11 October 2026 (the query is stored in each file); Depakote set_id 08a65cf4-7749-4ceb-6895-8f4805e2b01f, effective 2026-03-31; Phenytoin Sodium Injection set_id 035a8d4e-2063-4240-83cb-d7eebcabe301, effective 2024-04-09 | `ddddd55f984d389c` (Depakote), `71300476c1318504` (phenytoin) |
| `../public_pages/9513d9ed909c7063.html` | https://www.fda.gov/about-fda/about-website/website-policies, content current as of 01/19/2024, read 11 October 2026: "Unless otherwise noted, the contents of the FDA website (www.fda.gov) — both text and graphics — are not copyrighted. They are in the public domain and may be republished, reprinted and otherwise used freely by anyone without the need to obtain permission from FDA." | `9513d9ed909c7063` |

The outside sources were read by three research passes whose notes hold every quote with its source: the
repository inventory (333 quotes checked against the files), the combined-perturbation sources and the
single-perturbation sources. The notes are outside the repository. Limits they record: the ClinPGx counts are
keyword matches; no page read contained text addressed to an AI assistant. The passes also recorded several
licences and two statistics pages as not retrieved. Section 11 replaces that: every one of them was requested again
on 11 October 2026. Every licence named there is now quoted from the provider; the ClinGen terms page declines
automated readers and was read by the user.

Sources named above and not linked there: CRESCENDDI, Kontsioti et al. 2022
([doi:10.1038/s41597-022-01159-y](https://doi.org/10.1038/s41597-022-01159-y)); TWOSIDES, Tatonetti et al. 2012
([doi:10.1126/scitranslmed.3003377](https://doi.org/10.1126/scitranslmed.3003377)); ClinGen dosage curation,
[doi:10.1002/humu.24291](https://doi.org/10.1002/humu.24291); IMPC,
[doi:10.1093/nar/gkaf1148](https://doi.org/10.1093/nar/gkaf1148).
