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
| (12 October) "let's do the control version and the equal loss per share version" | the three controls built, equal loss shares turned on in the amendment draft | 13.1 |
| (12 October) "Is there any other sources / ways to deal with ligand binding competition?" | ten sources surveyed and measured against this study; one route remains | 13.4 |
| (12 October) "shouldn't better evidenced absence via placebo control should count more" | a loss weight for those negatives alone | 13.2 |
| (12 October) the placeholder for C "wouldn't work" | a concentration per target, measured, replacing the anchor | 13.5 |
| (12 October) the threshold grid is "a poorly correlated replacement for dispersion level" | accepted; what tissue labels the graph has, and what can be built without a new model variant | 13.8 |
| (12 October) "If you think there's enough enzyme statements to train on, include them" | not enough, and the checks on them test nothing; what is kept instead | 13.7 |
| (12 October) "Can you please do a deeper check for the Lombardo clearance table license?" | permission is required; two routes remain, both the user's | 13.9 |
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

**Decided by the user, 12 October 2026: "let's do the control version and the equal loss per share version."** So the
input stays 1 per target, the three controls are built, and `--equal-drug-shares` is turned on. What each is:

| control | what it is | where |
| --- | --- | --- |
| the input-size baseline | a model given a perturbation's number of seeds and its total input magnitude and nothing else: no graph, no node identity, no target identity. Per symptom a logistic regression on those two numbers. Whatever it scores is the part of a score that input size alone explains | `models/baselines/input_size_baseline.py`; `run_baselines.py --with-input-size` |
| results by number of targets | macro AUPRC inside strata of how many mechanism targets a drug has (1, 2, 3, 4 or more; gene perturbations their own stratum, since they outnumber the drugs six to one). An advantage that sits only in the many-target stratum is a reading of input size | `evaluation/input_size_strata.py`; reported by `run_baselines.py` and `aggregate_main_model_runs.py` without a flag |
| the permuted-target control | every drug's seeds exchanged with those of another drug that has the same number of them, inside groups of equal count, as a derangement. Each drug keeps its input size and loses which targets it acts on, so a run that scores as well here is reading size | `run_main_model.py --permute-drug-targets`; a negative control, never a result of the study |
| equal loss shares per drug | each drug's positives and, apart, its negatives rescaled to the mean counts over the training drugs, so a drug with many labels does not take a larger share of the gradient | `run_main_model.py --equal-drug-shares`, which existed and was off |

Two notes on the baseline, because they decide how to read it. Under the input rule as built the two features are the
same number (total magnitude = number of targets for every drug but one), so the fit carries a small ridge on the
slopes; it is kept at two features so the same baseline reads correctly under `--normalise-drug-input` (total 1, count
varying) and under an occupancy input (both varying, not proportional). And a gene perturbation has one seed of
magnitude 1, so for the genes this baseline is the popularity baseline: its reading is about the drugs.

`--equal-drug-shares` changes training, so turning it on in the confirmatory configurations is a registered change;
it is written into `docs/preregistration_amendment_draft.md` and not applied here. Nothing was trained
("No new pilots yet").

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

**And they should count for more than silence.** The user, 12 October 2026: "shouldn't better evidenced absence via
placebo control should count more (just in terms of loss share)?" Yes, and the asymmetry was in the build. An
unlabelled pair is a negative at `--negative-weight` 0.2, on no evidence at all; a relabelled placebo pair was taking
the same 0.2 although a controlled comparison stands behind it. `--placebo-negative-weight` now sets the weight of
those pairs alone (the loader marks them in `evidenced_negative`, which only a selection built with
`--never-above-placebo-as-negative` fills). Unset, they weigh what any negative weighs, so nothing changes until it is
set. The weight is in the loss only: a negative is a negative in every metric, which have no weights.

What the weight should be is the user's. The case for 1.0, the weight of a well-evidenced positive: a placebo-arm
comparison is evidence of the same kind as the label statement that makes a positive. The case against: 48 of the 61
pairs rest on one label entry with both arms, SIDER gives percentages without the sizes of the arms, so the comparison
can be two small numbers; and the seven pairs that are the drug's own indication are the strongest of the 61 as
evidence of absence but also the most confounded by the treated condition. A middle value (0.5) says the absence is
better evidenced than silence and weaker than a reported event.

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

**Corrected on 12 October 2026:** only the two-drug kind is a test of a model. The other three go through a relation
that no label of the study trains, so they are reported without a verdict. Section 13.7 gives the reason and what is
kept in their place.

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
  24 and others) and, for 41 rows over 20 drugs, a mode of binding (type "Allosteric modulator" 31, action "Pore
  blocker" 7, action "Irreversible inhibition" 3). Recounted on 11 October with that definition; an earlier count in
  this session gave 40 rows over 19 drugs.
- Of the 80 positive two-drug triples, 2 have opposing actions at one target (aripiprazole with loxapine at the D2
  receptor, a partial agonist with an antagonist).
- The displacement pairs are competition at a carrier, not at a target, and 13.3 says what they can carry.

Competition at a target depends on the two concentrations and affinities, θ_A = (A/K_A) / (1 + A/K_A + B/K_B). The
table of 13.5 holds those numbers for about half the drugs, so the quantity could be computed for a pair of drugs
at one target. That would be a model input built from pharmacology, not a label, and no outcome data exist to test
it on beyond the two triples.

**Other sources, surveyed on 12 October 2026** after the user asked whether any exist. A research pass measured each
against this study's 222 drugs and 148 mechanism targets; the counts below are that pass's, and any the study adopts
should be re-derived in code here before use. The family split of the targets decides how far a family-specific
resource reaches: voltage-gated channels 44, G-protein-coupled receptors 42, Cys-loop ligand-gated channels 20,
ionotropic glutamate receptors 10, SLC6 transporters 4, other 28, protein kinases 0.

| source | what it would give | reach into this study | licence, as read | verdict |
| --- | --- | --- | --- | --- |
| ChEMBL `drug_mechanism.site_id` and `binding_site_comment` | the name of the pocket a drug's mechanism acts at | **the best label found**: 23 rows over 21 drugs, separating the benzodiazepine site (13 drugs), the barbiturate site (4), the GABA site and the α2δ-1 site. Whole database: 251 of 7,561 mechanism rows carry a site | "The ChEMBL data is made available on a Creative Commons Attribution-Share Alike 3.0 Unported License." | **usable, for GABA-A.** Already on disk (`data/raw/chembl/mechanisms.json`). Three of its 24 comments are statements of uncertainty and must stay uncertain |
| DrugCentral `ACT_COMMENT` | the displaced radioligand of an assay, so two drugs displacing the same tracer at one gene occupy one site | **the widest reach, and already on disk**: a tracer for 278 drug-gene records over 54 genes and 59 tracers, giving 1,021 of 1,307 same-gene drug pairs a shared tracer | CC BY-SA 4.0 (already in `docs/data_sources.md`) | **usable as evidence, not as a label.** A negative allosteric modulator can also displace a tracer |
| ChEMBL assay descriptions | the same, as free text, at a larger scale | 11,027 of 56,741 assays on the study's 126 ChEMBL targets contain "Displacement of"; 675 of 1,345 drug pairs share a tracer at the 14 targets sampled | as ChEMBL above | usable, with tracer-name normalisation; the text also needs filtering, since "competitive" in an assay description usually names the assay format and not the mechanism |
| PDBe-KB ligand binding sites | the residues a ligand contacts, as UniProt positions, so two ligands' pockets can be compared | pockets definable on 107 of 148 targets; only 22 of 222 drugs are co-crystallised with one of their own targets | the PDB archive is CC0 1.0 | **usable to define the pockets** that the two sources above assign drugs to. The pass demonstrated it: on GABRA1 it separates the benzodiazepine-site ligands from propofol's transmembrane site, and on the dopamine transporter it puts bupropion in dopamine's own site |
| UniProt binding features | where the endogenous ligand binds | 85 of 147 targets, naming glutamate, dopamine, serotonin, GABA, histamine, noradrenaline, adenosine and glycine | CC BY 4.0 | **usable for the drug-against-endogenous-ligand half** of the question |
| GPCRdb | residue contacts with a generic numbering that is comparable across receptors | 42 of 148 targets | "The data in GPCRdb is freely available under the Creative Commons Attribution 4.0 International license." | usable for the aminergic receptors |
| pA2, pKB and KB (the Schild output, a typed field in ChEMBL, GtoPdb and DrugCentral) | the measurement that settles competition experimentally | **4 pA2 activities in all of ChEMBL for this study's molecules**, 5 pKB rows in GtoPdb, 0 in DrugCentral | as above | **empty.** Schild analysis belongs to isolated-tissue pharmacology; this drug set was characterised by radioligand binding |
| Allosteric Database (ASD) | curated allosteric sites, which is the other half | not measured | "Copyright © 2009-2023 Molecular Design Laboratory, SJTU All Rights Reserved." | **unusable.** Named here so it is not re-explored |
| KLIFS | separates orthosteric from allosteric ligands outright, which is exactly what this needs | **0 of 148 targets**: it covers kinases | no licence sentence found on its pages | not applicable |
| DrugBank, PDSP Ki database, Drug Target Commons | affinities and action vocabularies | not measured | **unverified**: DrugBank's terms page answered HTTP 403, the other two could not be reached. Addresses for the user to open: `go.drugbank.com/legal/terms_of_use`, `pdsp.unc.edu/databases/kidb.php`, `drugtargetcommons.fimm.fi` | PDSP's values already reach the study through DrugCentral (356 rows over 57 drugs) |

So the answer to the user's question is: **no source labels competition for these drugs, and one route computes it.**
The route in order of reach: take ChEMBL's site name where it exists (21 drugs, and it is the GABA-A split, which this
study's benzodiazepines, barbiturates and Z-drugs need); treat a shared displaced tracer as evidence of a shared site
(about 1,000 drug pairs, no download); use structures to define the pockets that the first two assign drugs to; and
take the mode of binding already counted above as the switch for the pairs that do not compete.

Once the site is shared, the strength of the competition needs nothing new. With both ligands at one site the loss of
A's occupancy caused by B is θ_A(alone) / θ_A(with B) = 1 + ([B]/K_B) / (1 + [A]/K_A), so B costs A a factor that is
entirely [B]/K_B, and 13.5 holds both numbers. For a pair at different sites the ternary-complex form needs a
cooperativity factor α, which no resource tabulates at scale; α ∈ {0, 1} chosen from the mode of binding is the
defensible simplification, and it is an assumption to register, not a label derived from data.

One caution for the Cheng-Prusoff conversion this study already makes: the correct correction from IC50 to K_i depends
on knowing whether the binding is competitive. So the absence of a competition label is not an alternative to fixing
the pooling of IC50 with K_i (section "Rules" of `docs/off_target_scoping.md`) but a prerequisite for it.

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

**A concentration per target, which is what the user asked for.** The user, 12 October 2026, on the scoping's
placeholder: it "is biased towards the drugs best-bound mechanism, which precludes some dose-dependent changes to
concentration modifying the effect via i.e. saturation of highest affinity target. Thus, this wouldn't work." That is
right, and it is not a limitation of the placeholder's accuracy but of its form: the placeholder *defines* C as nine
times the best target's dissociation constant, so the best target is 90 percent occupied by construction, for every
drug and at every dose. Nothing about saturation can be read from it, because saturation is what it assumes.

The measured concentration has no such anchor, so each target's occupancy is its own. `measure_target_occupancy.py`
now also reads each target at its own site: the unbound brain concentration for a target the graph records as
expressed in a brain tissue, the unbound plasma concentration otherwise and for a drug with no brain value. On the v4
table that is 305 mechanism pairs over 115 drugs, 122 of them read in brain. What the anchor hid:

| | placeholder | measured, at each target's own site |
| --- | --- | --- |
| drugs covered | 183 (an affinity at a mechanism gene is enough) | 115 with a mechanism pair, 121 with any pair |
| occupancy of the best-bound mechanism target | 0.90 for every drug, by construction | median 0.90, quartiles 0.65 and 0.98; at or above 0.9 (saturated) for 59 drugs, below 0.5 for 21 |
| off-target pairs at occupancy ≥ 0.25 | 540 over 117 drugs | 374 over 81 drugs |
| off-target pairs at ≥ 0.1 / ≥ 0.5 | 743 / 357 | 510 / 252 |

So the measured route costs 68 drugs of coverage and buys the thing the placeholder cannot express: at a therapeutic
dose a third of the drugs do not saturate their best target, and a fifth are below half occupancy there. The two can be
combined, measured where there is a measurement and the placeholder otherwise, at the cost of two scales in one
column; recommended instead is the measured concentration alone, with the drugs that lack one reported as their own
stratum, since that is what the measurement-coverage strata already do (section "Stratified evaluation by measurement
coverage" of `docs/off_target_scoping.md`). The user's decision 3 of that document.

Limits: active metabolites are not added to the parent's concentration; an EC50 or IC50 is pooled with Ki and Kd;
the concentrations are at steady state, not at a peak; local anaesthetics act at the site of injection, so their
plasma occupancy says little. The joined table stays out of git because the tables of Lombardo et al. and Benet et
al. are not licensed for redistribution; the licence of the copy of Lombardo's table that was read (the US EPA's
`httk` repository, MIT for the package) does not settle whether the table may be redistributed. **Checked again on
12 October 2026 at the user's request, and the answer is now negative rather than unknown** (13.9). The consensus guideline of Hiemke et al. 2018 would be the better source of therapeutic ranges for
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

### 13.7 What a check can test, and the enzyme statements

The user, 12 October 2026: "If you think there's enough enzyme statements to train on, include them. Otherwise checks
are fine I guess but they're almost certain to fail without having trained any enzymes so I don't see the point in
checking them."

The objection is right, and it is sharper than the enzyme statements. It applies to the carrier edge too, which was the
original design (`docs/drug_entry_nodes.md`). No label of this study varies a carrier or an enzyme, so neither
relation is trained, and a check that goes through one reads something other than carriage or clearance:

- **in the message-passing encoder** the relation has its own weight matrix per layer, 32 x 32 x 2 = 2,048 parameters
  at the registered state dimension, and no gradient ever reaches them. The check reads the random initialisation;
- **in the linear-response encoder** a signed relation's gain is `sigmoid(...)`, strictly positive, and the edge carries
  the sign, so the direction in which a lowered carrier or enzyme moves the drug's node is built into the graph. What
  the check then reads is whether the head is monotone in that node's field, which is a property of the head.

So the three kinds that go through `sequesters` or `clears` are not tests of a model, and
`experiments/run_comparative_checks.py` no longer reports a verdict for them: it reports the difference and says why
there is no verdict. Only the two-drug kind is judged, because a second drug enters through its mechanism edges, the
relation every drug perturbation trains.

**Is there enough to train on? No, and the reason is the same arithmetic.** The 23 enzyme statements cover 22 drugs
and 133 kept positive pairs, of which 12 statements, 12 drugs and 75 pairs are on the development side. Against 2,048
untrained parameters in the message-passing encoder, 75 ordering constraints fit nothing. In the linear-response
encoder the relation has four parameters, one gain per channel, which 75 constraints could fit — but the direction is
already built in, so training would set a magnitude against constraints that any positive magnitude satisfies. Neither
encoder learns anything about clearance from these statements. Training them would need a comparative loss and a
constrained gain, which is a new model variant, and the user's instruction stands: "I don't want to add extra
confirmatory model variations because that starts increasing probability of winning by luck."

What is kept, then:

- the **edges** stay. They are how the design represents a carrier and an enzyme, and the drug-node ablation the user
  asked for on 10 October reads them;
- the **direction** of each edge is checked once, against the literature, not per run: the sign is -1, the sources say
  more carrier or more enzyme means less free drug, and `tests/test_clearance_edges.py` checks that the field moves
  that way. That is a design check, and it is done;
- the **statements** stay in `configs/comparative_checks.csv` with their quotations, because they are the evidence the
  edge's direction rests on and the file is where a reader finds it;
- the **two-drug statements** stay as the only per-run check: 160 rows over 80 triples, 22 rows on the development side,
  and 21 negative controls.

The honest summary of what this buys: one real check (the two-drug comparison, mostly on held-out drugs), one design
check already passed, and a documented reason not to read the other 33 rows as tests.

### 13.8 Dispersion, and whether the graph has tissues

The user, on the occupancy threshold grid of `docs/off_target_scoping.md`: "I'm not convinced on the occupancy
threshold grid based on the number of nodes touched. It's a poorly correlated replacement for dispersion level. What
about replacing it with the reading at different points of maximum spread based on node tissue type".

The first half is accepted: the number of off-target nodes a threshold admits (median 1, maximum 31 per drug) is a
count of what was measured about a drug, not of how far its effect spreads, and the stratification by measurement
coverage exists precisely because that count tracks research attention. Choosing a threshold by it is choosing by
coverage.

The second half needs a label the graph does not carry, and the measurement is worth stating plainly:

- **there is no tissue label.** The confirmatory graph's 49,574 nodes carry 43 *subcellular* compartment strings (10
  single compartments over 18,928 nodes: cytosol 8,098, extracellular 3,713, mitochondrion 2,104, endoplasmic
  reticulum 1,909 and six smaller), and no gene or protein node carries one at all;
- **the nearest thing to a tissue is a brain cell class**, and it is a weight, not a partition: the cell-class table
  gives 11 classes (astrocyte, excitatory neuron, dopaminergic neuron, medium spiny neuron, oligodendrocyte lineage,
  vascular and five others) with a weight per node. 32,234 nodes have weights that differ between classes; a node is
  active in a median of 11 of the 11 classes, and only 311 in exactly one. The classes do discriminate (median
  pairwise correlation of the class weightings 0.36, none above 0.9; the largest class takes a median 23 percent of a
  node's weight against 9 percent for an even split), but they cover genes, proteins and 7,681 of the 13,793
  reactions, and no metabolite or protein entity;
- so "peak distribution per tissue type" resolves either to **11 weighted spread curves over one shared node set** (the
  cell classes, covering the neuronal layer) or to **10 compartment partitions** (covering the metabolic layer). Neither
  is a tissue, and the first is 11 readings of overlapping subsets.

The memory cost the user anticipated is not the obstacle. The linear-response encoder propagates 8 steps over 4
channels, so storing every step is 49,574 x batch x 4 x 8 floats, about 200 MB at a batch of 32 — affordable. The
obstacle is that reading the field at a per-class peak is a new encoder variant, and the instruction against new
confirmatory variants stands.

**What can be done without a new variant, and is the useful half of the idea:** compute the spread per step per cell
class as a *measurement*, not a readout, and use it to choose the occupancy threshold. That gives a dispersion
statistic defined on the propagation itself rather than on the off-target count: for each candidate threshold, how far
and into which cell classes a drug's field reaches, and at which step it peaks. It needs no change to any model, it
answers the question the threshold grid was standing in for, and it can be read on the development perturbations
before any threshold is fixed. Not built; recommended as the next measurement if the user wants the threshold decided
on evidence rather than by default.

### 13.9 The licence of the clearance table, checked again

The user, 12 October 2026: "Can you please do a deeper check for the Lombardo clearance table license?" A research pass
tried the routes `docs/exposure_sources/lombardo2018.md` had not: the publisher's author-rights and permissions pages,
whether supplementary material is licensed apart from the article, author and repository deposits, the licence fields
of Crossref and OpenAlex, and whether anyone has reused the table with a stated permission. **The position changed from
"no licence found" to "the stated terms require permission".**

- The copyright holder is the journal's society, confirmed in three records (Crossref assertion, PubMed, OpenAlex).
  OpenAlex reports the article as open access of the "bronze" kind with `license: null`: **free to read at the
  publisher, under no licence.**
- The only licence Crossref records for the article is a text-and-data-mining licence, and it says: "Substantially or
  systematically reproduce, retain or redistribute the Dataset" is not allowed, and text mining is "provided to
  subscribers for non-commercial research purposes" (`elsevier.com/tdm/userlicense/1.0/`, read 12 October 2026). The
  publisher's copyright policy says the same in its own words: a user may mine subscription articles "without sharing
  any adaptation of the original content with others".
- The journal is not in the publisher's open-archive programme, so **no post-embargo licence will attach** to a 2018
  article. That closes the "wait and it becomes reusable" route permanently.
- The data-as-facts route is narrower than it looks. The publisher's permissions policy allows "Creating an original
  figure or table from data or factual information that was not previously in figure or table format" — and Table S1
  was already a table.
- **The US EPA route does not carry.** A whole-repository search of the `httk` package finds the table named in eight
  files, every one a citation or build code, and no permission, provenance or rights statement anywhere. The EPA's own
  data licence covers data "produced by the U.S EPA", which this is not. The one EPA deposit that is licensed CC0
  (doi:10.23645/epacomptox.6062791) contains a single 0-byte `index.html`: a metadata stub, not a copy of the table.
  Recorded so that a later pass does not read that CC0 as covering it.
- **Two routes remain, and both are requests the user would make.** The publisher's permissions route is the Copyright
  Clearance Center, `marketplace.copyright.com/rs-ui-web/mp/search/journal/1521009X`, and the society states: "ASPET
  handles all permissions requests through the Copyright Clearance Center's Republication Service." The shorter route
  is the authors: the publisher's own copyright policy says authors retain "the right to use their research data freely
  and without restriction", so a deposit of Table S1 by its authors under CC BY or CC0 would settle it for everyone. No
  such deposit exists today (Zenodo, figshare and Dryad all searched).

**What this changes here: nothing, which is the point.** The joined exposure table was already kept out of git, the
transcription script and the source address are committed, and that posture is exactly what the terms permit. What
changes is the record: `docs/exposure_sources/lombardo2018.md` and `docs/data_sources.md` now say that redistribution
needs permission rather than that a licence was not found, and they carry the sentence that says so.

One correction the pass found in that file, and it matters for every future retrieval: the journal host is not blocking
the supplementary path in particular. Its robots file disallows the **whole host** for the user agents an AI
assistant's reader announces itself under, so every address there is refused, including the free full-text PDF. No
other path on that host is worth trying, and the user has to open them personally.

## 14. Answers of 12 October, 14:24

Five points the owner pressed on after the answers of 13.1 to 13.9. Each is measured by a script in
`experiments/label_review/` and nothing here changes a label, a table, a graph or a model input.

### 14.1 The binding-site comments of ChEMBL, and a correction

`experiments/label_review/chembl_binding_sites.py`, on the 222 drugs of `data/processed/evidence_full_v4`: 326
mechanism rows over 179 drugs, 23 of them with a curated `site_id` over 21 drugs and 24 with a
`binding_site_comment` over 22 drugs.

**Correction.** The pass of 12 October reported four comments carrying explicit uncertainty. There are three. The
fourth, gabapentin's, is not a hedge but a citation:

- acamprosate, GABA-A, POSITIVE MODULATOR, no `site_id`: "Likely allosteric but unclear whether benzodiazepine or
  barbiturate site involved"
- propofol, GABA-A, POSITIVE ALLOSTERIC MODULATOR, no `site_id`: the same sentence
- zolpidem, GABA-A alpha-1/beta-1/gamma-2, POSITIVE ALLOSTERIC MODULATOR, `site_id` 2619: "Not clear if allosteric,
  but looks like benzodiazepine"
- gabapentin, voltage-gated calcium channel, MODULATOR, `site_id` 2657: "alpha2delta-1 is binding site: 23836974,
  23642658" — a site stated plainly with two PubMed identifiers, not an uncertainty

What the uncertainty costs is unequal. Acamprosate and propofol carry no `site_id` at all, so they are already
outside any rule that pairs drugs by site and their hedge changes nothing. Zolpidem's does cost something, and not
because of the hedge:

**`site_id` is scoped to a target record, so it misses a site two target records share.** The 20 rows on the record
"GABA-A receptor; anion channel" split into the GABA site (2616, meprobamate), the benzodiazepine site (2617, twelve
drugs) and the barbiturate site (2618, four drugs). Zolpidem's row sits on a different record, "GABA A receptor
alpha-1/beta-1/gamma-2", with `site_id` 2619, its own number for what the comment calls the benzodiazepine site. A
rule that pairs drugs on equal `site_id` therefore does not pair zolpidem with diazepam, although flumazenil
reverses both. Guide to PHARMACOLOGY's `Receptor Site` column does not have this defect, because it holds a site
*name*: 325 human rows carry one, over 8 names, and `[3H]zolpidem` is listed under "benzodiazepine site" beside
alprazolam, clonazepam, diazepam, flunitrazepam and flumazenil. Its coverage of the study's drugs is thinner than
ChEMBL's, so the two are complementary: ChEMBL for reach, Guide to PHARMACOLOGY for a site label that crosses
records.

### 14.2 The cooperativity factor: the binary is the wrong shape, and a measured value exists for no study drug

The owner: "I don't like the 0 or 1 simplification for coopertivity factor alpha. Is there a better way that actually
uses calculated or measured values?"

The objection is right twice over, and the second reason is the stronger one.

**The binary does not span the quantity.** In the allosteric ternary complex model alpha multiplies the orthosteric
ligand's affinity while the modulator is bound: below 1 the modulator weakens the ligand, above 1 it strengthens it,
and 0 is the competitive limit where the two cannot bind together. So alpha is a positive real number and its
interesting values lie on both sides of 1. `{0, 1}` offers the competitive limit and "no interaction" and nothing
else, which leaves out every case the parameter exists for. The measured values bear this out. ChEMBL holds 350
activities of `standard_type` "Cooperativity", and they run from 0.2 to 3000 with a median of 0.995, 175 above 1 and
175 below (`experiments/label_review/cooperativity_sources.py`). A parameter whose measurements are spread
symmetrically over three orders of magnitude around 1 cannot be replaced by a choice between 0 and 1.

**But no study drug has one.** Every cooperativity-style activity in ChEMBL was fetched, five `standard_type` values
in all, 535 activities over 215 molecules, cached in `data/raw/label_review/chembl_cooperativity/`:

| `standard_type` | activities | where |
| --- | --- | --- |
| Cooperativity | 350 | human muscarinic M1 to M5 only |
| Alpha | 45 | human histamine H1 (41), rat NMDA GluN1 (4) |
| Shift | 71 | rat GABA-A |
| Dose ratio | 59 | guinea pig, mouse and rat opioid receptors |
| Fold shift | 10 | ADMET assays (9), one human |

"Schild slope", "Fold-shift" and "alpha factor" hold 0 activities in all of ChEMBL. Eight of the 535 records touch
any of the study's 223 ChEMBL molecule ids, and not one of the eight is a cooperativity:

- diazepam at rat GABA-A, "Shift" = 1.0, from inhibition of [35S]-TBPS binding to rat brain membranes;
- naltrexone, five "Dose ratio" records on mu, kappa and delta receptors in guinea pig ileum and mouse vas deferens;
- morphine, two "Dose ratio" records in guinea pig ileum after an irreversible antagonist.

The dose ratios are Gaddum relations in animal tissue, which is evidence of competition rather than a value of alpha,
and the study is restricted to human measurements. Guide to PHARMACOLOGY adds nothing: one human row in 19,594 has
"cooperativity" in its `Action comment`, for an unnamed research compound at the beta-2 adrenoceptor, and the
`Ligand Context` column, which looked like a field for assay conditions, is topological (extracellular, cytosolic,
luminal) and holds 78 values in all.

Two routes that were considered and are not recommended. Deriving alpha from the overlap of PDBe-KB residue sets
would require a function from a shared-site fraction to a cooperativity, which no measurement calibrates; it is the
same error as the occupancy threshold grid the owner rejected in 13.8, a plausible number standing in for an
unmeasured one. Fitting alpha as a parameter of the model needs a label in which two drugs at one target vary with
dose, and the study has none.

**What is recommended instead: do not carry alpha at all.** A competition statement is made only for a pair that
shares a site, where alpha is 0 by the definition of a shared site, and the strength is then the Gaddum factor
1 + [B]/K_B, which the study already has from the affinity and exposure tables. A pair at different sites on one
target gets no competition term and is reported as such. That removes the parameter rather than guessing it, and the
cost is coverage of the different-site pairs, which no source quantifies for these drugs anyway. The scoping document
should say this in place of the binary.

### 14.3 What the measured concentration costs, and whether a mixed input is sound

The owner: "Why does the concentration per target cost 67 drugs? We can still keep the drugs that don't specify and
train them the old way. Or are you saying because this allows the off targets to include additional mechanisms, which
grows the leakage group?"

Measured by `experiments/label_review/measured_occupancy_leakage.py`.

**The number is 68, not 67, and it is pure coverage.** Of the 222 drugs, 183 have an affinity at a mechanism gene,
which is all the placeholder needs. 128 have a measured unbound concentration. 115 have both. So 68 drugs have an
affinity but no concentration and lose their reading. (My message of 12 October said 67; the decision log's 68 was
right.) Nothing about off-targets or leakage is in that number.

**The leakage hypothesis is wrong, and the truth is worse.** Admitting off-targets merges the leakage groups
catastrophically, and it does so under every scheme alike:

| scheme, threshold 0.25 | off-target pairs | drugs | largest leakage group | drugs in it |
| --- | --- | --- | --- | --- |
| no off-targets (registered rule) | 0 | 0 | 241 perturbations, 14.7% | 71 of 222 |
| placeholder C = 9 K_primary | 540 | 117 | 437 perturbations, 26.7% | 209 of 222 |
| measured C at each target's site | 374 | 81 | 434 perturbations, 26.5% | 207 of 222 |
| mixed, measured where there is one | 540 | 125 | 438 perturbations, 26.8% | 209 of 222 |

So the measured concentration does not grow the leakage group: it admits fewer off-targets (374 pairs against 540)
and leaves a group two drugs smaller. What grows the group is admitting off-targets at all, which takes 203 to 209 of
the 222 drugs into one group at every threshold from 0.1 to 0.5 and under every scheme. That is the argument against
off-targets entering a registered input, and it is independent of how C is computed.

**The mixed input is not free, and the reason is the one the owner already asked to control for.** "Train them the
old way" means a unit input per seed for the drugs without a concentration, so the two groups carry input on
different scales. Over the same seeds: the 115 drugs with a measured occupancy carry a total input with median 1.00
(quartiles 0.87 and 1.97, maximum 12.7), and the 107 drugs that would keep the unit rule carry a median of 3.0
(quartiles 1 and 10, maximum 23). Three times the input at the median, and which group a drug falls in is "does this
drug have a published therapeutic range and fraction unbound", which follows how long and how widely it has been
studied. That is the input-size confound of 13.1 reintroduced through the back door, and it correlates with label
count by the same route. The input-size baseline and the strata of 13.1 would measure it rather than remove it.

There is also a hole inside the covered drugs: among the 115, the number of seeds whose gene has no affinity, and so
no occupancy, has a median of 0 and a third quartile of 1, reaching 21 for one drug. A mixed rule has to say what
those seeds carry too.

**Recommendation.** If an occupancy input is wanted, use it for every drug or for none, and take the coverage loss
openly: 115 drugs with a reading, 107 without. A mixed column trades a stated 68-drug gap for an unstated confound,
which is the worse trade for a preregistered study. The question stays the owner's; nothing is built either way.

### 14.4 Options for representing tissue type

The owner: "The tissue type of a particular mechanism matters and if that isn't represented in the graph, I can't
calculate when the tissue type is fully saturated regardless of measuring per-class spread per propagation so I would
be relegated to only the steady state measurement. Can you please look into options for tissue type representation?"

Measured by `experiments/label_review/tissue_representation_options.py` on
`data/processed/graph_full_neuronal_split_binders` (49,574 nodes: 13,793 reaction, 12,810 gene, 12,709 protein, 8,580
metabolite, 1,681 protein entity, 1 membrane potential).

**What the graph has today.** One number per gene, `brain_median_tpm_max`, on 12,460 nodes, and the boolean
`brain_expressed` derived from it. No node carries a tissue. The `compartment` column is subcellular and sits on
24,055 nodes, none of which has a gene symbol.

**What is already on disk, and what it covers.** All three sources were downloaded on 10 October and only the first
is used, collapsed to one scalar:

| source | granularity | graph gene and protein nodes covered |
| --- | --- | --- |
| GTEx v10 gene median TPM | 68 tissues, 13 of them brain regions; 59,033 genes | 25,249 of 25,519 (98.9%) |
| HPA `rna_brain_region_hpa` | 13 brain regions; 20,162 genes | 25,248 (98.9%) |
| HPA `rna_single_nuclei_cluster_type` | 34 brain cell types; 19,580 genes | 25,104 (98.4%) |
| the node descriptor table's brain cell classes | 11 classes, continuous weights | as measured in 13.8 |

Matching is by gene symbol for almost all of it: only 5,667 nodes match GTEx by Ensembl id, where 25,247 match by
symbol. So a 68-dimensional tissue vector on nearly every gene and protein node needs no new download.

**The metabolic layer can inherit one.** A reaction node has no gene of its own but 7,759 of 13,793 carry a gene
reaction rule naming their catalysing genes, and 6,996 of 8,580 metabolite nodes (81.5%) take part in at least one
reaction that has such a rule. So a tissue weight can be propagated genes to reactions to metabolites, which is what
tissue-specific extraction from a Human-GEM model does. 6,034 reaction nodes and 1,584 metabolite nodes would still
have none.

**The finding that decides the saturation question.** A tissue partition of the nodes does not exist. Of the 25,271
matched genes, the median is expressed at 1 TPM or above in **67 of the 68 tissues**; 11,376 are expressed in all 68,
73.7% in more than half and only 615 in exactly one. The largest tissue takes a median of 7.4% of a gene's total
expression, against 1.5% for an even spread, so there is specificity in degree but not in kind. This is the same
shape as the 11 brain cell classes of 13.8, now at 68 tissues: "how much of tissue T has moved" is well defined as an
expression-weighted average over all nodes, and is not well defined as a fraction of a set of nodes belonging to T.
A weighted reading is buildable; a per-tissue saturation point in the sense of "every node of this tissue is at its
ceiling" is not.

**The part that is worth building anyway, and it is not the saturation reading.** The mechanism targets are far more
tissue-specific than genes at large. Over the 148 mechanism target genes, all 148 found in GTEx, the largest tissue
takes a median of 17.0% of the gene's expression, above 25% for 38 genes and above 50% for 11; the 13 brain regions
together take a median of 47%, above half for 69 genes and below a tenth for 26. Some are nearly single-tissue:
CACNA1S 95% skeletal muscle, CACNG1 93% skeletal muscle, SLC6A3 68% substantia nigra, DRD3 61% nucleus accumbens,
GABRA6 58% cerebellar hemisphere with a brain share of 1.00. Others are spread so widely that no tissue label would
say much: COMT 3%, SIGMAR1 4%, MAOA 5%.

This bears directly on the only place a tissue label touches a number today, the site a target's concentration is
read at (13.5). That choice is the boolean `brain_expressed`, a threshold on the maximum brain median TPM. Of the
148 mechanism targets, 131 are flagged brain-expressed, and **42 of those hold under a quarter of their expression in
brain**: ADORA3, the five alpha-adrenoceptors and ADRA2C, ALDH2, BCHE, five calcium channel subunits, CACNG1, CALCRL,
CHRM3, COMT, DRD4, FAAH, GABRE, HCRTR1, HRH1, HTR1B, HTR3A, HTR4, HTR7, four potassium channel subunits, MAOA, MAOB,
PTGS1, PTGS2, SCN7A, SCN9A, SIGMAR1, SLC18A2, TRPV1. One target runs the other way, HCRTR2, not flagged yet over half
in brain. A threshold answers "is any of it in the brain" where the share answers "is this mainly a brain target",
and those differ for 42 of 148 targets. Replacing the boolean with the GTEx share would change the site for them.

**So, three options, ordered by what they cost.**

1. **A tissue share on the mechanism targets, for the site a concentration is read at.** 148 genes, data on disk,
   changes no encoder and no input; it corrects the site for up to 42 targets. The limit is that a concentration can
   only be read where a tissue has a partition coefficient, and the study has Kp,uu for brain alone, on 89 drugs. So
   the honest form is a three-way choice, brain, plasma and unknown, rather than 68 concentrations.
2. **A 68-dimensional tissue weight on every gene and protein node, propagated to reactions and metabolites.** A
   measurement, reported per propagation step, which is what 13.8 recommended at 11 classes. It is a new readout and
   not a new model, so it does not add a confirmatory variant.
3. **A tissue-specific extraction of the metabolic layer**, a different graph per tissue. This is a different study:
   it changes the node set, so every registered count, the held-out set and the leakage groups would be redrawn.
   Not recommended before the confirmatory run.

### 14.5 Enzymes: the claim is already supported, and a curated table would not be what supports it

The owner: "I can't realistically claim physiological mechanistic pathways without enzymes and without carriers. I
could probably live without carriers but a lot of very important meds act on enzymes themselves … keep thinking on
that. Even if we have to curate our own enzymes data table, I think it would be worth it but I'm willing to hear a
dissenting opinion if you have one."

Here is the dissent, and it rests on a distinction that the exchange of 12 October ran together. "Enzymes" named four
different things, and three of them are already in the study. `experiments/label_review/enzymes_in_the_graph.py`
counts each, taking an enzyme to be a gene whose reviewed UniProt entry carries an EC number (4,613 gene symbols).

1. **Enzymes as the graph's own structure.** 3,825 of the 12,810 gene nodes are enzymes (29.9%), as are 3,815 of the
   12,709 protein nodes; 2,831 gene nodes are named in a Human-GEM gene reaction rule, so they catalyse a reaction
   the graph holds; and 7,759 of the 13,793 reaction nodes name their catalysing enzymes. The metabolic layer is
   made of enzyme-catalysed reactions. There is no version of this graph without enzymes.
2. **Enzymes as drug mechanism targets.** 11 of the 148 mechanism target genes are enzymes: ABAT, ACHE, ALDH2,
   ALDH5A1, BCHE, COMT, FAAH, MAOA, MAOB, PTGS1, PTGS2. 20 of the 222 drugs act on one: acetaminophen, aspirin,
   deprenyl, diflunisal, disulfiram, donepezil, entacapone, galantamine, isocarboxazid, moclobemide, neostigmine,
   opicapone, phenelzine, pyridostigmine, rasagiline, rivastigmine, tolcapone, tranylcypromine, valproate,
   vigabatrin. The owner's "a lot of very important meds act on enzymes themselves" is right, and those drugs train
   the mechanism edges like every other drug perturbation.
3. **Enzymes as labelled gene perturbations**, which is what trains the graph around an enzyme directly: **603 of
   the 1,663 labelled gene symbols are enzymes (36.3%), and 345 of those catalyse a reaction the graph holds.** Over
   a third of the study's gene labels are enzyme perturbations.
4. **The `clears` relation**, an enzyme that metabolises a drug. `configs/drug_clearing_enzymes.csv` gives 23 edges
   over 22 drugs through four enzymes: CYP2C19, CYP2C9, CYP2D6, NAT2. All four are gene nodes and all four carry an
   EC number. **None of the four has a labelled gene perturbation of its own**, which is why 13.7 found that nothing
   trains the relation.

A curated enzyme table would add to the fourth only. The first three need nothing: 345 labelled enzyme perturbations
on the metabolic layer, 20 drugs acting on enzyme targets and 7,759 enzyme-catalysed reactions already support the
sentence "mechanistic pathways including enzymes", and they can be stated with these counts. The fourth is drug
metabolism, that is pharmacokinetics rather than pathway physiology, and curating more rows of it would not make it
trainable, because what is missing is not the enzyme list but a perturbation that moves the enzyme and a symptom
recorded against it.

**The constructive form, which costs a curation the study can do.** The fourth becomes trainable if the four
clearing genes get labelled perturbations of their own, and such a label exists in kind: a poor-metaboliser genotype
is a loss of function of CYP2D6 or CYP2C19, so a symptom reported in poor metabolisers is a gene perturbation of
sign -1 with an outcome. That is exactly what the FDA Table of Pharmacogenetic Associations records, and it is thin
(23 statements, 12 of them on the development side, 75 development pairs against the 2,048 weights of the relation in
the message-passing encoder). ClinPGx and the CPIC guidelines are the larger source of the same kind and were not
retrieved. So the recommendation is: if the owner wants the `clears` relation trained, curate **CYP and NAT gene
perturbation labels**, not a drug-clearing enzyme table. The table already exists and is not the limit.

The carrier edge, which the owner said he could live without, is in the same position for the same reason, and 13.7
holds the measurement.

## 15. Answers of 12 October, 15:30

Seven points, each measured by a script in `experiments/label_review/`. Nothing here changes a label, a table, a
graph or a model input.

### 15.1 Binding sites matched by name, over both sources

The owner: "Why can't you match on name for ChEMBL? Would it be better to replace ChEMBL with Guide to
PHARMACOLOGY's receptor site?"

Nothing stopped it. `experiments/label_review/binding_site_names.py` reads the site out of ChEMBL's
`binding_site_comment` and out of Guide to PHARMACOLOGY's `Receptor Site`, normalises both into one vocabulary of
nine site names and matches on the name. That removes the defect of 14.1 at once: zolpidem joins the fourteen other
drugs at the benzodiazepine site, where its `site_id` 2619 had kept it apart.

The second question answers itself on coverage:

| source | drug-gene-site rows | drugs | sites |
| --- | --- | --- | --- |
| ChEMBL `binding_site_comment` | 317 | 20 | 4 |
| Guide to PHARMACOLOGY `Receptor Site` | 15 | 5 | 1 |
| either | | 21 | 4 |

ChEMBL reaches four times as many of the study's drugs, so it is not replaced. The two are used together as a union,
and Guide to PHARMACOLOGY earns its place twice over: it corroborates, since 11 drug-gene-site assignments are named
by both sources and every one agrees, and it supplies a label of its own, since flunitrazepam is named at the
benzodiazepine site there and nowhere in ChEMBL's comments, which is 1 drug and 4 rows the union would otherwise
lose. The owner's point stands and the script already takes the union; my earlier description of the second source as
corroboration alone was wrong. The sites, with the study's drugs at each: benzodiazepine (alprazolam,
clobazam, clonazepam, diazepam, estazolam, flunitrazepam, flurazepam, lorazepam, midazolam, oxazepam, quazepam,
temazepam, triazolam, zolpidem, zopiclone), barbiturate (methohexital, pentobarbital, primidone, secobarbital), GABA
(meprobamate) and alpha2delta-1 (gabapentin).

Two comments yield no site, and that is the right answer rather than a gap: acamprosate's and propofol's "Likely
allosteric but unclear whether benzodiazepine or barbiturate site involved" names two sites and settles neither, so
pairing on it would be a guess. A hedge about the *mechanism* is different and keeps its site, which is why
zolpidem's "Not clear if allosteric, but looks like benzodiazepine" still counts.

The table is `configs/binding_sites.csv`.

### 15.2 Competition without a cooperativity factor

The owner, on the recommendation of 14.2: "Do that then". So the rule is settled and written here for the
amendment: **a competition statement is made only for a pair of drugs named at the same site on one target**, where
the shift in occupancy is the Gaddum factor 1 + [B]/K_B that the affinity and exposure tables already give. A pair
named at different sites gets no competition term and is reported as such. No cooperativity factor is carried,
estimated or fitted.

What that covers, over the drugs that share a mechanism gene:

| | triples (drug, drug, gene) | drug pairs |
| --- | --- | --- |
| share a mechanism gene | 9,908 | 2,241 |
| both named at the same site: competition applies | 1,184 | 111 |
| named at different sites: no competition term | 932 | 79 |
| at least one site unknown | 7,792 | 2,070 |

### 15.3 What a refusal of the Lombardo permission would cost, which is less than it looks

The owner: "Have to allow off-targets. Unless we can't without the Lombardo clearances table and they deny
permission."

**The premise is wrong in the owner's favour.** Permission governs redistribution, not use. The repository already
downloads the workbook in the build and keeps the joined table out of git, which is what the stated terms permit
(13.9); a refusal would stop the table and any derived copy being published, and would leave every occupancy the
study computes exactly where it is. Off-targets are not waiting on the authors' answer.

The dependency measured anyway, from `data/processed/exposure/drug_exposure_v4.csv`:

| | with Lombardo | with the drug labels alone |
| --- | --- | --- |
| drugs with an unbound plasma concentration | 128 | 76 |
| drugs with an unbound brain concentration | 69 | 46 |

So the table is worth 52 drugs in plasma and 23 in brain. The cheapest way to narrow that further: 33 drugs have a
therapeutic range from Schulz but no fraction unbound from any source, and the 38 that have one from a label were
read by hand from the label sentences (`experiments/read_label_protein_binding.py`), so most of those 33 are
probably readable the same way.

### 15.4 What a per-tissue saturation reading would be

The owner: "what do you mean my saturation reading doesn't survive in the strict form?"

A correction of my own wording. What 14.4 rules out is the literal reading only: "tissue T is fully saturated"
taken as "every node belonging to T is at its ceiling". There is no set of nodes that is T, because the median gene
is expressed at 1 TPM or above in 67 of GTEx's 68 tissues, so almost every node belongs to almost every tissue and
68 such curves would be near-copies of the global one.

The weighted reading is well defined and is what the owner should build:

    saturation of tissue T at step t  =  sum over nodes v of w(v, T) * x_v(t)  /  sum over nodes v of w(v, T)

with x_v(t) how far node v has moved from its unperturbed state at propagation step t, and w(v, T) the share of
v's expression that sits in T. Those curves do separate, because the weights are concentrated where membership is
not: the largest tissue takes a median 7.4 percent of a gene's expression against 1.5 percent for an even spread.
The peak per tissue is then a well-defined step and the owner's "work backwards to find the peak distribution point
for each tissue type" applies to it unchanged. How far the 68 curves separate is a measurement to make before
relying on it, the same measurement recommended in 13.8 for the 11 brain cell classes, where the median pairwise
correlation was 0.36.

The owner also accepted the narrower option as a fallback: "The only splitting out brain tissue option is a
reasonable fallback." That is option 1 of 14.4, the GTEx brain share in place of the `brain_expressed` boolean,
which would change the site for up to 42 of the 148 mechanism targets.

### 15.5 Enzymes: a correction, and where a CYP and NAT table would come from

**Correction of the owner's reading, not of the record.** "I don't see anything acting on AChE in that list" — ACHE
is in it. ACHE and BCHE are two of the 11 enzyme mechanism targets of 14.5, and five study drugs act on them:
donepezil, galantamine, neostigmine and pyridostigmine on ACHE, rivastigmine on ACHE and BCHE. Each of the five has
ACHE as its only mechanism target, so nothing dilutes them. The four-name list that holds no cholinesterase is the
`clears` relation (CYP2C19, CYP2C9, CYP2D6, NAT2), which is drug metabolism.

**A table for CYP and NAT: the edges are easy and the labels are not.** The Clinical Pharmacogenetics
Implementation Consortium publishes its database through an API at `api.cpicpgx.org`, which its robots file permits
(the web pages at `www.clinpgx.org` do not; see "Refused" below). Queried on 12 October 2026:

- `v1/pair` holds 138 gene-drug pairs for the four genes that are not marked removed: CYP2D6 79, CYP2C19 27,
  CYP2C9 25, NAT2 7, with a CPIC evidence level, a ClinPGx level, a testing category and PubMed citations.
- **56 of those pairs name a study drug, over 45 of the study's 222 drugs**, against the 23 edges over 22 drugs that
  the FDA Table of Pharmacogenetic Associations gives. A curated clearance table from this source would more than
  double the relation's edges, and with a graded evidence level per edge.
- But it carries no outcome. `v1/recommendation` holds 2,115 rows, 1,640 of them naming one of the four genes and
  1,150 of those on a study drug over 15 drugs; **none of the 1,150 names any symptom term of the study's crosswalk**.
  CPIC recommendations are dosing advice ("use standard dosing", "avoid", "consider an alternative"), not symptom
  frequencies.

So the position is sharper than it was. Curating a CYP and NAT table is worth doing and is nearly free, and it would
grow the `clears` relation from 23 edges to about 56 with evidence levels; it would not make the relation trainable,
because what is missing is a symptom recorded against a metaboliser phenotype, and the FDA table's 23 statements
remain the only source of that kind found. Training would need either primary-literature curation of
poor-metaboliser outcomes or a cohort carrying genotypes and symptoms together.

**The CPIC licence is verified and it is permissive**, which the Lombardo table is not: "All curated content
published by CPIC is available free of restriction under the CC0 1.0 Universal (CC0 1.0) Public Domain Dedication",
read in `LICENSE.md` of `github.com/cpicpgx/cpic-data` at commit `82aa2518f1051bb2db9c42c1023a9f08db5a5f80` on
12 October 2026, because the licence page on the website declines automated readers. The same file asks for
attribution to CPIC, the primary source `www.clinpgx.org`, the relevant publications, and for the URL, date accessed
and version number when content comes from the API. One thing the owner has to open personally: item 6 binds the
database and API to a further data-usage policy at `www.clinpgx.org/page/dataUsagePolicy`, which is behind the same
robots rule and was not read.

**Refused, and respected.** `www.clinpgx.org` disallows automated readers for its pages, so the CPIC licensing page
(`cpicpgx.org/?p=1245`, which redirects there), the API documentation page (`?p=2871`), the data-usage policy and the
ClinPGx downloads page were all left unread. `api.pharmgkb.org` is unreachable from this workspace (the proxy returns
502). No cache, mirror or other user agent was tried. The repository of the data model,
`github.com/cpicpgx/cpic-data`, carries agent-instruction files of its own; they were not read, since a third party's
instructions to an assistant have no standing here.

### 15.6 Placebo arm sizes, and the weight they give

The owner: "can you please look for placebo arm sizes for those specific studies and then calculate negative weight
based on that then use .5 as default for nothing found?"

Where the sizes are. SIDER gives frequencies with no arm sizes and OnSIDES v3.1.1 gives neither, so they are read
from each drug's own label, whose adverse-reaction tables head their columns "Placebo (n=385)".
`experiments/label_review/placebo_arm_sizes.py` queries openFDA, parses the table headers, falls back to the
adverse-reactions text and caches every label under `data/raw/label_review/placebo_arm_sizes/`.

The weight. With n_drug and n_placebo the effective size of a two-arm comparison is
n_effective = 1 / (1/n_drug + 1/n_placebo), which is the quantity the standard error of a difference in proportions
turns on. Then

    weight = n_effective / (n_effective + n_reference),   capped at 1.0,   and 0.5 where no size was found

with n_reference the median effective size over the pairs that have one, 214. That makes the owner's 0.5 both the
default and the weight of a comparison of typical size, so a pair with no sizes found is not an outlier: a larger
comparison earns more and a smaller one less, and nothing reaches the weight of a grade A positive.

What it gives over the 61 pairs:

| | pairs | drugs |
| --- | --- | --- |
| total | 61 | 39 |
| arm sizes from an adverse-reactions table | 13 | 6 |
| arm sizes from the adverse-reactions text | 28 | 19 |
| no arm size stated in the label | 17 | 11 |
| no United States label found | 3 | 3 (pergolide, reboxetine, rotigotine) |

41 of the 61 pairs get a measured weight and 20 take the default. The weights run from 0.057 to 0.772, median 0.500,
mean 0.467. The extremes are worth seeing, because they are what the flat 0.5 was hiding: pimozide's pairs rest on 20
placebo and 36 drug patients and earn 0.057, while fluoxetine's rest on 975 and 2,869 and earn 0.772. The table is
`configs/placebo_negative_weights.csv`, one row per pair with its arm sizes, its source and the caption of the table
the sizes came from.

This needs a flag that reads a weight per pair rather than one number, which `--placebo-negative-weight` does not
have yet; it is not built, because the owner's instruction against new pilots stands and the amendment is not
applied.

### 15.7 The target limit: one drug, not a count

The owner: "I really don't like having a target limit. Can you present another option if leakage changes are
needed?"

**What the current rule actually is, which matters.** `has_dominant_target` is an admission test: a drug with more
than three distinct ChEMBL mechanism targets is dropped from the study, not truncated to three. Its own docstring
says "Affinity margins are not yet used", so the rule is named after a dominance it does not measure. Eight drugs
are excluded by it: desflurane, enflurane, isoflurane, nefazodone, orphenadrine, sevoflurane, topiramate,
trimipramine.

`experiments/label_review/target_rule_options.py` compares the alternatives on the approved table with no limit (230
drugs), reporting the largest leakage group under the registered chained rule:

| rule | drugs admitted | largest group |
| --- | --- | --- |
| no limit | 230 | 194 drugs (84%) |
| three targets (the current rule) | 222 | 71 drugs (31%) |
| targets within 0.5 to 2.0 log units of the drug's best-bound gene | 187 | 59 to 60 drugs (26%) |
| the same pruning, every drug admitted | 230 | 171 to 190 drugs |
| only targets a curator flags as the mechanism (DrugCentral MOA, GtoPdb primary) | 213 | 170 drugs (74%) |

A measured dominance rule does beat the count on leakage, 59 drugs against 71, but it admits 35 fewer drugs, because
only 187 of the 230 have a measured affinity at a seeded gene. Pruning without an admission test keeps all 230 and
loses the split. The curator rule is the worst of both.

**But none of that is the answer, because the problem is one drug.** Removing one drug at a time from the no-limit
table and asking how large the biggest group is without it:

| drug removed | drugs the largest group loses |
| --- | --- |
| orphenadrine | 117 |
| tapentadol | 19 |
| agomelatine | 4 |
| topiramate, levetiracetam, fenfluramine, betahistine | 2 each |

With **orphenadrine alone set aside, and no limit on targets at all, the largest group holds 77 drugs** — against 71
under the three-target rule, which excludes eight drugs to get there. Setting aside tapentadol and agomelatine as
well changes nothing further; a fourth drug takes it to 75.

Why orphenadrine. It carries 18 target genes: the NMDA receptor subunits GRIN1, GRIN2A to GRIN2D and GRIN3A, ten
sodium channel subunits, HRH1 and SLC6A2. The NMDA and sodium-channel block and the monoamine block are otherwise
separate, and orphenadrine is in both. It holds 6 evidence rows over 6 symptoms, so setting it aside costs 6 rows
where the count rule costs eight drugs.

**The option, then, stated so that it can be registered before anything is run.** No limit on a drug's targets. A
drug is set aside when its target set joins two leakage blocks that are otherwise separate, which is measured on the
evidence table by removing each drug in turn, and the drugs set aside are named in the registration with their
measured effect. On the current table that rule names orphenadrine and nothing else. It is a rule about the
*structure of the split* rather than about pharmacology, which is the honest way to describe what the three-target
count was doing anyway, and it keeps 229 drugs instead of 222.

Two things the owner should weigh against it. The rule is defined by its effect on the split, so it has to be fixed
and published before any model is run, or it becomes a choice made on results. And 77 is worse than 71: the question
is whether seven more drugs are worth six more drugs in the largest group.

### 15.8 Enzyme descriptor columns, and the enzyme-targeting drugs

The owner: "Are there enzyme descriptor columns? If so, how many? I'll probably want the 11 enzyme targeting drugs in
the training actually even though it's not much training data."

**There are eight.** The node descriptor table `full_neuronal_split_descriptors_brain_expression.parquet` holds 138
columns over the 49,574 nodes, and eight of them are enzymatic: `reaction_has_ec` and `reaction_ec_class_1` to
`reaction_ec_class_7`, the seven top-level Enzyme Commission classes. 5,780 of the 13,793 reaction nodes carry an EC
number, split oxidoreductases 1,796, transferases 1,990, hydrolases 1,225, lyases 371, isomerases 145, ligases 341
and translocases 155. The other descriptor groups are 66 protein columns, 35 reaction columns in all, 27 gene and 10
metabolite.

So the model already sees what kind of enzyme catalyses each reaction, as one of seven indicators plus a presence
flag. What it does not see is an enzyme's identity beyond its protein descriptor, and the EC class sits on the
reaction rather than on the gene.

**The enzyme-targeting drugs are already in the training set.** The 11 of 14.5 are target genes, not drugs; the drugs
that act on them are 20, and nothing excludes any of them. They are ordinary drug perturbations: acetaminophen,
aspirin, deprenyl, diflunisal, disulfiram, donepezil, entacapone, galantamine, isocarboxazid, moclobemide,
neostigmine, opicapone, phenelzine, pyridostigmine, rasagiline, rivastigmine, tolcapone, tranylcypromine, valproate
and vigabatrin, and they train the mechanism edges like every other drug. Nothing has to be added for them. If the
owner meant the 23 enzyme *statements* of the FDA pharmacogenetic table, that is the separate decision of 13.7,
where they are reported as checks and do not train.

### 15.9 The per-tissue saturation reading, built

The owner: "For saturation, sure that works. Do it." Built as `experiments/measure_tissue_saturation.py`, reported in
`docs/tissue_saturation.md`. It is a structural measurement with every relation's gain at 1, because the gains are
learned and no model has been trained on this table, so peak steps read off untrained gains would measure the
initialisation.

Weights, which is where the coverage sits: a gene or protein node takes its gene's share of its own GTEx median
expression per tissue (25,180 nodes), a reaction node the mean over the genes of its Human-GEM reaction rule (7,755),
a metabolite node the mean over the reactions it touches (7,006). 9,633 nodes carry no weight and enter no tissue's
mean.

**Which reading has a peak.** Three were computed and only one is informative. Under the registered recurrence the
saturation *level* rises monotonically to a fixed point, so every tissue "peaks" at the last step computed and the
curves correlate at 0.998; that is the dynamics, not a finding. Under an impulse with a mean-aggregating adjacency
the field contracts from the first step, so 67 of 68 tissues peak at step 0. The quantity that does have a peak under
the dynamics the study uses is the **increment**, how much of a tissue fills at each step.

**What it shows**, over the 217 drug perturbations with 24 steps at the registered damping of 0.5:

- the step of fastest filling runs from 1 to 11 with a median of 2;
- the normalised curves separate, though not sharply: median pairwise correlation 0.901, against 0.985 for an
  impulse and effectively 1 for the level curves;
- the order is physiological. The brain tissues fill first (amygdala, frontal cortex BA9, cortex and cerebellum at
  step 1), then kidney cortex, liver and skin at 6 and 7, and whole blood at 10 with EBV-transformed lymphocytes at
  11;
- 66 of the 68 tissues fill fastest within the registered 8 propagation steps; only whole blood and lymphocytes do
  not. **But the level is not settled there**: at step 8 a tissue has reached a median of 0.727 of its settled level,
  the lowest being whole blood at 0.393 and the highest 0.953. So the window is not short of the fastest step, it is
  short of the steady state, which is the concern recorded in the architecture review.

**A correction of the first run of this measurement, made the same day.** It gave every relation a gain of 1. With an
adjacency whose rows sum to one that makes the transition non-contracting, so the iteration has no fixed point and
grows without limit, and the peak steps it produced (1 to 20, 26 tissues past step 8, correlation 0.183) were partly
reading the divergence. The gain is now the encoder's contraction bound of 0.9, which is the fastest spread a
registered model can have and does converge. The ordering survived the correction; the spread of peak steps and the
separation of the curves did not.

Nothing is changed by this: it is a measurement, and whether the propagation window moves is a registered decision
the owner has not made.

## 16. Answers of 12 October, 16:00

### 16.1 The propagation window: what the record says

The owner: "I also thought we did away with only having 8 propagation steps, opting instead for convergence as steady
state reading."

The record does not carry that decision. The decision log for 10 October 2026 holds it as **"Open with the user, not
built"**, from the exchange "Why aren't we doing message passing in a way that carries the signal until all sink
nodes are hit?" and "wouldn't convergence be preferred for physiological correctness?". I agreed there that
convergence is the better target, because the labels describe conditions that stay on, and the probes recorded what
it would need: the perturbation held on at every round, since the encoder gives it once and a stable shared layer
then forgets it (norm 0.15 at round 8 against 2e-6 at round 90, while a held one settles at 3.85 with all 16
perturbations settled by a median of round 33), and stability by construction, since at the registered scale both
encoders grow without limit. It was proposed as a **third encoder** — held perturbation, shared relation laws,
time-scale groups, a damped and bounded update, iteration to tolerance — and the row ends "the user's decision
whether to build it and whether beside or in place of message passing". Nothing was built, and both registered
encoders still run a fixed number of steps: message passing 3 edges, the linear response 8.

Two numbers that look contradictory and are not. That same row recorded the linear response at 8 steps as "within
0.14 percent of its fixed point at the median", and 15.9 reports a tissue at a median of 0.727 of its settled level
at step 8. The first is with the encoder's own gains, which sit well below the bound; the second is at the
contraction bound of 0.9, the fastest spread a registered model may have. How short 8 steps is depends entirely on
the gains, which is why 16.3 matters.

### 16.2 Reopening 13.7: what training the enzyme statements would take

The owner: "Yes, let's reopen 13.7."

**A correction of my own analysis first.** 13.7 said there was "not enough to train on either" encoder. That is wrong
for one of them, and the difference is large enough to change the answer.

| encoder | parameters the `clears` relation owns | development pairs per parameter |
| --- | --- | --- |
| message passing (`node_state_dim` 32, **3 layers**, separate weights per layer) | **3,072** | 0.024 |
| linear response (`propagation_channels` 4, one gain per channel) | **4** | **18.8** |

(Corrected on the owner's challenge: the registered encoder has three message-passing layers, not two. The weight
tensor is `[layers, relations, 32, 32]`, so a relation owns 3 x 32 x 32 = 3,072 weights. lockbox_v1 had two layers
and lockbox_v2 has three; I had carried the older figure. The conclusion is unchanged and the ratio is worse.)

There are 23 enzyme statements over 22 drugs, 12 of them on the development side with 75 development pairs. Against
3,072 weights that is hopeless and the earlier conclusion stands. Against four bounded gains it is thin but ordinary:
the linear-response encoder's gain for a signed relation is `sigmoid(...)` per channel, four numbers in (0, 0.9), and
75 labelled pairs can identify four bounded parameters.

So reopening has a definite shape rather than a yes or no:

1. **What must be built.** The combined perturbation type, a drug and a gene perturbation carried together in one
   example, which rows 4 to 6 of section 8 named and which nothing in the trainer has. Without it the statements
   cannot enter the loss at all, whatever the parameter count.
2. **Where it can train.** The linear-response encoder only. Training the same statements in message passing would
   fit 2,048 weights on 75 pairs.
3. **What that costs.** The two confirmatory models would then differ in what they are trained on, which is a design
   asymmetry the owner has been avoiding ("I don't want to add extra confirmatory model variations because that
   starts increasing probability of winning by luck"). It is also a change to a registered model, so it belongs in
   the amendment and not in a commit.
4. **What would make it stronger, and is nearly free.** The relation currently has 23 edges from the FDA table. CPIC
   gives 56 pairs on study drugs over 45 drugs for the same four genes (15.5), each with an evidence level and
   citations, under a verified CC0 licence. More edges do not add labels, but they give the same four gains more
   structure to act through, and the edge table is the part that can be built now without touching a model.
5. **What stays true either way.** As *checks*, the carrier, displacement and enzyme kinds still test nothing while
   no label varies those relations, because the check reads an untrained relation. Training the enzyme statements
   would make the enzyme checks meaningful and would leave the carrier checks exactly as they are, since no label
   varies a carrier.

Recommended, and the owner's to decide: build the combined perturbation type and the CPIC edge table, hold the
training decision until the owner has seen what a development run on 12 statements does, and keep the two
confirmatory models symmetric unless the owner accepts the asymmetry knowingly.

### 16.3 Gains: no source has them, and what exists instead

The owner: "Can you find gain numbers to train in?"

**No source tabulates a per-relation gain, and the reason is structural rather than a gap in the databases.** The
model's gain is one number per relation per channel over 20 relations, not a number per edge; it is a property of the
model's own parameterisation, not a measured constant of biology. Nothing measures "how much an activation edge
transmits on average", because that quantity only exists inside this model.

What does exist falls into three kinds, and only the second is the right kind.

**1. Per-edge confidence, which is the wrong kind.** SIGNOR carries a score per interaction, 0.2 to 0.981, over
43,616 human rows; 21,180 are protein to protein and 15,608 of those are signed. Against the confirmatory graph it
covers **7,640 of the 69,740 `activates` edges (11.0 percent) and 4,029 of the 18,157 `inhibits` edges (22.2
percent)** with a score of the same sign. OmniPath, already on disk, carries `curation_effort` and consensus flags,
and the CollecTRI table as downloaded carries no weight at all. Every one of these is a confidence: how much evidence
the curators found, not how much signal passes. Using a confidence as a magnitude is the error the owner rejected
twice already, in the occupancy placeholder and in the threshold grid, and I am not going to recommend it a third
time.

**2. Per-edge measured effect sizes, which is the right kind and is out of reach today.** A transcription factor
knocked down, with the log fold change of each target, is a measured transmission magnitude for
`regulates_transcription_of`; the same for a perturbation compendium over signalling edges. The two sources are
KnockTF 2.0 and LINCS L1000, and all three addresses tried refuse automated readers:
`bio.liclab.net/KnockTFv2` serves no robots file at all, which this project's reader takes as a refusal, and
`clue.io/data/CMap2020` and the GEO series directory for GSE92742 both disallow it. Recorded, respected and left for
the owner to open. If one of them is retrieved, the per-relation gain it would justify is an aggregate over edges,
and an aggregate is what the model wants, so this is the route worth the owner's time.

**3. Theory, which gives a constraint rather than a value.** For the metabolic relations the quantity that
corresponds to a gain is a flux control coefficient, and metabolic control analysis gives the summation theorem:
over an enzyme set the control coefficients sum to one. That is a constraint on the gains rather than a number to
look up, and it is the kind of thing the study already does elsewhere through the contraction bound. For `encodes`
the measured analogue is the genome-wide ratio of protein to messenger RNA, which is published and would give that
one relation an initialisation.

**What I would do instead of looking for gains.** The gains are four numbers per relation in the linear-response
encoder, 80 in all, and they are what the development runs are for. If the worry is that the development labels
under-determine them, the useful moves are to constrain rather than to import: keep the contraction bound, add the
summation constraint on the metabolic relations if it is wanted, and note in the amendment that the gains are fitted
and not measured. If per-edge weights are wanted, that is a different model — an edge-weighted propagation rather
than a relation-gain one — and a registered change rather than an initialisation.

**Licence, verified.** "SIGNOR 4.0 is licensed under Creative Commons Attribution 4.0 International"
(`signor.uniroma2.it/documentation`, read 12 October 2026). Nothing from SIGNOR is used by any build; it was read to
measure the coverage above.

## 17. Answers of 12 October, 16:26

### 17.1 Two corrections the owner caught

**The message-passing encoder has three layers, not two.** `num_message_passing_layers` is 3 and the weight tensor is
`[layers, relations, node_state_dim, node_state_dim]`, so a relation owns 3 x 32 x 32 = 3,072 weights, not the 2,048
I wrote in 16.2. The architecture review already recorded it as "Message passing stops at 3 edges. Three layers,
separate weights per layer", and `configs` moved from two layers under lockbox_v1 to three under lockbox_v2. I was
carrying the v1 figure. 16.2 is corrected; its conclusion is unchanged and the ratio against 75 development pairs is
worse, 0.024 pairs per parameter rather than 0.04.

**The propagation channels are not the descriptors.** The owner: "for the clears relation, shouldn't there only be 8
channels since there are only 8 enzyme descriptors?" They are different things that never meet.

- `propagation_channels` is 4 by default and is a property of the linear-response encoder, shared by every relation.
  The channels are time scales: each carries its own per-relation gain so the response settles at different rates,
  and a linear map expands them to the `node_state_dim` channels the heads read. The number 4 has nothing to do with
  any descriptor and does not vary by relation, so `clears` owns four gains because the encoder propagates four
  channels.
- The eight enzyme descriptors are node features. They enter through the descriptor projection into the node state,
  not through the propagation at all.
- And they are not even present at a `clears` edge. Measured: the EC descriptors are non-zero on 5,780 nodes and
  **every one of them is a reaction node**. A `clears` edge runs from an enzyme's gene node to a drug entry node, and
  none of CYP2D6, CYP2C19, CYP2C9 or NAT2 carries a non-zero EC descriptor on its gene or protein node.

### 17.2 Convergence in place of 8 steps

The owner: "Don't make it a third encoder. Use the convergence instead of 8 steps. 8 steps is arbitrary, not
physiological."

**For the linear response this is a change of stopping rule, not a new encoder**, because the encoder already has the
two properties convergence needs. It holds the perturbation on at every step, since the recurrence adds the sustained
input each round, and its gains are bounded by the contraction, so the transition is a contraction and a fixed point
exists. Iterating to a tolerance is then a one-line change to the loop.

What it costs, measured by `experiments/measure_propagation_convergence.py` over the 217 drug perturbations at the
registered damping of 0.5, with every relation given the same gain:

| gain | steps to a relative change below 1e-2 | below 1e-3 | below 1e-4 | field at 8 steps, by norm, over the settled field |
| --- | --- | --- | --- | --- |
| 0.9, the contraction bound | 10 | 23 | 45 | 0.993 |
| 0.7 | 7 | 14 | 22 | 0.999 |
| 0.5 | 5 | 11 | 16 | 1.000 |
| 0.3 | 3 | 8 | 12 | 1.000 |

So the honest cost is three to six times the propagation of a forward pass at the bound, and less at the gains a
trained model is likely to hold.

**The owner's objection is right, and the norm is what hides it.** At 8 steps the field is 99.3 percent of the
settled field *by norm*, which is the figure recorded on 10 October and which makes 8 steps look sufficient. By
distribution it is not: against a field run to 200 steps, a tissue has reached a median of **0.631** of its settled
saturation at step 8, the lowest being whole blood at 0.287 and the highest 0.934 (15.9, rerun against a settled
reference; the earlier 0.727 used a 24-step reference that was not itself settled). The norm is carried by the large
values near the seeds, which settle immediately; the far field, which is where a claim about pathway reach lives, is
between a quarter and two thirds of the way there. Both numbers are correct and they answer different questions,
and the question this study asks is the second one.

**For message passing it is not a stopping rule, and this is where "don't make it a third encoder" has to be
spent.** That encoder writes the perturbation once into the initial state and gives every layer its own weights, so
there is nothing to iterate: more rounds means more parameters. Each layer is `relations x 32 x 32`, which at 20
relations is 20,480 weights, so running to a depth that reaches the far field would multiply the parameter count
against 2,975 labelled pairs. Convergence there needs the perturbation held on and the weights shared across rounds,
which is what the owner proposed on 10 October himself: "something more creative like making a schema that reuses
some params across layers and not others or having multiple layer param sets that get collated".

So the instruction can be met without a third encoder, by changing the two that exist rather than adding one:

1. **Linear response**: iterate to a tolerance in place of 8 steps. The tolerance is the registered quantity and the
   step count follows from the learned gains. 1e-3 is the defensible choice, costing 23 steps at the bound.
2. **Message passing**: hold the perturbation on at every round, share the relation weights across rounds, and
   iterate to the same tolerance with the same bounded update. This keeps the parameter count at one layer's worth
   and makes depth free, and it replaces the registered message-passing encoder rather than joining it, so there are
   still two confirmatory models.

Both are changes to registered models, so both belong in the amendment and neither is applied. Two things the owner
should weigh. The second changes what message passing *is*: shared weights across rounds is a different hypothesis
class from three independent layers, and the head choice and stopping rule were selected on the current one, so the
development runs that chose them would have to be redone. And a convergent model has no "3 edges" story, so the
disclosure in the preregistration about reach would be rewritten rather than amended.

## 18. The message-passing encoder's parameter count, and whether to keep it

The owner, 12 October 2026: "there's gotta be a way to parameterize that up-projection so that it doesn't require
that many params. Otherwise, I think I'm gonna drop the MP encoder variant because its just getting too parameter
hungry."

**The instinct is right and the number is 14 to 1.** `experiments/count_encoder_parameters.py`, on the confirmatory
graph's 20 relations with the registered node state of 32, three layers and four propagation channels, against the
2,975 kept pairs:

| encoder | parameters | per labelled pair | where they are |
| --- | --- | --- | --- |
| message passing | 69,184 | 23.3 | relation transforms 61,440 (89 percent) |
| linear response | 4,816 | 1.6 | the shared feature projection 4,448 (92 percent); relation gains are 80 |

Almost all of the difference is one tensor, `[layers, relations, 32, 32]`. The linear response spends nearly nothing
on relations: 80 numbers, four per relation.

**There is a standard way to parameterise it down, and it is the method the architecture was named for.**
Schlichtkrull et al. 2018 introduced basis decomposition for exactly this problem, a relational graph network whose
parameters grow with the number of relations. Each relation's matrix becomes a mixture of B shared bases,
W_r = sum_b a_rb V_b, so the bases carry the capacity and a relation owns only its B coefficients. Priced on this
graph, with the matrices also shared across layers, which convergence needs anyway (17.2):

| scheme | relation transforms | one relation owns | per labelled pair |
| --- | --- | --- | --- |
| registered: a full matrix per relation per layer | 61,440 | 3,072 | 20.65 |
| the same matrices shared across layers | 20,480 | 1,024 | 6.88 |
| **basis decomposition, 4 bases, shared across layers** | **4,176** | **4** | **1.40** |
| basis decomposition, 8 bases | 8,352 | 8 | 2.81 |
| block diagonal, 4 blocks | 5,120 | 256 | 1.72 |
| diagonal, shared across layers | 640 | 32 | 0.22 |

Basis decomposition at four bases cuts the relation transforms by a factor of 15 and brings the whole encoder to
about 9,800 parameters, 3.3 per labelled pair, against the linear response's 1.6. Two further things fall out of it:

- **It removes the asymmetry of 16.2.** Under four bases the `clears` relation owns four coefficients, the same count
  as the linear-response gain. The enzyme statements would then be as trainable in message passing as in the linear
  response, and the reason for training them in one model only disappears.
- **It is defensible in a preregistration**, being a published method for this exact problem rather than a choice
  made to fit this dataset.

**Where the spectrum ends matters for whether two models are worth running.** The last row is not a reduction of
message passing, it is a change of identity: a diagonal relation weight at 32 channels *is* one gain per relation per
channel, which is what the linear-response encoder already does. Shrink that far and the two confirmatory models
differ only in the nonlinearity and in the sign treatment. That is not necessarily bad — a pair of models that differ
in one thing is a sharper comparison than a pair that differ in a dozen, and the registered choice was "Linear + one
nonlinear" — but it should be chosen knowingly rather than arrived at by shrinking.

**If the owner drops message passing instead.** What is lost is the only nonlinear model, so nothing in the study
tests whether nonlinearity matters, and the owner's own choice of 10 October of one linear and one nonlinear model is
abandoned. What is gained is half the compute and no architecture-selection surface; the multiplicity itself is
already controlled by the one-sided sequence of H1 then H2, so two models were never needed for rigour. Dropping it
is a defensible study, and a smaller one.

**Recommended, and the owner's to decide.** Basis decomposition at four bases with the matrices shared across layers.
It answers the parameter objection, it is what convergence needs anyway, it removes the enzyme asymmetry and it keeps
a nonlinear arm. Two costs to weigh: it is a registered change, so the head choice and the stopping rule would have
to be selected again on development runs, which the instruction "No new pilots yet" currently blocks; and basis
decomposition assumes the relations share a low-dimensional space of transforms, which is an assumption about
capacity that the development runs have to test. (The reason first given here, that `activates` and `inhibits`
plausibly differ mainly by sign, was withdrawn on the owner's challenge; see 19.3, which gives the argument that does
hold.)

A fourth option, not recommended on its own: halve the node state from 32 to 16, which quarters the relation
transforms to 15,360. It costs capacity everywhere, including in the heads, where the parameters are not the problem.

## 19. Answers of 12 October, 16:51

### 19.1 Two things I said would need redoing, and neither does

**The head choice is settled and I should not have reopened it.** The owner: "Head choice is never recoded. Only the
noisy or is defensible for mechanistic pathways due to the stoichiometric requirements of the graph." The decision
log for 8 October 2026 records exactly that: "Sigmoid head dropped from the confirmatory family: the noisy-OR head of
each encoder is tested", with the reason given as "only the noisy-OR head has modules to read as pathways". So both
confirmatory models already use the noisy-OR head, the four-model family of the preregistration is down to two, and
there is no head choice left for a change of encoder to invalidate. My statement in 17.2 and 18 that it would have to
be selected again was wrong.

**The training stopping rule is a rule, not a fitted value.** The preregistration fixes it: validation loss with
patience 8, at most 60 epochs. Those numbers were set in advance and are not read off a development run, so a change
of encoder does not invalidate them. A smaller model will stop at a different epoch, which is an outcome and not a
parameter. The owner's reading — "Oh you're just saying the max number of steps needs recalculated?" — is nearly
right, with one distinction: the 60-epoch cap does not need recalculating either. What changes is the **propagation**
budget, which under 17.2 stops being a step count at all and becomes a tolerance, and the step count then follows
from the learned gains rather than being registered.

So the cost of changing the encoder is development runs for their own sake, not for re-selecting a head or a stopping
rule.

### 19.2 One layer, repeated: yes, with one requirement

The owner: "use the 4 bases and I think just the 1 layer, repeated. Does the one repeated layer make sense to do?"

It does, and it is the same object as the convergence of 17.2 rather than a second change. A single set of relation
transforms applied repeatedly is a fixed-point iteration, h <- f(h, u) until the change falls below a tolerance,
which is what a steady-state reading means. Three things follow.

**The parameter count stops growing with reach.** At four bases shared across the rounds the relation transforms are
4,176 whatever the depth, so there is no longer a trade between how far a perturbation travels and how many
parameters the model holds. That is the whole of the owner's objection answered: depth becomes free.

**It needs an explicit contraction bound, and this is the one requirement.** A linear iteration converges when the
transition's spectral radius is below one, which is what the linear-response encoder's contraction of 0.9 already
buys it. With a nonlinearity the same argument works only if the activation is non-expansive, which tanh, ReLU and
the sigmoid all are, and if the linear part is a contraction. Under basis decomposition the linear part is
sum_b a_rb V_b, so the bound has to be placed on the bases and the coefficients together — normalising each basis and
bounding the coefficients is the direct way, and it is the same device the codebase already uses for the gains. Built
without that bound, a weight-tied iteration can diverge or oscillate, and the 10 October probe recorded that at the
registered scale both encoders grow without limit, which is what that failure looks like. This belongs in the
amendment as a stated constraint, not as an implementation detail.

**What is lost is depth-specialised transforms, and it is a real loss.** With three independent layers the model can
use the first round for one kind of step and later rounds for another, and the architecture review records that it
does: "A knockout spends one layer on gene to protein, so it reaches 2 edges past its protein". A repeated layer
applies the same transform at every round, so that specialisation is not available; gene-to-protein and
signalling-to-signalling steps must share one transform per relation. Whether it matters is empirical and would show
in a development run. It is worth naming, because it is the one thing the registered encoder can do that the
repeated one cannot.

The family is a recognised one — weight-tied layers run to a fixed point, with the gradient taken through the last
rounds or through the fixed point itself. The canonical references could not be verified from this workspace: the
arXiv API is unreachable through the proxy and Crossref does not hold the conference records, so none is cited here
rather than cited from memory. The owner should add them.

### 19.3 The activates and inhibits assumption, withdrawn

The owner: "I don't think I agree that activates and inhibits should mainly differ by sign. Double check your
assumptions there."

**The objection is right and the sentence should not have been written.** Three reasons, in increasing order of how
badly it was wrong.

1. Mechanistically they are not negatives of one another. Activation and inhibition proceed by different mechanisms
   — phosphorylation at an activating against an inhibitory site, proteolysis, sequestration, competitive
   displacement, degradation — and their kinetics are asymmetric, since inhibition is bounded below by the absence of
   the target and saturates there while activation is not bounded in the same way.
2. The model cannot express the claim anyway. The architecture review records that the message-passing encoder
   "reads no edge sign and treats the perturbation's sign as a feature", so `activates` and `inhibits` are two
   relations with two independent matrices and nothing in that encoder relates them by negation.
3. The claim is not what basis decomposition assumes. Basis decomposition says the 20 relation matrices lie near a
   four-dimensional subspace of the 1,024-dimensional space of 32 by 32 matrices. Whether any two of them are
   negatives of each other has no bearing on whether four bases span the set. I offered a physiological-sounding
   reason for a purely numerical constraint, which is the kind of dressing this project exists to avoid.

**The argument that does hold is about data support, not physiology.** A full matrix gives every relation 1,024
parameters per layer whatever its size, and the relations are nothing like equal:

| | relations | edges |
| --- | --- | --- |
| relations with more edges than 1,024 | 9 | 278,386 |
| relations with fewer edges than 1,024 | **11** | **1,398** |

Eleven of the twenty relations have fewer edges in the entire graph than the parameters a full matrix would give
them: `changes_membrane_potential` 405, `changes_ion_pool` 395, `changes_redox_pool` 324, `regulates_release` 89,
`cofactor_of` 83, `oxidatively_modifies` 41, `produces_oxidant` 28, `redox_capacity` 19, `voltage_gates` 7,
`driving_force` 4 and `catalyzes_reverse_transport` 3. That last one owns 3,072 parameters under the registered
encoder and has three edges. Three relations carry 53.2 percent of all edges.

So the case for basis decomposition is that the registered parameterisation spends capacity in inverse proportion to
the evidence: the relations that could support a free matrix are a minority, and the rest are given one anyway. Under
four bases every relation owns four coefficients, which is proportionate to what the small ones can identify and is a
constraint on the large ones that the development runs will show the cost of. **B is a capacity knob and should be
chosen like one**: compare 2, 4, 8 and the full matrix on development folds and read the underfitting. If the
relations genuinely need distinct transforms, too small a B shows as underfitting, which is testable rather than
assumed in either direction.

## 20. Answers of 12 October, 17:09

### 20.1 Depth specialisation: the owner is right, and I over-weighted it

"Why does it need to specialize per depth? I thought there was separate matrices per relation type even if the bases
are shared."

The second sentence is correct and the two axes are independent. Basis decomposition acts on the **relation** axis:
with four shared bases there are still 20 distinct relation matrices, each a different mixture. Weight tying acts on
the **depth** axis: the same 20 matrices are applied at every round rather than three different sets of 20. Repeating
the layer therefore costs nothing per relation; it costs only the ability for round 1 to use a different `encodes`
matrix from round 2.

And on the first question: it does not need to. I listed depth specialisation as a loss without asking whether it is
a capability worth having, and it mostly is not. The law of a relation should not depend on how many edges a signal
has already crossed — `encodes` means the same thing at round 1 and round 5 — so a transform that varies with hop
count is a fitting convenience rather than a physiological one. The argument for keeping it is raw capacity, and at
2,975 labelled pairs capacity is the thing in short supply. The example I gave, that a knockout spends its first
layer on gene to protein, is not an argument either: with a repeated layer the `encodes` relation still has its own
matrix and still fires at round 1, it simply fires with the same matrix again later, which is what a physiological
law does. **Withdrawn as a cost.**

### 20.2 The message-passing encoder does not read edge signs, and this is a real defect

"Why would MP not read the edge sign?! That seems like a bug..."

**Verified in the code, not from the review.** `RelationalMessagePassingEncoder.build_relation_adjacencies` sets
every entry of every relation's adjacency to `1.0 / in_degree[destination]`, a positive number, and the constructor
takes no edge sign for ordinary relations: `entry_edge_sign` exists and applies to the drug entry edges alone. The
linear-response encoder is passed `data.edge_sign` and uses it. So the two registered encoders differ here, and the
message-passing one discards the sign of every ordinary graph edge.

**How much that costs, measured.** For a relation whose sign never varies, nothing is lost: the relation index
carries the sign and the relation's own matrix can learn to negate. That covers 12 relations and 236,881 edges,
`activates` at +1 and `inhibits` at -1 among them. The loss is confined to the 8 relations whose sign varies within
the relation:

| relation | edges | of which negative |
| --- | --- | --- |
| `regulates_transcription_of` | 41,638 | **7,087** |
| `changes_membrane_potential` | 405 | 52 |
| `changes_ion_pool` | 395 | 36 |
| `changes_redox_pool` | 324 | 49 |
| `regulates_release` | 89 | 65 |
| `oxidatively_modifies` | 41 | 28 |
| `voltage_gates` | 7 | 2 |
| `driving_force` | 4 | 2 |
| total | 42,903 (15.3 percent of all edges) | **7,321** |

So the encoder models 7,087 transcriptional repressions identically to activations of transcription. For a study
whose claim is mechanistic, a repressor read as an activator is not a small thing, and it is concentrated in the
second largest relation in the graph.

**It is not newly found, and it was not acted on.** The architecture review of 10 October listed it among the points
most likely to differ from the owner's expectation: "message passing reads no edge sign and treats the perturbation's
sign as a feature". It was reported and left, and I repeated it in 19.3 as though it were a settled design choice
rather than an open defect. It is the owner's call whether it was ever a choice.

**The fix is one line and costs no parameters**: make the adjacency value `sign / in_degree[destination]` instead of
`1 / in_degree[destination]`, which is what the linear-response encoder effectively does. One consequence to state in
the amendment rather than discover later: with signs in the values, the in-degree mean can cancel, so a node with one
activator and one repressor of equal weight receives zero. That is arguably correct and it is a change in behaviour,
so it belongs in the amendment with the rest. It is a registered change either way.

### 20.3 Shared bases assume too much, and per-relation rank does not

"Why wouldn't we just use 4 unique bases per matrix? Assuming that all relation types are close to the same 4d
manifold is insane..."

The objection holds. Basis decomposition asserts that all 20 relation matrices lie near one four-dimensional
subspace, which is an assumption across relations that nothing here tests, and I offered it as a default without
saying so plainly. What the owner describes — each relation with its own four bases — is a per-relation low-rank
factorisation, W_r = U_r V_r transposed with U_r and V_r of shape [32, 4]. It assumes only that each relation's
transform is low rank, which is an assumption inside a relation and not across them, and it costs almost the same:

| scheme | total | one relation owns | per labelled pair | what it assumes |
| --- | --- | --- | --- | --- |
| full matrices, shared across layers | 20,480 | 1,024 | 6.88 | nothing |
| **per-relation low rank 4** | **5,120** | **256** | **1.72** | each transform is low rank |
| per-relation low rank 2 | 2,560 | 128 | 0.86 | the same, more strongly |
| per-relation low rank 8 | 10,240 | 512 | 3.44 | the same, more weakly |
| shared basis, 4 bases | 4,176 | 4 | 1.40 | all 20 lie near one 4-d subspace |

Per-relation rank 4 costs 5,120 against the shared basis's 4,176. For a thousand parameters it drops an assumption
across 20 relations. That is the better default and the owner should have it.

**One thing the shared basis did better, and how to keep it.** Under a shared basis a relation with 3 edges owns 4
parameters, which is proportionate; under a flat per-relation rank 4 it owns 256, which is the over-parameterisation
of 19.3 in a milder form. The fix is to let the rank follow the evidence rather than be flat.
`experiments/count_encoder_parameters.py` prices one such scheme, rank 8 above 10,000 edges, 4 above 1,000, 2 above
100 and 1 below:

| tier | relations | rank | parameters each |
| --- | --- | --- | --- |
| above 10,000 edges | 8 (`activates`, `regulates_transcription_of`, `product_of`, `substrate_of`, `binds`, `catalyzed_by`, `inhibits`, `encodes`) | 8 | 512 |
| 1,000 to 10,000 | 1 (`member_of`) | 4 | 256 |
| 100 to 1,000 | 3 (the three pool relations) | 2 | 128 |
| below 100 | 8 | 1 | 64 |

**5,248 parameters in all, 1.76 per labelled pair** — the same budget as a flat rank 4, spent where the evidence is.
Even this is generous at the bottom: 64 parameters for `catalyzes_reverse_transport`'s 3 edges is still more than
those edges can identify, so the lowest tier should probably be diagonal (32 parameters) or a single gain per
channel, which is what the linear-response encoder gives every relation.

**Recommended, replacing the four shared bases in the amendment draft**: a per-relation low-rank factorisation with
the rank set by the relation's edge count, the tiers named in the registration, and the bottom tier diagonal. It
answers the parameter objection as well as the shared basis does, it assumes nothing across relations, and the ranks
are a stated rule rather than a fitted choice. The amendment still records the owner's "use the 4 bases" until he
says otherwise.

## 21. The rank each relation's evidence supports (12 October, 17:18)

The owner: "I was actually thinking that the rank of each relation type should be estimated based on the number of
training examples that pass through the relation type and how much that can support."

That is the better criterion and it replaces the edge-count tiers of 20.3. An edge is structure; what identifies a
parameter is gradient, and gradient reaches a relation only through labelled examples whose perturbation actually
sends signal along its edges. `experiments/measure_relation_flow.py` measures it: the 896 development perturbations
that hold a kept pair, 2,468 kept pairs behind them, propagated to the settled state at the contraction bound, with
each relation's flow taken as the encoder takes it, the sum over that relation's edges of the source's magnitude
divided by the target's in-degree **within** the relation.

### 21.1 Edge count is a bad proxy, which is the owner's point measured

| relation | edges | training pairs reaching it | flow share |
| --- | --- | --- | --- |
| `encodes` | 12,736 | 2,468 | **33.69 percent** |
| `catalyzed_by` | 24,266 | 2,357 | 13.73 |
| `substrate_of` | 36,562 | 2,333 | 10.63 |
| `activates` | **69,740** | 2,144 | **8.69** |
| `product_of` | 37,343 | 2,331 | 8.10 |
| `binds` | 31,177 | 2,113 | 8.02 |
| `regulates_transcription_of` | 41,638 | 2,100 | 7.26 |
| `inhibits` | 18,157 | 2,110 | 4.72 |
| `member_of` | 6,767 | 2,204 | 4.66 |
| the other 11 | 1,398 | 0 to 451 | 0.51 in total |

`activates` has five and a half times the edges of `encodes` and carries a quarter of its flow. The edge-count tiers
of 20.3 would have given both rank 8. On this measure `encodes` deserves four times what `activates` does.

Why `encodes` dominates is worth stating, because it is arithmetic rather than biology: every gene perturbation is
seeded on a gene node and `encodes` is the gene-to-protein edge, so all 896 development perturbations send their
whole signal through it at the first step, and a protein has exactly one gene, so the in-degree divisor is 1 and
nothing is averaged away. `activates` targets have many incoming activations, so each message is divided down. That
is faithful to what the encoder does, and it is a reason to read the measure as "how much signal this relation
carries" rather than as "how well this relation is identified". The two differ, and the second needs a trained model.

The bottom of the table is the firmer half of the finding: `changes_membrane_potential` is reached by **no** training
pair at all, `driving_force` by one and `changes_ion_pool` by eleven, while 9 relations carry 99.5 percent of the
flow. Whatever budget is chosen, almost all of it belongs to those nine.

### 21.2 What that evidence can support, and the honest limit of the calculation

Parameters have to come out of one budget rather than be allocated per relation independently, because the same 2,468
pairs flow through every relation at once. Taking the budget as kept pairs divided by a stated number of pairs per
parameter, split by flow share, and converting to a rank at 2 x 32 parameters per rank:

| pairs per parameter | budget | ranks it buys |
| --- | --- | --- |
| 5 | 494 parameters | `encodes` 2, `catalyzed_by` 1, everything else 0 |
| 10 | 247 | `encodes` 1, everything else 0 |
| 20 | 123 | nothing |

**This is austere and the austerity is the point, but the rule is cruder than the numbers look.** Pairs per parameter
is a rule of thumb from regression with independent predictors; here the parameters are shared across every example
rather than fitted per predictor, early stopping and weight decay regularise on top, and the heads and the feature
projection are consumers of the same evidence. So the calculation should be read as a bound of the right order rather
than a precise allocation: it says the registered 61,440 relation parameters are two to three orders of magnitude
beyond what 2,468 pairs identify, and it does not establish that rank 0 is correct.

The firm conclusions are the two that do not depend on the budget constant: the flow is concentrated in nine
relations, and the ordering within those nine is not the ordering by edge count.

### 21.3 What to build from it

A floor and a budget above it, which is the shape the measurement supports:

- **Floor: a diagonal gain for every relation**, 32 parameters each, 640 in all. That is one gain per relation per
  channel, which is what the linear-response encoder gives every relation, and it is the least that lets a relation
  scale and sign its messages. Below this a relation is not a relation.
- **Above the floor, rank bought by flow share** out of a budget the owner names. At 5 pairs per parameter that is
  rank 2 on `encodes` and rank 1 on `catalyzed_by`, for 832 parameters in all against the registered 61,440.
- **The ranks are fixed before any run and recorded in the registration**, with the flow measured on development
  perturbations only, which this script does. Measuring them on the held-out side would be a leak.

Two things the owner should weigh. The measure is the encoder's own aggregation, so a relation whose targets have a
low in-degree scores high for that reason alone; if the intent is identifiability rather than throughput, the
rigorous version is the gradient variance per relation, which needs a development run and so is blocked by "No new
pilots yet". And a relation that no training pair reaches, `changes_membrane_potential` among them, cannot be
identified at any rank, so its parameters are decoration whatever the scheme: the honest options for those eleven
relations are a diagonal gain held at its initialisation, or dropping them from the encoder and recording that the
study does not model them.

## 22. The electrical layer: why no training evidence reaches it (12 October, 17:29)

The owner, on the three relations of 21.3 that no training pair identifies: "How do the node descriptors look for the
node types involved in those three edge types?"

They look like the explanation. The three relations all belong to the electrical layer, and between the descriptors
and the aggregation there are two structural reasons the layer is silent that have nothing to do with how many
labels the study holds.

### 22.1 What the three relations connect

| relation | edges | sources | targets |
| --- | --- | --- | --- |
| `changes_membrane_potential` | 405 | 405 reaction nodes | **1**: `VM_c` |
| `driving_force` | 4 | 4 metabolite nodes | **1**: `VM_c` |
| `changes_ion_pool` | 395 | 356 reaction nodes | 5 ion metabolites |

`VM_c`, the graph's single `membrane_potential` node, is the target of every edge of the first two. It is not a dead
end: it feeds seven voltage-gated channel reactions through `voltage_gates`, with the signs the physiology asks for
(potassium and calcium channels at +1, the two HCN currents at -1). So the loop is built as intended, currents into a
potential and a potential into channels.

### 22.2 The hub carries no descriptors at all

The descriptor table has 138 columns in four groups: 10 metabolite, 35 reaction, 66 protein, 27 gene. A node's state
before anything propagates is the projection of its row.

| nodes | non-zero descriptor columns |
| --- | --- |
| every node, graph-wide | median 27 (quartiles 9 and 65) |
| the 405 reaction sources of `changes_membrane_potential` | median 27, and better filled than reactions generally: 23.3 of the 35 reaction columns active against 16.0 graph-wide |
| the 5 ion metabolites of `changes_ion_pool` | 9 of the 10 metabolite columns, all five |
| the 4 metabolite sources of `driving_force` | 9 of 10, all four |
| **`VM_c`** | **0 of 138** |

The sources are well described, several of them better than their node type at large, because they are annotated
transport reactions. The sink is blank. `VM_c` belongs to none of the four groups the descriptor table has — it is
not a metabolite, a reaction, a protein or a gene — so every column is zero and its initial state is the projection
bias alone, the same vector every featureless node gets. The one quantity the electrical layer exists to represent
enters the model with no identity of its own and can be distinguished only by what arrives along its edges.

### 22.3 And almost nothing arrives, because the mean divides by 405

The encoder averages a relation's messages over the target's in-degree **within** that relation. For these targets
that is not a mild normalisation:

| relation | edges | targets | in-degree of the target |
| --- | --- | --- | --- |
| `changes_membrane_potential` | 405 | 1 | **405** |
| `changes_ion_pool` | 395 | 5 | 338 for the largest, then 32, 9, 9, 7 |
| `driving_force` | 4 | 1 | 4 |
| `voltage_gates` | 7 | 7 | 1 each |
| `encodes`, for comparison | 12,736 | 12,709 | median 1 |

So one perturbed transport reaction moves the membrane potential by one four-hundred-and-fifth of its own magnitude,
while `VM_c`'s output to each channel passes undivided. That asymmetry is the opposite of the physiology: a membrane
potential is the integral of the currents crossing the membrane, which is a sum, and the channels then each read the
one potential. Averaging the currents and not the readout turns the layer into a near-silent node with a loud output.

This is why 21 found `changes_membrane_potential` reached by no training pair at the 0.1 percent threshold and
carrying 0.001 percent of the flow. It is not that perturbations fail to reach those reactions; it is that whatever
reaches them is divided by 405 before it is heard. The relation would be nearly silent at any label count.

### 22.4 What follows

Three things are now tangled together and they should be decided apart.

1. **The aggregation is wrong for an additive quantity.** A sum, or a sum scaled to keep the transition a
   contraction, is what a membrane potential wants. Summing 405 messages raw would take the row sum to 405 and break
   the contraction that 17.2 depends on, so the scaled forms are the real options: divide by the square root of the
   in-degree, or give the relation a learned scalar inside the contraction bound. This is a registered change and
   it is independent of every parameter question in 18 to 21.
2. **The hub has no features.** `VM_c` could carry a small group of its own columns, which costs almost nothing, or
   it could be left featureless and the amendment could say that the model distinguishes it only by position.
3. **Whether the electrical layer should be in the confirmatory graph at all.** It is in it today: these relations
   are part of `graph_full_neuronal_split_binders`, which is registered. If the aggregation and the descriptors are
   not fixed, the honest alternative is to disclose that the layer is present in the graph and carries no fitted
   behaviour, rather than to leave a reader to infer that a membrane potential was modelled.

None of this is a parameter-count question, and none of it is fixed by a rank. A relation whose messages are divided
by 405 into a featureless node will not be identified at rank 8 any more than at rank 1.

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
| `run_baselines.py --with-input-size`, `run_main_model.py --permute-drug-targets` | the first and third controls of 13.1; the second is reported without a flag by `run_baselines.py` and `aggregate_main_model_runs.py` |

Scripts of section 14, under `experiments/label_review/`:

| script | what it does |
| --- | --- |
| `chembl_binding_sites.py` | every ChEMBL mechanism row of a study drug with its `site_id` and `binding_site_comment`, and the comments marked where they hedge (14.1) |
| `cooperativity_sources.py` | every cooperativity-style activity in ChEMBL and every site and cooperativity field of Guide to PHARMACOLOGY, against the study's drugs (14.2); `--fetch` writes the ChEMBL cache |
| `measured_occupancy_leakage.py` | the coverage of each concentration scheme, the two scales of a mixed input and the leakage groups each scheme produces (14.3) |
| `tissue_representation_options.py` | the tissue sources on disk, what each covers of the graph's nodes, how many tissues a gene reaches and the tissue spread of the mechanism targets (14.4) |
| `enzymes_in_the_graph.py` | enzymes as graph structure, as drug mechanism targets, as labelled gene perturbations and in the `clears` relation (14.5) |

Scripts of section 15:

| script | what it does |
| --- | --- |
| `label_review/binding_site_names.py` | binding sites matched by name over ChEMBL and Guide to PHARMACOLOGY, and the drug pairs a site settles (`configs/binding_sites.csv`; 15.1, 15.2) |
| `label_review/measured_occupancy_leakage.py` | also the Lombardo dependency of 15.3, through `data/processed/exposure/drug_exposure_v4.csv` |
| `label_review/placebo_arm_sizes.py` | the arm sizes of the 61 placebo-controlled pairs from openFDA and the weight each gives (`configs/placebo_negative_weights.csv`; 15.6) |
| `label_review/target_rule_options.py` | the alternatives to the target limit, the leakage group each leaves and the drugs that bridge the blocks (15.7) |
| `label_review/enzymes_in_the_graph.py` | also the enzyme descriptor columns of 15.8 |
| `measure_tissue_saturation.py` (under `experiments/`) | the per-tissue saturation curves and the step each tissue fills fastest (`docs/tissue_saturation.md`; 15.9) |

Scripts of sections 16 to 21, under `experiments/`:

| script | what it does |
| --- | --- |
| `label_review/cooperativity_sources.py` | also the SIGNOR and OmniPath survey of 16.3 |
| `measure_propagation_convergence.py` | the steps the linear response needs to settle, at several gains (17.2) |
| `count_encoder_parameters.py` | where each encoder's parameters go, and what each reduction scheme would cost (18, 20.3) |
| `measure_relation_flow.py` | the training evidence that passes through each relation and the rank it supports (21) |

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
| `chembl_cooperativity/cooperativity_activities.json` | the ChEMBL activity API, every activity of standard_type Cooperativity, Alpha, Fold shift, Shift and Dose ratio, fetched 12 October 2026 (535 activities) | `17dda554a91e3985` |
| `placebo_arm_sizes/*.json` | openFDA label queries of 12 October 2026, one file per drug, for the arm sizes of 15.6 | one file per drug |
| CPIC API, `api.cpicpgx.org/v1/{pair,recommendation,drug}` | read 12 October 2026 for 15.5; CC0, quoted in docs/data_sources.md | not cached in the repository |
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

Two research passes ran on 12 October 2026, for sections 13.4 and 13.9. Both followed the source rules of `CLAUDE.md`:
no permission request was made or left waiting, every page was read with `experiments/read_public_page.py` or through a
programmatic service, and a refusal was recorded rather than worked around. Refused and still unread, for the user to
open: `go.drugbank.com/legal/terms_of_use` (HTTP 403), `pdsp.unc.edu/databases/kidb.php` and
`drugtargetcommons.fimm.fi` (both unreachable through the proxy's certificate chain), every address on
`dmd.aspetjournals.org` (the host disallows an AI assistant's user agents), and
`marketplace.copyright.com/rs-ui-web/mp/search/journal/1521009X` (served, but a script application with no text).
Text addressed to an AI assistant, reported as the rules require and not acted on: BindingDB's download pages advise AI
data users to cite an archived release when training a model, which is advice about reproducible versioning and not an
instruction to an assistant. Nothing else of the kind was seen.

Sources named above and not linked there: CRESCENDDI, Kontsioti et al. 2022
([doi:10.1038/s41597-022-01159-y](https://doi.org/10.1038/s41597-022-01159-y)); TWOSIDES, Tatonetti et al. 2012
([doi:10.1126/scitranslmed.3003377](https://doi.org/10.1126/scitranslmed.3003377)); ClinGen dosage curation,
[doi:10.1002/humu.24291](https://doi.org/10.1002/humu.24291); IMPC,
[doi:10.1093/nar/gkaf1148](https://doi.org/10.1093/nar/gkaf1148).
