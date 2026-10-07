# Literature appraisal, staging file

**Status: UNREVIEWED.** Nothing here has been checked by the author. No row may be
copied into `docs/experiment_design.md`, `docs/references.bib` or any other document
until it has been read against the cited source. Rows are written by the paced
appraisal task described at the end of this file.

## What each row must carry

A row is admissible only with a PMID or DOI and a quoted sentence from the source.
A paraphrase without a quote is not a row.

The appraisal records three fields rather than one, because "supports" alone does not
say how much weight the paper can bear.

### Direction

| value | meaning |
|---|---|
| `supports` | the paper's result points the same way as the claim under test |
| `contradicts` | it points the other way |
| `mixed` | it points both ways across conditions the paper reports |
| `no_bearing` | it does not test the claim, whatever its topic |

### Bearing: how close the paper's test is to our claim

| value | meaning |
|---|---|
| `direct` | the same comparison, on the same kind of object and outcome |
| `same_construct_other_domain` | the same comparison, a different network or a different outcome |
| `analogous_mechanism` | an argument about mechanism, without running the comparison |
| `background` | supplies a definition, a dataset or a method, with no bearing on the comparison |

### Weight: how much the paper's design licenses its own conclusion

| value | meaning |
|---|---|
| `strong` | held-out evaluation, the baseline at issue is run, the effect exceeds the reported fold or seed variation |
| `moderate` | held-out evaluation, but the baseline is weaker than ours, or variation is not reported |
| `weak` | in-sample fit, one split, no variability, or the comparison appears only in passing |
| `assertion` | stated without a measurement (a review, an introduction, a discussion aside) |

Every row also states, in one line each:

- **what it would take to settle the question**: the measurement the paper would have had
  to report for the claim to be closed either way
- **what it leaves open**: the part of our claim the paper does not reach

## Questions under appraisal

Q1 is live: it bears on the result recorded in `docs/graph_content_null_results.md`, where
the graph's wiring carries signal that an unsigned random walk already extracts and the
signed mechanistic encoder adds nothing on top.

- **Q1.** Has signed or relation-typed propagation on a biological network been shown to
  beat unsigned diffusion for predicting a phenotype or a clinical label?
- **Q2.** Where edge signs are reported to help, what fraction of nodes receive messages of
  both signs? (Our graph: 7,811 of 21,337 destinations, 36.6%, and all of them annihilate under the
  encoder's cross-relation mean.)
- **Q3.** Is noisy-OR pooling over pathway modules reported against a linear or sum-pooled
  alternative anywhere, with a held-out comparison?

## Rows

| PMID / DOI | question | direction | bearing | weight | quoted sentence | to settle | leaves open |
|---|---|---|---|---|---|---|---|
| [10.1038/s41591-024-03233-x](https://doi.org/10.1038/s41591-024-03233-x) (PMID 39322717) | Q1 | `supports` | `same_construct_other_domain` | `strong` | "However, TxGNN outperformed all methods across all nine disease areas for both tasks, demonstrating its broad generalizability and accuracy in zero-shot drug repurposing." | the per-baseline AUPRC under the disease-area split with its spread (diffusion state distance, network proximity and the plain relational network beside TxGNN), a run with all relation types merged into one, and a diffusion baseline given a trained readout on the same labels, so that typing is separated from supervision; for Q1 proper, a graph that does not contain the label relation | the graph has 29 undirected edge types and no signs, so only the relation-typed half of Q1 is touched; indication and contraindication are edge types of the same graph and the training diseases' label edges stay in it; the untyped baselines are unsupervised statistics, so typing is confounded with supervision, pretraining and the disease-similarity module; per-baseline numbers are in figures the served text omits, and the text names a language model as the best baseline for held-out indications in seven of nine areas and the plain relational network for contraindications in eight of nine, so the diffusion baselines were not the bar to clear; gains over the next best method run from 0.5% to 59.3% across areas; the abstract's 49.2% and 35.1% were not found in the body as served; the negative sampling that sets the AUPRC base rate is in a supplementary note that was not read |
| [10.1186/1471-2105-15-238](https://doi.org/10.1186/1471-2105-15-238) (PMID 25015298) | Q1, Q2 | `mixed` | `analogous_mechanism` | `moderate` | "As in the previous example, performance results showed a more robust behavior of network-based features as compared to gene selection, for the majority of the learning algorithms used." | the same classifiers on backbone values computed three ways (the signed Laplacian as published, all edges made positive, signs permuted with the topology kept), scored across cohorts with intervals sized for test cohorts of about 24 patients; and a plain diffusion of the gene fold-changes over the same backbone as the unsigned baseline | no run removes or permutes the signs for prediction, so Q1 is not reached: the baselines are gene-level signatures and gene sets, not an unsigned propagation; the K statistic permutes backbone edges together with their signs and tests the perturbation score, not the classifier; for infliximab response the backbone signature is not the best performer (mean G-performance across cohorts 0.83, against 0.85 for a classifier on the genes downstream of one node and 0.80 for nearest shrunken centroids on all genes), with 47 pre-treatment patients over both cohorts; for smoking status all five backbone classifiers score 0.83 to 0.91 across data sets while gene-level ones range from 0.00 to 0.88; edge direction is discarded by the signed Laplacian; the inference runs backward from expression to upstream activity in a patient, not forward from a perturbation to a phenotype; standard errors are in a supplementary table that was not read; the fraction of nodes receiving both signs is not reported |
| [10.1038/s41467-019-10887-6](https://doi.org/10.1038/s41467-019-10887-6) (PMID 31289271) | Q1 | `supports` | `same_construct_other_domain` | `moderate` | "We found that the orientation-based scores outperform the ones obtained when using an unoriented network (Fig.), with a mean rank for a true drug target of 6045 (out of 15,501 network genes) compared to a mean rank of 8837 obtained by a completely unoriented network" | the same diffusion on three networks (unoriented, curated directions only, curated plus inferred directions), scored by a top-k or area measure with the chance level beside it and the spread over drugs, and a fourth network whose edges also carry signs; for Q1 proper, a phenotype outcome | the attribute is edge direction, with no signs anywhere, and the outcomes are drug-target and driver-gene ranking, not a phenotype; a random ranking of 15,501 genes has an expected mean rank of 7,751, so the unoriented baseline (8,837) is below chance and the oriented network (6,045) puts the true target at the 39th percentile on average; "completely unoriented" appears to strip the 33,756 curated directed interactions as well, so the gain may come from known directions and not the inferred ones; the driver-gene result prints numbers for the oriented network only, with the comparison in a figure; the p-value exponents are missing from the PMC text as served; checks 3 and 6 hold as described, since the drug-target test orients from cancer data alone and the driver test leaves out the scored disease |
| [10.64898/2026.04.29.721775](https://doi.org/10.64898/2026.04.29.721775) (bioRxiv preprint, version 1, no PMID) | Q1, Q2 | `mixed` | `same_construct_other_domain` | `moderate` | "Without task-specific fine-tuning, FLASH consistently outperforms or matches nine state-of-the-art unsigned, relational, and signed graph baselines across drug mode-of-action prediction, clinical response modeling, and drug-drug interaction prediction, while substantially improving computational efficiency." | the same model with signs removed or permuted and nothing else changed, on a split that holds out whole drugs or diseases, with the majority-class score beside each accuracy; for Q1 proper, an outcome that is not itself an edge of the graph | the "clinical response" label is the sign of a compound-disease edge of the same graph (indication against contraindication), split by random edges, so this is the association-label design of 39976387; the unsigned baselines lose edge types as well as signs and differ in architecture, and the signed models receive the training labels as input edge signs, so the gain cannot be credited to sign propagation; contraindications are 92.6% of the labels and no majority-class score is given beside the binary accuracy of 0.9893; in the one drug-disjoint setting, on an unsigned outcome, the lead over the best baseline is not significant (p = 0.20 as relayed, three runs); balance is imposed by pruning edges in unbalanced cycles, which our graph could not survive; a preprint with unresolved table references; the body was read through a summarising fetch tool, so every figure beyond the abstract is relayed |
| [10.1093/bioinformatics/btm170](https://doi.org/10.1093/bioinformatics/btm170) (PMID 17646318) | Q1 | `supports` | `analogous_mechanism` | `weak` | "In cross-validation tests, SPINE obtains very high accuracy in predicting knockout effects (99%)." | for Q1, a readout that an unsigned propagation over the same physical network can also be scored on (whether a knockout changes a gene at all, or a phenotype), since direction gives it nothing to predict; for the 99%, the counts of up and down effects among the held-out pairs and a split that holds out whole knockouts instead of single pairs | a direction task on expression, so Q1 is not reached and the readout is not a phenotype; the signs are fitted to the same knockout compendium they are scored on, with one or five pairs hidden per iteration, so the other pairs of the same knockout stay in the fit; the only baseline is another signed method (Yeang et al. 2004), and by the table as relayed the leave-one-out gap is 101 against 100 correct of 103 held-out pairs; no class balance or majority-direction rate in the part read; the figure is for the edge variant on a mating subnetwork of 58 interactions and 149 knockout pairs; the full text is not in PMC and was read through a publisher-page fetch that stopped in section 3.2 and relays only short quotes, so every figure beyond the abstract needs checking against the PDF |
| [10.1093/bioinformatics/btaf674](https://doi.org/10.1093/bioinformatics/btaf674) (PMID 41429577) | Q1 | `supports` | `analogous_mechanism` | `weak` | "Evaluation by precision-recall analysis gave an area under the curve (AUPRC) of 0.8, indicating that SIGNAL accurately recapitulates TLM knockout phenotypes." | the telomere evaluation with its positive class and base rate stated, against an all-positive network (which predicts every path positive, the majority-direction baseline) and with the telomere genes' own knockouts excluded from the sign features | a direction task, so an unsigned propagation has no prediction to compare and Q1 is not reached; the AUPRC of 0.8 comes with no positive class or base rate, and with 456 'long' against 37 'short' pairs a constant predictor would score about 0.925 if 'long' were positive; those pair counts also contradict the gene counts (63 short-phenotype genes, 38 long), and the discussion calls the same figure an AUC; the knockout-signature reconstruction (precision-recall area 0.87) excludes the knockout's own features in training but also reports no base rate; the method rests on negatives being rare, the fourth row developed in that regime |
| [10.1038/s41467-022-30684-y](https://doi.org/10.1038/s41467-022-30684-y) (PMID 35654811) | Q1 | `no_bearing` | `background` | `moderate` | "Here we develop a recurrent neural network framework constrained by prior knowledge of the signaling network with ligand-concentrations as input and transcription factor-activity as output." | the same cross-validation with the network constraint removed (a dense or randomly sparse recurrent network with as many parameters) and with the prior edge signs removed or permuted, r reported per fold | the abstract reports no model without the prior network and no unsigned variant, so whether the topology or its signs contribute is not shown; the readout is transcription factor activity, not a phenotype; whether the cross-validation (r = 0.8) holds out whole ligands is not stated; the knockout result (r = 0.8) is on synthetic data whose generating network is not stated. The full text is in PMC (PMC9163072) and would decide whether a network ablation exists; it was not fetched, since even a positive ablation on transcription factor activity would not reach Q1's phenotype comparison |
| [10.1073/pnas.1720589115](https://doi.org/10.1073/pnas.1720589115) (PMID 29925605) | Q1 | `supports` | `analogous_mechanism` | `moderate` | "We find that simple dynamically agnostic models are sufficient to recover the strength and sign of the biochemical perturbation patterns observed in 87 biological models for which the underlying kinetics are known." | on the same 87 models, a sign-blind propagation beside the signed one for strength, and for direction the accuracy set against each model's majority-direction rate; for Q1 proper, a phenotype readout with unsigned diffusion run on the same network | the comparison it runs is topology against full kinetics, not signed against unsigned, so Q1 is not reached; the 65-80 percent and the ~80 percent chemotaxis figures come with no majority-direction baseline in the abstract or significance statement, and our own sign check shows why that matters (signed response 0.615 against 0.806 for always guessing the majority direction); the full text is not served through PMC, so the methods, the network sizes and whether the models contain frustrated loops are unread |
| [10.3390/e26020161](https://doi.org/10.3390/e26020161) (PMID 38392416) | Q2 | `no_bearing` | `background` | `moderate` | "negative cross-correlations have become pyramidally rare in the past three decades" | nothing: the network is US foreign exchange rates, inferred from detrended cross-correlation, with no biology and no phenotype | everything about Q1; it is kept only as a third instance of signed-network method development on a negative-sparse, balance-satisfying network |
| [10.1093/bib/bbae573](https://doi.org/10.1093/bib/bbae573) (PMID 39523622) | Q1, Q2 | `supports` | `same_construct_other_domain` | `weak` | "combined with balance theory and state theory, molecular features are extracted from the point of social relations through the propagation and aggregation of signed graph attention network" | an ablation isolating the signed graph attention from the topological embedding, the denoising autoencoder and the GBDT classifier, and an unsigned baseline on the same network; the abstract reports neither | everything about attribution, and whether any of it holds on a network that violates balance, which ours does at 26.9% of edges |
| [10.1093/bib/bbaf062](https://doi.org/10.1093/bib/bbaf062) (PMID 39976387) | Q1 | `supports` | `same_construct_other_domain` | `moderate` | "CSGDN uses a signed graph diffusion method to uncover the underlying regulatory associations between genes and phenotypes" | an unsigned diffusion baseline on the same graph and split, and an ablation separating the signed diffusion from the contrastive objective; they report neither | whether the gain comes from the signs or from the contrastive views, and whether anything survives when the signed graph is a mechanistic substrate rather than the label graph itself |
| [10.1016/j.neunet.2021.04.007](https://doi.org/10.1016/j.neunet.2021.04.007) (PMID 33906083) | Q1, Q2 | `no_bearing` | `background` | `strong` | "we obtain some sufficient conditions which can guarantee that networks with signed graph topologies realize bipartite synchronization under any initial conditions and arbitrary switching signals" | nothing: it is a control-theory result about coupled reaction-diffusion systems, with no biological network, no phenotype and no unsigned comparison | everything about Q1, but it supplies the condition under which signed dynamics reduce to unsigned, which is testable on our graph |
| [10.1186/s12915-024-01968-0](https://doi.org/10.1186/s12915-024-01968-0) (PMID 39148051) | Q1 | `supports` | `same_construct_other_domain` | `weak` | "A novel strategy for propagating signed message in signed networks addresses heterogeneity and consistency among nodes connected by signed edges" | an ablation replacing their signed propagation with an unsigned one, same features and same split; the abstract reports neither | whether the reported AUROC of 0.9742 owes anything to the signed strategy, and whether any of it transfers off a similarity-feature bipartite graph |
| [10.1186/1756-0500-7-516](https://doi.org/10.1186/1756-0500-7-516) (PMID 25113603) | Q1, Q2 | `supports` | `analogous_mechanism` | `weak` | "the aggregation method we used is limited to a specific class of causal network models called \"causally consistent\", which is equivalent to the notion of balance of a signed graph used in graph theory" | a held-out comparison of their signed aggregation against an unsigned random walk on the same networks and the same readout; they report neither | whether an unsigned walk would have scored the same, which is our whole question |

## How the task is paced

One paper per firing, each firing scheduled about five minutes after the last, through a
one-shot self-directed reminder rather than a background shell. A firing does nothing but
read one source, write one row, push the file and schedule the next. If any registered job
under `runs/jobs/` needs attention, the firing handles that and skips the paper.

The limit is enforced by `scripts/literature_pacing.sh`, not by intention. The firing must
call it with `claim` before reading anything, which checks and consumes the slot in one
step, so a firing cannot read a second source by calling the gate twice. The gate refuses
inside 300 seconds of the last appraisal and refuses outright after 19 appraisals, which is
the number of PubMed records matching the Q1 query. A refusal means no source is read that
firing. The chain also stops on a message from the author or when a question is closed.

### Handoff to a separate session (7 October 2026, 17:55 UTC)

The paced chain was stopped in the experiment session at the user's request after slot 10 of 19 (nine rows and one
replacement search), so that it can continue in a session started after the user connected a Consensus Pro account.
From here on that literature session owns this file; the experiment session reads it but does not edit it.

The gate's state lives in runs/literature_appraisal/pacing_state, which git ignores, so a fresh clone starts the count
at zero. Before the first claim in the new session, seed it with the count spent here:

    mkdir -p runs/literature_appraisal
    printf 'last_appraisal_epoch=0\nappraisals_consumed=10\n' > runs/literature_appraisal/pacing_state

Then each firing claims one slot with `bash scripts/literature_pacing.sh claim` and spends it on one source or one
search, under the rules in this file. The next source is the one entry left in the replacement queue (Ourfali 2007);
after it, the next slot goes to a new mechanism-named search. Commit only this file, and pull before every push,
since the experiment session pushes to main as well. The literature session runs no experiment jobs: its container
has no runs/ directory, so scripts/resume_jobs.sh has nothing to relaunch there, and it should not register any.

### Slot log, literature session

The gate's state is not in git, so the count is repeated here. A fresh clone seeds `appraisals_consumed` from the
last line of this list.

- Slot 11, 7 October 2026, 18:00 UTC: Ourfali 2007 appraised as PMID 17646318. Slots spent: 11 of 19.
- Slot 12, 7 October 2026, 18:29 UTC: one Consensus search, 20 records, five queued. Slots spent: 12 of 19.
- Slot 13, 7 October 2026, 18:38 UTC: Mottaqi 2026 appraised as DOI 10.64898/2026.04.29.721775. Slots spent: 13 of 19.
- Slot 14, 7 October 2026, 19:00 UTC: Silverbush 2019 appraised as PMID 31289271. Slots spent: 14 of 19.
- Slot 15, 7 October 2026, 19:08 UTC: Martin 2014 appraised as PMID 25015298. Slots spent: 15 of 19.
- Slot 16, 7 October 2026, 19:19 UTC: one Consensus search on the relation-typed half of Q1, 20 records. Slots spent: 16 of 19.
- Slot 17, 7 October 2026, 19:27 UTC: Huang 2024 (TxGNN) appraised as PMID 39322717. Slots spent: 17 of 19.

### Triage checks: the original four and three additions

The replacement-queue section refers to "the four checks recorded at the top of the rows". They were never in this
file. They come from the firing prompt that set the experiment session's check-ins (16:43 UTC, 7 October 2026), and
the author supplied the wording at 19:06 UTC:

"Four checks have decided every row so far: whether the sign is a mechanistic edge property or the association label
being predicted; whether the propagation runs over a substrate carrying no label information or over the label graph
itself; whether the split holds out whole entities or is transductive; and whether the method assumes structural
balance, which a stoichiometric metabolic graph violates by construction at 26.9% of edges."

The 26.9% is quoted as the prompt gave it and was not re-measured for this entry. Numbered as the notes use them:

1. Is the sign a mechanistic edge property, or the association label being predicted?
2. Does the propagation run over a substrate carrying no label information, or over the label graph itself?
3. Does the split hold out whole entities, or is it transductive?
4. Does the method assume structural balance?

Checks 5 to 7 are additions by the literature session, in use from slot 11. None of them was in the original.

5. Is the signed component isolated: the same pipeline run without signs or with signs permuted, on the same
   network and split?
6. Where signs are fitted, are they kept apart from the outcome they are scored on?
7. Is each direction or classification figure reported beside its base rate, the score of always predicting the
   majority class?

Correction, 7 October 2026, 19:06 UTC. Slots 11 to 14 were triaged against a reconstruction that this list
replaces. The reconstruction had check 1 right, merged original checks 2 and 3 into one, left out original check 4
as a check and numbered the three additions 3, 4 and 5. Check numbers in the slot 12 search section and in the notes
on 10.64898/2026.04.29.721775 and 31289271 have been renumbered to the list above. The three rows written in those
slots, against the original four:

- 17646318: check 1 holds, since the signs sit on physical interactions, though they are fitted (check 6 fails).
  Check 2 fails in part, because the fitted signs carry the knockout effects that are scored. Check 3 fails: one or
  five pairs are hidden and the rest of each knockout stays in the fit. Check 4 is not established: a pair is scored
  only when its paths agree (relayed) and no balance measure is reported.
- 10.64898/2026.04.29.721775: checks 1 and 2 fail. Check 3 fails for the sign tasks, which split random edges, and
  holds for the drug-disjoint setting. Check 4 fails, since balance is imposed by pruning.
- 31289271: checks 1, 2 and 3 hold. Check 4 does not arise, because the network carries no signs.

No direction, bearing or weight changes as a result. The queue order from slot 12 also stands, with one caution
added: if Martin 2014 uses the framework of the first row, it is restricted to the balanced case and fails check 4.

### Note on the first row

The paper is a methods note, not a benchmark, so its weight is `weak` and its bearing is
`analogous_mechanism`: it compares two aggregation algorithms on causal biological networks
and never runs an unsigned walk as a baseline. What it supplies is the standard name for the
property this project measured by hand. Signed-graph balance, which the authors call causal
consistency, is the condition under which the signs carry no structure of their own: a
balanced signed graph is switching-equivalent to the all-positive graph, so propagation over
it differs from unsigned propagation by a node-wise sign that a linear readout absorbs. If
our graph is close to balanced, the null result in `docs/graph_content_null_results.md` has
a stated cause rather than an open puzzle.

**Where this does not apply.** Balance is a condition on the sign product around each cycle,
so a sign-homogeneity count at destinations is not a balance measurement either way. The test
is the frustration index, or its cheap exact check for balance: two-colour the nodes so that
every positive edge joins same-coloured nodes and every negative edge joins differently
coloured ones, by breadth-first search per component, and count the edges that fail.
`experiments/measure_signed_graph_balance.py` does this. Run on the stored edge table of
`data/processed/graph` it returns zero frustrated edges, but that result is vacuous: the
stored table has no negative edges at all, since the negative signs are derived inside the
encoder. The script therefore has to be run on the derived relation stack before it says
anything, which is recorded as open work rather than a result.

## Queue for Q1

Nineteen PubMed records match `("signed network" OR "signed graph" OR "signed networks") AND
(propagation OR diffusion OR "random walk")` in title or abstract. Appraised: 25113603, 39148051, 33906083, 39976387, 39523622, 38392416, 29925605, 35654811, 41429577, 17646318.
Set aside as out of scope on the abstract, with the reason:

- 28297881, 31661504 — epidemic and information spread on social signed networks; no
  biological network and no phenotype readout
- 28129187 — antisynchronization control of coupled reaction-diffusion neural networks; a
  multilayer signed graph used as a coupling topology, not a predictor
- 35131567 — signed network representation scored on sign prediction and link prediction,
  not on a phenotype label; keep for Q2 if the node-proximity metric reports sign mixing

Superseded unread on 7 October 2026 by the replacement query below, not discarded: 35364845, 40103114,
40680519, 36049951, 42481621 and the four records beyond the first fifteen. They matched only the generic phrase,
and the slots are better spent on the candidates the mechanism-named query found; read them only if the cap allows
after the replacement queue is done.

### Replacement queue (slot 7, one Consensus search, 10 records)

Query, run through Consensus because it matches meaning rather than phrase: "propagating gene knockout or
inhibition effects through a signed signaling or metabolic network outperforms unsigned network diffusion for
predicting phenotype or disease genes". The free tier returns no PMID or DOI, so candidates are listed by first
author, year and journal, and the firing that reads one resolves its PMID in the same single metadata call. Triage
used only the search records, against the four checks recorded at the top of the rows.

Remaining to read: none. The four queued records are appraised: 29925605 (slot 8), 35654811 (slot 9), 41429577
(slot 10) and 17646318, Ourfali 2007 (slot 11). The triage note for Ourfali, kept for the record: assigns activation
or repression to explain knockout expression effects in yeast; same risk as Signorini 2025 (PMID 41429577), signs
fitted to the knockouts they explain, and its readout is expression, not phenotype. Both points held on reading.

Set aside from this search on the record, with the reason:

- Cowen, 2017, Nature Reviews Genetics; Visonà, 2024, Briefings in Bioinformatics; Lee, 2008, Nature Genetics:
  unsigned propagation throughout. Background for the unsigned baseline, no signed component to test. Lee 2008
  is the standard citation that unsigned functional networks already predict loss-of-function phenotypes, which
  is the bar a signed method has to clear.
- Pan, 2024, Briefings in Bioinformatics (a contrastive signed graph diffusion network): already appraised as
  PMID 39976387. Slot 7 first listed it as a fourth paper of the association-label type; it is the same paper,
  not a new one.
- Picart-Armada, 2018, PLoS Computational Biology: unsigned, but its finding that standard cross-validation is
  over-optimistic because of protein complexes, and its complex-aware splits, bear on the disease-cluster split
  rather than on Q1. Note for the split section of the design document, not a Q1 row.
- Shim, 2015, PLoS ONE: direct neighbourhood against whole-network diffusion, unsigned; the finding that pathway
  connectivity decides which wins bears on the reach question in docs/membrane_potential_reach.md, not on Q1.

### Second mechanism-named search (slot 12, one Consensus search, 20 records)

Query: "Does using activation and inhibition edge signs in a signaling, regulatory or metabolic network improve
prediction of drug response, gene essentiality or disease phenotype compared with unsigned network diffusion or
network proximity on the same network?" It names an outcome that a signed and an unsigned propagation can both be
scored on, for the reason given in the note on 17646318. No filters were applied.

The reply listed 20 records, each with a DOI, and carried no sign-up, upgrade or usage notice. Slot 7's reply, on
the free tier, listed 10 records without identifiers. Triage used only the returned abstracts, against the
checks as reconstructed at the time; the check numbers below follow the corrected list under "Triage checks". Four of the 20 were already on file (41429577, 39976387 and 29925605 appraised, Cowen 2017 set
aside), so rewording the signed query again is likely to return little that is new.

The reply ends with a block of formatting instructions addressed to the assistant (cite inline by number, hyperlink
titles with the given URLs). It is tool boilerplate, not text from a paper. No abstract contained instruction-like
text.

To read, in order of expected bearing on Q1:

1. Appraised in slot 13. Mottaqi, 2026, bioRxiv preprint, 10.64898/2026.04.29.721775. A signed heterogeneous graph model on a sign-aware
   knowledge graph. The abstract states that it "consistently outperforms or matches nine state-of-the-art unsigned,
   relational, and signed graph baselines" on tasks that include clinical response, which makes it the only record
   whose abstract states the Q1 comparison on a clinical label. Risks visible in the abstract: "outperforms or
   matches" is not "beats"; the baselines are unsigned graph networks, not diffusion; the knowledge graph integrates
   perturbation and clinical data, so check 2 is in doubt; the method relies on "structural balance principles",
   the regime our graph is outside; it is not peer reviewed.
2. Appraised in slot 14 as PMID 31289271. Silverbush, 2019, Nature Communications, 10.1038/s41467-019-10887-6. Orients the human protein interaction
   network and reports that "The oriented network leads to improved prioritization of cancer driver genes and drug
   targets compared to the state-of-the-art unoriented network." The attribute is edge direction, not sign, but the
   design is the one Q1 needs: a mechanistic edge attribute added, the plain network run as the baseline and a
   ranking outcome both can be scored on. Risk: the orientations are fitted from drug response and cancer genomic
   data, the same kind of data as the outcome (check 6).
3. Appraised in slot 15 as PMID 25015298. Martin, 2014, BMC Bioinformatics, 10.1186/1471-2105-15-238. Scores signed cause-and-effect networks from
   transcriptomics and builds a signature "to predict the response to the treatment" in ulcerative colitis, a
   clinical label. The abstract says the value of the method's components is "substantiated", which may include a
   test with the signed structure removed or permuted (check 5). Probably the same framework as the first row
   (25113603); to confirm on reading.
4. Trinh, 2016, Bioinformatics, 10.1093/bioinformatics/btw464. Boolean dynamics on signalling networks: "the
   edgetic sensitivity predicted drug-targets better than the node-based sensitivity". Drug-target status is an
   outcome an unsigned measure can be scored on, and the abstract mentions connectivity and betweenness, so a
   structural baseline may be reported.
5. Gates, 2021, PNAS, 10.1073/pnas.2022598118. Compares a logic-weighted effective graph with the plain
   interaction graph for how perturbations propagate in 78 Boolean models. Mechanistic detail against topology
   alone, scored on the models' own dynamics, so the limit met in 29925605 is likely to apply.

Set aside from this search on the abstract, with the reason:

- Chen, 2025, Methods (10.1016/j.ymeth.2025.05.005); Nakis, 2025, Bioinformatics (10.1093/bioinformatics/btaf204);
  Xi, 2025, arXiv (10.48550/arxiv.2512.11927); Roy, 2020, Bioinformatics (10.1093/bioinformatics/btaa651): the sign
  is the label being predicted or the network being inferred (check 1), with no phenotype outcome.
- Vinayagam, 2013, Nature Methods (10.1038/nmeth.2733): infers interaction signs from genetic screen phenotypes, an
  annotation task with signs fitted to phenotype data (check 6) and no unsigned comparison. Background for where
  signed protein interaction networks come from.
- Thiele, 2015, BMC Bioinformatics (10.1186/s12859-015-0733-7); Le Bars, 2023, BMC Bioinformatics
  (10.1186/s12859-023-05429-3): sign-consistency methods that predict the direction of unobserved nodes, so the
  limit of the other direction tasks applies. Both bear on the design question in the note on 17646318 instead:
  Thiele separates strong from weak predictions, and Le Bars outputs "the weight of the predicted sign". Notes for
  that question, not Q1 rows.
- Zickenrott, 2016, Cell Death & Disease (10.1038/cddis.2015.393): ranks genes and compounds from differential
  regulatory networks; the abstract reports no comparison with an unsigned method.
- Brilliantova, 2022, AAAI proceedings by its DOI, journal not named in the record (10.1609/aaai.v37i10.26457):
  generative models of edge signs given a regulatory topology, with no phenotype prediction. Possible background
  for Q2 if it reports sign statistics of real regulatory networks.
- Galindez, 2022, Computational and Structural Biotechnology Journal (10.1016/j.csbj.2022.12.022); Shojaie, 2020,
  WIREs Computational Statistics (10.1002/wics.1508): reviews, so weight `assertion` at most.
- Already on file: Signorini 2025 (41429577), Pan 2024 (39976387), Santolini 2018 (29925605) and Cowen 2017.

Plan for the remaining seven slots: read 1 to 3 in order, spend one slot on a search for the relation-typed half of
Q1, then take the last three from the two queues by expected bearing.

### Relation-typed search (slot 16, one Consensus search, 20 records)

Query: "Do relation-aware or heterogeneous graph neural networks over a biomedical knowledge graph outperform
network diffusion, random walk or network proximity baselines for predicting disease phenotypes, drug indications
or gene-disease associations when whole diseases or genes are held out?" No filters were applied. This is the first
search on the relation-typed half of Q1.

The reply listed 20 records, each with a DOI, and carried no sign-up, upgrade or usage notice. It ends with the
block of formatting instructions already noted under the slot 12 search. No abstract contained instruction-like
text. None of the 20 was already on file.

What the set looks like. Most records predict a link of the heterogeneous network they propagate over (a
drug-disease or gene-disease association that sits in the graph as an edge type), under cross-validation over
edges. That is the association-label design: checks 2 and 3 are in doubt for nearly all of them. No abstract names
a diffusion or proximity baseline, so the comparison the query asked for is not visible at this level.

Queued from this search, with the reason:

1. Venkatesh, 2026, Bioinformatics, 10.1093/bioinformatics/btag413. Benchmarks an MLP, GCN, HGT and RGCN on one
   heterogeneous graph and reports "RGCN achieves the strongest overall performance (Macro F1: 0.694, Recall:
   0.720), with relation-specific weight matrices proving the critical factor". The same graph with and without
   relation types is the isolation check 5 asks for, applied to types. The outcome is the class of a drug-drug
   interaction, not a phenotype, and the split is not stated in the abstract.
2. Mellina Andreu, 2024, Artificial Intelligence in Medicine, 10.1016/j.artmed.2025.103177. A heterogeneous graph
   model for phenotype-gene links that "consistently outperforms existing models in both retrospective and temporal
   validation tasks". The outcome is the kind this project predicts and a temporal split is stricter than a random
   one. The abstract does not say whether any baseline is a diffusion.
3. Hu, 2026, Nature Biomedical Engineering, 10.1038/s41551-025-01598-z. Path-based reasoning that learns "by
   considering all relations along paths" and "leverages a background regulatory graph for enhanced message
   passing"; it "outperforms or matches existing methods" on tasks that include drug-disease indication. A typed
   propagation over a background graph kept apart from the predicted links is close to our design. The baselines
   are not named in the abstract.
4. Macaulay, 2024, Proceedings of Machine Learning Research, 10.1101/2024.09.24.614782. A multi-relational graph
   model that "matches and often outperforms traditional single-relation approaches", tested on external data
   sets. Node features come from text embeddings, which changes more than the relation typing.
5. Hu, 2025, Bioinformatics Advances, 10.1093/bioadv/vbaf187. A benchmark of 42 techniques on eight biomedical
   link-prediction data sets. It may hold typed against untyped comparisons at scale, all on transductive tasks.

Added from the appraiser's recall, not returned by any search. The identifier and every claim are unverified until
the slot that reads it:

- Appraised in slot 17 as PMID 39322717; the recalled DOI resolved to this paper and the recalled baselines were
  confirmed in the text. Huang, 2024, Nature Medicine, "A foundation model for clinician-centered drug repurposing" (TxGNN), DOI recalled
  as 10.1038/s41591-024-03233-x. Recalled as a relation-typed graph model on a medical knowledge graph, scored on
  indications and contraindications for held-out disease areas against baselines that include diffusion state
  distance, network proximity and untyped or typed graph networks. If the recall is right it is the most direct
  test of the relation-typed half of Q1 on a clinical label, with check 3 met and check 2 in doubt, since the
  graph holds the indication edges of the training diseases.

Set aside from this search on the abstract, with the reason:

- Gu, 2022 (10.1016/j.compbiomed.2022.106127); Zhao, 2025 (10.1016/j.ins.2024.121360); Zhao, 2021
  (10.1093/bib/bbab515); He, 2024 (10.1186/s12859-024-05705-w); Mary, 2026 (10.25258/ijddt.16.21s.106); Jia, 2024
  (10.1186/s12859-024-05841-3); Zeng, 2026 (10.1109/jbhi.2026.3679534); Borah, 2025 (10.64898/2025.12.25.696543);
  Li, 2023 (10.1093/bib/bbac578); Keichin, 2025 (10.1093/bib/bbaf369): association or interaction prediction on a
  network that contains the predicted edge type, compared with other methods of that kind (checks 2 and 3), with
  no untyped or diffusion variant in the abstract.
- Zhao, 2025, Briefings in Bioinformatics (10.1093/bib/bbaf555): tests "unseen drugs and diseases", so check 3
  holds, but the abstract names no untyped or diffusion baseline.
- Bang, 2023, Nature Communications (10.1038/s41467-023-39301-y): a random walk guided by semantic information
  against link-prediction models, up to 16.8% better on drug-disease association. A typed walk against untyped
  ones, on a transductive task. Next in line if a slot frees up.
- Mastropietro, 2023, Bioinformatics (10.1093/bioinformatics/btad482): a graph network for disease-gene ranking on
  an untyped interaction network. No relation types, so it does not bear on Q1. It bears on our null result from the
  other side: a learned model reported ahead of existing gene discovery methods on an untyped graph.
- Ye, 2025, Advanced Science (10.1002/advs.202412402): drug screening with a heterogeneous graph and expression
  data; several inputs change at once.
- Yue, 2019, Bioinformatics (10.1093/bioinformatics/btz718): a benchmark of untyped graph embeddings. Background
  for the untyped baseline.

Plan for slots 17 to 19, the best three unread entries across this search and the slot 12 queue:

- Slot 17: Huang 2024 (TxGNN), because it is the only candidate expected to run diffusion and proximity baselines
  on a clinical label with whole disease areas held out. If the recalled identifier does not resolve to that paper,
  the slot goes to entry 1 above.
- Slot 18: Venkatesh 2026, the typed against untyped comparison on one graph.
- Slot 19: Mellina Andreu 2024, the phenotype outcome with a temporal split.

Left unread when the cap is spent: Hu 2026, Macaulay 2024, Hu 2025 and Bang 2023 from this search, and Trinh 2016
and Gates 2021 from slot 12. The two from slot 12 dropped below these because both score Boolean models' own
dynamics or drug-target status under a signed logic with no unsigned variant visible in the abstract.

### Note on 39148051

The paper names our problem and does not measure it. Its signed propagation exists precisely
because nodes joined by signed edges disagree in sign, which is the condition measured on our
graph at 36% of destinations, so it is prior art for the aggregator work in task #24 rather
than evidence about it. Three things keep the weight at `weak`. The sign is a property of the
association label, the direction of an abundance change, not of an edge in a causal network.
The outcome is association and sign classification, not a phenotype from a perturbation. And
the abstract reports no unsigned-propagation ablation, so none of the 0.9742 AUROC can be
attributed to the signed strategy; a similarity-feature bipartite graph with XGBoost on top
can reach that figure without the propagation contributing anything.

What it does settle: a dedicated strategy for signed message propagation, in place of a plain
mean over signed neighbours, is published and not novel. Task #24 is therefore a measurement
on a mechanistic network rather than a new idea, which is the weaker and more defensible claim.

### Note on 33906083, and the balance question it closes

The row is `no_bearing` on Q1 and `strong` on weight, which is the point of keeping the two axes apart:
the paper proves its own claim properly and its claim is not about anything we are asking. It is
Lyapunov analysis of pinning control for coupled reaction-diffusion systems whose coupling topology is
a signed graph.

What it contributes is the term for the condition that would make our signed encoder redundant.
Bipartite synchronisation holds when a signed network is structurally balanced, and a balanced signed
graph is switching-equivalent to the all-positive one, so propagation over it differs from unsigned
propagation by a node-wise sign that a linear readout absorbs. Reading that against our construction
settles it on paper, without a run.

Our graph is frustrated at every stoichiometric edge, by construction. Balance requires a node
labelling in plus and minus ones with node_sign[source] * node_sign[target] equal to the edge sign for
every edge. Each substrate_of edge from metabolite m to reaction r carries +1 and so demands
node_sign[m] * node_sign[r] = +1, while the depletes_substrate edge the encoder derives from that same
edge runs r to m with -1 and demands the same product be -1. Every substrate_of edge therefore forms a
two-cycle of negative sign product with its own derived reverse edge: 34,791 such pairs on the
metabolic graph and 36,562 on the neuronal graph. No labelling satisfies them, so the graph is as far
from balanced as a signed graph gets, and switching equivalence does not apply.

That closes the balance thread opened by the first row. The signs are not reducible to an unsigned
graph, so the null result cannot be explained by the signs carrying no structure. It is explained, if
at all, by the encoder discarding the structure they carry, which is what the row-sum measurement in
docs/membrane_potential_reach.md shows.

**Measurement that now has a point:** run `experiments/measure_signed_graph_balance.py` on the derived
relation stack rather than the stored edge table, to report the frustration per relation and confirm
the count above rather than deriving it.

### Note on 39976387, the nearest published analogue so far

This is the closest hit in the queue: signed graph diffusion, gene-phenotype associations, a held-out
link-prediction comparison and a reported gain of up to 9.28% AUC on Gossypium hirsutum. It still does
not answer Q1, for three reasons that are worth separating because each one is a difference between
their setting and ours rather than a flaw in their work.

The sign is a label, not a mechanism. Their edges are positive or negative gene-phenotype
associations, the thing being predicted. Ours are stoichiometric and regulatory signs on a metabolic
graph, which is the substrate the prediction propagates through.

They diffuse over the label graph; we diffuse over a mechanistic graph to produce labels. Theirs is
transductive link prediction on the association network itself, so part of the performance can come
from the structure of the labels. Our split holds out whole genes and whole disease clusters, and the
graph carries no label information at all.

The signed component is not isolated. The gain is reported for the whole method, signed diffusion plus
stochastic perturbation views plus a multiview contrastive loss, against other methods rather than
against its own unsigned variant. Nothing in the abstract attributes the 9.28% to the signs.

Taken with the rows above, Q1 is so far unanswered in the strict form: no paper in this queue reports
signed or relation-typed propagation on a mechanistic biological network beating unsigned diffusion
for a phenotype, with the baseline run and the signed component ablated. Two papers do the comparison
on association or similarity graphs without the ablation, one supplies the balance condition and one
has no bearing. If that holds through the rest of the queue, the strict form of Q1 is an open question
rather than a settled negative, which makes the unsigned random-walk baseline in
docs/graph_content_null_results.md more unusual as a reported result than as a failure.

### The pattern across the five rows, which is more informative than any single one

Three of the five appraised papers (39148051, 39976387, 39523622) share one shape: a signed graph
neural network over a biomedical **association** network, scored on link or sign prediction, with a
headline figure reported for the whole pipeline and no unsigned ablation. None is a weak paper for its
own question; all three are weakly informative for ours, for the same three reasons each time, which
is why the firing instructions now name those three distinctions explicitly.

A fourth observation cuts deeper and is worth stating as a claim rather than a row. Two of those three
build on **balance theory**, 39523622 naming it outright, and 25113603 is about the aggregation
method's restriction to the balanced, or causally consistent, case. Balance is the assumption that
sign products around cycles are positive, and the published signed-propagation machinery for
biological networks leans on it or on its near-miss.

A stoichiometric metabolic graph violates balance by construction, not by accident. Every substrate_of
edge forms a negative two-cycle with the depletes_substrate edge the encoder derives from it, giving
34,791 frustrated edges of 129,577 (26.9%) on the metabolic graph and 38,958 of 143,079 (27.2%) on the
neuronal graph, measured by `experiments/measure_signed_graph_balance.py`. So the balance-theoretic
tools are unavailable here, and the papers that work best on association networks do so partly because
those networks are close to balanced.

That reframes the project's negative result. The question is not only whether signed propagation beats
unsigned diffusion, but whether it can when the signed graph is maximally frustrated, where the
switching trick that makes balanced signed propagation tractable does not exist. Nothing in the queue
so far addresses that case, which makes it the open question rather than a settled negative.

### Note on 39322717, typed message passing ahead of diffusion on held-out diseases

The full text was read from PMC (PMC11645266): main text and methods. Figures, equations and the supplementary notes
are not in the served text, and those hold the per-baseline numbers, the split details and the negative sampling.

What it is. A heterogeneous graph network with "a relationship-specific weight matrix" per edge type, trained on a
medical knowledge graph that is "heterogeneous, with 10 types of nodes and 29 types of undirected edges" (123,527
nodes, 8,063,026 edges), with a disease-similarity module added for diseases that have few neighbours. It predicts
indications and contraindications, of which the graph holds 9,388 and 30,675.

The comparison Q1 asks for is run, for relation types. "We compared TxGNN with eight methods": two divergence
statistics, "graph-theoretical network proximity approach, diffusion state distance (DSD)", three graph networks
(RGCN, HGT, HAN) and a language model (BioBERT). Three splits, in the order the paper presents them:

- Random drug-disease pairs: "three of eight existing methods achieving AUPRC > 0.80 and HAN as the best at 0.873
  AUPRC", against 0.913 for TxGNN.
- Random held-out diseases with all their drugs removed: gains over the next best method of 19.0% AUPRC for
  indications and 23.9% for contraindications.
- Nine held-out disease areas: "relative AUPRC gains of 0.5–59.3% (average 25.72%) across 9 disease areas for
  indications and 11.8–35.6% (average 18.67%) for contraindications."

What the text does not let a reader credit to relation typing.

- The untyped comparators are not trained. Diffusion state distance, proximity and the divergence statistics are
  unsupervised scores, while TxGNN is pretrained on every relation and fine-tuned on the labels. Typing changes
  together with supervision, pretraining and the disease-similarity module.
- The plain typed network is itself a baseline, and it is not always the runner-up: "BioBERT performed best for
  indication prediction in seven of nine disease areas, whereas RGCN was the best baseline for contraindications in
  eight of nine." For held-out indications the strongest competitor reads text and uses no graph. Where the
  diffusion baselines rank is shown only in figures.
- No run merges the relation types. "Ablation analyses confirmed that each component of the TxGNN Predictor is
  essential", with the components and numbers in a supplementary figure.
- The labels are edges of the graph. For a held-out disease area "all drug indications and contraindications were
  removed from the training dataset, along with a fraction of relationships between drug and other medical
  concepts", but the label edges of the training diseases remain, and during fine-tuning "Other relationship types
  remained in the KG to facilitate indirect information flow."

Against the checks. Check 1: there are no signs, and the typed substrate is curated biomedical relations, but the
two predicted relations are themselves edge types. Check 2 fails in part, as above. Check 3 holds: whole diseases
and whole disease areas are held out. Check 4 does not arise. Check 5 fails for typing. Check 6 does not arise.
Check 7 cannot be judged, since the negative sampling is in a supplementary note.

Two details that touch our encoder. The aggregation is the one our encoder uses: "we aggregated on the incoming
messages from neighboring nodes of each relation ... by taking the average of these messages", then combined the
relations. Without signs that average cannot cancel, so the paper is no evidence about the cancellation at
mixed-sign destinations. And the edges are undirected, so this row and 31289271 are evidence for two different
attributes, type and direction, neither of them sign.

The abstract's headline figures, 49.2% for indications and 35.1% for contraindications, were not found in the body
text as served, whose figures are the three sets quoted above. They may come from a figure.

No instruction-like text was found in the article.

Where Q1 stands after slot 17. In its literal wording the relation-typed half now has one supporting row of strong
weight: a typed model ahead of diffusion and proximity baselines on a clinical label, with whole disease areas held
out. In the strict form this file has used (a mechanistic substrate without the labels, the plain baseline run
through the same pipeline and the typed or signed component isolated) Q1 is still open, and for signs nothing
appraised so far reaches it in either form. The chain continues on the strict reading. It also bears on how the
null result is described: in this paper a trained typed model beats untrained diffusion by wide margins, which is
the opposite of our pattern. One paper with the confounds listed above does not settle why, but the null result
should not be written up as a general rule that an unsigned walk extracts all a graph has to give. It may reflect
our graph, labels or split.

### Note on 25015298, a signed method built for networks that are not balanced

The full text was read from PMC (PMC4227138), all sections. Reference numbers, figure numbers and most formula
symbols are missing from the text as served, and the result tables arrive as flat lists of cells, so the table
figures below rest on reading each method name with the seven numbers that follow it.

Check 4, and a correction to this file's running claim. The notes above say that the published signed machinery
leans on balance and that nothing in the queue addresses a frustrated graph. This paper does. Its predecessor was
"restricted to causally consistent networks (e.g., no negative feedback loops are allowed)", and two of the networks
used here are stated to violate that: "the xenobiotic metabolism network is not causally consistent", which
"further motivated our method that does not assume causal consistency of the backbone", and "The cell cycle model
contains many negative feedback loops, hence is not causally consistent". The paper equates the two terms: "a graph
is balanced if and only if all its cycles are positive. This property is called “causally consistent”". The
restriction named in the quoted sentence of the first row (25113603) is the one this paper removes. Whether the
first row belongs to the same framework could not be confirmed, because the reference list is not in the served
text.

How it handles signs. Node values on the signed backbone are fitted by least squares to the gene fold-changes: the
"smoothest" vector under the signed Laplacian with the measured genes as a boundary condition. The network-level
score is a sum of per-edge quadratic terms, chosen because "we should avoid canceling out (“destructive
interference”) cumulative signed edge scores". Two limits, the first stated in the paper and the second an
inference made here from the objective as described:

- Direction is lost: "the directionality of the edges has no more importance at this stage."
- The quadratic score prevents cancellation in the network summary, not at a node. A backbone node whose signed
  neighbours pull with equal strength in opposite directions still receives a fitted value near zero. So this is
  not a ready answer to the encoder's mixed-sign destinations; it is a precedent for carrying a second,
  sign-invariant quantity beside the signed value.

The paper also ties balance to the spectrum of the signed Laplacian (the symbols are stripped, so the exact
statement needs the PDF). The standard result is that the smallest eigenvalue of the signed Laplacian of a connected
graph is zero exactly when the graph is balanced. That eigenvalue, computed on our derived relation stack, would be
a continuous measure of frustration to set beside the count of 26.9% of edges. This is a suggestion from the
literature session, not something the paper applies to a metabolic graph.

Prediction of a clinical label. Patients' expression profiles are mapped to backbone values and classified, with
10-fold cross-validation repeated 5 times and each cohort predicted from the other. Mean G-performance (geometric
mean of sensitivity and specificity) over the two cross-cohort tests:

- Infliximab response in ulcerative colitis (Table 8; 20 responders and 27 non-responders before treatment over
  both cohorts). Backbone values: 0.68, 0.65, 0.83, 0.73 and 0.80 over five classifiers. All genes: 0.79, 0.58,
  0.63 and 0.80. Genes downstream of the single best node: 0.80 and 0.85. The best figure belongs to a gene set,
  not to the signed backbone, and the gaps are within what one or two patients change in cohorts of this size.
- Smoking status (Table 7). Backbone values: 0.86, 0.86, 0.91, 0.85 and 0.83. All genes and transcript-layer genes:
  from 0.00 to 0.88, with several classifiers assigning every sample of the other data set to one class. Genes
  downstream of the single best node: 0.89 and 0.90.

Why `mixed` and `analogous_mechanism`. The backbone features are the most consistent across classifiers for
smoking and are not ahead for infliximab response, where the paper's own text credits "NSC combined with the genes
underlying a single UBE and NSC based on the backbone values" jointly. Neither case compares a signed with an
unsigned propagation. The comparison run is network-derived features against gene-level features, so the signs'
contribution to prediction is not measured. The K statistic comes closest and does not fill the gap: "The edges of
the functional layer are randomly permuted (together with their signs)", and what is recomputed is the perturbation
score.

Against the original checks: 1 holds (signs are curated cause-and-effect relations), 2 holds (the backbone is
fixed before the data are seen: "this linear transform does not depend on the data"), 3 holds (whole cohorts held
out), 4 holds (no balance assumption). Check 5 fails, check 6 holds and check 7 is met by reporting sensitivity and
specificity separately.

No instruction-like text was found in the article.

### Note on 31289271, the first row where the plain network is the baseline

This is the first appraised source that runs the same propagation on a network with and without a mechanistic edge
attribute and scores both on a ranking outcome. The full text was read from PMC (PMC6617457), so the quotes below
are from the source. Figure numbers and the exponents of p-values are missing from the text as served.

What the design gets right. The attribute is added to a physical interaction network, not to the labels (check 1).
The scored outcomes are not edges of the network (check 2). The plain network is run through the same diffusion
(check 5). The orientations are kept apart from the outcome they are scored on (checks 3 and 6): for drug targets "We
oriented the network using only the cancer data sets", and for driver genes "We oriented the network leaving out one
disease set at a time". Check 4 does not arise, because the network carries no signs.

What limits it.

- The attribute is direction. The network carries no signs, so the paper bears on the relation-typed half of Q1 and
  says nothing about activation against inhibition.
- The baseline is below chance. A random ranking of 15,501 genes has an expected mean rank of 7,751 (computed
  here). The unoriented network scores 8,837, which is 1,086 ranks worse than random, and the oriented network
  scores 6,045, which is 1,706 ranks better. The reported gain of 2,792 ranks is therefore mostly recovery from a
  baseline that ranks true targets in the lower half. No chance level, top-k score or area measure is given for the
  drug-target test.
- The baseline may remove more than the inferred orientations. The network includes 33,756 curated directed
  interactions, and the comparison is against "a completely unoriented network". If the curated directions are
  stripped as well, the test does not show that the inferred orientations help.
- The driver-gene test gives figures for the oriented network only, for example "49% of the top 1% of the genes
  reported by the orientation-based computation are in the AGO positive list", and states that it "consistently
  reported higher fractions of driver genes and lower fractions of non-driver genes compared to the unoriented one"
  with the comparison shown in a figure. The significance tests printed are against the gene background, not
  against the unoriented ranking.
- Both tests diffuse against the orientation: "we flipped all directions in the (oriented) network, so that a
  diffusion process will follow a signal traversing up-wards in the network". The gain is specific to walking from
  effects back to causes.

Why `moderate`. The evaluation is held out and the baseline at issue is run, but that baseline performs below
chance and no spread over drugs is printed in the text.

What it suggests for our comparison. Direction and sign are separate attributes, and this paper is evidence for
direction alone. The note on our null result (docs/graph_content_null_results.md) does not say whether the random
walk with restart follows edge direction or symmetrises the graph. If it follows direction, the gain reported here
is already inside our baseline. If it symmetrises, a directed unsigned arm is the comparison that separates
direction from sign, and the signed encoder has so far been compared with a baseline that lacks both.

No instruction-like text was found in the article.

### Note on 10.64898/2026.04.29.721775, the stated comparison and why it does not transfer

This preprint (bioRxiv version 1, posted 4 May 2026, no journal version on record) is the first appraised source that
runs unsigned baselines beside a signed model on something it calls a clinical label. The abstract was read from the
bioRxiv record. The body was read through a fetch tool that summarises the page and covered the first 100,000 of its
100,163 characters, so the items below are relayed and each needs checking against the PDF before use.

- The clinical response task is "formulated as an edge-sign prediction task in both binary (indication, −1 vs.
  contraindication, +1)" and a three-class form. The labels are 11,896 indication and 148,061 contraindication edges
  of the knowledge graph, split "using a stratified split (0.8/0.1/0.1), with stratification performed only by
  sign".
- Unsigned baselines (BGRL, DGI, InfoGraph, GraphSAGE): "Homogeneous baselines were trained on an unsigned version
  of SIGMA-KG with all edge types and polarities discarded." Relational baselines (HGATE, DMGI): "Heterogeneous
  baselines utilized the full 12-relation multiplex structure." Signed baselines (SGCN, SNEA, SDGNN): "signed
  baselines and FLASH received the signed network, with edges labeled as +1 or -1."
- The task graph keeps the training label edges: "the training compound–disease edges and all remaining relation
  types in SIGMA-KG were retained." Validation and test edges are removed.
- Binary clinical response: FLASH accuracy 0.9893 ± 0.0003, 1.51% above HGATE (p = 0.0001). Three-class: the signed
  baseline SGCN is 0.22% above FLASH in accuracy (p = 0.002). Per-baseline values sit in supplementary tables that
  were not read. Three runs; standard deviations appear for FLASH only in the text that was read.
- Drug-drug interaction, the one setting with held-out entities ("a drug-disjoint split") and an outcome that is
  not a sign: on DrugBank FLASH has AUROC 0.7960 ± 0.0157 with p = 0.20 against the best baseline, and on MUDI the
  signed baseline SDGNN has accuracy 0.6989 against 0.6934 for FLASH (p = 0.445).
- No run with signs removed, made all positive or permuted was found, and none with the balance component removed.
- "Structural balance is enforced by pruning low-confidence edges involved in unbalanced signed cycles with sign
  product < 0." After pruning, "Across closed triads in the atlas, 99.6% were balanced".
- The abstract's 69.6% is 16 of 23 newly prioritised drug-disease pairs found in approval records afterwards, not a
  held-out score.

Why the direction is `mixed`. Where the label is itself a sign and the split is over random edges, the signed and
relational models are ahead of the unsigned ones. Where the outcome is unsigned and whole drugs are held out, the
lead over the best baseline is not distinguishable from zero. Three runs cannot separate no gain from a small one,
so that null is weak evidence as well, and the direction says only that the paper's own conditions disagree.

Why none of it reaches Q1. Checks 1 and 2 fail in the way they did for 39976387: the sign is the label, and the
model propagates over the graph that holds the training labels. A test edge between a drug and a disease sits among
training edges of that drug and that disease whose signs the signed models see and the unsigned baselines do not.
Check 3 fails for the sign tasks, which split random edges, and check 4 fails because balance is imposed by pruning.
Check 5 fails because the unsigned baselines drop edge types along with signs and use other architectures, so three
things change at once. Check 7: always predicting contraindication scores 0.9256 (148,061 of 159,957, computed
here), and no such baseline is printed beside the 0.9893.

What it adds for Q2 and for our graph. This is one more appraised source whose signed method lives in the balanced
regime, after 25113603, 39523622 and 38392416, and the first to produce that regime by deleting the edges that
violate it. On our metabolic graph the
violating edges are the derived depletes_substrate reverses, 26.9% of all edges, each in a negative two-cycle with
its own substrate_of edge. Pruning to balance would remove the stoichiometry the encoder exists to represent. The
fraction of nodes that receive both signs is not reported.

Odd text in the source, none of it instruction-like: unresolved "Table ??" references, a "[?]" citation, a
truncated reference and an empty section headed "Signed network balance semantics". The figure captions describe
paired t-tests and the methods describe unpaired Welch tests.

### Note on 17646318, the third direction task and a precedent for abstaining

SPINE is the third of the four mechanism-named candidates whose test is the direction of a knockout's effect, after
29925605 and 41429577. The same limit applies: a network without signs has no direction to predict, so the
comparison Q1 asks for cannot be run on this readout. The row is `supports` for the premise only, that signs placed
on a physical interaction network carry the direction of a perturbation's effect on expression.

The abstract was read from the PubMed record. Everything else comes from the publisher page through a fetch tool
that summarises the page, relays quotes of at most about 125 characters and stopped at "Interestingly, the" in
section 3.2. The items below are therefore relayed, not read, and each needs checking against the PDF before use.

- Data: the mating subnetwork defined by Yeang et al. (2004), "33 protein–DNA interactions, and 25 PPIs, covering a
  set of 149 knockout pairs", with expression from Hughes et al. (2000).
- Design: "at each iteration one or five knockout pairs were hidden", and "The accuracy represents the percentage
  of correct predictions among predictions made."
- Table 1 as relayed, in percent correct, incorrect and undecided. Leave-one-out, 103 trials: Yeang et al. 97.1,
  2.9 and 0; SPINE edge variant 98, 1 and 1. Leave-five-out, 200 trials: Yeang et al. 96.5, 3.5 and 0; SPINE edge
  variant 96.9, 0.4 and 2.7.

Those percentages correspond to whole counts, checked by running the arithmetic. Leave-one-out, 103 held-out pairs:
100 correct and 3 incorrect for Yeang et al.; 101 correct, 1 incorrect and 1 undecided for SPINE. Leave-five-out,
1,000 held-out pairs: 965 correct and 35 incorrect against 969 correct, 4 incorrect and 27 undecided. Two points
follow, both conditional on the relayed table being right.

The abstract's 99% is accuracy among predictions made (101 of 102, and 969 of 973). Counted over all held-out pairs
the two signed methods stand at 98.1% against 97.1% and at 96.9% against 96.5%, a difference of one pair in 103 and
four in 1,000. In the leave-five-out test SPINE has 31 fewer incorrect pairs than the baseline, and 27 of that
difference appears as undecided, not as correct. The counts are not paired, so they do not show which pairs moved.

That is a published precedent for the design question raised on 7 October, what to do at a destination whose
incoming signs conflict. A held-out pair is scored as predicted only when all of its explanatory paths agree
(relayed in paraphrase, not as a quote), and the rest are reported as undecided instead of being resolved toward one
sign. On the mating subnetwork that withholds 1 to 2.7% of pairs. On our metabolic graph the same rule would
withhold a prediction at the 36.6% of destinations that receive both signs. The precedent supports reporting
coverage beside accuracy; it does not show that abstaining is affordable here. The low undecided rate also fits the
earlier pattern, signed methods developed where sign conflicts are rare, though the paper as read reports no
balance or frustration measure.

Two cautions on the 99%. No count of up-regulated against down-regulated effects and no majority-direction baseline
was found in the part of the page that was read, and the hold-out unit is a pair, so the other pairs of the same
knockout and of the same affected gene stay in the fit. The cross-validated figure also belongs to the edge variant,
run for comparison with Yeang et al.; for the whole-network application (861 paths, 183 genes) the abstract reports
"high agreement with current biological knowledge", not a held-out accuracy.

No instruction-like or filler text was found in the abstract or reported from the page.

Where Q1 stands after ten rows. The replacement queue is finished and the strict form of Q1 is still unanswered.
Three of its four candidates scored direction, where unsigned diffusion has nothing to predict. The unsigned
propagation papers set aside from the same search (Cowen 2017, Lee 2008) score ranking, where a sign is not needed.
The two bodies of work are evaluated on different tasks, so a query asking whether signed beats unsigned keeps
returning one or the other. The next search therefore names an outcome both can be scored on (gene essentiality,
drug response, disease-gene or phenotype ranking) together with both comparators.

A second gap. Every slot so far has gone to the signed half of Q1. No search has covered relation-typed
propagation, for example message passing over a typed biomedical knowledge graph compared with diffusion or
proximity baselines on a held-out disease split. One of the remaining slots should.

### Note on 41429577, and the rarity of negative edges

The full text was read because its telomere-length validation is a phenotype readout. It does not reach Q1 for the
same reason as 29925605: the reconstruction predicts the direction of a knockout's effect from the fraction of
negative shortest paths, and a network without signs has no direction to give. The headline AUPRC of 0.8 cannot be
credited without the positive class and its base rate, and the paper's own pair counts make the base rate decisive
(456 'long' against 37 'short').

It bears on one design question raised on 7 October: whether exact cancellation at a node should be resolved toward
the negative message because negatives are rare. The paper's premise is the same rarity ("we hypothesize that negative
interactions are less prevalent than positive ones in physical interaction networks"), supported by curated counts in
signalling databases (yeast KEGG 136 positive against 53 negative; human kinase-phosphatase 2,506 against 923). Those
counts describe curated signalling edges, where a database records what was studied. In our graph the negatives are
structural, the derived depletes_substrate reverse of every substrate edge, so their frequency is set by stoichiometry
and not by biology's preference, and the rarity argument does not transfer.

The methods section of the PMC text contains a stray filler pangram ("The quick brown fox jumps over the lazy dog.",
twice). It is an editorial leftover with no bearing on the row, recorded only because the user asked that any odd or
instruction-like text found in sources be flagged.

### Note on 29925605, the first mechanistic-network row

This is the first row whose network is mechanistic rather than an association graph: each of the 87 models has
known kinetics, so the truth the topology-only prediction is scored against is the model's own response, not a
label that leaks through the graph. It passes checks 1 and 2 cleanly, and check 3 does not arise because nothing
is fitted. It still does not reach Q1, for a reason worth keeping: on a direction task an unsigned diffusion has no
prediction to make, so the comparison Q1 asks for only exists for the strength of a response or for a phenotype
ranking, and the abstract reports neither against a sign-blind model. Its headline accuracies also need the
majority-direction rate beside them before they can be credited. Our sign check on the HPO laboratory labels
(docs/linear_response_sign_check.md) has the majority direction at 0.806, above the untrained signed response at
0.615, so an accuracy near 80 percent is uninformative until the base rate is known. The paper is the best support so
far for the premise behind the linear-response encoder, that topology without kinetics carries the sign of a
perturbation's effect, and no support yet for the claim that this sign improves a phenotype prediction.

### Note on 38392416, and a problem with the queue

The paper is US foreign exchange rates and has no bearing on anything here; PubMed indexes it because
it appeared in Entropy. It is kept for one line in its abstract, that negative cross-correlations have
become rare in that network, and for identifying balanced strong triads. That is the third appraised
paper whose method is developed on a signed network which is negative-sparse and close to balanced,
after the two association-network papers that invoke balance theory. The signed-network toolkit was
built in that regime, social and financial networks among them, and a stoichiometric metabolic graph
is in the opposite one: 26.9% of its edges are frustrated and every stoichiometric edge is frustrated
by construction.

The queue itself is the problem now. `("signed network" OR "signed graph" OR "signed networks") AND
(propagation OR diffusion OR "random walk")` in title or abstract returns 19 records, of which six are
appraised and three of the remaining were already set aside on their abstracts. The phrase is generic,
so the recall is poor in the direction that matters: a paper that propagates inhibition through a
signalling network and never calls its graph signed will not match. Before the cap is spent on the
remainder, the query should be replaced with ones naming the mechanism rather than the formalism, for
example inhibition or repression together with propagation and a phenotype, or metabolic or signalling
network together with perturbation and prediction, and Consensus should be used beside PubMed since it
searches semantically rather than by phrase.
