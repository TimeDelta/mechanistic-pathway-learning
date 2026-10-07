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

### Triage checks, reconstructed in the literature session

The replacement-queue section refers to "the four checks recorded at the top of the rows" and the note on 29925605
cites checks 1 to 3 by number, but the checks are not in this file or in the git history: they were part of the
experiment session's firing instructions. The list below is a reconstruction from the notes in this file, used for
triage from slot 11 on. Checks 1 to 3 follow the three distinctions in the note on 39976387 and agree with how the
note on 29925605 uses the numbers. Check 4 is inferred from the slot-7 triage remark that signs were "fitted to the
knockouts they explain" and is the least certain. Check 5 is added here from the base-rate remarks in the notes on
29925605 and 41429577. The author should replace the list with the original wording where it differs.

1. The sign belongs to a mechanistic edge (activation, inhibition, stoichiometry), not to the association label
   being predicted.
2. The propagation runs over a substrate that carries no label information and the scored units are held out, so
   the result cannot come from the structure of the label graph.
3. The signed component is isolated from the other fitted parts of the pipeline: the same pipeline is run without
   signs or with signs permuted, on the same network and split.
4. Where signs are fitted, they are not fitted to the outcome they are scored on.
5. A direction or classification figure is reported beside its base rate, the score of always predicting the
   majority class.

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
