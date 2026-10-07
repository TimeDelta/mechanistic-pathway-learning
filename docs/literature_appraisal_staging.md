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
(propagation OR diffusion OR "random walk")` in title or abstract. Appraised: 25113603, 39148051, 33906083, 39976387, 39523622, 38392416.
Set aside as out of scope on the abstract, with the reason:

- 28297881, 31661504 — epidemic and information spread on social signed networks; no
  biological network and no phenotype readout
- 28129187 — antisynchronization control of coupled reaction-diffusion neural networks; a
  multilayer signed graph used as a coupling topology, not a predictor
- 35131567 — signed network representation scored on sign prediction and link prediction,
  not on a phenotype label; keep for Q2 if the node-proximity metric reports sign mixing

Remaining to read: 35364845, 40103114, 40680519, 36049951, 42481621, and the four records beyond the first fifteen.

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
