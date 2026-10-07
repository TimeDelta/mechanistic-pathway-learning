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
(propagation OR diffusion OR "random walk")` in title or abstract. Appraised: 25113603.
Set aside as out of scope on the abstract, with the reason:

- 28297881, 31661504 — epidemic and information spread on social signed networks; no
  biological network and no phenotype readout
- 28129187 — antisynchronization control of coupled reaction-diffusion neural networks; a
  multilayer signed graph used as a coupling topology, not a predictor
- 35131567 — signed network representation scored on sign prediction and link prediction,
  not on a phenotype label; keep for Q2 if the node-proximity metric reports sign mixing

Remaining to read: 39148051, 39976387, 39523622, 38392416, 35364845, 40103114, 33906083,
40680519, 36049951, 42481621, and the four records beyond the first fifteen.
