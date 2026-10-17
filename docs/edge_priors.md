# Edge priors: evidence strength as a prior on whether an edge is there

Record, written 16 October 2026. The design of decision 9.57 as built: `mechanistic_pathway_learning/graph/edge_evidence.py`
(what each edge's evidence is), `mechanistic_pathway_learning/models/edge_gates.py` (the gate and its penalty), how the two
encoders read a gate, the wiring into `experiments/run_main_model.py`, the controls, what a run records and what is
provisional. `docs/open_decisions.md` 9.57 is the log entry; this page does not say what is current.

## 1. What the owner asked, and what stays from the earlier rulings

The owner, 16 October 2026: "That was not MY ruling. Must have been misunderstood. I've ALWAYS wanted differently
weighted prior observations based on evidence strength (not just labels) from the very beginning of this project.
Please find a way." The three rulings the record had generalised against this say something narrower (9.57): the
occupancy placeholder (a constant standing in for a measurement), the threshold grid ("a poorly correlated replacement
for dispersion level") and the overexpression coupling (9.20: a non-physiological context is excluded, not
down-weighted). What stays from them: **a confidence is not a transmission magnitude.** An edge prior never sets how
much signal an edge carries. It states how sure the sources are that the edge is there, and how far the fit may move
that belief.

## 2. The evidence of an edge

`experiments/build_edge_evidence.py` writes `<graph>/edge_evidence.parquet`, one row per row of the graph's
`edges.parquet`: the edge's family (read from its `evidence_source`), the kind of its evidence, its count of distinct
records, whether IntAct curators recorded the pair as not interacting (`contested`), and the edge's own `source_id`,
`target_id` and `relation_type`, which the trainer checks against the graph before it reads a row. Its summary records
the SHA-256 of the `edges.parquet` it was built from. The graph files are not touched.

| family | kinds | records |
| --- | --- | --- |
| signalling (OmniPath) | curated, text_mined_only (references only from ProtMapper, which assembles statements read by the REACH and Sparser text-mining readers), unreferenced, unmatched | distinct PubMed ids, counted apart for ProtMapper |
| transcription (CollecTRI) | curated, text_mined_only (supported only by the ExTRI layers), unreferenced, unmatched | distinct PubMed ids |
| enzyme_regulation (BRENDA) | in_vitro | the reference count the build wrote into `evidence_source` |
| curated (this project's tables) | curated | the PubMed id `evidence_source` names |
| small_molecule, other | none | not read |

On `graph_full_neuronal_split_binders_drugs_v4` (268,954 edges), from its `edge_evidence_summary.json`: 84,108
signalling curated, 32,817 signalling unreferenced, 491 signalling text-mined only, 700 signalling unmatched (their
endpoints were renamed by the build and match no raw row), 27,875 transcription curated, 13,654 transcription
text-mined only, 109 transcription unreferenced, 4,134 BRENDA, 88 curated, 440 small-molecule and 104,538 other (the
Human-GEM metabolic edges, whose reaction confidence score is 0 for 12,828 of 12,877 reactions and so discriminates
nothing, and the definitional and derived layers: encodes, the membrane-potential and ion-pool layers, drug mechanism
edges with their own magnitude, plasma binders, cofactors, Reactome's neuronal layer). 54 edges are contested.

## 3. The prior, the gate and the penalty

`configs/edge_prior_tiers.csv` gives each (family, kind) a prior probability that the edge is present and a unit
strength. Every value in it is provisional, to be chosen on the development folds as the label tiers are (9.50):

| family | kind | prior probability | unit strength |
| --- | --- | --- | --- |
| signalling | curated | 0.9 | 1.0 |
| signalling | unreferenced | 0.7 | 0.5 |
| signalling | text_mined_only | 0.6 | 0.5 |
| transcription | curated | 0.85 | 1.0 |
| transcription | unreferenced | 0.7 | 0.5 |
| transcription | text_mined_only | 0.55 | 0.5 |
| enzyme_regulation | in_vitro | 0.6 | 0.5 |
| curated | curated | 0.9 | 1.0 |

The `transcription, unreferenced` row was added on 16 October with the trainer wiring: without it those 109 edges were
ungated, which reads an edge with less evidence as more certain (gate 1) than a curated one (gate 0.85). It takes the
values of the signalling row of the same kind. A (family, kind) without a row is ungated: its edges carry exactly what
they carry today.

For a gated edge `e` with prior probability `p` and `r` records, the pseudo-count is

    n_e = s * unit_strength * min(max(r, 1), cap) ** tempering

with `s` the global `--edge-prior-strength`, `cap` 2 and `tempering` 0.5 as on the label side (a record count of 0
gives the tier's unit strength once, so an unreferenced tier is not a zero prior). A curated IntAct absence is a Beta
update: alpha = `p n`, beta = `(1 - p) n` + `s`, and the edge's prior becomes alpha / (alpha + beta) at strength
alpha + beta.

Each gated edge carries a gate `g_e = sigmoid(theta_e)`, initialised at `p_e`, which multiplies the edge's weight
wherever an encoder reads it. The training loss adds

    penalty = sum_e n_e * KL(Bernoulli(p_e) || Bernoulli(g_e)) / W

**Why this form.** `n KL(p || g) = -n p log g - n (1 - p) log(1 - g) + constant`, which is the negative logarithm of a
Beta(`n p + 1`, `n (1 - p) + 1`) density on the gate: the posterior of a uniform prior after `n p` pseudo-observations of
the edge present and `n (1 - p)` of it absent. Its mode is `p`, so without data the gate stays at its prior, and `n` is
how many observations the data must outweigh to move it. A gate never reaches 0 or 1 (the prior is held inside
[1e-4, 1 - 1e-4]), so no edge is deleted and the sparsity of every adjacency is fixed.

**Why `W`.** The binary cross-entropy is a weighted mean over a batch (`sum w l / sum w`), so the objective it estimates
is the full weighted sum divided by `W`, the total loss weight of the pairs the fit trains on. Dividing the penalty by
the same `W` makes one pseudo-count weigh as one unit of label weight. `W` is computed per phase (the fit, then the
refit on training and validation rows) from the weights the phase trains on, soft priors included.

**The optimiser.** The gates' logits are a parameter group of their own, at `--edge-gate-learning-rate` (the main
rate when unset) and without weight decay: decay pulls a logit toward 0 and every gate toward 0.5, which is a prior of
its own that no source states. A gate the data does not reach has no data gradient and stays at its prior (checked
under AdamW at rates 1e-3 and 1e-2: no movement in 3,000 steps).

## 4. How the two encoders read a gate

**Linear response.** Since decision 4.11 an edge's weight enters as a share at its target,
`sign * w / sum(w at the target)`; the gate multiplies the weight on both sides, `sign * w g / sum(w g)`. Equal gates on
a node's inputs cancel, so only the relative certainty among a node's inputs matters, and a node with one input is not
affected by its gate at all. The gated share is the ratio of expected weights, a first-order (mean-field) reading of
the expected share of an uncertain input, not its exact expectation. The adjacency's values are recomputed from the
gates each forward over entries laid out once at construction; the spectral normalisation is refused with gates,
because its power iteration is not repeated per step.

**Message passing.** The adjacency divides by the count of incoming edges, so the gate scales its edge's message: with
the divisor fixed at the graph's in-degree, `g` times the message is the expected message of an edge present with
probability `g`. The conjunction relation (`member_of`, a minimum over members), the drug entry edges and the
change-only edges (sequestration, clearance) are aggregated apart from the stacked relations and are not gated. The
gated values are rebuilt from the gates in every propagation over the compressed layout read once from the static
adjacencies; the encoder refuses adjacencies other than those its gates were attached to, so a rewired graph needs a
fresh encoder, as `experiments/evaluate_on_rewired_graph.py` builds one.

The two encoders therefore differ: uniform uncertainty about a node's inputs does nothing in the linear-response
encoder and shrinks the node's input in the message-passing encoder. That follows from how each normalises its weights,
which predates the priors. The run records, per encoder, how many gated edges it reads
(`gated_edges_the_encoder_reads`).

**The product.** Both encoders multiply by their gated adjacency through `edge_gates.SparseValuesProduct`. PyTorch's
backward of `torch.sparse.mm` with respect to a sparse matrix forms the dense product `grad @ dense.T`, [rows, columns],
before masking it to the matrix's pattern. While the values are fixed, in every run without gates, no gradient is asked
of them and nothing dense is formed; once the values carried the gates, one linear-response step on the v4 graph asked
for 15.8 GB, a [K, N] float matrix per propagation round, which the toy graphs of the tests were too small to show. The
product gathers each stored value's gradient directly, `sum_m grad[row, m] * dense[column, m]`, in chunks of at most
2^25 floats. It is tested against the dense product with `gradcheck` and on a 100,000 by 100,000 pattern whose dense
gradient would need 40 GB.

**Cost**, measured 16 October 2026 on `graph_full_neuronal_split_binders_drugs_v4` with two CPU threads
(`--timing-batches 4`, median step): the linear-response encoder at batch 32 takes 3.15 s per step without gates and
4.65 s with them (an estimated 2.0 and 2.9 minutes per epoch); the message-passing encoder, at two fixed layers and batch
4, 0.89 s and 1.25 s. The registered message-passing configuration ran out of memory in that 7 GB container at its
default batch size and at batch 4, with gates and without, so it was timed in the lighter configuration.

**With drug entry nodes**, the v4 default, the linear-response encoder delivers the mechanism edges apart and lays out
its adjacency on the remaining edges; its gates index the graph's edges, not the remaining ones, so a gate cannot shift
onto its neighbour when the mechanism edges are taken out. `tests/test_edge_gates.py` holds a case with the mechanism
edge first and its gate nearly closed.

## 5. Off, on and the controls

`--edge-prior-strength 0`, the default, builds no gate: the run is the registered model bit for bit and its
configuration fingerprint is unchanged, with every other edge-prior flag set (tested end to end). Above 0:

- **The permuted control** (`--permute-edge-priors`): the priors of the gated edges shuffled among the gated edges of each
  relation, each prior's probability and strength moving together. A relation keeps its distribution of priors and loses
  which edge holds which; the real arm must beat it for the evidence to carry edge-specific information.
- **The rewiring control** keeps each prior at its edge position, as it keeps the signs: the rewired arm has the same
  priors per relation and per source node, and only the wiring is randomised. The priors are read for the graph's own
  edges before the rewiring.
- **The arms** in `experiments/run_main_model_batch.py`: `v4_linear_response_edge_priors`, `v4_message_passing_edge_priors`
  and their `_permuted` controls, at strength 1.0. They are secondary arms; nothing about them is a primary reading.

The tier values, the global strength and the gate learning rate are to be chosen on the development folds once they are
finalised, beside the soft-prior grid (the owner, 16 October: the grid "should run after finalization of the dev folds").

## 6. What a run records

`results.json` gains `edge_priors`: the table's path and SHA-256, the tier table's path and SHA-256, the strength,
tempering and cap, the counts by family and kind, the contested count, the permutation (seed and a digest of the
shuffled priors) or null, the penalty denominator per phase, and a summary of the scored model's gates over the edges its
encoder reads (mean prior, mean gate, mean and largest move, how many moved by more than 0.1, how many closed below 0.5
against how many started there, and the penalty). Each epoch of the history and of the refit history records the
penalty's value at its end (`edge_prior_penalty`, part of `train_penalty`). `edge_gate_values.npy` beside the
predictions holds the scored model's gate per edge, NaN where an edge has none. The fingerprint adds the strength and
the settings and the digests of the two tables, only while the strength is above 0.

## 7. What a closed gate does and does not say

A gate the fit closes says that, at this strength, the fit predicted the labels better with less of that edge's signal.
It is not evidence that the edge is absent: an edge can be real and silent for these symptoms, and an encoder that
normalises by shares cannot move a node's only input. Reading closed gates as a finding about the graph needs the
permuted control, agreement across folds and seeds, and agreement between the two encoders; until then they are a
diagnostic of the fit.

## 8. Not built, and why

- **The 440 small-molecule edges and the 700 unmatched signalling edges** stay ungated: the build does not read the
  small-molecule rows' references, and the unmatched edges' endpoints were renamed by the graph build so they match no
  raw row. Both are listed in the summary of every table.
- **The metabolic layer**: Human-GEM states no usable per-reaction confidence (above), so its edges have no prior to
  carry. A prior there would need another source, such as the enzyme's measured presence in the tissue, which is a
  question for the owner and not a tier.
- **A prior on edge sign** is not part of this: the sign is the graph's, and a gate never flips it.

## 9. Correction, 16 October 2026: why 700 signalling edges are unmatched

Sections 2 and 8 say the 700 unmatched signalling edges "were renamed by the graph build so they match no raw row".
Measured against the raw OmniPath table and HGNC's complete set, that is true of 7 at most: only 7 have an endpoint
whose HGNC previous symbol or alias appears in the raw table. 81 have both symbols in the raw table under the same
names, so another matching rule failed for them (a direction or a kind). The examples are gene families the build
expands, the histone cluster genes (`H4C14`, `H3C11`) and a node naming two genes (`CCL3L1; CCL3L3`). Matching them
needs the build's own expansion of a raw row into graph nodes, read edge by edge; it is not a renaming table.

## 10. Addition, 16 October 2026: the tiers under hard gates, and what planted worlds say about them

Decision 9.62 (e) turned each tier's (p, n) into a binary prior π = sigmoid(n logit p) and made logit π the cost of
closing an edge: 2.20 nats for a curated edge (p 0.9, n 1), 0.42 for an unreferenced one (0.7, 0.5) and 0.10 for a
text-mined one (0.55, 0.5). The simulations of `docs/edge_gate_simulation.md`, in worlds with these tiers in the v5
proportions, bear on the tier values in three ways. First, the cost sets where closures fall: with calibrated priors
165 of the 234 closures in five worlds were of text-mined edges, where 0.10 nats is a bar that noise clears often, and
95 of them were false (0.58, against 0.47 of the tier's edges). Second, the tier values weigh as much as the data:
when the truth had more real edges than the tiers say (curated 0.97, unreferenced 0.85, text-mined 0.75 real), 33 of
the 47 closures were of real edges, while the closures with at least 1 nat of evidence held up better (9 of 12
false). Third, taking the curated tier as known (decision 9.63) changed nothing measurable in the other tiers. The
tier values are therefore a first-order input to what the gates close, which decision 9.63 (d), tier values estimated
from data as a control, addresses; until then every value in `configs/edge_prior_tiers.csv` stays provisional, as
section 3 says.

## 11. Addition, 17 October 2026: curated edges taken as known (decision 9.63 (b))

With edge priors on, the trainer now takes as known, with no gate, every edge of a curated kind that has at least two
independent records (distinct PubMed ids), at least one of them human evidence, and no IntAct absence
(`curated_edge_labels` in `graph/edge_evidence.py`; `--curated-edge-labels`, on by default). A paper is human evidence
when NIH iCite's human score, the fraction of its MeSH terms in the human category, is above 0
(`experiments/fetch_icite_species.py`). Such an edge carries no prior and no pseudo-count, as an ungated family does
(section 3), and the permutation control of section 5 leaves it in place. On the v5 graph 38,435 edges qualify, 23,359
signalling and 15,076 transcription, and 124,841 edges stay gated. `edge_evidence.parquet` now also holds each edge's
PubMed ids, the resources that state it and each paper's species scores, which the review list of the closures reads
(decision 9.63 (c)). The species is the paper's, not the interaction's, and a transcription edge's ids can include
ids of the ExTRI text-mining layer, since CollecTRI attributes most of its ids to itself; both are put to the owner in
`docs/open_decisions.md`, 9.63.

## 12. Addition, 17 October 2026: the evidence behind a sign, apart from the pair's (decision 9.70)

Section 11's rule counts the papers behind a pair of nodes. OmniPath records the support for an interaction's positive
sign, its negative sign and its direction apart (`fields=evidences`), and CollecTRI states for each regulation how it
decided the sign: from a paper ("PMID"), by a default when no paper stated one ("default activation") or from the
factor's role across its regulon ("regulon"). `graph/sign_evidence.py` reads the two evidences downloads
(`scripts/fetch_raw_sources.sh`; SHA-256 in `docs/data_sources.md`) and gives each signed signalling and transcription
edge the PubMed ids behind its stored sign (`sign_pubmed_ids`, `sign_records`; the text-mining layers left out as in
section 2) and, for a regulation, CollecTRI's decision (`sign_basis`). `experiments/build_edge_evidence.py` adds the
three columns when the files are present. The label rule reads them only when asked:
`--curated-label-minimum-sign-records N` (off by default) keeps a signed signalling label only with N ids behind its
sign, and a regulation label only when CollecTRI decided its sign from a paper and N of CollecTRI's ids stand behind
it. A binding has no sign and the rule does not concern it; a signed edge the downloads do not hold fails the rule.

On v5, 40,847 of the 129,040 signed signalling and transcription edges have no PubMed id behind their sign. Of the
14,680 signed signalling labels, 3,942 have none, 6,124 have one and 4,614 have two or more. Of the 15,076 regulation
labels, CollecTRI decided the sign from a paper for 10,573, by default for 3,673 and from the regulon for 617, and 213
carry more than one decision. With a minimum of one record the rule keeps 29,551 of the 38,435 labels (19,417
signalling, 8,679 of them bindings, and 10,134 transcription); with two it keeps 20,007 (13,293 and 6,714).

Two limits. The rule counts the ids a source puts behind a sign and cannot tell whether a paper states that sign:
CollecTRI2 cites Shiio et al. 1992 behind TP53 activating RB1, and the paper's abstract reads "This study shows that
wild-type p53 suppresses transcription of the retinoblastoma (Rb) gene" (decision 9.72). And for a regulation the rule
reads only CollecTRI's 2023 release, the one entry that states a sign decision. OmniPath also holds a later release,
CollecTRI2, which states none; the two state opposite signs for 3,795 of the 41,638 regulation edges, and one of them
states both signs for 92 more. The sign the graph stores for those is OmniPath's consensus, which follows CollecTRI2 in
2,445 of the 3,795 (decision 9.72). The table marks the 3,887 (`sign_disputed`), and
`--curated-label-exclude-disputed-signs` (off by default) leaves them gated: 2,019 of the 38,435 labels are among
them, and with a sign minimum of one record the rule keeps 29,306 labels instead of 29,551.

## 13. Addition, 17 October 2026: the three label options measured together (decisions 9.70 and 9.72)

The trainer now carries three label options, each off by default and each recorded in the run's labels summary:
`--curated-label-minimum-sign-records N` (section 12), `--curated-label-exclude-disputed-signs` (a regulation whose
two CollecTRI releases dispute its sign is left gated; 9.72) and `--curated-label-bulk-papers` (9.70 (3)): a PubMed
id behind more than 100 signalling or transcription edges, of which v5 has 114 behind 48,887 edges, counts one record
each as before (`count`), all of an edge's together as one (`one`) or not at all (`none`;
`records_with_bulk_papers` in `graph/edge_evidence.py`). The labels each combination leaves on v5:

| bulk papers | sign minimum | disputes left gated | labels | signalling | transcription |
| --- | --- | --- | --- | --- | --- |
| count | - | no | 38,435 | 23,359 | 15,076 |
| count | - | yes | 36,416 | 23,359 | 13,057 |
| count | 1 | no | 29,551 | 19,417 | 10,134 |
| count | 1 | yes | 29,306 | 19,417 | 9,889 |
| count | 2 | no | 20,007 | 13,293 | 6,714 |
| count | 2 | yes | 19,831 | 13,293 | 6,538 |
| one | - | no | 36,867 | 21,793 | 15,074 |
| one | - | yes | 34,848 | 21,793 | 13,055 |
| one | 1 | no | 28,123 | 17,990 | 10,133 |
| one | 1 | yes | 27,878 | 17,990 | 9,888 |
| one | 2 | no | 18,690 | 11,977 | 6,713 |
| one | 2 | yes | 18,514 | 11,977 | 6,537 |
| none | - | no | 33,796 | 18,950 | 14,846 |
| none | - | yes | 31,794 | 18,950 | 12,844 |
| none | 1 | no | 25,751 | 15,679 | 10,072 |
| none | 1 | yes | 25,508 | 15,679 | 9,829 |
| none | 2 | no | 17,016 | 10,339 | 6,677 |
| none | 2 | yes | 16,840 | 10,339 | 6,501 |

The "count" rows with no sign minimum and no dispute rule are the rule of section 11; the "none" row matches the
practices review's count (4,639 labels with fewer than two records without the bulk papers). Bulk papers fall mostly
on signalling edges; disputed signs fall only on regulations, and most of them already fail a sign minimum of one.
