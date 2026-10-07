# On the slice the wiring carries signal, and the mechanistic encoders have not been shown to add to diffusion (7 October 2026)

Scope: one metabolic slice graph, 451 monogenic perturbations and ten symptoms in five disease-cluster folds. It is
not a general rule about signed or typed propagation against unsigned walks; published work on other graphs and
labels reports typed models ahead of diffusion (docs/literature_appraisal_staging.md, state of Q1 at the cap; unreviewed).

This note was first written with the opposite conclusion, from a rewiring control that turned out to be too weak. The
correction is the point of the note, so both versions are stated.

## The control, at two strengths

A degree-preserving rewiring keeps every node's in- and out-degree per relation type and destroys which node joins
which. Grouped disease-cluster split, 451 perturbations, random walk with restart, one run and one seed:

| graph | macro AUPRC per fold | pooled | macro AUROC |
|---|---|---|---|
| real | **0.293 ± 0.014** | 0.245 | 0.530 |
| rewired, 2 swaps per edge | 0.272 ± 0.014 | 0.227 | 0.523 |
| rewired, 50 swaps per edge | **0.250 ± 0.033** | 0.211 | **0.485** |
| labels permuted within degree strata, a different control | 0.269 ± 0.028 | 0.224 | 0.529 |
| degree-scaled popularity, for scale | 0.259 ± 0.024 | 0.211 | 0.500 |

A second error, found after the first correction and worth recording because it caused it. The figure 0.269 in the
fourth row was read as the real graph's random-walk score for most of the day, from a grep whose context window
caught the tail of the label-permutation table in docs/phase2_baselines_disease_cluster.md rather than the
grouped-split table above it. The real graph has scored 0.293 in that document since 3 October and in both runs made
today, so the committed baseline was never stale; the comparison was simply against the wrong row, which is how a
rewired graph came to look as good as the real one. Read a figure's own section heading before comparing to it.

At two swaps per edge the rewiring barely moved the score, and that reading supported the conclusion that the wiring
carried nothing beyond the degree sequence. Two swaps per edge is below the rewiring function's own default of ten
and far below what the configuration-model literature treats as mixed, so enough local structure survived that the
control had not destroyed what it was meant to destroy. At fifty swaps the real graph beats the rewired one by 0.043
per fold, and the rewired graph falls below chance on AUROC and below the degree baseline on AUPRC.

**The wiring carries real information.** The first conclusion here was an artefact of an under-mixed control, and
anything resting on it is withdrawn.

## What the corrected picture is

The graph matters, and the random walk already extracts what it has to give:

| | per fold | pooled, paired against the random walk |
|---|---|---|
| b3_typed_nodes, highest per fold | 0.293 ± 0.031 | -0.031 [-0.055, -0.011] |
| b6_mechanistic_gate_time_scales, closest pooled | 0.292 ± 0.024 | -0.020 [-0.043, +0.006] |
| random walk with restart on the real graph | 0.293 ± 0.014 | |

Eighteen trained configurations in docs/phase3_main_model.md span 0.241 to 0.293 per fold. None beats the random walk,
and most fall below it: on the pooled paired difference, the pre-registered endpoint, 14 of the 18 have an interval
entirely below zero, and the four that reach zero are b6_default, b6_mechanistic_gate_time_scales,
b6_mechanistic_expected_gate and b6_linear_response_gate_time_scales_cofactors (upper bound -0.000). "Best" depends on
the metric: b3_typed_nodes ties the random walk per fold and is reliably below it pooled. This note said earlier that
sixteen configurations span 0.241 to 0.292 with none separating from the random walk; that undercounted the table and
missed that the separations run against the trained models. So the open question is not whether the graph holds
signal, nor whether the labels are learnable: it is why a signed, relation-typed, stoichiometric propagation with a
noisy-OR head extracts less of that signal than unsigned diffusion does, or at best as much.

One candidate explanation stood here and has been withdrawn, because the measurement under it was wrong. It read
that at 99.7% of destinations every incoming message carries the same sign, so a signed propagation on this graph is
close to an unsigned one. That count omitted the encoder's derived depletes_substrate relation, which supplies nearly
every negative sign in the model. On the relation stack the encoder builds, 7,811 of 21,337 destinations (36.6%) on
this graph and 8,712 of 24,121 (36.1%) on graph_neuronal receive both signs.

The corrected measurement points the other way. At equal gains the aggregate row sum is zero at 7,811 of the 21,337
nodes with any input, and every one is a metabolite, 7,811 of 8,460 (92.3%; recorded first as 7,437 and 88%, an
undercount from a tolerance below float32 rounding), while no reaction has one, because a
metabolite that is both produced and consumed receives product_of with mean +1 and depletes_substrate with mean -1.
Those metabolites are blind to whatever their production and consumption inputs share and pass on only the
difference. The field there is not zero, since the two inputs carry different upstream states; what is lost is the
shared part, and under attenuation the difference of two small similar numbers is where precision goes.
Learned gains do not undo it. relation_gain() keeps signed relations positive and shares every gain across all
nodes, so one global ratio of the product_of gain to the depletes_substrate gain fixes the net direction at every
both-produced-and-consumed metabolite at once. The candidate cause of the null is therefore not that the signs carry
too little structure, but that the encoder's node-wise averaging destroys the structure they carry at a third of its
destinations. docs/membrane_potential_reach.md records the measurement and the two fixes it implies. This is still a
hypothesis about the null rather than a demonstrated cause: it predicts that the spectral normalisation arm, which
keeps no per-node divisor, and an aggregator that survives sign mixing should both move the score, and neither has
been run against symptoms yet.

## The wiring does carry signal: 20 degree-preserving rewirings (7 October 2026)

The single rewiring reported earlier left the question open, since one draw cannot say whether a gap is
real. Twenty degree-preserving rewirings at 50 double-edge swaps per edge, with the graph, the
evidence, the fold assignment and the restart probability all held fixed and only the rewiring seed
varied, give the null distribution in docs/rewiring_null_distribution.md:

| | per-fold macro AUPRC |
|---|---|
| real graph | **0.293** |
| rewirings, mean of 20 | 0.268 ± 0.010 |
| rewirings, best of 20 | 0.287 |
| rewirings, worst of 20 | 0.250 |

No rewiring reached the real graph, so p = 0.048, and the real score sits 2.4 null standard deviations
above the null mean. The p-value is held at 0.048 by the number of rewirings, not by the size of the
effect: with 20 draws and none at least as good, (0 + 1) / (20 + 1) is the smallest value the test can
return. More seeds would lower it.

Two things keep this from being a strong result. The margin over the best single rewiring is 0.006,
which is smaller than the real graph's own fold-to-fold spread of ±0.014, so a reader comparing the
real graph with one rewiring rather than with the distribution would see nothing. And a rewiring
preserves in-degree and out-degree per relation but not the bipartite metabolite-reaction structure,
so part of the gap may be the plausibility of the graph as a graph rather than the correctness of its
biology. What the control does establish is the claim it was run for: the specific wiring matters, and
the random walk is reading the graph rather than the degree sequence.

## Which earlier interventions were and were not null

Rescored 7 October 2026 by experiments/compare_twin_runs.py (docs/twin_comparisons.md), which takes every two slice
configurations differing in one argument with five folds each, 20 pairs, and reads each three ways: paired per fold
with a t interval on four degrees of freedom, pooled out-of-fold with a paired bootstrap over the 451 perturbations,
and pooled within degree strata. The table this replaces gave single numbers and, where it gave an interval, a
percentile bootstrap over the five fold differences, which runs about a third narrower than the t interval (the
manganese row read [-0.005, +0.002]). Macro AUPRC, A minus its twin; 95 percent intervals, not corrected:

| intervention (twin) | per fold, t | pooled | within degree strata |
|---|---|---|---|
| manganese cofactor layer (b3 cofactors laboratory) | -0.002 [-0.008, +0.005] | +0.004 [-0.002, +0.014] | -0.004 [-0.015, +0.009] |
| neuronal and oxidative variant (same) | +0.012 [-0.028, +0.053] | +0.004 [-0.009, +0.015] | +0.014 [-0.007, +0.037] |
| that variant's curated layers (variant without them) | +0.007 [-0.029, +0.043] | +0.005 [-0.005, +0.013] | +0.015 [-0.004, +0.031] |
| signed logarithm (b3 cofactors) | +0.002 [-0.012, +0.015] | +0.009 [-0.006, +0.029] | -0.004 [-0.023, +0.015] |
| signed logarithm (b6 time scales) | -0.001 [-0.027, +0.024] | +0.000 [-0.017, +0.014] | -0.005 [-0.037, +0.019] |
| auxiliary laboratory loss (b3 cofactors) | -0.002 [-0.008, +0.005] | +0.000 [-0.001, +0.001] | +0.000 [-0.002, +0.002] |
| auxiliary laboratory loss (b3 typed nodes) | -0.007 [-0.036, +0.022] | -0.000 [-0.013, +0.012] | +0.002 [-0.015, +0.020] |
| cofactor relations (b3 linear response) | +0.004 [-0.009, +0.017] | -0.003 [-0.011, +0.004] | +0.006 [-0.002, +0.016] |
| protein descriptors (b3 typed nodes) | -0.030 [-0.056, -0.004] | -0.001 [-0.017, +0.018] | -0.018 [-0.039, -0.002] |

Read across the 20 pairs, three per-fold intervals exclude zero (protein descriptors, p = 0.034, and the
local-descriptor encoder against the plain message-passing encoder and against the linear response, p = 0.009 and
0.015), one within-strata interval does (protein descriptors) and no pooled interval does. Twenty comparisons at the
0.05 level should produce about one such interval by chance, and after Holm's correction across the 20 per-fold tests
none survives (smallest adjusted p 0.18). So the protein descriptors are the only addition with any sign of an effect,
in the harmful direction on two of the three readings, and calling it established was too strong. Individual
additions to the graph's content have not moved the score, while removing the wiring wholesale does; the signal lives
in the gross connectivity the random walk reads, and the layers added on top of it are small relative to the noise of
a five-fold comparison on 451 perturbations.

## The baseline is undirected and the encoders are not (7 October 2026)

The comparison behind this document sets directed models against an undirected walk. Both trained encoders follow the
stored edge direction: the linear-response encoder puts each edge at row relation x N + target and column source, and
adds one reverse edge, the depletion edge from a reaction back to its substrate
(mechanistic_pathway_learning/models/linear_response_encoder.py, signed_stacked_adjacency); the message-passing
encoder takes edge_index as [edge_source, edge_target] and adds none. The random walk with restart symmetrises the
graph ("The walk runs on the undirected graph with currency metabolites removed",
mechanistic_pathway_learning/models/baselines/random_walk_with_restart_baseline.py), and so does every rewiring in
the null distribution above. The null result therefore reads: a directed, signed model has not beaten an undirected,
unsigned walk.

Silverbush and Sharan (2019, doi:10.1038/s41467-019-10887-6) report that diffusion over an oriented network ranks drug targets better than diffusion over
an unoriented one (docs/literature_appraisal_staging.md, slot 14). On this graph a walk that follows direction from a
gene cannot return to another gene, since genes have outgoing edges only, so the guess this section was written to
test is that direction costs the encoders here. experiments/score_directed_random_walk.py (docs/directed_random_walk.md)
runs three unsigned walks with one readout fixed in advance (cosine similarity of the held-out gene's walk distribution
over reactions and metabolites with the mean distribution of the symptom's training positives), differing only in the
edges they follow. Macro AUPRC, minus the undirected walk with the same readout; 95 percent intervals, not corrected:

| walk | per fold, t | pooled | within degree strata |
|---|---|---|---|
| the encoder's own edge directions | +0.005 [-0.003, +0.013] | +0.006 [+0.001, +0.011] | -0.002 [-0.011, +0.009] |
| downstream only (substrate to reaction to product) | +0.006 [-0.015, +0.027] | +0.008 [+0.003, +0.014] | +0.012 [-0.002, +0.028] |

Following direction does not cost the walk anything, so the guess is not supported. Nor is a gain established: one
reading of three excludes zero for each arm, the pooled one, and for the encoder's edges the within-strata reading
puts the point estimate below zero, so that arm's pooled gain is consistent with ordering genes by degree. The
committed walk (B1, which scores the mass landing on the training genes rather than a profile similarity) is within
noise of all three (+0.004 [-0.007, +0.015] per fold against the undirected profile walk). On this slice direction
alone neither explains the encoders' deficit nor adds to the walk's score. Signs are the other half of the
difference: b3_linear_response_cofactors_unsigned (every edge +1, every gain positive, otherwise identical to
b3_linear_response_cofactors) is running, and it removes the cancellation at mixed-sign nodes together with the sign
information.

## A trained readout on the walk's own distributions does not beat the walk either (7 October 2026)

Every trained configuration fits a readout to the labels and the random walk does not, so the null compares a change
of encoder and the addition of supervision at once. Huang et al. (2024, doi:10.1038/s41591-024-03233-x) has the same
confound in the opposite direction, a trained typed model against untrained diffusion statistics
(docs/literature_appraisal_staging.md, slot 17). experiments/score_trained_diffusion_readout.py
(docs/trained_diffusion_readout.md) keeps the diffusion and adds only supervision: each gene's walk distribution over
reactions and metabolites, square-rooted, reduced to 32 components fitted on the training genes, read by one
L2-regularised logistic regression per symptom with its strength chosen by inner cross-validation. Macro AUPRC minus
the untrained walk, 95 percent intervals, not corrected:

| trained readout on | per fold, t | pooled | within degree strata |
|---|---|---|---|
| the undirected walk | -0.012 [-0.031, +0.006] | -0.026 [-0.044, -0.011] | -0.004 [-0.026, +0.015] |
| the downstream-only walk | +0.001 [-0.036, +0.037] | -0.020 [-0.037, -0.002] | +0.005 [-0.018, +0.023] |

Supervision on top of the walk's own information ties the walk per fold and within degree strata and trails it
pooled, the pattern of the trained encoders in research summary row 5. So on 451 genes in five folds a trained
readout has not been shown to beat the untrained walk even when it is given everything the walk computes, and the
encoders' failure to beat it is weak evidence against their mechanism specifically. That reading has two limits: it
rests on one readout design, fixed before scoring, and a richer or better-regularised readout might do better; and
the pooled deficit of trained readouts may partly reflect score scales that follow each fold's training base rate,
which pooling penalises and the other two readings remove.

## What is still not established

This section was written before the 20-rewiring null existed and said there was no distribution to place the real
score against; the section above supersedes that. What the null leaves open:

- The p-value of 0.048 is the floor that 20 draws allow, not a measure of the effect's size, and the margin over the
  best single rewiring (0.006) is below the real graph's fold-to-fold spread (±0.014).
- Rewiring keeps in-degree and out-degree per relation but not the bipartite metabolite-reaction structure, so part of
  the gap may reflect the graph's plausibility as a graph rather than the correctness of its biology.
- The full-data rewiring control in docs/phase2_baselines_full_disease_cluster.md ran at 2 swaps per edge, the strength
  shown here to be under-mixed, so its reading that the full-data signal is degree and not topology is withdrawn until
  it is rerun. experiments/run_baselines.py now defaults to 50 swaps per edge.
- Every number here is the random walk. Whether the mechanistic encoder loses more than the random walk under rewiring
  is the follow-up below.

## The decisive follow-up

The same rewiring control against a mechanistic encoder rather than the random walk. The encoder uses edge signs,
relation types and stoichiometry that a random walk discards, so if those contribute anything the encoder should lose
more from rewiring than the random walk's 0.043. If it loses the same amount, the mechanistic layer is decorative.
That is a five-fold run and it bears on the experiment's central claim, so it waits for the user.
