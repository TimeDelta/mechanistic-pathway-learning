# The wiring carries signal; the mechanistic machinery adds nothing on top of diffusion (7 October 2026)

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

These paired differences stand, since they do not depend on the rewiring control. Each intervention against its own
twin, paired over five folds:

| intervention | difference | interval |
|---|---|---|
| manganese cofactor layer | -0.002 | [-0.005, +0.002] |
| signed logarithm (b3 cofactors) | +0.001 | within noise |
| signed logarithm (b6 time scales) | -0.001 | within noise |
| auxiliary laboratory loss (b3 cofactors) | -0.002 | within noise |
| protein descriptors | -0.030 | the one clear effect, and it is harmful |

So individual additions to the graph's content have not moved the score, while removing the wiring wholesale does.
Those are consistent: the signal lives in the gross connectivity the random walk reads, and the specific layers added
on top of it are small relative to the noise of a five-fold comparison.

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
