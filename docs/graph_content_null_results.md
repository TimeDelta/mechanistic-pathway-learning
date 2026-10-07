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
| degree-scaled popularity, for scale | 0.259 ± 0.024 | 0.211 | 0.500 |

At two swaps per edge the rewiring barely moved the score, and that reading supported the conclusion that the wiring
carried nothing beyond the degree sequence. Two swaps per edge is below the rewiring function's own default of ten
and far below what the configuration-model literature treats as mixed, so enough local structure survived that the
control had not destroyed what it was meant to destroy. At fifty swaps the real graph beats the rewired one by 0.043
per fold, and the rewired graph falls below chance on AUROC and below the degree baseline on AUPRC.

**The wiring carries real information.** The first conclusion here was an artefact of an under-mixed control, and
anything resting on it is withdrawn.

## What the corrected picture is

The graph matters, and the random walk already extracts what it has to give:

| | per fold |
|---|---|
| b6_mechanistic_gate_time_scales, the best trained configuration | 0.292 ± 0.027 |
| random walk with restart on the real graph | 0.293 ± 0.014 |

Sixteen trained configurations span 0.241 to 0.292 and none separates from the random walk by a paired interval
excluding zero. So the open question is not whether the graph holds signal, nor whether the labels are learnable: it
is why a signed, relation-typed, stoichiometric propagation with a noisy-OR head extracts no more of that signal than
unsigned diffusion does.

One candidate explanation ties this to a measurement in docs/membrane_potential_reach.md: at 99.7% of destinations
every incoming message carries the same sign (65 of 23,709 nodes see both signs across every relation feeding them).
A signed propagation on a graph with almost no sign mixing is close to an unsigned one, so the edge signs can
discriminate very little, and the extra machinery is acting where it has almost nothing to act on. That is a
hypothesis the numbers suggest rather than a demonstrated cause.

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

## The decisive follow-up

The same rewiring control against a mechanistic encoder rather than the random walk. The encoder uses edge signs,
relation types and stoichiometry that a random walk discards, so if those contribute anything the encoder should lose
more from rewiring than the random walk's 0.043. If it loses the same amount, the mechanistic layer is decorative.
That is a five-fold run and it bears on the experiment's central claim, so it waits for the user.
