# Sign and relation-typing ablations of the linear-response encoder (slice pilot, 7 October 2026)

The monogenic slice is a pilot (amendment of 7 October in docs/preregistration.md): one layer, 451 genes, no brain
expression. These results say what the signs and relation types do on that graph; they decide nothing about the
full graph, where the signaling and regulatory layers carry most of the signed edges.

Each arm is one argument away from b3_linear_response_cofactors (experiments/run_main_model_batch.py), scored with
experiments/compare_twin_runs.py: per fold with a t interval (4 degrees of freedom), pooled with a paired bootstrap
over perturbations (1000 resamples), and within degree strata. Intervals are 95 percent and uncorrected.

| arm (A) | change | A per fold | twin per fold | A minus twin, per fold | pooled | within degree strata |
|---|---|---|---|---|---|---|
| b3_linear_response_cofactors_unsigned | every edge sign +1 (`--edge-signs all_positive`) | 0.266 ± 0.022 | 0.266 ± 0.025 | +0.000 [-0.005, +0.006] | +0.003 [-0.005, +0.010] | +0.002 [-0.011, +0.014] |

Per-fold differences for the unsigned arm: +0.004, -0.004, -0.004, +0.005 and +0.001. Removing every sign leaves the
slice score where it was, so on the metabolic slice the encoder does not use the signs to predict symptoms. Every
stored slice edge is +1 (catalyzed_by, product_of and substrate_of); the only negative edges are the ones the encoder
derives, depletes_substrate (one per main substrate edge; carrier cofactors get none), so on the slice a sign carries
stoichiometry and nothing else. On graph_full_neuronal, OmniPath and CollecTRI add measured inhibition and repression
(18,149 inhibits edges and 7,087 repressive regulon edges), so the sign ablation has to be repeated there.
The shared-gain arm (one gain for every relation) and the sign permutation (signs shuffled among edges) are still to
come.
