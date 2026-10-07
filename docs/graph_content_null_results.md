# Every graph-content change lands within noise, and a rewired graph scores the same (7 October 2026)

This note collects the slice results that bear on one question: does the content of the graph, as opposed to its
degree structure, contribute anything measurable to symptom prediction? As of today the answer appears to be no, with
one confirmatory measurement outstanding.

## The table

Per-fold macro AUPRC, five folds, disease-cluster grouped split, 451 perturbations.

| configuration | per fold |
|---|---|
| b6_mechanistic_gate_time_scales | 0.292 ± 0.027 |
| b6_linear_response_gate_time_scales_cofactors | 0.288 ± 0.044 |
| b3_typed_nodes_laboratory | 0.285 ± 0.035 |
| **random walk, DEGREE-PRESERVING REWIRED graph** | **0.272 ± 0.014** |
| random walk, real graph | 0.269 ± 0.028 |
| b3_linear_response_cofactors_log | 0.267 ± 0.029 |
| b3_linear_response_cofactors | 0.266 ± 0.025 |
| b3_linear_response_cofactors_laboratory | 0.264 ± 0.024 |
| b3_linear_response_cofactors_laboratory_manganese | 0.263 ± 0.020 |
| b3_typed_nodes_descriptors | 0.263 ± 0.029 |
| b3_linear_response | 0.262 ± 0.019 |
| b3_descriptors_only | 0.260 ± 0.047 |
| b6_linear_response_time_scales_cofactors | 0.260 ± 0.017 |
| b6_linear_response_time_scales | 0.259 ± 0.023 |
| degree-scaled popularity | 0.259 |
| b6_linear_response_time_scales_cofactors_log | 0.259 ± 0.037 |
| b6_linear_response | 0.258 ± 0.024 |
| b3_degree_only | 0.246 ± 0.029 |
| b3_local_structural | 0.241 ± 0.030 |

A random walk on a graph whose wiring has been destroyed, keeping only each node's in- and out-degree per relation
type, scores 0.272 ± 0.014. That is above ten of the sixteen trained configurations and above the same random walk on
the real graph.

## The paired differences, which is the honest form of the comparison

Each intervention against its own twin, paired over the five folds with a 2000-draw bootstrap:

| intervention | difference | interval |
|---|---|---|
| manganese cofactor layer | -0.002 | [-0.005, +0.002] |
| signed logarithm (b3 cofactors) | +0.001 | within noise |
| signed logarithm (b6 time scales) | -0.001 | within noise |
| auxiliary laboratory loss (b3 cofactors) | -0.002 | within noise |
| protein descriptors | -0.030 | the one clear effect, and it is harmful |

Every interval includes zero except the descriptors, which hurt.

## Why this explains the rest of the project's findings

- The four propagation variants sit at 0.288 to 0.293, with paired intervals against the random walk that include
  zero. If they were reading mechanism rather than connectivity, they should separate from it.
- The signed logarithm exists to make distant changes readable despite the field decaying to about 1e-10 four edges
  out. It moves nothing. If the distant field carried task-relevant information, making it readable should have
  helped, so the decay is real but is not what limits accuracy.
- The gated noisy-OR modules are hollow, with activation standard deviation 0.000 to 0.004 against a mean of 0.071 to
  0.076. There is no perturbation-specific structure for a module to lock onto.
- The protein descriptors hurt rather than help, which is what a second route to the same degree-and-annotation
  shortcut would do.

Three normalisation schemes for the propagation were tried and all failed (docs/membrane_potential_reach.md). That
work was aimed at transmitting a signal more faithfully; this table says the signal is not carrying anything.

## The outstanding measurement

The committed rewiring control used two degree-preserving double-edge swaps per edge. The rewiring function's own
default is ten and the configuration-model literature generally wants considerably more before treating a graph as
mixed, so at two swaps substantial local structure survives and the control may simply be too weak to have destroyed
what it was meant to destroy. The registered job `rewiring_check` repeats it at fifty swaps per edge, writing to
runs/rewiring_check/ so the committed baselines are untouched. Until that lands the conclusion here is provisional.

If the null survives at fifty swaps, the decisive follow-up is the same control against a mechanistic encoder rather
than the random walk, because the encoder uses edge signs, relation types and stoichiometry that a random walk
discards, and that is the one remaining way the wiring could be earning its place. That is a five-fold run and a
question about the experiment's central claim, so it is the user's call rather than a maintenance decision.
