# Gene and protein split: slice results (8 October 2026)

Development readings. The split graph build, the user's decisions and the open questions are in
docs/gene_protein_split.md; nothing there is changed by this document, and no confirmatory configuration moved.

The `split_slice` job (runs/module_fix/split_slice.sh, 20 of 20 folds finished) trained four arms on the slice, five
disease-cluster folds each, seed 0:

| arm | graph | descriptors |
|---|---|---|
| `b6_mechanistic_gate_time_scales_descriptors_split_seed_masked` | `graph_split` | split table, protein block on the protein nodes |
| `b6_mechanistic_gate_time_scales_descriptors_resolved_seed_masked` | `graph` (merged) | resolved table, protein block on the gene nodes |
| `b6_linear_response_gate_time_scales_cofactors_descriptors_split` | `graph_split` | split table |
| `b6_linear_response_gate_time_scales_cofactors_descriptors_resolved` | `graph` (merged) | resolved table |

Each split arm is paired with the merged twin that carries the same descriptor values resolved per entry, so the
comparison is the graph change alone.

## Status: not a verdict (user, 8 October 2026)

These readings do not decide whether the split helps. The user's instruction: the slice result is not accepted as
evidence against the split until the final two confirmatory models are ready and tested on it. Until then the split
stays in, and this document is a development reading of the metabolic slice, not a decision about the full graph. The
three intervals below that exclude zero are reported as they came out and are not read as a null result.

## Macro and micro AUPRC

`experiments/compare_twin_runs.py`, A = split arm, B = merged twin, 451 perturbations scored by both. Readings: paired
per fold with a t interval on 4 degrees of freedom; pooled out-of-fold with a paired bootstrap over perturbations
(1,000 resamples); and the pooled difference after ranking each score inside one degree stratum of one test fold.

| encoder | A per fold | B per fold | per fold, t | pooled, bootstrap | within degree strata |
|---|---|---|---|---|---|
| message passing (seed masked) | 0.266 ± 0.035 | 0.275 ± 0.026 | −0.010 [−0.035, +0.015] | −0.008 [−0.028, +0.013] | −0.015 [−0.033, −0.001] |
| linear response (cofactors) | 0.264 ± 0.031 | 0.275 ± 0.042 | −0.011 [−0.034, +0.012] | +0.008 [−0.012, +0.026] | −0.008 [−0.032, +0.010] |

Micro AUPRC, the same three readings:

| encoder | per fold, t | pooled, bootstrap | within degree strata |
|---|---|---|---|
| message passing (seed masked) | −0.016 [−0.038, +0.005] | −0.013 [−0.030, +0.002] | −0.015 [−0.029, −0.004] |
| linear response (cofactors) | −0.012 [−0.023, −0.002] | −0.005 [−0.016, +0.007] | −0.007 [−0.025, +0.010] |

The degree strata above come from the merged slice (`data/processed/graph`), the registered strata. Recomputing them
from the split graph gives:

| encoder | macro within split-graph strata | micro within split-graph strata |
|---|---|---|
| message passing (seed masked) | −0.012 [−0.028, +0.003] | −0.012 [−0.024, −0.001] |
| linear response (cofactors) | −0.007 [−0.025, +0.005] | −0.010 [−0.024, +0.001] |

See "Degrees collapse" below for why that second reading carries no information on the split graph.

Reading: the split does not raise either encoder's symptom ranking on the slice. Every pooled interval covers zero,
in both metrics and for both encoders. Three of the twelve intervals exclude zero, and all three point against the
split: the message-passing within-strata difference in macro AUPRC (−0.015, upper bound −0.001) and in micro AUPRC
(−0.015, upper bound −0.004), and the linear-response per-fold micro difference (−0.012, upper bound −0.002). At 95
percent about 0.6 of twelve intervals would exclude zero by chance, so three is more than chance alone would give,
and the direction is consistent. The split costs about 0.01 AUPRC or nothing; it does not pay for itself in accuracy
on the slice.

## Degrees collapse on the split graph

Every one of the 451 slice perturbations seeds a gene node, and on `graph_split` every one of those gene nodes has
degree exactly 1: its single `encodes` edge to its protein. (Gene nodes in general reach degree 21 there, through
several reviewed entries and the regulatory edges that stay on genes; the perturbed ones do not.) Protein nodes carry
the rest: median degree 3, maximum 258, against median 2 and maximum 145 for the merged slice's gene nodes.

`ExperimentData.perturbation_degrees` sums the degree over a perturbation's seed nodes, so on the split graph it is 1
for every perturbation. Three parts of the design read it:

- the H1 reading "macro and micro AUPRC within degree strata" (docs/preregistration.md), which exists to show the
  model is not ordering perturbations by degree. With one stratum it equals the pooled reading and shows nothing;
- the degree-stratified label permutation control (`permute_symptom_labels_within_degree_strata`), which becomes a
  free permutation;
- `macro_auprc_by_degree_bin`, whose equal-count bins become arbitrary.

The fix, which the user proposed ("can't you just add a 1-hop to the degree for the strata and normalization purposes
only?"): `ExperimentData.perturbation_degrees_for_strata` returns the seed degree on a merged graph and, on a graph
with protein nodes, `perturbation_degrees_with_encoded_proteins` — the degree of the seeds and of the proteins they
encode, less the `encodes` edges themselves, which count in both degrees and did not exist in the merged graph. It is
read for strata, normalisation and the degree-scaled baselines only; no model reads it as a feature.

Which hop matters. Reading *every* node one edge away was tried first and is wrong on the full graph: a gene node
keeps its regulatory edges there, so the sum picks up its neighbours' degrees as well.

| degree on the split graph | median | max | Spearman against the merged seed degree | same tercile | equal value |
|---|---|---|---|---|---|
| full graph, seed only | 2 | 493 | 0.648 | 0.545 | 0.086 |
| full graph, every node one hop away | 424 | 37,585 | 0.733 | 0.568 | 0.000 |
| full graph, the encodes hop | 16 | 1,953 | 0.998 | 0.998 | 0.992 |
| slice, seed only | 1 | 1 | undefined (constant) | 0.379 | 0.386 |
| slice, every node one hop away | 4 | 147 | 1.000 | 1.000 | 0.000 |
| slice, the encodes hop | 2 | 145 | 1.000 | 1.000 | 1.000 |

The merged graph's seed degree has median 16 and maximum 1,953 on the full graph (1,539 perturbations,
`evidence_full_v2` with the better_v2 selection) and median 2, maximum 145 on the slice. So the encodes hop does not
approximate the registered stratification, it reproduces it: the same degree value for 1,527 of 1,539 full-graph
perturbations and for all 451 slice perturbations. Reading every neighbour agrees on the slice only because there
every perturbed gene node had exactly one neighbour; on the full graph it moves 43 percent of the strata assignments
and inflates the degree by a factor of 26 at the median.

Twelve full-graph perturbations keep a different value. Eleven differ by 1 to 3 edges, from edges the merged build
deduplicated and the split keeps apart. The twelfth is the H3-3A knockout: 3 in the merged graph against 118 under the
encodes hop, because `PROTEIN:H3-3A` is the reviewed entry of both H3-3A and H3-3B (one protein, two genes), so the
protein node carries the product-level edges of both while the merged graph gave each gene its own. That is the split
working as specified, not an artifact of the degree rule.

What that costs the registered reading: of 1,539 perturbations, six change degree quintile, the stratification the H1
within-strata reading uses (CASP2, FGF8, GBA1, H3-3A, OTX2, WT1), and three change tercile, the bins of
`macro_auprc_by_degree_bin` (CNTNAP2, H3-3A, PARK7). CNTNAP2 and PARK7 keep their degree value exactly; they move
because the bins hold equal counts and a tie at a boundary falls either side depending on the order of the whole
vector. So the stratification is preserved to within six of 1,539 assignments, and no re-registration of the strata is
needed for the move to the split graph.

## Module health

`experiments/module_health.py`, five-fold means, split arm minus merged twin (rule in docs/module_health.md):

| reading | message passing: split, merged, difference | linear response: split, merged, difference |
|---|---|---|
| responding modules | 7.2, 6.0, +1.2 (3 of 5 folds higher) | 0.0, 0.0, 0.0 |
| floor share | 0.545, 0.635, −0.090 | 0.998, 0.998, −0.000 |
| median link correlation | 0.961, 0.523, +0.437 (5 of 5) | 0.258, 0.629, −0.371 |
| nonempty modules | 7.4, 8.0, −0.6 | 8.0, 7.8, +0.2 |
| median support size | 2.1, 19.9, −17.8 (0 of 5) | 30.7, 18.2, +12.5 (4 of 5) |

Readings:

- **Message passing.** One more module responds, but the hard support per module falls from about 20 nodes to 2, and
  the responding modules' link vectors become nearly identical (median correlation 0.96 against 0.52). Many modules
  carrying the same link pattern over a handful of nodes each is the degenerate case docs/module_health.md names, not
  a repaired one. The lower floor share (0.545 against 0.635) is the one reading that improves: the modules add less
  uniformly across perturbations.
- **Linear response.** No module responds in either arm (standard deviation of activation below 0.05 in every fold),
  as in every earlier linear-response run, and the floor share stays at 0.998: the modules add the same amount to
  every perturbation. The split does not change this. The larger support (30.7 against 18.2 nodes) is support of
  modules that do not respond, so it is not a gain.

The split therefore does not fix the unused modules of the linear-response encoder, and it makes the message-passing
modules narrower and more redundant.

## What this does not settle

- The slice is the metabolic slice. The split's point is that a gene and its protein are different nodes where data
  tell them apart, which should matter most where protein-level descriptors and protein complexes carry signal. The
  full neuronal graph has 1,681 `protein_entity` complexes and 12,709 protein nodes; the slice has far less of that
  structure, so a null here is weak evidence about the full graph.
- Only one treatment per encoder ran (seed masking for message passing, none for linear response), as the job was
  registered. The `zero_init_slow` arm on the split still waits for the user.
- No confirmatory configuration was moved by this document. The user has since asked that the confirmatory
  configuration be the split graph with the protein descriptors, brain expression and compartments; that change needs
  a dated amendment in the preregistration before `configs/confirmatory_v2` is touched.
- The two confirmatory models themselves have not been run on the split. Whatever the slice shows, the comparison the
  user accepts is those two models, so this document cannot close the question either way.

## Reproducing

```
PYTHONPATH=. OMP_NUM_THREADS=1 python experiments/compare_twin_runs.py --run-root runs/module_fix \
  --cache-dir runs/twin_comparisons_split_merged \
  --pairs b6_mechanistic_gate_time_scales_descriptors_split_seed_masked:b6_mechanistic_gate_time_scales_descriptors_resolved_seed_masked \
          b6_linear_response_gate_time_scales_cofactors_descriptors_split:b6_linear_response_gate_time_scales_cofactors_descriptors_resolved
PYTHONPATH=. OMP_NUM_THREADS=1 python experiments/module_health.py --run-root runs/module_fix --pairs <the same two pairs>
```

`--pairs` takes a list, so both pairs go after one flag; a second `--pairs` replaces the first.
