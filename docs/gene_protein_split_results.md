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
| message passing (seed masked) | −0.016 [−0.038, +0.005] | −0.013 [−0.030, +0.002] | −0.015 [−0.034, +0.001] |
| linear response (cofactors) | −0.012 [−0.023, −0.002] | −0.005 [−0.016, +0.007] | −0.007 [−0.025, +0.005] |

The degree strata above come from the merged slice (`data/processed/graph`), the registered strata. Recomputing them
from the split graph moves the within-strata readings to −0.012 [−0.028, +0.003] and −0.007 [−0.025, +0.005]; see
"Degrees collapse" below for why that reading carries no information on the split graph.

Reading: the split does not raise either encoder's symptom ranking on the slice. Every macro interval covers zero.
The two intervals that exclude zero both point against the split and both are small: the message-passing
within-strata macro difference (−0.015, upper bound −0.001) and the linear-response per-fold micro difference (−0.012,
upper bound −0.002). Each arm is one of four comparisons here, so at 95 percent about 0.4 of these intervals would
exclude zero by chance. The split costs about 0.01 AUPRC or nothing; it does not pay for itself in accuracy on the
slice.

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

This needs a decision before any confirmatory run moves to a split graph. The options are to read a perturbation's
degree from the protein nodes its seed gene encodes, to keep the merged graph's degree as the stratifying variable, or
to drop the degree readings on split graphs and say so in the amendment. The first keeps the reading's meaning: what
the perturbation reaches is the protein's neighbourhood.

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
- No confirmatory configuration was moved. `configs/confirmatory_v2` still names `graph_full_neuronal`.

## Reproducing

```
PYTHONPATH=. OMP_NUM_THREADS=1 python experiments/compare_twin_runs.py --run-root runs/module_fix \
  --cache-dir runs/twin_comparisons_split_merged \
  --pairs b6_mechanistic_gate_time_scales_descriptors_split_seed_masked:b6_mechanistic_gate_time_scales_descriptors_resolved_seed_masked \
          b6_linear_response_gate_time_scales_cofactors_descriptors_split:b6_linear_response_gate_time_scales_cofactors_descriptors_resolved
PYTHONPATH=. OMP_NUM_THREADS=1 python experiments/module_health.py --run-root runs/module_fix --pairs <the same two pairs>
```

`--pairs` takes a list, so both pairs go after one flag; a second `--pairs` replaces the first.
