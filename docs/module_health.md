# Module health: the rule for the module fix

Written on 8 October 2026, before any run it reads existed. The user's decisions of that day: drop the sigmoid head
from the confirmatory family, fix the modules before any second-family run, and run nothing of the second family
until then.

## Why

On the first family's full-graph pilots the noisy-OR leaks stayed at the raw base rate (0.061 against a loss-optimal
constant of 0.204), and the modules supplied the difference (docs/best_epoch_zero.md). The slice shows the same at the
confirmatory leak learning rate of 0.0002: leaks of 0.193 to 0.194 against a raw rate of 0.191 and an optimum of 0.434.
In b6_linear_response_gate_time_scales_cofactors no module varies across perturbations. Modules that add the same
amount to every perturbation cannot be read as pathways.

## Readings (experiments/module_health.py, per fold, then five-fold means)

- **responding_modules:** modules whose activation over the test perturbations has a standard deviation above 0.05.
- **floor_share:** the modules' part of P is m = 1 - (1 - P) / (1 - leak). For each symptom, take the 10th percentile
  of m over the test perturbations divided by its mean, then the median over the symptoms whose mean m is above 0.01.
  Near 1, the modules add about the same to every perturbation; near 0, they add to some perturbations only.
- **median_link_correlation:** the median correlation between the link vectors of the pairs of modules with a nonempty
  hard support.
- **nonempty_modules** and **median_support_size:** both on the hard support (gate above 0.5).
- **macro AUPRC and micro AUPRC.**

Readings of runs that predate the fix (smoke test of the script, 5 folds each):

| run | responding_modules | floor_share | median_link_correlation |
|---|---|---|---|
| b6_linear_response_gate_time_scales_cofactors | 0.0 | 0.997 | 0.382 |
| b6_linear_response_time_scales | 8.0 | 0.707 | 0.908 |
| b6_mechanistic_gate_time_scales | 4.4 | 0.709 | 0.467 |

## Runs

runs/module_fix/module_fix_slice.sh has two arms, five disease-cluster folds each, seed 0, run fold by fold:
- b6_mechanistic_gate_time_scales and b6_linear_response_gate_time_scales_cofactors, each rerun under the current code;
- each of them with --start-at-weighted-optimum (the `_weighted_start` configurations, which differ from their twins
  by that flag only).

## Rule (per encoder, fix minus twin, five-fold means)

The fix repairs the modules if all three hold:
1. floor_share falls by at least 0.2, or responding_modules rises by at least 2;
2. median_link_correlation rises by no more than 0.05;
3. macro and micro AUPRC each fall by no more than 0.01.

If 1 holds but 2 or 3 does not, the reading is partial.

The targets for a module set read as pathways are floor_share at most 0.5, at least 4 of 8 responding modules and a
median link correlation of at most 0.7. If the leak start alone falls short of them, the next candidates are tried on
the slice in turn: a faster leak learning rate, mean pooling over the module support, and fewer modules.

The reading goes to the user, who decides whether the confirmatory_v2 configurations take the fix. A pass on the slice
is a development reading and confirms nothing.
