# Why early stopping picks epoch 0 on the full graph

Investigation of 8 October 2026. Material: the first family's four development pilots (fold 0 of lockbox_v1,
seed 0; runs/full/confirmatory_*_disease_cluster_and_targets_development/fold0_seed0), their stored best and final
states, and three slice runs. Each statement is marked as measured or inferred. The analysis scripts stayed outside the
repository. No lockbox output and no score of the second family was read.

## Findings

1. **Calibration sets the stopping epoch (measured).** On the pilots' validation set, a constant per-symptom
   prediction fitted on the training pool has a validation loss of 0.488. The best validation losses of the four
   pilots are:

   | pilot | best validation loss | best epoch |
   |---|---|---|
   | message passing, sigmoid | 0.534 | 0 |
   | linear response, sigmoid | 0.494 | 4 |
   | message passing, noisy-OR | 0.519 | 8 |
   | linear response, noisy-OR | 0.497 | 8 |

   For comparison, one global constant scores 0.647 and P = 0.5 scores 0.693. For the message-passing sigmoid pilot,
   the stored states were centred on the constant and their per-perturbation deviations scaled by a. The validation
   loss is 0.488 at a = 0, 0.505 at a = 0.5 and 0.548 at a = 1 (the epoch-0 state). The training fit does not carry
   over to validation: on training perturbations macro AUROC rises from 0.59 to 0.85, on validation from 0.515 to
   0.573. Validation AUPRC is flat between 0.19 and 0.23 in all four pilots.

2. **The sigmoid head starts far from the optimum on the full graph (measured).**
   - The loss weighs negatives at 0.2, so the loss-optimal constant is about three times the raw base rate. Over the
     pilots' training pool it is 0.204 against a raw base rate of 0.064 (means over 22 symptoms).
   - The head starts at P = 0.5 and reaches the constant within the first full-graph epoch (53 steps).
   - On the slice, positives are 19 percent of training pairs and the optimum is near 0.5. The same head starts near
     its optimum there, and no finished slice run has best epoch 0.
   - Counted in steps, the message-passing sigmoid head reaches its minimum after 40 to 80 steps on the slice and on
     the full graph alike.

3. **The noisy-OR leak stays at the raw base rate, and the modules make up the gap (measured).**
   - `--init-leak-from-base-rate` starts each leak at the raw base rate. At `--leak-learning-rate 0.0002` the
     leak's logit moved by a median of +0.01 by the best epoch.
   - Means over 22 symptoms:

     | pilot | raw base rate | loss-optimal constant | leak at the best epoch | mean test prediction | part of P from the modules, 1 - (1 - P) / (1 - leak) |
     |---|---|---|---|---|---|
     | linear response, noisy-OR | 0.064 | 0.204 | 0.061 | 0.179 | 0.136 |
     | message passing, noisy-OR | 0.064 | 0.204 | 0.061 | 0.160 | 0.114 |

   - Across symptoms, the part of P from the modules correlates with the gap between the optimum and the leak (0.97
     for linear response, 0.96 for message passing).
   - In the linear-response pilot the module activations are nearly constant (standard deviation at most 0.008). A
     constant activation adds the same amount to every perturbation, so those modules work as a second leak.
   - The correlation alone does not separate a module that calibrates from one that explains the common symptoms,
     since the gap grows with the base rate. The near-constant activations do separate them for linear response.
   - Inferred: this is one cause of the degenerate modules recorded for these pilots.

4. **A shifted validation set made the rise steep (measured).**
   - The fold-0 validation set holds cluster:ARNT2: 58 of its 140 perturbations and 49.6 percent of the validation
     weight.
   - Positives are 14.4 percent of ARNT2's pairs, against 5.7 percent in training. ARNT2 accounts for +0.125 of the
     +0.133 rise in validation loss from epoch 0 to epoch 8. The rest of the validation set is flat (0.370 to 0.385).
   - `--validation-draw rotated_stratified` (second family) keeps groups of more than about 37 perturbations in
     training, so ARNT2 no longer enters validation (inferred from the draw's limit).

5. **Sum pooling gives some perturbations extreme fields (measured).**
   - GLI2's pooled vector has a norm of 117 to 167, against a median of 7, and accounts for 24 percent of the rise.
     Its apathy probability falls from 1.8e-5 at epoch 0 to 2.8e-12 at epoch 8.
   - The pooled norm does not correlate with each perturbation's rise (Spearman 0.05).

6. **Development fold 0 asks the model to extrapolate in seed count (measured, lockbox_v1 and lockbox_v2 alike).**
   - The 52 development drugs with more than four seed genes all belong to cluster:ABCA3, which falls in fold 0. The
     training pool of fold 0 holds no perturbation with more than four seeds; the other folds' pools hold all 52.
   - The 28 lockbox_v2 drugs have one to four seeds, so the lockbox does not ask for this extrapolation.
   - Fold 0 enters the five-fold means of experiments/choose_heads.py.

7. **The probability clamp hides part of the rise; it does not cause it (measured).** At epoch 8, 124 validation
   pairs sit beyond the clamp at 1e-6. In log space the validation loss would be 0.704 instead of 0.667. On the test
   set, 3.7 percent of pairs are beyond the clamp, 13 of them positive.

8. **No defect found in the stopping code (measured).**
   - The validation loss uses the training loss's weights, targets and clamp.
   - best_epoch is the 0-based index of the end-of-epoch evaluation, and the refit trains best_epoch + 1 epochs.
   - The training loss in each run's history is a running mean over the epoch, so the epoch-0 value includes the
     start near P = 0.5.

## The flag `--start-at-weighted-optimum` (default off)

experiments/run_main_model.py starts each symptom's rate at the loss-optimal constant over the fitted perturbations
(`weighted_constant_optimum`: the weighted mean of the loss's targets with the loss's weights):
- for the noisy-OR head, the leak, in place of the raw base rate;
- for the sigmoid head, the output bias.

It applies to both the early-stopped run and the refit. Left off, it changes neither the model nor the configuration
fingerprint, so runs resumed from older checkpoints report no change (tests/test_weighted_optimum_start.py). No
configuration carries it. Whether the confirmatory_v2 configurations take it is the user's decision.

Expected effects (inferred; a full-graph test would be training outside the registered jobs):
- The sigmoid head's first epoch no longer goes to calibration.
- The noisy-OR modules no longer stand in for the leak, so with the module bias at -3 a module stays off unless the
  perturbation drives it.

What it does not change (finding 1): the per-perturbation fit does not lower the validation loss below the constant,
so stopping on the validation loss can still pick an early epoch. On the slice, a sigmoid output bias started at the
base rates kept the validation loss in a flat band of 0.559 to 0.587 around the constant's 0.562.

Other remedies considered:
- Evaluate for stopping every 20 or so steps instead of once per epoch.
- Normalise the pooled field by the number of seeds, or use mean pooling.
- Lower the learning rate. On the slice this moves the minimum later (step 45 at 0.002, step 190 at 0.0005) without
  lowering it.
- Stop on validation AUPRC, which is flat and noisy.

None is applied.
