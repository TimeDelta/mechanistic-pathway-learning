# First-family development pilots against the development baselines

Steps 1 and 1b of docs/confirmatory_runbook.md, read on 8 October 2026 when confirmatory_pilots and refit_pilots had
ended. The first family (confirmatory_* configurations, configs/lockbox_v1.json, better_v1 labels on evidence_full_v2)
is not tested: these readings check that each configuration trains. Development fold 0, seed 0, the 317 lockbox
perturbations removed: 837 training, 140 validation and 245 test perturbations; the refit trains on the 977 training and
validation perturbations for best epoch + 1 epochs. Sources: results_early_stopped.json and results.json in
runs/full/confirmatory_<configuration>_disease_cluster_and_targets_development/fold0_seed0, and the fold 0 entries of
`per_fold` in runs/full/baselines_v2_development/results.json (docs/phase2_baselines_full_v2_development.md shows pooled
scores and per-fold means only). Both score the same 245 test perturbations and average over the symptoms with at least
five positives among them (18 of 22).

## Did each pilot train

| configuration | epochs | best epoch | training loss, first to lowest | validation loss, first to best | refit epochs | refit training loss, first to last | NaN |
|---|---|---|---|---|---|---|---|
| confirmatory_message_passing_sigmoid | 9 | 0 | 0.578 to 0.304 | 0.534 (epoch 0) | 1 | 0.596 (one epoch) | none |
| confirmatory_message_passing_noisy_or | 17 | 8 | 1.616 to 0.357 | 0.663 to 0.519 | 9 | 1.523 to 0.389 | none |
| confirmatory_linear_response_sigmoid | 13 | 4 | 0.567 to 0.328 | 0.539 to 0.494 | 5 | 0.564 to 0.397 | none |
| confirmatory_linear_response_noisy_or | 17 | 8 | 1.613 to 0.383 | 0.676 to 0.497 | 9 | 1.524 to 0.416 | none |

All four trained by the runbook's checks: the training loss falls, nothing is NaN, early stopping ends each run, and
each refit ran. The message-passing sigmoid pilot stops at epoch 0. Its validation loss never improves on the first
epoch, so its scores are those of a model trained for one epoch. Its one-epoch refit has no second point to show a
falling loss. The research summary row of 8 October ("Best epoch 0 explained") records why: calibration sets the
stopping epoch. The user dropped the sigmoid head the same day, so this configuration is not among the tested models.

## Fold 0 development scores

| model | macro AUPRC, early-stopped | macro AUPRC, refit | micro AUPRC, early-stopped | micro AUPRC, refit |
|---|---|---|---|---|
| confirmatory_message_passing_sigmoid | 0.131 | 0.131 | 0.257 | 0.251 |
| confirmatory_message_passing_noisy_or | 0.162 | 0.175 | 0.280 | 0.276 |
| confirmatory_linear_response_sigmoid | 0.237 | 0.241 | 0.312 | 0.318 |
| confirmatory_linear_response_noisy_or | 0.143 | 0.152 | 0.298 | 0.288 |
| popularity | 0.130 | | 0.281 | |
| degree_popularity | 0.230 | | 0.316 | |
| random_walk_with_restart | 0.163 | | 0.130 | |
| knowledge_graph_embedding_transe | 0.160 | | 0.124 | |

One fold and one seed, with no interval; these numbers cannot carry a comparison between configurations. Fold 0 is
atypical: it holds the leakage group of 232 perturbations (63 drugs) and most positives of seven symptoms (amendment of
8 October on the early-stopping validation set). Degree-scaled popularity scores 0.230 macro there against a per-fold
mean of 0.159, which suggests that ordering perturbations by degree is unusually effective on this fold.

What the scores show:
- Only the linear-response sigmoid pilot reaches degree-scaled popularity, by 0.007 to 0.011 macro and -0.004 to +0.002
  micro. Those differences are smaller than the fold-to-fold spread of any baseline.
- Both noisy-OR pilots, the head the user chose for the tested family, sit 0.055 to 0.087 macro below degree-scaled
  popularity. Their micro AUPRC is close to plain popularity (0.276 to 0.298 against 0.281). The refit adds 0.009 to 0.013
  macro to each noisy-OR pilot and lowers micro by 0.004 to 0.010.
- Every pilot beats the random walk and TransE on micro AUPRC, by 0.12 to 0.19. Only the linear-response sigmoid pilot
  and the message-passing noisy-OR refit (0.175 against the random walk's 0.163) beat them on macro.

Graph paring (docs/preregistration.md, "Graph paring"): pared graph variants are tested on the development set only if
the development pilots do not beat the development baselines. On fold 0 three of the four pilots do not beat
degree-scaled popularity, and the fourth matches it within a margin one fold cannot resolve. The rule does not say
whether it reads the first family's pilots or the second family's development runs (step 1c), which carry the module fix,
lockbox_v2 and better_v2 labels. That is for the user to decide; nothing has been started on it.
