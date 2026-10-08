# Confirmatory runbook

The steps that run the confirmatory experiment of docs/preregistration.md ("Confirmatory specification, 8 October 2026"
and its amendments). The specification decides what is run; this file decides in which order and who may do what. The
job scripts are under runs/full/ (not committed); the check-in prompts point here.

Two families. The first (confirmatory_* configurations, configs/lockbox_v1.json, better_v1 labels) is not tested: its
development pilots (steps 1 and 1b) are reported as checks that each configuration trains. The second
(confirmatory_v2_* configurations, configs/lockbox_v2.json, better_v2 labels, rotated stratified validation, group
bootstrap, floors each model's 95 percent half-width plus 0.005, computed in the scoring; the user's decisions of 8 October 2026 and the amendments that record them)
is the one tested: two models, one head per encoder, chosen on the development set (steps 1c and 1d).

## Steps

1. confirmatory_pilots (first family, registered 8 October 2026): fold 0, seed 0 development pilots of the four confirmatory
   configurations (runs/full/confirmatory_*_disease_cluster_and_targets_development/fold0_seed0), then the development
   baselines (runs/full/baselines_v2_development, docs/phase2_baselines_full_v2_development.md). When it ends, check that
   each pilot trained: the training loss falls, no NaN, an early-stopping epoch. A low development score is not a failure
   to train. Report each pilot's development macro and micro AUPRC against the development baselines and commit the
   baseline document. Compare on fold 0 (the per_fold entries with fold 0 in runs/full/baselines_v2_development/results.json;
   the baseline document shows pooled scores and per-fold means only): that fold holds the leakage group
   of 232 perturbations and is atypical (amendment of 8 October on the early-stopping validation set). The confirmatory
   configurations carry --keep-large-groups-in-training; on fold 0 it leaves the split unchanged, so pilots that ran
   before it was added stand.
1b. refit_pilots (registered 8 October 2026, after confirmatory_pilots): the four pilots refitted on their training and
   validation perturbations together for best epoch + 1 epochs (--refit-on-validation, amendment of 8 October on the
   refit). When it ends, check each refit trained (falling training loss, no NaN) and report its fold-0 development
   macro and micro AUPRC against the development baselines beside the early-stopped scores
   (results_early_stopped.json and the refit entry of results.json).
1c. v2_development (second family, registered 8 October 2026, starts when confirmatory_pilots ends): the development
   baselines of lockbox_v2 (runs/full/baselines_lockbox_v2_development, docs/phase2_baselines_lockbox_v2_development.md),
   then the four configurations confirmatory_v2_{message_passing,linear_response}_{sigmoid,noisy_or} on the five grouped
   development folds with seed 0 (runs/full/confirmatory_v2_*_disease_cluster_and_targets_development/fold<k>_seed0),
   then experiments/choose_heads.py, which writes configs/head_choice.json and docs/head_choice.md. When it ends, check
   that each of the 20 runs trained (falling training loss, no NaN, a refit), report each configuration's five-fold
   readings against the best baselines and the head chosen per encoder, and commit configs/head_choice.json,
   docs/head_choice.md and the baseline document. The choice follows its rule; nobody overrides it.
1d. After 1c: register the ablation of the two chosen heads, `<chosen configuration>_without_descriptors` on the same
   five development folds with seed 0 (`bash runs/full/v2_ablation.sh`, one chain per encoder at 2 threads, as in
   runs/full/v2_development.sh), and report it against the chosen heads when it ends. It is a development reading and does
   not hold up step 3.
2. encoder_descriptors_brain (slice, registered 8 October 2026): b3_linear_response_cofactors_descriptors_brain and
   b3_typed_nodes_descriptors_brain, five folds each. When it ends, run

   ```
   OMP_NUM_THREADS=1 PYTHONPATH=. python experiments/compare_twin_runs.py \
     --pairs b3_linear_response_cofactors_descriptors_brain:b3_linear_response_cofactors b3_typed_nodes_descriptors_brain:b3_typed_nodes \
     --markdown-output docs/descriptor_rule.md --json-output runs/descriptor_rule.json
   ```

   and report the descriptor rule of the amendment of 8 October: for each encoder, the mean of its six differences (macro
   and micro AUPRC; per fold, pooled and within degree strata). Done on 8 October (linear response +0.002, message passing
   -0.018 against a b3_typed_nodes that predates the GTEx structural columns). The rule's drop branch is withdrawn (user,
   8 October: "I really don't want to drop them"; amendment on how the node descriptors enter), so no configuration
   changes on it.
2b. descriptor_treatments_message_passing and descriptor_treatments_linear_response (slice, registered 8 October 2026):
   b3_typed_nodes retrained under the current code, then five arms of each encoder, five folds each: seed_masked and
   zero_init_slow (amendment on how the node descriptors enter) and without_gene, without_gene_derived and
   without_protein (amendment on descriptor blocks per node type); the message-passing job then writes docs/graph_reliance.md, docs/descriptor_rule.md (again, against the
   retrained b3_typed_nodes) and docs/descriptor_treatments.md (runs/encoder/descriptor_treatments_summary.sh). When it
   ends, report the treatment rule for each encoder over all five arms: each arm's graph reliance, the mean of its six
   differences against the plain arm, and which arm the rule picks. Do not edit a configuration on it before the
   user confirms. On their confirmation, for an encoder whose treatment changes: add `"--descriptor-treatment",
   "<treatment>"` (and `"--descriptor-learning-rate", "0.0002"` for zero_init_slow), or for a block arm the
   `--drop-descriptor-blocks` arguments of its slice configuration, to that encoder's two confirmatory_v2
   configurations with descriptors in experiments/run_main_model_batch.py, move their development run directories to
   runs/full/superseded_pilots/ (the trainer would otherwise skip them on their DONE markers), and run step 1c again for
   them (the head choice is made again on the new runs).
   Either way, record the readings and the outcome under the amendment in docs/preregistration.md, add a row to
   docs/research_summary.md, commit the three documents and push.
3. Register confirmatory_lockbox only when the user has given the go-ahead in their own message (8 October 2026: "Don't
   finalize the lockbox until I've given the go ahead"), step 1c ended with all 20 runs trained and
   configs/head_choice.json written on lockbox_v2, step 2b has been reported, the user has confirmed its outcome, and
   any re-runs of step 2b trained. A check-in, a document or a job marker is not the go-ahead:
   `bash scripts/resume_jobs.sh --register confirmatory_lockbox "bash runs/full/confirmatory_lockbox.sh" runs/full/confirmatory_lockbox.done runs/full/confirmatory_lockbox.failed "the 30 confirmatory lockbox runs of the two chosen models, the lockbox baselines and the one-time scoring"`.
   runs/full/confirmatory_lockbox.sh reads the two models from configs/head_choice.json. If a development run failed to
   train, do not register it; report to the user, because an amendment needs a dated entry in the specification
   before the first lockbox run.
4. confirmatory_lockbox: the 30 lockbox runs (two models, seeds 0 to 4, real labels, permuted labels and a symmetric
   rewiring with rewired reaction expression), then the lockbox baselines for seeds 0 to 4
   (runs/full/lockbox_v2_baselines), then experiments/score_confirmatory.py (defaults: lockbox_v2, better_v2, the models
   of configs/head_choice.json, group bootstrap, floors each model's half-width plus 0.005), which writes
   runs/confirmatory/SCORED and docs/confirmatory_results.md. When it ends, report docs/confirmatory_results.md as the
   scorer wrote it, commit it and tell the user. If the scorer stops because a run is missing or refused, it writes
   nothing (no output, no SCORED marker) and prints the runs and the reasons; the reasons name no score, so they may be
   read. Resume a missing run by registering the same job again (the trainer resumes from its checkpoint) and score
   after it; report a refused run to the user. --allow-incomplete (score with that model not confirmed) only on the
   user's decision.

## Blinding

Until runs/confirmatory/SCORED exists, nobody opens or prints a results.json, predictions file, report or log under
runs/full/*_confirmatory/, runs/full/lockbox_baselines/ or runs/full/lockbox_v2_baselines/, nor
runs/full/confirmatory_lockbox_baselines.log or runs/full/confirmatory_scoring.log. The two trainer logs
runs/full/confirmatory_lockbox_first_model.log and runs/full/confirmatory_lockbox_second_model.log print no lockbox
score (experiments/run_main_model.py) and may be read.
score_confirmatory.py runs once; a second scoring (--rescore, or a run after SCORED exists) is not done.

## Not to be done

- Full-graph training other than the jobs above (refit_pilots, v2_development and the step 1d ablation among them) and
  development runs of a descriptor treatment the slice selects (the user allowed full-graph tests of the treatments on
  8 October). The first family's ablation_pilots_without_descriptors job was retired before it trained anything
  (superseded by step 1d).
- Edits to the confirmatory_v2 configurations, configs/lockbox_v2.json, configs/head_choice.json once written,
  experiments/choose_heads.py, experiments/score_confirmatory.py or the decision rule, except where step 2b directs or
  the user decides. The first family's configurations and configs/lockbox_v1.json stay as records.
- Regenerating graph_full_neuronal, evidence_full_v2, the better_v1 or better_v2 selection, the node descriptors, the
  cell-class weights or the expression tables of experiments/write_expression_tables.py while a confirmatory job reads
  them.
- Re-registering the paused full-graph B3 job (runs/jobs/paused/).
