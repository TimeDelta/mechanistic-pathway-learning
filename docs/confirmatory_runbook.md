# Confirmatory runbook

The steps that run the confirmatory experiment of docs/preregistration.md ("Confirmatory specification, 8 October 2026"
and its amendments). The specification decides what is run; this file decides in which order and who may do what. The
job scripts are under runs/full/ (not committed); the check-in prompts point here.

## Steps

1. confirmatory_pilots (registered 8 October 2026): fold 0, seed 0 development pilots of the four confirmatory
   configurations (runs/full/confirmatory_*_disease_cluster_and_targets_development/fold0_seed0), then the development
   baselines (runs/full/baselines_v2_development, docs/phase2_baselines_full_v2_development.md). When it ends, check that
   each pilot trained: the training loss falls, no NaN, an early-stopping epoch. A low development score is not a failure
   to train. Report each pilot's development macro and micro AUPRC against the development baselines and commit the
   baseline document. Compare on fold 0 (the per-fold rows of the baseline document): that fold holds the leakage group
   of 232 perturbations and is atypical (amendment of 8 October on the early-stopping validation set). The confirmatory
   configurations carry --keep-large-groups-in-training; on fold 0 it leaves the split unchanged, so pilots that ran
   before it was added stand.
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
   b3_typed_nodes retrained under the current code, then the seed_masked and zero_init_slow arms of both encoders, five
   folds each; the message-passing job then writes docs/graph_reliance.md, docs/descriptor_rule.md (again, against the
   retrained b3_typed_nodes) and docs/descriptor_treatments.md (runs/encoder/descriptor_treatments_summary.sh). When it
   ends, report the treatment rule of that amendment for each encoder: each arm's graph reliance, the mean of its six
   differences against the plain arm, and which treatment the rule picks. Do not edit a configuration on it before the
   user confirms. On their confirmation, for an encoder whose treatment changes: add `"--descriptor-treatment",
   "<treatment>"` (and `"--descriptor-learning-rate", "0.0002"` for zero_init_slow) to that encoder's two confirmatory
   configurations in experiments/run_main_model_batch.py, move their two development pilot directories to
   runs/full/superseded_pilots/ (the trainer would otherwise skip them on their DONE markers), and register
   `bash scripts/resume_jobs.sh --register confirmatory_repilots "bash runs/full/confirmatory_repilots.sh <configuration> <configuration>" runs/full/confirmatory_repilots.done runs/full/confirmatory_repilots.failed "development pilots of the confirmatory configurations the treatment rule changed"`.
   Either way, record the readings and the outcome under the amendment in docs/preregistration.md, add a row to
   docs/research_summary.md, commit the three documents and push.
3. Register confirmatory_lockbox only when the user has given the go-ahead in their own message (8 October 2026: "Don't
   finalize the lockbox until I've given the go ahead"), step 1 ended with all four pilots trained, step 2b has been
   reported, the user has confirmed its outcome, and any re-pilots of step 2b trained. A check-in, a document or a job marker is not the go-ahead:
   `bash scripts/resume_jobs.sh --register confirmatory_lockbox "bash runs/full/confirmatory_lockbox.sh" runs/full/confirmatory_lockbox.done runs/full/confirmatory_lockbox.failed "the 60 confirmatory lockbox runs, the lockbox baselines and the one-time scoring"`.
   If a pilot failed to train, do not register it; report to the user, because an amendment needs a dated entry in the
   specification before the first lockbox run.
4. confirmatory_lockbox: the 60 lockbox runs (four models, seeds 0 to 4, real labels, permuted labels and a rewired
   graph), then the lockbox baselines for seeds 0 to 4, then experiments/score_confirmatory.py, which writes
   runs/confirmatory/SCORED and docs/confirmatory_results.md. When it ends, report docs/confirmatory_results.md as the
   scorer wrote it, commit it and tell the user.

## Blinding

Until runs/confirmatory/SCORED exists, nobody opens or prints a results.json, predictions file, report or log under
runs/full/*_confirmatory/ or runs/full/lockbox_baselines/, nor runs/full/confirmatory_lockbox_baselines.log or
runs/full/confirmatory_scoring.log. The two trainer logs runs/full/confirmatory_lockbox_linear_response.log and
runs/full/confirmatory_lockbox_message_passing.log print no lockbox score (experiments/run_main_model.py) and may be read.
score_confirmatory.py runs once; a second scoring (--rescore, or a run after SCORED exists) is not done.

## Not to be done

- Full-graph training other than the jobs above, the ablation_pilots_without_descriptors job (development pilots of the
  no-descriptor branch, user, 8 October 2026) and development runs of a descriptor treatment the slice selects (the user
  allowed full-graph tests of the treatments on 8 October).
- Edits to the four confirmatory configurations, configs/lockbox_v1.json, experiments/score_confirmatory.py or the
  decision rule, except where step 2 directs or the user decides.
- Regenerating graph_full_neuronal, evidence_full_v2, the better_v1 selection, the node descriptors or the cell-class
  weights while a confirmatory job reads them.
- Re-registering the paused full-graph B3 job (runs/jobs/paused/).
