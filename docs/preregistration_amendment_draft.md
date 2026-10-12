# Draft of an amendment to the pre-registration: the evidence table v4 and the checks

**Superseded on 12 October 2026 by the amendment in `docs/preregistration.md`** (the owner: "Let's lock in everything
you have no objections to"). That amendment registers eleven items of this draft and names what it leaves open:
off-targets in a model input, the enzyme statements as training labels, the weight of the Open Targets class, the
rows of the admitted-groups list, the electrical layer and the held-out set itself. This file is kept as the working
record of how the draft was built, and the registered text is the preregistration's own.

Still true of both: no held-out set has been drawn on the new table and no model has been trained on it. The code for
the propagation and encoder items is not written, so those are registered and have been read on no run. Written on
11 October 2026, superseded 12 October 2026.

The registered state it would amend: the evidence table `evidence_full_v2`, the `better_v2` label selection,
`configs/lockbox_v2.json`, and the graph `data/processed/graph_full_neuronal_split_binders` (amendment of 9 October
2026, fourth).

## 1. What the owner has decided

Each line gives the owner's words and where they are recorded (`docs/research_summary.md`, decision log).

1. **Three label changes** (10 October: "Up to three targets", "Fix gain-of-function sign", "Mask never-above-placebo";
   11 October: "I accept all the suggestions under section 8 of the label review."; "Please do the recommended label
   fixes").
2. **Where the added drugs go** (11 October: "I accept put the three mechanism target pairs being mostly in the held
   out set"): an added drug that acts on a node of a held-out perturbation is held out, as its own stratum outside the
   primary score (`docs/label_source_review.md`, section 9.2, option a).
3. **Combined statements are checks** (10 October: "Build, report as checks"; 11 October: "Add the valproate, etc
   displacements as checks then since those are important and survived"; "Add the enzyme labels as the carrier twin
   instead then").
4. **A reading at the peak of the unbound tissue concentration** (11 October: "Peak tissue unbound concentration can
   be relegated to confirmation only").
5. **The input of a multi-target drug, and the loss shares** (12 October: "let's do the control version and the equal
   loss per share version"). The input stays 1 per target; `--equal-drug-shares` is turned on; three controls are
   reported beside the score.
6. **The occupancy anchor** (12 October): the placeholder that fixes the best-bound target at 90 percent occupancy is
   rejected; where a drug has a measured unbound concentration, occupancy is read at each target's own site.
7. **What the comparative checks test** (12 October): only the two-drug kind is a test of a model; the three kinds that
   go through a lowered carrier or enzyme are reported without a verdict, and the enzyme statements do not train.

## 2. What would change

- **The evidence table.** `data/processed/evidence_full_v4` in place of `evidence_full_v2`, built by
  `scripts/rebuild_data_layer.sh` with `--max-drug-targets 3` and `--seed-gain-of-function-genes-positive`; the
  selection `label_selection/better_v2_full_v4.parquet` with `--mask-grades C --mask-never-above-placebo`. SHA-256,
  first 16: records `e0206dbe181747e8`, reports `50d847913fa8ef9f`, selection `1d0cde3b1292217c`. 1,636 perturbations
  (1,414 genes and 222 drugs), 24 symptoms, 2,975 kept pairs of which 1,076 are drug pairs. Against `evidence_full_v2`
  the table also carries two earlier changes that are not registered either: drugs with a target of any type and the
  salt-form fix (the table called v3, `docs/drug_targets_any_type.md`), and parkinsonism as the 24th symptom (decision
  of 9 October). The amendment would register all of them together. 316 of the 317 perturbations of
  `configs/lockbox_v2.json` are in the table.
- **The sign of 17 genes.** A gene whose every causal association is a gain of function is seeded with sign +1. The
  registered text discloses the old sign and leaves it ("so the signed encoder reads the wrong sign for them"); the
  amendment would replace that sentence.
- **The held-out set.** The membership of `configs/lockbox_v2.json`, carried to v4 by its leakage groups (320
  perturbations, 31 drugs), stays the set the primary score is read on. The 29 added drugs that act on a held-out
  node (`configs/held_out_additions_v4_draft.csv`; 160 kept pairs) are held out and reported as their own stratum.
  They do not train, do not validate and do not enter H1 or H2. The 39 other added drugs train. Under the rule 349 of
  1,636 perturbations, 60 of 222 drugs and 307 of 1,076 kept drug pairs are held out. The registered group rule chains
  and is not used to place the added drugs, because it would pull 149 further perturbations into the held-out set;
  the amendment would say so and state the one-step rule in its place for added perturbations.
- **Checks.** `configs/comparative_checks.csv` (214 statements) is read by `experiments/run_comparative_checks.py`
  after scoring. The checks are named in advance, test no hypothesis of the study, change no decision rule and are
  reported apart from every score. A finished development run is asked the development checks only.
- **The loss and its shares.** The four confirmatory configurations and their four no-descriptor ablations gain
  `--equal-drug-shares`: each drug's positives and, apart, its negatives are rescaled to the mean counts over the
  training drugs, so a drug with many labels does not take a larger share of the gradient. Gene rows are unchanged and
  the evidence weights still differ between drugs. This changes training, which is why it belongs in the amendment and
  not in a commit; it is the one item here that does.
- **What is reported beside the score.** Three controls for the input-size objection, none of which changes a model:
  the input-size baseline (`run_baselines.py --with-input-size`), macro AUPRC inside strata of a drug's number of
  targets (reported without a flag), and the permuted-target control (`run_main_model.py --permute-drug-targets`), a
  negative control whose run is never a result of the study. The amendment should say that the controls were named
  before any v4 run, which they were.
- **Optionally, the weight of a placebo-controlled negative.** If the owner takes the never-above-placebo pairs as
  negatives (point 2 of section 4 below), `--placebo-negative-weight` sets their loss weight above that of an
  unlabelled pair. The amendment has to name the number, because it is a free parameter of the loss.
- **The occupancy anchor, if off-targets ever enter the input.** The measured unbound concentration at each target's
  own site, not the 90-percent placeholder. No off-target enters any registered input today, so this binds nothing
  yet; it is written here so the anchor is not re-chosen later by default. Measured on 12 October and worth stating in
  the amendment as the reason none does: admitting off-targets merges the leakage groups whatever the concentration
  scheme, taking the largest group from 241 perturbations and 71 drugs to between 429 and 449 perturbations and
  between 203 and 209 of the 222 drugs, at every threshold from 0.1 to 0.5 (`docs/label_source_review.md`, 14.3). The
  measured concentration is not what grows the group; it admits fewer off-targets than the placeholder and leaves a
  group two drugs smaller. A mixed input, measured where there is one and a unit per target otherwise, is not
  recommended: the two groups then carry a median total input of 1.00 against 3.0, and which group a drug falls in
  follows whether its pharmacokinetics have been published.
- **The comparative checks.** Named in advance, reported after scoring, testing no hypothesis. The amendment should
  state which of them can test a model at all: the two-drug kind, because a second drug enters through the mechanism
  edges that every drug perturbation trains. The carrier, displacement and enzyme kinds go through relations no label
  varies and are reported without a verdict.
- **The peak reading.** Not built. The owner's sentence is read here as: a reading at the peak concentration would be
  reported at the confirmation stage, if at all, and is not trained on or selected on. The owner should correct this
  reading if it is wrong.

## 3. What would not change

The two tested models and their heads (`configs/head_choice.json`), the decision rule, the one-sided 0.025 sequence of
H1 then H2, the floors, the degree strata, the confirmatory graph and the node descriptor table. No drug-node graph and
no clearance edge is registered: `graph_full_neuronal_split_binders_drugs_v4` and its variant with clearance edges
exist for the ablation the owner asked for on 10 October and for the checks, and whether a confirmatory configuration
reads either is a separate decision.

- **The propagation window, if the owner applies his decision of 12 October.** "Use the convergence instead of 8
  steps. 8 steps is arbitrary, not physiological." For the linear-response encoder this is a stopping rule and not a
  new architecture: the encoder holds the perturbation on at every step and its gains are bounded by the contraction,
  so a fixed point exists and the loop can run to a tolerance. The amendment has to name the tolerance, because it
  replaces a fixed count as the registered quantity; 1e-3 costs 23 steps at the contraction bound, 14 at gain 0.7 and
  11 at 0.5. For the message-passing encoder it is not a stopping rule, since each layer owns its weights and more
  rounds means more parameters, so convergence there means holding the perturbation on and sharing the relation
  weights across rounds, which replaces the registered encoder rather than joining it. That second change would make
  the registered head choice and stopping rule invalid, because both were selected on three independent layers, so it
  needs its own development runs before any held-out run (`docs/label_source_review.md` 17.2).

- **The message-passing encoder, if the owner applies his decision of 12 October.** "use the 4 bases and I think
  just the 1 layer, repeated." One set of relation transforms, each a mixture of four shared bases, applied
  repeatedly until the change falls below the tolerance of the item above. The relation transforms fall from 61,440
  to 4,176 and stay there at any depth, and the `clears` relation owns four coefficients, the same as the
  linear-response gain. Three things the amendment has to state, not leave to the implementation: the number of
  bases, because it is a capacity choice and should be read off development folds against 2, 4, 8 and the full
  matrix; the norm bound on the bases and their coefficients, because a weight-tied iteration converges only if the
  update is a contraction and a non-expansive activation alone does not give that; and the loss of depth-specialised
  transforms, since a repeated layer applies one transform per relation at every round where three independent
  layers do not. Neither the noisy-OR head nor the training stopping rule is reopened by this: the head was settled
  on 8 October and the stopping rule is preregistered as validation loss with patience 8 and at most 60 epochs
  (`docs/label_source_review.md` 18, 19). **Open on the owner's challenge of 17:09**: four shared bases assume all 20
  relation matrices lie near one four-dimensional subspace, which nothing tests. A per-relation low-rank
  factorisation assumes only that each transform is low rank and costs 5,120 against 4,176. The owner's later
  criterion supersedes the edge count: the rank follows the training evidence that passes through each relation, and
  measured on the 896 development perturbations with a kept pair, nine relations carry 99.5 percent of the flow and
  the ordering is not the ordering by edge count (`encodes` 33.7 percent on 12,736 edges against `activates` 8.7
  percent on 69,740). The shape the amendment would register is a diagonal gain for every relation as the floor, 640
  parameters, with rank above it bought by flow share out of a budget the owner names, and the ranks fixed before any
  run and measured on development perturbations only (21). The draft keeps "use the 4 bases" until the owner
  chooses (20.3, 21.3). The loss of depth-specialised transforms is withdrawn as a cost (20.1).

- **The message-passing encoder's edge signs.** Verified on 12 October: `build_relation_adjacencies` gives every
  entry of every relation's adjacency the positive value `1 / in_degree`, and the encoder is passed no edge sign for
  ordinary relations, where the linear-response encoder is passed `data.edge_sign` and uses it. For the 12 relations
  of constant sign nothing is lost, since the relation index carries it; for the 8 whose sign varies within the
  relation, 7,321 negative edges are read as positive, 7,087 of them transcriptional repressions among
  `regulates_transcription_of`'s 41,638 edges. The fix is to make the adjacency value `sign / in_degree`, which costs
  no parameters, and the amendment has to state one consequence: with signs in the values the in-degree mean can
  cancel, so a node with one activator and one repressor of equal weight receives zero. Reported in the architecture
  review of 10 October and not acted on (20.2).

- **The electrical layer, which is a separate decision from every parameter question.** Measured on 12 October:
  `VM_c`, the single membrane-potential node, carries none of the 138 descriptor columns, and the in-degree mean
  divides each of its 405 incoming `changes_membrane_potential` messages by 405 while its seven `voltage_gates`
  outputs pass undivided. A membrane potential is the integral of its currents, so the mean is the wrong
  aggregation for that relation, and the layer is nearly silent whatever the label count. The amendment has to
  settle three things apart: the aggregation for `changes_membrane_potential` and `changes_ion_pool` (a sum scaled
  to keep the transition a contraction, since a raw sum would take the row sum to 405 and break the convergence item
  above); whether `VM_c` is given descriptor columns of its own; and, if neither is changed, a disclosure that the
  electrical layer is present in the registered graph and carries no fitted behaviour
  (`docs/label_source_review.md` 22). **One change would settle all three**: compute `VM_c` from the ion pools
  through a linearised Goldman-Hodgkin-Katz voltage equation instead of learning it. The 405 direct edges and their
  mean go, since 348 of those reactions already change an ion pool; `VM_c` takes eight edges with weights fixed by
  the voltage equation and no free parameters; and the featureless hub stops mattering because its state is computed.
  It needs three missing neuronal extracellular pools, the 57 Reactome potassium and sodium channels wired to the
  pools, the constants recorded with a citable source, and one statement in the registration: all seven reactions
  `VM_c` gates feed back into it, an excitable membrane is not contractive, so the claim is the subthreshold steady
  state and not an action potential (23).

## 4. What is still missing before this can be applied

1. **The limit on a drug's targets.** Three stays the rule until the owner decides (`docs/label_source_review.md`,
   section 12). Without a limit eight drugs enter and one of them joins the two large drug groups into one unless a
   written rule sets it aside.
2. **The placebo pairs**, and if they become negatives, their weight. Set aside, as decided on 10 October, or counted
   as negatives, as the owner's question of 11 October suggests ("Why would working no better than the placebo not be
   counted as evidence against a label?"). Both selections are built; the second is
   `better_v2_full_v4_placebo_negative.parquet`, SHA-256 first 16 `3f4decacd783ab08`. On 12 October the owner asked
   that a placebo-controlled absence count for more than silence in the loss share; `--placebo-negative-weight` does
   that and needs a number. **The number is read as 0.5**, from the owner's "I think -.5 is best unless you
   disagree" (12 October). The leading minus is read as a typing slip, because the weight multiplies a per-pair
   cross-entropy: at -0.5 the optimiser would maximise the loss on those pairs and train the model to predict the
   symptom for exactly the pairs a placebo-controlled comparison found no excess for. `run_main_model.py` now refuses
   a weight below zero rather than run it (`refuse_a_negative_loss_weight`, `tests/test_input_size_controls.py`).
   0.5 sits where the owner put it, above the 0.2 of a pair that is a negative on no evidence and below the weight of
   a reported event. One caveat belongs in the amendment beside the number: 48 of the 61 pairs rest on a single label
   entry carrying both arms and SIDER gives no arm sizes, so the evidence behind a 0.5 is thinner than the evidence
   behind a positive at its evidence weight. The owner confirms the value.
3. **The draft list of admitted drug groups** (`configs/admitted_atc_prefixes_draft.csv`): its rows need the owner's
   sign-off, and one drug of the list joins the two large drug groups in the same way.
4. **The other questions of `docs/confirmatory_architecture_review.md`, section 8.** Question 6, the input of a
   multi-target drug, was answered on 12 October and is item 5 of section 1 above; the other eleven are open.
5. **The Open Targets statements**: their evidence grade, and whether held-out statements enter the primary score
   (`docs/label_source_review.md`, section 9.3).
6. **Development runs on v4.** The registered head choice and the stopping rule were made on `evidence_full_v2`. A new
   table needs its development runs before a held-out run, and the owner's instructions "No new pilots yet" and "don't
   rerun the same pilots until we've completely nailed down the design options" stand until the points above are
   settled.
7. **The owner's go-ahead for the held-out set** ("Don't finalize the lockbox until I've given the go ahead").

## 5. What the amendment should state plainly

- The label changes were chosen from properties of the evidence (the number of curated targets, the type of a gene's
  association, the frequency in a placebo arm) and never from a model's predictions. No model was trained on v4
  before the changes were fixed.
- The membership of the held-out set was read to place the added drugs; no held-out prediction or score was.
- The 29 held-out added drugs share targets outside the held-out set with 43 development drugs. That is why they are
  their own stratum and do not enter the primary score.
- Catatonia's one kept pair is in the held-out set, so catatonia has no training positive.
