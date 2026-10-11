# Confirmatory architectures: what the code does, set against what you have said

Written 10 October 2026 at the user's request: "I'd also like a full review of the confirmation architectures so that
I can be sure you're doing everything the way I expect before finalizing the lockbox because there have been plenty
of miscommunications between us".

Nothing here changes a model, a configuration or a table. Nothing was trained, and no lockbox output was opened.

## How this was made, and how far to trust it

Four readers worked apart from one another and from me. Three read code only (the linear-response encoder; the
message-passing encoder; data loading, heads, loss, training and scoring) and were told to report what the code
does, not what its comments or the documents say. The fourth read documents only and listed every statement of
yours that bears on the models, word for word, with the line that records it. I then checked the claims this review
leans on most against the source myself (the perturbation injection, the edge signs, the conjunction default, the
step loop, the leak start, the drug-entry default, the configured paths) and recomputed the path-weight table of
section 3.

Code references are to commit `2a278dc`; the files cited have not changed since. "Measured" means a reader built
the encoder as the trainer does, untrained, and inspected it. Seed statistics are from `evidence_full_v3`, because
the configured `evidence_full_v2` is not on disk.

Limits. The job scripts under `runs/full/` are outside git and were not opened, so the exact launch flags are not
verified. No run history was read, so nothing here says how far training moves a parameter from its start.

## 1. The points most likely to differ from what you expect

Ranked by how much each bears on what you have said you want. Sections 3 to 7 give the detail and the code lines.

| # | What the code does | What you have said | State |
|---|---|---|---|
| 1 | **Message passing stops at 3 edges.** Three layers, separate weights per layer, the perturbation written once into the initial state. A knockout spends one layer on gene to protein, so it reaches 2 edges past its protein: a median of 176 of 49,574 nodes, and 107 of 1,397 knockouts reach only their gene and protein. | "Why aren't we doing message passing in a way that carries the signal until all sink nodes are hit?" and "wouldn't convergence be preferred for physiological correctness?" | Not built. The registered encoder is unchanged. |
| 2 | **The linear response is not at a steady state.** It runs exactly 8 steps with damping 0.5 and stops. A path of 5 edges enters at 0.36 of its steady-state weight, 7 edges at 0.035, 9 or more at zero. A knockout also spends its first step on gene to protein. | The same two messages. I had described this encoder as a held input under a contraction, which reads as if it settled. It is stable, and it is cut off before it settles. | Not what I led you to believe. Needs your decision (section 8, question 2). |
| 3 | **No drug node exists in the confirmatory graph.** The configuration reads `graph_full_neuronal_split_binders`. Drug nodes, mechanism edges, carrier edges, the presence rule and the change-carrying edge all live in a second directory that no configuration names. A drug is seeded on its target proteins. | "i definitely want drugs as nodes with an ablation on them"; "The whole point of including exogenous drug nodes was to gain the carrier edges for them." | Built, tested, not in the family. Needs an amendment. |
| 4 | **Message passing does not read edge signs.** Activation and inhibition are separate relation types, but eight relation types hold both signs under one weight matrix, among them `regulates_transcription_of` (34,551 activating and 7,087 repressing edges). | No statement of yours covers it. The linear-response module's own text says it keeps the signs "which the message passing encoder's single unsigned relation dropped". | As designed, and easy to miss. |
| 5 | **In message passing the perturbation's sign is a feature, not a direction.** The injected vector is `sign * a + magnitude * b + c` with learned `a`, `b`, `c`. Flipping the sign does not mirror the vector. All 1,397 knockouts inject the same vector. A seed with sign 0 still injects. | You approved sign +1 for the 17 gain-of-function genes. | In the linear response a gain is the exact negative of a loss. In message passing it is whatever the model learns from 17 examples. |
| 6 | **A perturbed node reads its own descriptors.** The treatment is `plain`: the perturbed gene's or protein's own 138 descriptor columns, its node type and its log degree reach the head with no edge used. The module that defines the treatments says of it (`descriptor_treatments.py:10`): `A readout can then learn "which protein is this" and skip the propagation`. | "Seed [masking] would then be included." | Not applied. The preregistration says the picks (seed masking for message passing, a slow zero start for the linear response) "wait for the user's confirmation". |
| 7 | **The registered family cannot start.** The configurations name `evidence_full_v2` and `better_v2_full_v2`, which are not on disk and cannot be rebuilt. `configs/lockbox_v2.json` pins their hashes and refuses the v3 tables. `configs/head_choice.json` does not exist. The scorer's default graph is the merged one. | "v3 is the usable one for now" | Paths, lockbox file and amendment all have to move together. |
| 8 | **A registered model changed when a default changed.** My commit `9df62ee` made the per-channel minimum over complex members the default, so both message-passing configurations now use it although neither names it. The scorer would accept a run made before the flag existed. | "go ahead and implement the 'multiple mins'" | The behaviour is what you asked for. The registration does not say so, and the configuration should name the flag. |
| 9 | **The models have 22 or 23 outputs, not 24.** A symptom is an output only if an evidence row names it. Anhedonia has none. Catatonia has no kept positive on the present table. Psychomotor retardation has 2. A symptom is scored only with 5 kept positives in the lockbox. | "add it as a symptom" (parkinsonism, the 24th) | Parkinsonism is in. Anhedonia and catatonia cannot be learned or scored on the present labels. |
| 10 | **Module membership is one free number per module and node.** 396,592 of the noisy-OR head's 397,270 parameters. It is the same for every perturbation, so it is a learned membership and not an attention. | The design document: parameters are "never per node, since a per-node parameter is memorisation that cannot transfer to a held-out gene". | A conflict between the design sentence and the head. Yours to rule on. |
| 11 | **Every unlabelled pair is a negative** at weight 0.2, all of them, none sampled. A third of the perturbations (517 of 1,568) have no kept positive and train as rows of negatives. | The design document says so too: "Every unobserved pair enters as a negative at reduced weight". | As documented. Listed because it shapes everything the labels do. |
| 12 | **Each lockbox training run writes its own lockbox scores to disk.** `results.json` holds macro and micro AUPRC on the lockbox. Only printing is suppressed. "Scored once" rests on the marker file and on nobody opening those files. | The runbook's blinding rule. | Works only while the rule is kept. A code change could stop the scores being written. |

Three statements of mine that were wrong or incomplete, so you can discount what followed from them:

- I wrote in `docs/label_source_review.md` that the objection to multi-target drugs ("carries total magnitude k")
  "no longer holds with drug nodes". Two errors. The confirmatory graph has no drug nodes (row 3). And with drug nodes
  each mechanism edge still delivers magnitude 1 per target. Measured on the three-target table: 46 drugs carry a
  total input of 2 and 23 carry 3. Section 8, question 6.
- I described the linear response as a held input under a contraction without saying that 8 steps stop it short
  (row 2). The help text of `--propagation-damping` says damping "sets the transient, not the fixed point"; at 8
  steps it changes the output.
- The correction already on record: "reaches 96 to 98 percent" counted any nonzero value. At 0.01 percent of the
  largest value an untrained linear response covers a median of 17 nodes for a knockout and 584 for a drug.

## 2. What the confirmatory family is today

| | Registered (preregistration and runbook) | Your later decisions | On disk or built |
|---|---|---|---|
| Models | two, `confirmatory_v2_message_passing_noisy_or` and `confirmatory_v2_linear_response_noisy_or` | "Linear + one nonlinear", the nonlinear one chosen by tests stated in advance | the two registered encoders; no candidate for the nonlinear slot is built |
| Graph | `graph_full_neuronal_split_binders`: 49,574 nodes, 279,784 edges, 21 relation types, no drug node | drugs as nodes with an ablation | `..._binders_drugs`: 154 drug nodes, 791 mechanism edges, 71 carrier edges; nothing trained |
| Evidence | `evidence_full_v2` | "v3 is the usable one for now"; three label changes approved today | v3 with parkinsonism; a counting build with the three changes |
| Labels | `better_v2`: gene pairs at frequency 0.30 or more, drug pairs at 1 percent or more or in both label sources; the rest masked; grade C pairs masked | the placebo mask, approved today | the flag exists, off by default |
| Lockbox | `lockbox_v2.json`: 317 perturbations, 289 genes, 28 drugs | "Don't finalize the lockbox until I've given the go ahead"; "Decide the lockbox first" | the file pins tables that no longer exist |
| Head | noisy-OR on both encoders; the sigmoid head dropped ("I guess we should drop the sigmoid head then") | | both heads exist in code |
| Training | seeds 0 to 4; early stopping on validation loss, patience 8, at most 60 epochs; a refit on training plus validation | "No new pilots yet"; "3 folds, capped epochs" for the selection tests | nothing has run on the second family |

## 3. The linear-response encoder

**Input.** For each seed, `sign * magnitude` times a learned weight per channel, added at the seed node. A knockout
is one gene node, sign -1, magnitude 1. A drug is the protein nodes of its targets (1 to 23 seeds). The input is
added again at every step. It is a source term: the seed's state is free to move, it is not held at a value.

**One step**, for node `n` and channel `c` (`linear_response_encoder.py:587-616`):

    g[r, c]  = 0.9 * sigmoid(theta[r, c])        (0.9 * tanh for `binds`, which may take either sign)
    m[n, c]  = sum over relations r of g[r, c] * (1 / (deg_r(n) * F(n))) * sum over edges s -> n of r of sign_e * h[s, c]
    h'[n, c] = (1 - d) * h[n, c] + d * w[n, c] * (m[n, c] + u[n, c])

`deg_r(n)` is the number of edges of relation `r` into `n`, `F(n)` the number of relation types that feed `n`,
`w[n, c]` the cell-class weight, `d` = 0.5. Run 8 times from `h(0) = w * u`, with the same parameters at every step.

| Property | What the code does | Line |
|---|---|---|
| Channels | 12, one per cell class (astrocyte, cortical interneuron, ependymal and choroid, excitatory neuron, immune, medium spiny neuron, oligodendrocyte lineage, other inhibitory neuron, other neuron, vascular, dopaminergic neuron, all cells). `--propagation-channels 4` is overridden and has no effect. | `:431` |
| Gains | one per relation and per class: 24 relations x 12 classes = 288, each below 0.9 by construction | `:455-459`, `:566-571` |
| Relations | the graph's 21, plus two carrier relations and `depletes_substrate`, a reverse edge reaction to substrate of sign -1 for every substrate edge that is not a carrier edge | `:521-535` |
| Edge signs | used as stored. `inhibits` is -1; `regulates_transcription_of` keeps +1 and -1 | `:540-548` |
| Averaging | a mean within each relation, then a mean across the relation types that feed the node. One producer and fifty consumers of a metabolite weigh the same | `:559-562` |
| Currency metabolites | 222 nodes send nothing through any relation and still receive | `:537-539` |
| Cell classes | the class weight in [0, 1] multiplies the input at the seed and every incoming message at every step | `:590-592`, `:610-611` |
| Extracellular pool | after every step the 1,684 extracellular metabolites take the plain mean over the 11 classes other than all cells | `:630-636` |
| Stability | guaranteed during training too: gains below 0.9, row sums of the averaging exactly 1, class weights at most 1, so one step multiplies any difference between two states by at most 0.95. No clamp, no test for convergence | `:422`, `:427-428` |
| Linearity | exactly linear in `sign * magnitude`. Sign +1 gives the exact negative of sign -1 (measured difference 0.0). No rectifier in the encoder | |
| Output | `[batch, nodes, 32]`: the 12 channels expanded by a learned map to 32, times a gate in (0, 1) that is a linear map of the node's 160 feature columns through a sigmoid | `:645`, `:658-660` |
| Parameters | 5,836: 288 gains, 12 input weights, 384 for the expansion, 5,152 for the gate | |

**What 8 steps compute.** The weight of a path of `k` edges from a seed, against 1 at a steady state (recomputed
here):

| edges from a seed | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 or more |
|---|---|---|---|---|---|---|---|---|---|---|
| 8 steps (registered) | 1.000 | 0.996 | 0.965 | 0.855 | 0.637 | 0.363 | 0.145 | 0.035 | 0.004 | 0 |
| 16 steps | 1.000 | 1.000 | 1.000 | 0.998 | 0.989 | 0.962 | 0.895 | 0.773 | 0.598 | 0.402 at 9, 0.038 at 12 |
| 40 steps | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0.997 at 12 |

Untrained, the gap between 8 steps and the settled state is 0.1 to 4 percent of the field's total size, because at
random gains most of the size sits on the seed (57 to 79 percent for a knockout). That says little about a trained
model, whose gains are free to grow to 0.9.

**Things in this encoder you may not expect.**

1. A knockout's signal at its own protein is divided by the number of relation types feeding that protein (1 for 296
   genes, 2 for 409, 3 for 363, 4 for 319, 5 for 10). A protein that also receives `binds`, `activates` and `inhibits`
   edges gets a quarter of what a protein with none gets. A drug is seeded on the protein and skips this.
2. The cell-class weight applies twice on the knockout path (at the gene and again into the protein) and once on the
   drug path.
3. The steps stayed at 8 when the gene and protein split added a hop. Message passing got a third layer for it.
4. The seed's own row of the cell-class table is read out directly: at a seed the state starts at
   `sign * magnitude * input_weight[class] * w[seed, class]`, so the head receives the perturbed gene's expression
   profile across classes with no edge used. The no-descriptor ablations keep the cell-class weights, so this path
   stays in them.
5. Each class has its own input weight of free sign (at seed 0, positive for 7 of the 11 pooled classes and negative
   for 4), and the extracellular pool then averages across those classes without aligning their signs.
6. Descriptors enter only through the output gate. They never touch the propagation.
7. Reversibility and transport of a reaction are gate features only. Every reaction propagates as substrate to
   reaction to product plus the depletion edge.
8. At equal gains 8,864 nodes have incoming signs that cancel exactly; 7,685 metabolites are fed by exactly
   `product_of` and `depletes_substrate`.
9. Four drugs have sign 0 on every seed and so give a field of exactly zero, the same as no perturbation. The trainer
   prints a warning and carries on.
10. A defect: the carrier rule's recurring-pair clause keys on the metabolite's base id, and every protein entity has
    the empty id, so all of them count as one species. 86 of the 2,293 carrier edges are carriers only for that
    reason (57 at protein entities, 29 at metabolites: H+ 10, Ca2+ 7, serotonin 5, dopamine 4, GABA 2). Those edges
    get the low-gain carrier relation and, on the substrate side, no depletion edge.
11. On a graph with drug nodes, `--sequestration-carries` is not read by this encoder. The carrier edge is an
    ordinary signed edge into the drug node. Since this encoder's state is already a change from rest, that is the
    change-carrying reading.

## 4. The message-passing encoder

**Starting state.** Every node starts from one shared linear map of its 160 fixed feature columns (22 structural, 138
descriptor). There is no learned vector per node and none per node type; "typed" means the node type is among the
columns. The perturbation adds `sign * a + magnitude * b + c` at each seed (`a`, `b`, `c` learned 32-vectors), once.

**One layer** (`relational_message_passing_encoder.py:491-507`; a reader's independent re-implementation of this
equation agrees with the code to 4.8e-7):

    x'[v] = ReLU( x[v] S_l
                  + sum over relations r other than member_of of (1 / in_r(v)) * sum over edges u -> v of r of x[u] W_{l,r}
                  + (per-channel minimum over the members u -> v of member_of of x[u]) W_{l,member_of}
                  + b_l )

| Property | What the code does | Line |
|---|---|---|
| Layers | 3, each with its own weights. The class has no option to share or repeat a layer | `:191-195` |
| Aggregation | a mean over incoming edges within each relation, then a plain sum across relations | `:257-258`, `:503-504` |
| Edge signs | not read | `:247-261` |
| Edge direction | one way, as stored. 56.6 percent of `substrate_of` edges have no stored edge back from the reaction, so a change at a reaction cannot reach that substrate. This encoder has no depletion edge | `:259` |
| Complexes | the 6,767 `member_of` edges are reduced by a hard minimum per channel, by the trainer's default | `run_main_model.py:830-833` |
| Normalisation, residual, dropout | none. The self term is a learned matrix | `:492` |
| Currency metabolites | they pass messages like any node (13,013 edges leave them), unlike in the linear response | |
| Cell classes | 11 input columns on gene and reaction nodes. No per-class propagation; passing `--cell-class-weights` to this encoder raises an error | `run_main_model.py:428-429` |
| Output | the 3-layer state of the perturbed pass minus the 3-layer state of one unperturbed pass, `[batch, nodes, 32]` | `:335-336` |
| Parameters | 72,928, of which 64,512 are the per-layer, per-relation matrices | |

**Reach** (nodes within `k` directed edges of the seeds, of 49,574):

| | 1 | 2 | 3 (registered) | 4 |
|---|---|---|---|---|
| knockouts, median | 2 | 7 | 176 | 4,950 |
| drugs, median | 25 | 1,523 | 13,822 | 23,823 |

At 3 layers the quartiles for knockouts are 13 to 1,034 nodes; 309 knockouts reach 10 nodes or fewer.

**Things in this encoder you may not expect.**

1. The perturbation is not held on. It is the initial state and nothing adds it again.
2. Sign and magnitude are two inputs to a map with a bias (table of section 1, row 5). The bias is added once per
   seed, so the number of seeds of a drug enters the field.
3. A single incoming edge of one relation weighs the same as the average of 2,048 incoming edges of another, because
   the mean is per relation and the relations are then summed.
4. The minimum over complex members acts on learned channels that the code ties to no quantity, and at the first
   layer on states that have passed no rectifier and can be negative.
5. Outside the reach the difference is not exactly zero when the batch holds two or more perturbations: the same 11
   nodes carry about 2.4e-7, because the perturbed pass and the reference pass run with different batch shapes.
6. The reference pass is recomputed for every batch and gradients flow through it.

**The drug-node options, which the confirmatory configuration does not use.** On the graph that has drug nodes:

| Option | What the code does |
|---|---|
| Seeding | the drug's own node, sign 1, magnitude 1 |
| Presence | a drug node's state is set to zero before the first layer and after every layer in every perturbation that does not seed it. You confirmed this: "The presence rule implementation is correctly inferred." |
| Mechanism edges | summed into their targets with no in-degree divisor, at every layer: `weight * (x U) + sign * weight * (x V)` with two extra matrices per layer. The sign multiplies the magnitude here, unlike at the injection |
| `--drug-mechanism-before-first-layer` | adds the mechanism messages once before the first layer, so a drug reaches 3 edges past its targets and not 2. Off by default |
| Carrier edges, `change` (the default on such a graph) | the edge sends the carrier's state minus the carrier's state in an unperturbed row. Its sign of -1 is not read |
| Carrier edges at 3 layers | silent for every one of the 49 drugs with a carrier when the drug is given alone: a carrier would have to lie within 2 directed edges of a seed to have changed, and none does. The edge speaks only in an example that changes the carrier, which is what the combined perturbations you approved are for |
| Default | the trainer reads an unset `--drug-entry` as `nodes` (`run_main_model.py:944`); the loader, the baselines and the scorer default to `targets`. Pointing a run at the drug-node graph would seed drugs one way for the model and another for the baselines unless each is told |

## 5. From field to prediction: the noisy-OR head

Both encoders hand the head a field `[batch, 49,574 nodes, 32]`. The head receives nothing else: no perturbation
type, no degree (`noisy_or_pathway_module_model.py:182-260`).

1. **Module membership.** 8 modules. Each has one free parameter per node, turned into a gate in [0, 1]. In training
   the gate is sampled once per step (at the start 36 percent of gates are exactly 0 and 8 percent exactly 1); in
   evaluation it is the expected gate.
2. **Gated sum.** For each module, the field is summed over nodes with the gate as weight. No division by the number
   of nodes or by the gate's total.
3. **Activation.** `a[k] = sigmoid(w[k] . pooled[k] + c[k])`, with `c[k]` starting at -3.
4. **Symptoms.** `P[s] = 1 - (1 - leak[s]) * product over k of (1 - link[k, s] * a[k])`.

| Consequence | Detail |
|---|---|
| The prediction never falls below the leak | a module can only add risk. The leak starts at the constant that minimises the loss and trains at 0.0002, the slowest rate in the model |
| A module is not silent without a signal | with a zero field its activation is `sigmoid(-3)` = 0.047 |
| The number of modules is fixed at 8 | the description-length selection function exists and is never called |
| The description-length term | about 2.7 at the start against a cross-entropy below 0.7; it is in the training loss and not in the validation loss that stops training |
| Half the link and leak parameters are unused | a second slice for the "relieves" relation is never read in the forward pass but is counted in the penalty |

What the model can use without the graph: the perturbed node's own features (row 6 of section 1), the perturbation's
type (a knockout is one gene node; a drug is 1 to 23 protein nodes with fractional magnitudes), base rates through
the leak, and node identity through the gates.

## 6. Loss and training

    weight[p, s] = mask[p, s] * (evidence weight of the pair if positive, else 0.2)
    target[p, s] = 0.99 if positive, else 0
    loss(batch)  = sum of weight * cross-entropy(target, P clipped to [1e-6, 1 - 1e-6]) / sum of weight
                   + 1e-6 * description length            (noisy-OR only)

| Item | What the code does |
|---|---|
| Positive weight | the evidence weight, 0.15 to 1.0 (grade A 0.25 to 1.0, grade B 0.15 to 0.6) |
| Masked pairs | weight 0: positives set aside by the selection, and pairs with grade C evidence only |
| Class balance | none beyond the two weights. The loss is normalised within each batch, so batches count equally whatever their positives |
| Optimiser | AdamW, 0.002, weight decay 1e-4; links 0.02, leaks 0.0002, gates 0.05 without decay. No schedule, no clipping, no dropout |
| Batch | 16 perturbations, each with all its symptoms |
| Stopping | validation loss, patience 8, at most 60 epochs |
| Validation set | one of 7 grouped folds (the seed modulo 7), drawn from groups no larger than a quarter of the expected set: about a seventh of the eligible perturbations, which is less than 15 percent |
| Refit | a fresh model from the same initial weights, trained on training plus validation for best epoch + 1 epochs. More steps than the run that chose the epoch count, and no held-out check |
| Flags with no effect | `--init-leak-from-base-rate` is overwritten by `--start-at-weighted-optimum` on the next line; `--keep-large-groups-in-training` is not read under the rotated validation draw, though the scorer requires it to be recorded |
| A failed run | the batch runner prints a line and goes on; its exit status stays zero |

## 7. Scoring

| Item | What the code does |
|---|---|
| Test set | the lockbox perturbations; predictions of the refit model |
| Metrics | macro AUPRC over symptoms with 5 or more kept positives in the lockbox; micro AUPRC over all pairs of symptoms with 5 or more kept positives in the whole data |
| Baselines | popularity; popularity times log degree; random walk with restart; TransE. The best one per reading is picked on the lockbox labels and held fixed |
| H1 | the largest p-value of four differences from the best baseline (macro, micro and the same two within degree strata), with both floors met |
| H2 | H1, and both differences from the same configuration trained on a rewired graph |
| Level | one-sided 0.025 per model, H1 then H2. Nothing is divided across the two models; the output states the bound of 0.05 |
| Bootstrap | 4,000 draws of whole leakage groups |
| Floors | the half-width of the same bootstrap interval plus 0.005, so each floor is computed from the result it gates |
| Rewiring | swaps within each relation, `encodes` held fixed ("Encodes should not be rewired in H2."), seeds and the seed nodes' own features untouched |

Three points to weigh:

- **None of the four baselines knows whether a perturbation is a gene or a drug, and the model does.** Type-specific
  base rates are something the model can express and the baselines cannot, and H1 does not separate that from graph
  content. H2 does. The preregistration lists per-type baselines as "Open, for the user's decision".
- **The within-strata readings remove credit for degree only between quintiles.** The top stratum spans degree 59 to
  1,953.
- **The label permutation mixes genes and drugs**, so a drug can receive a gene's label row.

## 8. What I need from you

Each question has the option I would take and why. None is acted on until you answer.

1. **Depth of message passing.** You asked for propagation to all reachable nodes, with convergence preferred. The
   registered encoder cannot do that (section 1, row 1). The two candidates I proposed and have not built: the
   registered 3-layer block repeated with shared weights, and a steady-state encoder with the perturbation held on.
   You have said two variants, linear response plus one nonlinear. I would build both candidates with tests and no
   training, then run the selection you approved once the held-out set is fixed. Is that still what you want?
2. **Steps of the linear response.** 8 steps, 16, or iterate to a tolerance? Iterating to a tolerance matches what
   you said about convergence, and the stability bound guarantees it ends. At a shrink factor of 0.95 per step a
   tolerance of 1e-4 needs at most about 180 steps; in practice far fewer at trained gains below the bound. Cost is
   linear in steps (4.9 seconds per training step at 8). I would iterate to a tolerance with a cap, and report the
   steps used.
3. **Drug nodes in the confirmatory graph.** Yes, with target seeding as the ablation, as you said? If yes: mechanism
   delivered before the first layer or not, and the carrier edge carrying change (both were proposed to you and are
   not confirmed).
4. **Descriptor treatment.** Apply seed masking to message passing and the slow zero start to the linear response,
   the picks the registered rule made? I would apply both. Without them H1 can pass on "which node is this" alone.
   Seed masking hides the 138 descriptor columns at the seed; the node type and the log degree stay visible, and the
   readings within degree strata are what answer for degree.
   The cell-class readout at the seed (section 3, item 4) is a second path of the same kind that neither treatment
   closes; I would close it by leaving the class weight off the input at the seed.
5. **Sign in message passing.** Leave the sign as a learned feature, or make the injected vector flip with the sign
   (inject `sign * magnitude * a`, no bias)? The second makes a gain of function the mirror of a loss at the input,
   which is what the sign fix assumes. The rectifiers after it still make the response differ, as they should.
6. **Input of a multi-target drug.** Each target at magnitude 1 (a three-target drug carries 3), or the drug's total
   at 1? Pharmacology favours 1 per target: occupancy at one target does not fall because the drug has others. The
   cost is that the size of the field then tracks the number of targets. I would keep 1 per target and report
   results by number of targets.
7. **Per-node module gates.** Keep them, or build membership from node features so that it transfers to a node no
   training example reached? Your design sentence excludes per-node parameters. This is the largest change on the
   list and would need its own pilots.
8. **Outputs with no labels.** Drop anhedonia and catatonia from the output list until they have labels, or keep
   them as unscored outputs? Under the three-target rule catatonia has one kept positive, and it falls in the
   held-out set.
9. **The carrier-rule defect** (section 3, item 10). Fix it? It changes 86 edges of the linear response.
10. **Baselines per perturbation type.** Add them to the comparison? I would, since without them H1 credits the model
    for telling a gene from a drug.
11. **Lockbox scores in `results.json`.** Stop the trainer writing them on a lockbox run? I would.
12. **Name every default that matters in the configuration**, so that a changed default cannot change a registered
    model again (the complex minimum, the drug entry, the steps, the treatment). I would.

## 9. Your statements, checked one by one

"Matches" means the code read here does what the statement asks in the confirmatory configuration.

| Your words | Recorded at | Code | Verdict |
|---|---|---|---|
| "the confirmatory config should definitely be the split graph and include the protein descriptors and brain expressions, compartments, everything that physiologically makes sense." | `docs/preregistration.md:842-843` | split graph; 64 protein components on protein nodes, 27 expression columns on gene nodes, compartments among the structural columns | matches |
| "yes, include binder edges in confirmatory graph" | `docs/preregistration.md:1019` | the configured graph is the binder variant, 27 `binds` edges | matches |
| "The genes shouldn't need descriptors but the proteins should." | `docs/gene_protein_split.md:12` | gene nodes carry the expression columns only | matches |
| "I would also require an equivalent several proteins from one gene." | `docs/gene_protein_split.md:19-21` | after a late correction CALCA is the only gene with several protein nodes | built to the rule; one live case. Worth your look |
| "Seed [masking] would then be included." | `docs/gene_protein_split.md:21` | treatment is `plain` | **differs** (question 4) |
| "the registered degree strata SHOULD change to the new split" | `docs/gene_protein_split.md:16-18` | the degree for the strata hops over `encodes`, which reproduces the merged graph's degree for 1,527 of 1,539 perturbations; the preregistration calls the strata "preserved, not redefined" | **two documents disagree, and the code follows the second**. Which did you mean? |
| "Encodes should not be rewired in H2." | `docs/preregistration.md:879` | `encodes` held fixed | matches |
| "Cell types enter as propagation channels with expression-weighted edges, one per cell class, not as graph copies; a dopaminergic class is required" (the log's paraphrase of your direction) | `docs/research_summary.md:157` | linear response: 12 class channels including a dopaminergic one. Message passing: classes are input columns only | **matches for one encoder**; the preregistration records the asymmetry |
| "go ahead and implement the 'multiple mins'" | `docs/research_summary.md:224` | per-channel minimum over complex members, message passing only, on by default | matches; not registered (section 1, row 8) |
| "i definitely want drugs as nodes with an ablation on them" | `docs/drug_entry_nodes.md:71` | built on a second graph directory; no configuration reads it | **built, not in the family** |
| "The presence rule implementation is correctly inferred." | `docs/drug_entry_nodes.md:138-139` | as described in section 4 | matches, on that graph |
| "The whole point of including exogenous drug nodes was to gain the carrier edges for them." | `docs/drug_entry_nodes.md:210` | 71 carrier edges exist on that graph; silent for a drug given alone | **built; has an effect only with carrier-varied examples** |
| "Why aren't we doing message passing in a way that carries the signal until all sink nodes are hit?" | `docs/drug_entry_nodes.md:349-350` | 3 layers | **not built** (question 1) |
| "something more creative like making a schema that reuses some params across layers and not others or having multiple layer param sets that get collated" | `docs/drug_entry_nodes.md:352-355` | each layer has its own weights; no sharing option exists | **not built** |
| "wouldn't convergence be preferred for physiological correctness?" | `docs/drug_entry_nodes.md:352-355` | neither encoder iterates to convergence | **not built** (questions 1 and 2) |
| "I don't want to add extra confirmatory model variations because that starts increasing probability of winning by luck." | `docs/research_summary.md:241` | two tested models in the registration; the second is still message passing as registered | the count matches; the nonlinear slot is open |
| "I guess we should drop the sigmoid head then" | `docs/preregistration.md:804` | the sigmoid configurations exist and are not to be run; `configs/head_choice.json`, which would name the two tested models, is missing | matches in the documents; the file has to be written |
| "Ya. Start at the optimum." | `docs/preregistration.md:817` | leaks start at the loss-optimal constant | matches |
| "I really don't want to drop them" (the descriptors) | `docs/preregistration.md:244` | descriptors in both encoders | matches |
| "Which means don't drop proteins" | `docs/preregistration.md:344` | no configuration drops a block | matches |
| "I guess 64 is fine for now." | `docs/preregistration.md:886-887` | 64 protein components | matches |
| "it would be very bad clinically if that side effect was not predicted so add it as a symptom" | `docs/research_summary.md:212` | parkinsonism is an output on the v3 table with parkinsonism | matches on that table; the registered v2 table lacks it |
| "v3 is the usable one for now" | `docs/research_summary.md:217` | the configuration names v2 | **differs** (section 1, row 7) |
| "Don't finalize the lockbox until I've given the go ahead" | `docs/confirmatory_runbook.md:97-98` | nothing drawn, registered or scored | kept |
| "No new pilots yet" | `docs/preregistration.md:1067` | nothing trained since | kept |

Decisions built without a statement from you, each reversible: mechanism edges summed with no divisor; leakage groups
and degree strata taken from the targets in both arms of the drug-node ablation; the rewiring control holds mechanism
and carrier edges fixed; seed masking for a drug covers its targets; a drug node has no descriptor; input
normalisation is refused with drug nodes. `docs/drug_entry_nodes.md:134-168` lists them.

## 10. Documents and help texts that are out of date

| Where | Says | Is |
|---|---|---|
| `docs/research_summary.md:50-56` | four tested models, a bootstrap over perturbations, floors of 0.041 and 0.028 | two models, a bootstrap over groups, floors from the half-width |
| `configs/evaluation.yaml`, `configs/model_noisy_or_pathway_modules.yaml` | lockbox_v1, two layers, the module count chosen by description length | read by no code; lockbox_v2, three layers, eight modules |
| `docs/confirmatory_runbook.md:136-138` | names the first family's graph | the binder variant of the split graph |
| help text of `--propagation-damping` | "sets the transient, not the fixed point" | changes the output at 8 steps |
| help text of `--normalisation in_degree` | divides "by the number of edges feeding the node" | divides by the in-degree under the relation times the number of relations feeding the node |
| help text of `--node-features typed` | "fixed structural features only" | two expression columns and, with descriptors, 138 more |
| module text of the message-passing encoder | a mechanism edge is "the edge form of the injection" | the injection is `sign * a + magnitude * b + c`; the edge is `magnitude * (x U) + sign * magnitude * (x V)` |
| `docs/preregistration.md:831-832` | 37 genes have several protein nodes | one does, after the late correction |
| `docs/label_source_review.md`, section 4 | the magnitude objection "no longer holds with drug nodes" | it holds (section 1) |

I have corrected the last one. The others wait for the amendment, since several change with your answers.

The readers' full notes (file and line for every statement, and the open questions each could not settle from code)
are outside the repository. I can add them under `docs/` if you want them kept with the project.
