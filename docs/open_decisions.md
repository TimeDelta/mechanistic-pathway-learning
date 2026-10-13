# Open decisions

Every decision waiting on the owner, in one place, as of 12 October 2026. The registered subset of this list is also
in the "What is still open" block of `docs/preregistration.md`, which is canonical for anything registration touches;
this file adds the build decisions and the actions only the owner can take, which that block does not carry.

Nothing on this list is decided by a session. A session may measure an option and recommend, and several entries record
a recommendation and an objection, but the choice is the owner's.

## 1. Gating the next run

| | decision | what hangs on it |
| --- | --- | --- |
| 1.1 | **The held-out set.** Drawn only on the owner's go-ahead ("Don't finalize the lockbox until I've given the go ahead") | every confirmatory run. Nothing downstream can start |
| 1.2 | **The v4 rebuild.** The table in use still carries the three-target limit; the amendment of 12 October registers no limit with orphenadrine set aside by the stated bridging rule, which needs the table rebuilt and the held-out set redrawn from its leakage groups | 1.1, and the drug counts in every later item |
| ~~1.2b~~ | **Closed 12 October: one pass, orphenadrine only** (amendment item 16). Tapentadol stays in; its noradrenaline reuptake inhibition is part of what the study represents, so removing it would remove a mechanism rather than an artefact. The cost is registered: the opioid and monoamine blocks stay joined, largest block 50.2 percent of drugs, no worse than the 51.4 percent of the table in use | closed |
| ~~1.2c~~ | **Closed 12 October: the rebuild is approved** (owner: "Rebuild package numbers seem acceptable") | closed |
| 1.4 | **The Open Targets class: whether an activation is admitted as a perturbation, and what the agreement is on the rebuilt table.** Amendment item 12a. The weight is built (0.45, with the registered 0.0 arm); the class is not. The crosswalk half is closed: the owner ruled out "weight gain" to increased appetite and "sexual dysfunction" to decreased libido on 13 October, which removes both symptoms from the class. Two things remain. 35 of the remaining 65 statements are activations, which are perturbation rows the study does not have rather than labels on rows it does. And the agreement the weight is derived from does not reproduce: 0.6418 on the capped table against the registered 0.7606, which would give 0.35 rather than 0.45 | the class adds perturbation rows, so it changes the leakage groups: it enters at the same rebuild as 1.2 or not at all, which puts it ahead of 1.1. The agreement is recomputed on that rebuild and the weight follows it |
| 1.3 | **Off-targets in a model input.** The owner asked for them on 12 October; my objection is recorded and stands: admitting them at any threshold from 0.1 to 0.5, under any concentration scheme, takes the largest leakage group from 71 drugs to between 203 and 209 of 222, which ends the drug split item 2 of the amendment depends on | the split, so 1.1 |

## 2. What trains

| | decision | state |
| --- | --- | --- |
| ~~2.1~~ | **Closed 12 October: they train, in both encoders** (owner: "Implement"; amendment item 17). Symmetric, because the asymmetry I warned of rested on `clears` being the only route and that was wrong — `catalyzed_by` carries them too |
| ~~2.2~~ | **Closed 12 October: they stay out** (owner: "That's fine") |
| ~~2.3~~ | **Closed 12 October: accepted** (owner: "Weight accepted") |

## 3. The drug set

| | decision | state |
| --- | --- | --- |
| 3.1 | **The four held drugs of amendment item 14.** Pindolol (penetrant at 0.50, no label sentence of its own), oxprenolol (penetrant at 0.20, no US label in openFDA), timolol (a central-action sentence, no measured Kp,uu,brain), carvedilol (neither). Each needs one part supplied, or to stay out | the rule is registered; these four are its edge |
| 3.2 | **Whether to read label sentences by hand for the positives the screen shortlists.** The screen's negatives are evidence; its positives include class boilerplate, so a reading decides. Five or six drugs | 3.1 |
| ~~3.3~~ | **Closed 12 October: not pursued** (owner: "Don't bother since it doesn't seem feasible") |

## 4. The model and the encoder

| | decision | state |
| --- | --- | --- |
| 4.1 | **The reach of the electrical layer.** Item 13 registers the computed potential; the magnitude arriving at `VM_c` is set by the ion pool's in-degree, and both normalisation alternatives the encoder carries were measured worse. So it is disclosed rather than fixed, unless the owner wants a new aggregation | disclosed |
| 4.2 | **The cell-class weight of a protein complex.** 413 reactions regulated only by a complex carry a flat class weight, adenylyl cyclase among them, because the weight reaches a reaction through its catalysing gene. 1,527 of the 1,681 complexes have a member that already carries one, so a complex could take its limiting member's weight | measured; a build, not yet proposed as a change |
| 4.3 | **Whether the cell-class conditioning should reach the propagation.** Today the weights gate the input at the seeds only, so each class's channel traverses the same structure and differs only in its gains | a registered-encoder change |
| 4.4 | **The weighted constant start** (`weighted_constant_optimum`): whether the confirmatory_v2 configurations take it | no configuration carries it |
| 4.7 | **The head refactor, in the two-level noisy-OR form the owner proposed** (13 October: "that needs a modification for the noisy or equivalent ... what about just putting the noisy or in between the modules and the lower rank matrix?"). `z_j(p) = 1 - prod_k (1 - u_kj a_k(p))` for each of d latent factors, then `P(s|p) = 1 - (1 - leak_s) prod_j (1 - v_js z_j(p))`, replacing the K x S link matrix with U (K x d) and v (d x S): 128 parameters per relation at d = 4 against 192, 96 at d = 3. Both matrices freeze for the symptom folds of item 26 d1, which is what pins the module basis. Preferred over factorising the link logit, which imposes a rank in a space the noisy-OR gives no meaning to; here every stage stays "any one cause suffices" and a latent factor is a mechanism group, which is the nosology object | three costs, measured or structural: the links already travel 90 percent of what Adam can move them in a run and still only reach a probability of 0.09 to 0.15 (docs/b6_module_diagnosis.md), so another multiplicative layer between the modules and the loss needs its own time scale or the modules move less still; a latent factor can die as a module can; and d is a free number, which the owner has objected to before. Decided before the first confirmatory run or not at all |
| 4.8 | **The same rank reduction with no extra layer: factor the link matrix itself as a noisy-OR.** `link_ks = 1 - prod_j (1 - u_kj v_js)`, which is in [0, 1] by construction with no clipping, gives the same U (K x d) and v (d x S) and the same frozen shared basis for item 26 d1, and keeps the head one layer deep. The difference from 4.7 is where the latent factor acts: in 4.7 it pools activations, so two modules of one factor saturate against each other and the symptom does not count them twice, which is the physiology; in 4.8 it only structures the parameters and each module reaches the symptom on its own | 4.7 is the better model and 4.8 the safer fit. How much the saturation is worth cannot be read off the runs in hand: at the measured link probabilities of 0.09 to 0.15 the noisy-OR is near its linear regime and the two forms nearly agree, but those links were travel-limited rather than converged (docs/b6_module_diagnosis.md), so that near-agreement is evidence about the optimiser and not about the model |
| 4.9 | **No refactor at all: freeze the link matrix the head already has.** What item 26 d1 needs from the head is one thing, a part of the output that is fixed across folds so module identity is pinned. The existing K x S link matrix does that: fit it on all 24 symptoms, freeze it with the leaks, and each fold refits only the encoder, the supports and the readout weights on the other 23. Module k's contribution runs through a fixed row, so it cannot permute, and the symptom profile stays in the basis of named modules with nothing anonymous between the pathways and the symptoms | what it gives up: with the readout frozen the symptom profiles are the same in every fold, so a symptom fold measures module stability and not profile stability. That is not a gap, because the two folds then divide the work - the symptom fold answers whether a module is just one symptom wearing a module's name, which is the circularity the nosology readout has to rule out, and the perturbation resampling of 26 d2 answers whether the grouping reproduces. 4.7 and 4.8 exist to make profile stability measurable under symptom folds too, by splitting the output into a frozen part and a free part; that is the only thing they buy, and it is what the latent layer costs readability for |
| 4.5 | **The probability clamp** at [1e-6, 1-1e-6] in the training loss: a positive beyond it has no gradient, and in the fold-0 pilots 3.7 percent (message passing) and 1.2 percent (linear response) of sigmoid-head test probabilities sat below 1e-6 | open in the preregistration |
| 4.6 | **A confirmatory baseline per perturbation type.** None is fitted; the development base rate is 0.159 for drugs against 0.063 for genes, so a predictor that knows only whether a perturbation is a drug earns a share of every reading | open in the preregistration, with the pseudo-lockbox measurement done |

## 5. The pathway test

| | decision | state |
| --- | --- | --- |
| 5.1 | **Multiplicity against H1 and H2.** Three ways are set out; the first two keep the stated bound and the third raises the chance of a false claim to 0.10 | the owner's, and it must be fixed before any confirmatory module exists |
| 5.2 | **Whether P1 is primary at all**, which needs the power simulation first: at a plausible effect with supports of 1 to 20 genes, a P1 failure may say nothing | power statement outstanding |

## 6. Actions only the owner can take

| | action | state |
| --- | --- | --- |
| 6.1 | **The CTD notice.** CTD's terms require notifying them of use and giving periodic access to publications. `docs/ctd_notice_draft.md` is written for the owner to send, with two of its own pages for the owner to read first | not sent |
| 6.2 | **The Lombardo permission request** for the clearance table, written and given to the owner on 12 October and deliberately not in the repository. A refusal costs less than it first looked | **sent, awaiting a reply** (owner: "Requested. Waiting") |
| ~~6.3~~ | **Closed 12 October: the owner set it aside** ("Don't worry about it"). Recorded and not raised again |
| 6.4 | **Hodgkin & Katz 1949** could not be read: the PMC copy is a page scan with no text layer and its PDF is behind a proof-of-work gate, which under the ruling of 11 October is not solved or routed around. Only a human reader can check the squid permeability ratio against the original | the address is in `docs/membrane_potential_constants.md` |
| 6.5 | **Raimondo 2015 Table 1's footnote markers** are in neither the Europe PMC nor the publisher HTML, so which cells the caption calls estimated is unknown from what was read | the typeset PDF would settle it |

## 7. Long-standing, from the preregistration trailer

Phase 1 counts per symptom and grade; the final symptom set after go/no-go; the B5 language model and prompt; the power
statement for the GWAS enrichment test; and the open questions 8 to 10 of design section 11 — frequency as a weight or
as a target, the subsystem against curated pathway split, and the deterioration-term exclusion.

## 8. Code registered and not yet written

Not decisions, but they gate 1.1 as surely as a decision would: the fixed-point propagation loop, the one repeated
layer with rank-factored relations and an explicit norm bound, and the signed adjacency for message passing are all
registered in the amendment of 12 October and none is implemented. Nothing has been run on them.
