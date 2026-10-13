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
| 1.3 | **Off-targets in a model input.** The owner asked for them on 12 October; my objection is recorded and stands: admitting them at any threshold from 0.1 to 0.5, under any concentration scheme, takes the largest leakage group from 71 drugs to between 203 and 209 of 222, which ends the drug split item 2 of the amendment depends on | the split, so 1.1 |

## 2. What trains

| | decision | state |
| --- | --- | --- |
| 2.1 | **The enzyme, carrier and two-drug statements as training labels.** The combined perturbation type they needed is built; 106 combinations and 317 positive pairs are expressible. Training them in the linear-response encoder but not in message passing makes the two confirmatory models differ in what they are trained on, which the owner has been avoiding | `--combined-perturbations` has no default, so no run loads them |
| 2.2 | **The 21 negative controls of the check table.** They are not written as labels at all. A labelled negative would assert the symptom does not occur, which no source states | recorded in the builder; not proposed |
| 2.3 | **The Open Targets class at 0.45 with a 0.0 secondary arm.** Registered as item 12 on the owner's "no preference", so this is open only in the sense that it may be overridden | registered |

## 3. The drug set

| | decision | state |
| --- | --- | --- |
| 3.1 | **The four held drugs of amendment item 14.** Pindolol (penetrant at 0.50, no label sentence of its own), oxprenolol (penetrant at 0.20, no US label in openFDA), timolol (a central-action sentence, no measured Kp,uu,brain), carvedilol (neither). Each needs one part supplied, or to stay out | the rule is registered; these four are its edge |
| 3.2 | **Whether to read label sentences by hand for the positives the screen shortlists.** The screen's negatives are evidence; its positives include class boilerplate, so a reading decides. Five or six drugs | 3.1 |
| 3.3 | **Whether to seek a Kp,uu,brain for the 492 candidates that have a brain-expressed target and no measurement.** The feasibility pass says there is no accurate route at scale | measured and reported; no recommendation to pursue |

## 4. The model and the encoder

| | decision | state |
| --- | --- | --- |
| 4.1 | **The reach of the electrical layer.** Item 13 registers the computed potential; the magnitude arriving at `VM_c` is set by the ion pool's in-degree, and both normalisation alternatives the encoder carries were measured worse. So it is disclosed rather than fixed, unless the owner wants a new aggregation | disclosed |
| 4.2 | **The cell-class weight of a protein complex.** 413 reactions regulated only by a complex carry a flat class weight, adenylyl cyclase among them, because the weight reaches a reaction through its catalysing gene. 1,527 of the 1,681 complexes have a member that already carries one, so a complex could take its limiting member's weight | measured; a build, not yet proposed as a change |
| 4.3 | **Whether the cell-class conditioning should reach the propagation.** Today the weights gate the input at the seeds only, so each class's channel traverses the same structure and differs only in its gains | a registered-encoder change |
| 4.4 | **The weighted constant start** (`weighted_constant_optimum`): whether the confirmatory_v2 configurations take it | no configuration carries it |
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
| 6.2 | **The Lombardo permission request** for the clearance table. Written and given to the owner on 12 October, deliberately not in the repository. A refusal costs less than it first looked | not sent |
| 6.3 | **The GitHub exposure.** An e-mail address was committed to the public repository and removed, and the history was rewritten, but GitHub still serves the orphaned commit and its public events feed exposes the SHA. A support request was drafted in conversation | not sent, and outside the repository by the owner's instruction |
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
