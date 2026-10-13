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
| 1.4 | **The Open Targets class: settled on 13 October except the weight.** The crosswalk wordings "weight gain" and "sexual dysfunction" are out; the inhibition statements are dropped because the gene knockout already represents that perturbation and is scored, and the labels are weak evidence that AUPRC would score at full strength; the activation statements become 23 `receptor_activation` rows that are never scored | what remains is the weight. 0.45 is provisional, the agreement it is derived from recomputes to 0.6418 on the capped table against the registered 0.7606, and if the rebuild leaves it near 0.64 the weight becomes 0.35. It now prices 35 training labels and nothing scored, so the choice moves little |
| 1.3 | **Off-targets in a model input.** The owner asked for them on 12 October; my objection is recorded and stands: admitting them at any threshold from 0.1 to 0.5, under any concentration scheme, takes the largest leakage group from 71 drugs to between 203 and 209 of 222, which ends the drug split item 2 of the amendment depends on | the split, so 1.1 |

## 2. What trains

| | decision | state |
| --- | --- | --- |
| ~~2.1~~ | **Closed 12 October: they train, in both encoders** (owner: "Implement"; amendment item 17). Symmetric, because the asymmetry I warned of rested on `clears` being the only route and that was wrong — `catalyzed_by` carries them too |
| ~~2.2~~ | **Closed 12 October: they stay out** (owner: "That's fine") |
| ~~2.3~~ | **Closed 12 October: accepted** (owner: "Weight accepted") |

| ~~2.4~~ | **Closed 13 October: the class enters the reliability fit** (the owner: "add it to the Dawid-Skene fit estimates"). `reliability_reports` emits a report only where another source already reports that item, so it creates no pair, reaches the fit and not the grade, and is estimated exactly where the hand-computed 0.761 was measured. 0.45 stays as the fixed-grades fallback and is labelled a convention | closed |

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
| 4.3 | **Whether the cell-class conditioning should reach the propagation.** Restated 13 October: under the registered `limiting_step` default the weights are out of the loop entirely, so they do not gate the input at the seeds either - the class's expression re-enters only as the post-hoc ceiling. The conclusion stands and the mechanism in the old wording was the product rule's. 34.0 percent of nodes are flat at 1 in every class, and under a max-min rule a route through a flat node costs nothing | unquantified, and `experiments/measure_limiting_step_rule.py` item 4 was written to measure it. Three documents still describe the product rule as current: this entry until today, `docs/label_source_review.md` 28.4 and the `brain_expression_weights.py` docstring |
| 4.4 | **The weighted constant start** (`weighted_constant_optimum`): whether the confirmatory_v2 configurations take it | no configuration carries it |
| ~~4.7~~ | **Set aside 13 October: the two-level noisy-OR has no job left.** `z_j(p) = 1 - prod_k (1 - u_kj a_k(p))` then `P(s|p) = 1 - (1 - leak_s) prod_j (1 - v_js z_j(p))`, the owner's form ("what about just putting the noisy or in between the modules and the lower rank matrix?"). It was the better of the two latent forms, because the factor pools activations and two modules of one factor saturate rather than counting twice. What it bought item 26 was profile stability under symptom folds, and the owner does not need that ("i guess i don't need profile stability under symptom folds so much as symptom profile stability in general"), which 26 d2 measures without it. Kept here with its costs - a second multiplicative layer between the modules and the loss in a head whose links already travel 90 percent of what Adam can move them, dead factors, and a free `d` - in case shrinkage for thin symptoms is later wanted on its own | set aside, not refused |
| ~~4.8~~ | **Set aside 13 October with 4.7.** `link_ks = 1 - prod_j (1 - u_kj v_js)`: the same rank reduction with no new layer, safer to fit and weaker as a model, since the factor structures the parameters only and every module reaches the symptom on its own. Falls with 4.7, for the same reason | set aside, not refused |
| ~~4.9~~ | **Closed 13 October: nothing to freeze.** Freezing the link matrix existed to pin module identity across symptom folds, and the owner dropped the folds ("i don't know that the joint reading refit is worth it, tbh"). Item 26 now asks for no change to the registered head at all: no latent layer, no frozen readout, no new parameter, no free rank. What survives of the question the folds were for is the link concentration `max_s link_ks / sum_s link_ks`, reported per module against the permuted-label null, which costs nothing beyond runs already registered | closed |
| 4.5 | **The probability clamp** at [1e-6, 1-1e-6] in the training loss: a positive beyond it has no gradient, and in the fold-0 pilots 3.7 percent (message passing) and 1.2 percent (linear response) of sigmoid-head test probabilities sat below 1e-6 | open in the preregistration |
| 4.6 | **A confirmatory baseline per perturbation type.** None is fitted; the development base rate is 0.159 for drugs against 0.063 for genes, so a predictor that knows only whether a perturbation is a drug earns a share of every reading | open in the preregistration, with the pseudo-lockbox measurement done |

| ~~4.10~~ | **Closed 13 October: a perturbation can now drive a reversible reaction backwards** (the owner: "I don't want to accept a perturbation not being able to reverse the reaction"). Amendment item 27: linearised mass action, dv = k_f dS - k_r dP, puts the product into the rate at -1, so the encoder derives product -> reaction at -1 (`product_opposes`) and the build writes the nominal direction only. Reversal turned out to be linear; I had framed it as switching edges rather than as the sign of the rate. Left open inside it: the resting direction from a thermodynamic or flux source rather than Human-GEM's nominal one, the relation in the message-passing encoder, and 16 Reactome pairs where a consumed-and-regenerated species (scaffolds, CoQ) needs a modifier edge | closed; reaches a run at the rebuild |
| 4.11 | **The linear-response encoder discards `edge_weight`**, which it reads only on drug entry edges, so item 13's six fixed conductance weights are inert there and K+, Na+ and Cl- reach `VM_c` with equal magnitude - the "weight of 1 on every edge" the membrane-potential module says it is replacing. The message-passing encoder does apply them | the two registered encoders disagree about the physics. Found 13 October; undisclosed before |

| 4.12 | **The cell-class bottleneck is computed on the unsigned structure** (the owner, 13 October: "Whats computed on the unsigned structure? It sounds like it probably shouldn't be"). `b(t) = max over routes of min over nodes` walks the edge structure with the comment "a bottleneck cares only whether a route exists, so signs, weights and the in-degree divisor play no part in it". For **existence** that holds - expression gates whether a cell can run a step whether the step activates or inhibits - but `b` then scales a **signed** sum, so the ceiling can come from a well-expressed inhibitory route while the response travels an activating one, and those are different routes with different expression | excluding the in-degree divisor is right, since it is a modelling artifact and not physiology. Excluding the sign is the questionable part and the fix is available where `b` already lives: the bottleneck is applied after the fixed point, so the settled sign at the target is known, and the max-min can be taken over the routes whose sign product matches it. No nonlinearity enters the loop. Accepted by the owner 13 October ("Sounds good"); unbuilt |

| 4.13 | **The aggregation across relations is a per-node mean, which is not physiology and the code says so** (the owner, 13 October: "setting mean as the aggregator doesn't sound very physiological though. Can you justify it?"). It cannot be justified as physiology. Mass balance at a metabolite is a signed **sum**; the per-node divisor is there so the row sums are at most one and the iteration contracts, which is a numerical requirement and not a biological claim. The encoder states the consequence itself: the divisor "is what discards the stoichiometric counts that decide direction at a metabolite", and 92 percent of metabolites end up with an aggregate row sum of zero | **the remedy is already written and is not the default.** `total_in_degree` divides by the destination's in-degree summed over every relation, so the unsigned row sums are one exactly as under `in_degree` - the same contraction, no global shrink - while the counts survive: five producers against one consumer give (5-1)/6 rather than 0, and the sign at a metabolite is set by its own stoichiometry rather than by the one gain ratio every metabolite shares. The `spectral` arm also keeps the counts but shrinks the field nearly everywhere, which makes it two changes at once. Given that the in-degree mean underlies most of the costed deviations in the adherence review, this is the single highest-value change available. The owner's, and it belongs with the rebuild |

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
