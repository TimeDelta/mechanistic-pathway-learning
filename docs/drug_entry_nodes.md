# Drugs as graph nodes: what the record says, and what adding them would take

Written 9 October 2026, at the user's question: "i thought we had decided to include drugs as nodes that aren't
connected to a gene (essentially new graph entry points), which should allow for more carriage edges".

## What the record holds

The intention is in the design document, once, as a forward reference. `docs/experiment_design.md`, in the
neuronal-variant paragraph, says of the reverse-transport relation that "the distinction bites once a drug that
reverses the transporter is representable, which needs the exogenous drug nodes of section 5.7". **Section 5.7 is the
ablation list and names no drug nodes**, so the reference has no referent: no specification was written, and nothing
was built. The graphs confirm it. `data/processed/graph_full_neuronal_split_binders`, the graph the confirmatory
family reads, holds six node types and no drug among them:

| node type | nodes |
|---|---|
| reaction | 13,793 |
| gene | 12,810 |
| protein | 12,709 |
| metabolite | 8,580 |
| protein_entity | 1,681 |
| membrane_potential | 1 |

A drug reaches the graph only through its target's nodes, or, for the one drug that is an endogenous metabolite, as
that metabolite (`configs/drugs_acting_as_graph_compounds.csv`, tryptophan). So the decision the question remembers
was never taken in a form anything acts on, and the limitation `docs/plasma_binder_graph.md` records stands as
written: "placing them needs drug nodes, a change to where a drug perturbation is seeded and a signed sequestration
relation".

## What it would buy

Carriage edges. The 27 carriage edges in the confirmatory graph run between an endogenous cargo and its binder's gene,
because a xenobiotic has no node to attach to. With drug nodes, every drug whose carrier is known gets one. Two
sources now name a carrier: the pinned fraction-unbound database for 27 of the 154 drugs
(`docs/fraction_unbound_coverage.md`) and the FDA labels for 31, 22 of them drugs the database does not hold
(`docs/label_plasma_binders.md`), so **49 of 154 drug perturbations would get a carriage edge** and 105 would not.

It would also make representable two things the graph currently cannot state: a drug that reverses a transporter
rather than blocking it, which is the amphetamine mechanism the neuronal variant's `catalyzes_reverse_transport`
relation was added for, and competition between two drugs for one binder.

## The three changes it needs, and what each costs

1. **A drug node and its edges.** One node per drug perturbation, with its mechanism edges to the target's nodes
   carrying the sign the mechanism table already assigns, and a `binds` edge to its carrier's gene or protein node.
2. **The seed moves.** A drug perturbation is currently seeded on its target's nodes with magnitude spread over them.
   With a drug node, it is seeded on that node and the mechanism edges carry it to the targets. This changes how every
   drug perturbation enters, not only the 49 with a carrier, and it costs one propagation hop: the message-passing
   configuration runs three layers, so a seed one hop further from the targets reaches a third less far through the
   metabolic layer. The linear-response encoder, which iterates one transition to a fixed point, pays an
   in-degree-averaged attenuation instead of a hard cut.
3. **A signed sequestration relation.** `binds` carries sign 0 in all 31,177 of its rows, which reads as moving a
   substrate without acting on it, so a carriage edge as it stands does not make bound drug unavailable. Representing
   that needs a signed depletion edge from binder to cargo, as `depletes_substrate` was added for stoichiometry. This
   applies to the 27 existing carriage edges too: albumin binding bilirubin does not lower free bilirubin in this
   model either.

## Two objections to weigh before deciding

- **A drug node is a place to memorise.** The design rule is that parameters are "never per node, since a per-node
  parameter is memorisation that cannot transfer to a held-out gene". A drug node under identity embeddings gives
  each drug its own vector, and a held-out drug's vector is then untrained, so the identity-embedding configurations
  would have to leave drug nodes featureless or structural. The inductive variant, whose node features are
  structural, has no such problem.
- **It changes the confirmatory graph and the perturbation model.** Both are preregistered, so this needs a dated
  amendment and the development runs repeated, which is the same cost as adopting v3. It is worth doing with that
  move rather than separately.

## The decision, 10 October 2026

The user: "sorry that i was unclear then but i definitely want drugs as nodes with an ablation on them". So drug
entry nodes are adopted with an ablation arm that removes them, and the two objections above are costs to carry
rather than reasons not to. Nothing is built yet; what follows is the specification that session should implement, and
the three choices in it were stated to the user and not contradicted.

1. **A drug node per drug perturbation, seeded on itself.** The mechanism edges carry the perturbation to the target
   genes, so the seed moves off the targets and onto the drug. This changes how all 154 drug perturbations enter, and
   costs the propagation hop item 2 above measures.
2. **A drug that is itself a graph compound attaches to that compound rather than duplicating it.** The user:
   "it just has to interact with the endogenous and exogenous drug nodes consistently". Today no drug perturbation is
   seeded on a metabolite (all 154 land on gene nodes), and `configs/drugs_acting_as_graph_compounds.csv` holds one
   row, tryptophan, so the case is latent rather than live; the rule has to exist before it is, or the same molecule
   becomes two nodes with separate states.
3. **A signed sequestration relation, not `binds`.** Item 3 above gives the reason: all 31,177 `binds` rows carry
   sign 0, so reusing it would place a carriage edge that does not lower free drug.

The ablation is the arm that drops the drug nodes and seeds on the targets as now, which is what makes the entry
point a measured claim rather than a modelling preference. It belongs in the ablation list of
`docs/experiment_design.md` section 5.7, which is where the original forward reference pointed and found nothing.

Two things this needs that are not code. A dated amendment to `docs/preregistration.md`, because the confirmatory
graph and the perturbation model are both preregistered; the amendment should be drafted for the user rather than
applied. And the development runs repeated, which is the same cost as the move to v3, so the two should move
together rather than separately.

## What was built and measured, 10 October 2026 (the session after the decision)

Nothing below was trained and no lockbox file was read. The measurements are in `docs/drug_entry_nodes_measured.md`,
which `experiments/measure_drug_entry_nodes.py` generates; the numbers quoted here are copied from it.

### A correction to "What the record holds"

The section above says the only mention of drug nodes in the design document is a forward reference. That is wrong.
`docs/experiment_design.md`, section 4.1, lists the node types as "metabolite (compartment-specific), reaction, gene,
protein (enzyme, transporter, receptor, transcription factor) and drug" and the edge types as ending in "targets (drug
to protein, with action type and affinity from ChEMBL [16])". The relation list of the release graphs
(`data/releases/v0.4`) and of every variant built from them holds `targets` at index 7, with no edge of that relation.
So the design named the node type and the relation, and no graph held a node or an edge of either. The build below uses that reserved relation for the mechanism edges and adds nothing to the
relation list except `sequesters`.

### The build

| piece | where |
| --- | --- |
| node and edge construction | `mechanistic_pathway_learning/graph/drug_entry_nodes.py` |
| graph variant | `experiments/build_drug_entry_node_variant.py`, written to `<graph>_drugs` |
| loader, both arms | `load_experiment_data(..., drug_entry="targets" or "nodes")` |
| encoders | `entry_node_mask` and the entry edges in both encoders |
| trainer | `--drug-entry {targets,nodes}`, `--drug-mechanism-before-first-layer` and `--sequestration-carries {change,presence}` in `experiments/run_main_model.py` |
| carriers | `configs/drug_plasma_carriers.csv`, written by `experiments/scope_label_plasma_binder_coverage.py --carrier-table` |
| tests | `tests/test_drug_entry_nodes.py` (18 tests) |

On `graph_full_neuronal_split_binders` with the 154 drugs of `evidence_full_v3_parkinsonism` the variant adds 154 nodes
of type `drug` (`DRUG:<perturbation id>`), 791 `targets` edges (338 of sign -1, 19 of -0.5, 46 of 0, 309 of 0.5, 79 of
1) and 71 `sequesters` edges of sign -1 from ALB, ORM1 or ORM2 to the 49 drugs with a known carrier: 49,728 nodes and
280,646 edges. Rows are appended and no existing row changes, so the source graph is the variant's prefix. No drug is a
graph compound in this table, so the rule of decision 2 is implemented and tested on a toy graph and has no live case.

The ablation is one flag on one graph. `--drug-entry targets` drops the drug nodes, their edges and the carriage edges
that were added with them and seeds each drug on its targets: 0 of the 40 fields the loader returns differ from the
source graph's. `--drug-entry nodes` seeds each drug on its own node with sign 1 and magnitude 1, and the mechanism
edges carry the sign and the magnitude the seed carried.

### Choices the specification left open, as built

Each is a choice made to get a working arm, and each can be reversed without touching the others.

1. **A drug node exists only where it is seeded.** The user, shown this rule on 10 October: "The presence rule
   implementation is correctly inferred." In every other perturbation its state is held at zero before the
   first layer and after every layer, and the reference pass of the difference field holds no drug node. The reason
   is measured: with the drug nodes always present and averaged like any neighbour, 175 of 240 gene knockouts change
   under message passing and 143 of 240 under linear response, and all 46 of the 46 knockouts at a drug target or a
   carrier change (median relative change 0.43 and 0.056). Under the rule 0 of 240 change, with a largest difference
   of exactly 0. The precedent is the multiscale interactome: "We then make the drug or disease of interest a source
   node (i.e., no in-edges) and all other drugs and diseases sink nodes (i.e., no out-edges)", because "drugs and
   diseases do not propagate their effect by using other drugs and diseases as intermediates" (Ruiz, Zitnik and
   Leskovec 2021, [doi:10.1038/s41467-021-21770-8](https://doi.org/10.1038/s41467-021-21770-8)).
2. **Entry edges are summed, with no in-degree divisor.** 82 of the 129 target nodes are shared by two or more study
   drugs (744 of the 791 edges), and three GABA-A subunits are each the target of 26. Under the mean aggregation of
   the other relations, a drug's input to its target would be divided by the number of other study drugs that share
   the target, which is a property of the drug table and not of the body.
3. **The mechanism edge carries the seed's sign and magnitude**, as an edge weight the encoders read. In message
   passing the message is `weight * (x W_unsigned) + sign * weight * (x W_signed)`, two matrices per layer that exist
   only when the graph has entry nodes, so an encoder without them has the parameters it always had.
4. **Leakage groups and degree strata come from the targets in both arms.** A drug node's own degree is its number
   of targets and carriers, which would put every drug in the lowest degree stratum; the strata and the groups are
   therefore computed from the seeds of the `targets` arm, and the perturbation ids, groups, outcomes, weights, label
   mask and strata are identical in the two arms.
5. **The rewiring control of H2 holds `targets` and `sequesters` fixed.** Rewiring a drug's mechanism edges would
   give the drug other targets, which is a different control from rewiring the body's edges.
6. **`seed_masked` hides the descriptors of a drug's targets in both arms.** The treatment hides the descriptors of
   each perturbation's perturbed nodes and their `encodes` partners. With the drug seeded on its own node, which has
   no descriptor, the mask would have covered nothing; it now covers the drug node, the targets of its mechanism
   edges and their `encodes` partners, one way, so a knockout of a target does not reach back to the drugs. This is
   the treatment the descriptor rule picks for message passing, so it had to be defined before that arm could run.
7. **A drug node has no descriptor and one more structural column.** The descriptor tables get a zero row per drug
   node, and the one-hot node type gains a `drug` column in the `nodes` arm only.
8. **`--normalise-drug-input` is refused with drug nodes**, because the seed on a drug node is one node of magnitude 1.

### What the measurements say

**The hop costs much more than item 2 of "The three changes" estimated.** That item says a seed one hop further
from the targets "reaches a third less far". At the registered three layers a drug seeded on its targets has a nonzero
field at a median of 13,823 nodes and a drug seeded on its node at 1,532, a share of 0.11, and all 154 drugs reach
fewer nodes. The estimate assumed reach grows evenly with depth; the third hop from a receptor is where the signalling
layer opens into the metabolic one. Two remedies were measured: delivering the mechanism once before the first layer
(`--drug-mechanism-before-first-layer`) gives a median of 13,823 nodes, with 150 of 154 drugs reaching at least as
many nodes as when seeded on their targets, and a fourth layer gives 13,839.

**The same section shows something about the registered design that does not depend on drug nodes.** A gene
knockout is seeded on a gene node, one `encodes` edge before its protein, and a drug is seeded on protein nodes, so
at three layers a drug is one edge ahead. Counted by breadth-first search over the edges (the encoder counts above
differ from these by a few nodes), the median gene knockout reaches 175 nodes within three hops (quartiles 13 to
1,033) and the median drug 13,822 (9,649 to 14,246). Decision 9 of `docs/gene_protein_split.md` added the third layer
for this reason and the gap is still a factor of 79 at the median. A drug on its own node sits where the gene
knockout sits, one edge before the protein, and reaches 1,523. So "the hop a drug node costs" can also be read as the
head start a drug has today.

**In the linear-response encoder the two arms are one model.** The response of a drug seeded on its node equals its
response seeded on its targets multiplied by one gain per channel: median cosine similarity 0.999886 at 8 steps and
1.000001 at 32, where the largest relative difference is 1.2e-7. The remaining difference at 8 steps is the one step
the node costs. The ablation therefore tests something only in the message-passing encoder; in the linear-response
encoder it can show a difference only through that gain and the missing step. Four drugs have a zero response in
both arms because every sign of their mechanism is 0: gabapentin, pregabalin, brivaracetam and pentazocine, the
zero-sign drugs the preregistration and `docs/drug_targets_any_type.md` already list.

**The sequestration edges do not do what decision 3 wanted.** In the linear-response encoder they change 21 drugs'
fields by at most 8.2e-7: the response is linear in the input, so in a drug's own perturbation a carrier holds only
whatever part of the drug's signal reaches it, and that is what the edge sends back. In message passing they
change all 49 drugs that have one, by a median relative 2.2, and none of the 105 that do not: the carrier's base
state arrives at the drug node whatever the drug is, so the edge is a constant input that says "this drug has a
known carrier". That flag is not neutral. A carrier is known for the drugs whose labels or measurements discuss
binding, and those drugs have more kept positives: 6.41 a drug against 4.19 (median 6 against 2; Mann-Whitney
two-sided p = 0.002). A model can use the flag to raise every prediction for the better-studied drugs, which is the
degree shortcut in another form ("performance attributable to factors other than degree is often only a small portion
of overall performance", Zietz et al. 2024, [doi:10.1093/gigascience/giae001](https://doi.org/10.1093/gigascience/giae001)).

**An edge that carries the carrier's change removes the flag and keeps the edge.** The flag arises because message
passing sends the carrier's resting state into the drug node. The user asked why the carrier is not varied in any
perturbation: "The whole point of including exogenous drug nodes was to gain the carrier edges for them." The answer
is that no label varies it (every example is one gene or one drug with its symptoms), and that the edge can be made to
act only when the carrier moves. `--sequestration-carries change`, now the trainer's default on a graph with the
relation, sends along the edge how far the carrier has moved from its unperturbed state; `presence` keeps the reading
measured above as an arm. With `change`, 0 of the 49 drugs with a carrier and 0 of the 105 without have a field that
differs from the graph without the edges (largest difference 3.6e-7). The linear-response encoder is a deviation from
rest already, so it needed no change. Simulated there with nothing trained, lowering a drug's carriers by one unit
beside the drug raises the drug's effect in every channel for all 49 drugs: the part of the response that passes
through the drug node is the drug's own response times 0.59 at 8 steps (0.38 to 0.64 over the drugs) and 0.30 at 32
(0.27 to 0.76), with a cosine of 0.9999 or more against the drug's own response. The sign of that factor is the
graph's; its size is the untrained gains' and means nothing yet. In message passing the same simulation would pass
through a weight matrix that no training example reaches, so its output is not reported.

**Message passing at three layers sees about 1 percent of what a gene knockout can reach.** The user asked why message
passing does not carry the signal until every reachable node is reached. By breadth-first search over every
perturbation, 39,801 of the graph's 49,574 nodes can be reached at all; a gene knockout needs a median of 18 hops to
reach all of its reachable nodes (26 at most) and 8 to reach 95 percent, and a drug seeded on its targets 16 (20 at
most) and 7. At the registered 3 layers the median gene knockout has reached 0.9 percent of its reachable set and the
median drug 35 percent; at the 8 steps of the linear-response default, 96 and 98 percent. The fixed depth is the
design's choice (section 5.2: "the difference field is exactly zero more than L edges from the perturbed nodes"), and
no record shows the depth was ever varied except by the one layer the gene and protein split added. The two
registered encoders therefore differ in reach, in form and in how sign enters, all at once.
`docs/message_passing_depth_probe.md` gives what reusing a layer for more rounds does before anything is built. One
untrained shared layer at the registered initialisation keeps the field's size from 3 to 8 rounds and then grows it
about a thousandfold by 26 (median norm 2.3, 1.9, 35 and 2,060 for gene knockouts at 3, 8, 18 and 26 rounds), and
with the weights halved shrinks it at every step, so a reused layer needs a damped, averaged or gated update before
perturbations can stop at different rounds. A round costs the same time and memory whichever weights it reads: 0.19 GB
a round at batch 16. And the linear response at its 8 steps is within 0.14 percent of its fixed point at the median
and 2.3 percent at most, so that encoder is already the converged case.

### What the literature says about the sign of a carrier edge

The edge states that more carrier means less free drug and so less effect. The pharmacokinetic literature supports
that for the free fraction and, for most drugs, not for the free concentration at steady state, which is what an
effect follows. The quotes below were read in the full text named; hyphens that the retrieval stripped are restored,
subscripts are written on the line and nothing else is changed.

- Oral dosing, hepatic elimination: "the area under the curve of unbound drug (AUC oral,u) is independent of protein
  binding", so "no adjustment of drug dosing should be required for any drug dosed orally for which the liver is the
  major organ of elimination" (Benet 2018, full text, [doi:10.1002/jcph.1105](https://doi.org/10.1002/jcph.1105)).
- Where binding does matter: "protein binding changes will only be relevant for high-extraction-ratio drugs dosed
  intravenously, and I estimate that fewer than 3% of drug dosings would require concerns about changes in response
  due to changes in binding caused by drug-drug interactions or disease states. And these drug dosings will almost
  exclusively be anesthetics and IV dosing of antiarrhythmic drugs, such as lidocaine" (the same paper; its Table 1
  reads "Orally dosed; hepatic elimination: No", "Orally dosed; renal elimination: Possible, but no known examples",
  "IV dosed; low ER: No", "IV dosed; high ER: Yes").
- The same table by another group, with one cell Benet's lacks: for the unbound AUC, intravenous dosing changes it
  at a high extraction ratio and not at a low one whatever the clearance route, and oral dosing changes it only at a
  high extraction ratio with nonhepatic clearance (Hochman et al. 2015, Table 2, full text,
  [doi:10.1002/jps.24306](https://doi.org/10.1002/jps.24306)). Peak and trough unbound concentrations change in more
  cells than the average does.
- For the brain: the unbound partition coefficient "is dependent on passive diffusion, uptake and efflux transport
  Cl, brain interstitial bulk flow, and brain metabolic Cl of a drug [...] and is independent of protein binding in
  blood and in brain" (Bohnert and Gan 2013, full text, [doi:10.1002/jps.23614](https://doi.org/10.1002/jps.23614)).

So the gate is the route and the extraction ratio, not oral against everything else: an oral drug cleared by the
liver is outside it at any extraction ratio, an intravenous drug of low extraction is outside it and an intravenous
drug of high extraction is inside it. Reading the non-intravenous routes that bypass the first pass (intramuscular,
transdermal, sublingual) as the intravenous row is my inference from the same equations; neither table has a row
for them.

`docs/drug_administration_routes.md` (generated by `experiments/tabulate_drug_routes.py` from openFDA's label routes)
gives the route half for the study's drugs. Of the 49 with a known carrier, 18 have a label that names the
intravenous route, 2 have another route past the first pass and no intravenous one, 26 are oral or local only and 3
have no label route. So at most 20 of the 49 edges could mean what they say even before the extraction ratio is
applied, and the adverse-effect reports do not say which product a patient took: caffeine and diphenhydramine count
as intravenous here because an injection exists. **The extraction-ratio half is still owed.** It needs a clearance per
drug, which no source in this repository holds; the candidate is the intravenous dataset of Lombardo, Berellini and
Obach 2018 ([doi:10.1124/dmd.118.082966](https://doi.org/10.1124/dmd.118.082966); 1,352 drugs), which
`docs/plasma_protein_binding.md` already names as the next source of the fraction unbound. Its supplement has not been
fetched or checked against the study's drugs.

What the literature does support is the other direction, for the body's own cargo: a drug that occupies a carrier
frees what the carrier held. That is drug to carrier to cargo, on the 27 carriage edges the graph already has, and
it needs a molar concentration on the drug and on the cargo that this model does not carry: a drug displaces cargo
only when its plasma concentration is of the order of the carrier's binding capacity. `docs/off_target_scoping.md`
names the source such a concentration would come from (therapeutic total plasma concentration times the fraction
unbound) and records that the route is not built.

### What the literature says about the design itself

- **No like-for-like comparison of a drug node against a seeded target set was found.** The nearest are benchmarks
  in which the representation changes together with other things: the multiscale interactome against target-overlap
  and proximity baselines, where the model also adds biological-function nodes (Ruiz et al. 2021). The ablation here,
  with the graph, the labels and the learner fixed, has no published effect size to expect. The search for one was a
  keyword scan of the papers read and two web searches, not a citation-graph search, so a comparison in a supplement
  could have been missed.
- **Depth.** "To allow a node to receive information from other nodes at a radius of K, the GNN needs to have at least
  K layers, or otherwise, it will suffer from under-reaching" and, with more layers, "information from the
  exponentially-growing receptive field is compressed into fixed-length node vectors" (Alon and Yahav 2021,
  [arXiv:2006.05205](https://arxiv.org/abs/2006.05205)). A layer added for the drugs is added for the gene knockouts
  too.

### Decisions that are the user's, with a recommendation for each

**A. Keep the sequestration edges in the graph, carrying the carrier's change** (`--sequestration-carries change`),
and run `presence` as an arm that shows what the flag alone is worth. The first recommendation written here was to
leave the edges out (`--carrier-table none`), because they were inert in one encoder and a label-density flag in the
other. The user's reply was that gaining the carrier edges was the point of the drug nodes, and the change-carrying
edge answers both: it is silent when a drug is taken alone, so the ablation's two arms differ only in how the drug
enters, and it acts when the carrier is perturbed with the drug. What it does not give is a trained or scored effect
of the carrier, because no label varies one. Three things would: simulated perturbations of a drug with its carrier
lowered or raised, read in the linear-response encoder only; the free fractions the pinned database holds for 27
study drugs in neonates, the elderly, hepatic and renal impairment and inflammatory disease, which can check the
direction of the drug node's response and not a symptom; and an edge weight equal to the bound fraction, since for
one binding site the free fraction moves with the carrier in proportion to the fraction bound. None of the three is
built. The pharmacokinetic limit above still holds: the simulated direction is a statement about the free fraction
and the peak free concentration, and for the average free exposure it holds for the intravenous high-extraction
drugs only.

**B. Deliver the mechanism before the first layer in the `nodes` arm** (`--drug-mechanism-before-first-layer`).
It is the only one of the three options under which the two arms of the ablation have the same reach at the same
depth, so a difference between them is a difference of representation. Keeping three layers without it makes the
`nodes` arm one hop short (a share of 0.11 of the nodes); a fourth layer changes every gene knockout's reach from a
median of 175 to 4,936 nodes and adds a layer's parameters, in both arms or in one. The cost is two more matrices in
the `nodes` arm. The alternative reading is also defensible: three layers and no early delivery puts drugs and gene
knockouts at the same number of edges from their proteins, which no configuration so far has done. That is a change
to how drugs are modelled and not a like-for-like ablation, so it would be a separate arm.

**C. State in the amendment that the ablation is informative for message passing only**, and report the
linear-response pair as the check that it comes out null.

**D. The lockbox on v3.** `docs/drug_targets_any_type.md` already records the membership: doxepin leaves the data,
and apomorphine, pergolide, pyridostigmine and phenelzine join groups that hold lockbox_v2 perturbations, with no
group holding lockbox and development perturbations together. The tables rebuilt in this session agree: 316 of the
317 lockbox_v2 perturbations are in the data, and the 266 groups that hold one hold 320 perturbations, the 316 and
those four drugs. That carried lockbox has 31 drugs, 22 of them with a kept positive, and 492 kept positive pairs,
against 28, 11 and 397 in lockbox_v2 on the v2 tables; development has 1,248 perturbations, 123 of them drugs, and
2,161 kept positive pairs. Carrying the groups over keeps every perturbation the pilots never read out of
development; a redraw would move perturbations the pilots have read into the lockbox. Carry-over is the
recommendation, and the file has to pin tables that travel with the repository (`docs/data_layer_rebuild.md`).

**E. Three things an amendment should state because they change a registered model without a configuration edit.**
Commit 9df62ee made `minimum` the default of `--conjunction-aggregation`, so the message-passing configurations
changed when the default did. `FULL_GRAPH_ARGUMENTS_V2` still names `evidence_full_v2` and `better_v2_full_v2`,
which cannot be rebuilt. And the descriptor treatment the rule picked for message passing, `seed_masked`, now has a
definition for drug nodes that the slice arms never ran.

**F. Whether message passing gets more rounds with reused weights.** Raised by the user on 10 October: "Can set it to
reuse parameters across layers (there are multiple ways this could be tried) and I meant until all reachable sink
nodes have been reached, which shouldn't need convergence just a test for whether the new round added any new nodes
reached and, if not, then stop." The stopping test is sound: the reached set only grows and is bounded, so it ends
whatever cycles the graph has, and each perturbation's stop round is a property of the graph that the search above
already gives. The options, set out for the user and not built: which weights are reused (one layer; the registered
three-layer block repeated, of which the registered encoder is the case of one repeat; an own first layer and a shared
rest; a gated update), how the repeated update is kept stable (damping, averaging across relations, a gate), and when
it stops (each perturbation's own saturation round, 13 to 26; one number for all, 26 for full reach or 9 for about 96
percent; or convergence, which the linear-response encoder already is). The recommendation given was the registered
block repeated with one round count for every perturbation, starting at three repeats, with the per-perturbation stop
as a switch. With enough rounds decision B stops mattering, because one hop is then a small share of the depth.

No amendment is drafted yet. A, B and F decide what it would say, and they were open with the user when this was
written; the draft follows those decisions and will be a file of its own, for the user, outside
`docs/preregistration.md`.
