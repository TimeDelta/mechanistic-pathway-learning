# Physiological adherence review

Every place the model departs from the physiology it represents, with the cost of each stated. 13 October 2026.

**How to read it.** A deviation is not a defect by itself: a model is a set of chosen simplifications, and several
here were chosen deliberately and are better than the alternative. What matters is whether the cost is known. Each
entry below is marked **measured** (a number in this repository sizes it), **open** (registered as a decision with
no number yet) or **undisclosed** (found in this review; no document named it before today). The undisclosed ones
are the finding of this review, and three of them contradict text the preregistration currently asserts.

Counts are taken on `data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes` (49,796 nodes, 280,966
edges) and `evidence_full_v4` unless another source is named.

---

## 1. Node-wise averaging erases the count of producing and consuming reactions — measured

**Biology.** Mass balance at a metabolite is a signed sum. Five producing reactions against one consuming reaction
is a different node from one against five.

**Model.** `LinearResponseEncoder.signed_stacked_adjacency` divides each entry by the in-degree under its relation
and by the number of relations feeding the node, so `product_of` arrives at mean +1 and the derived
`depletes_substrate` at mean -1. The count divides out.

**Cost, measured.** At equal gains the aggregate row sum is **exactly zero at 7,811 of the 21,337 destinations with
any input, and every one is a metabolite** - 7,811 of 8,460, 92.3 percent. No reaction is affected. Those
metabolites pass on only the difference of two small similar numbers. Because the relation gains are shared across
every node, the net direction at every both-produced-and-consumed metabolite is fixed by **one global ratio**, so
the model cannot raise one metabolite and lower another from the same relation pair
(`docs/membrane_potential_reach.md`).

**The price of the propagation law as a whole, and the most important number in this review.** Untrained, over
1,067 signed gene-metabolite pairs from HPO laboratory abnormalities, the linear response predicts the **direction**
of measured metabolite change **worse than always guessing "increased"**: sign agreement 0.615 at best against
0.806 for the majority direction, AUROC 0.574 (`docs/linear_response_sign_check.md`). Nothing yet separates "the
law is linear" from "the aggregation cancels" as the cause, and that separation is worth having.

**Stoichiometric coefficients proper are the smaller half of this.** The build reads each Human-GEM coefficient only
for its sign and writes every edge at 1.0. Of 6,172 metabolites both produced and consumed, the production and
consumption coefficient totals differ at only 92 of the 3,271 with equal reaction counts (2.8 percent), and 14.5
percent touch any non-unit coefficient. So it is the **counts**, not the coefficients, that the in-degree mean
destroys. The per-edge `bottleneck_edge_demand` that item 20 registers is `None`: every demand is 1, and the
owner's two-molecule acetylcholine example is registered and inert.

## 2. A reversible reaction cancels against itself — undisclosed

**Biology.** A reversible reaction has a net direction set by mass action; displacing a product drives it backwards.

**Model.** The build writes the mirror `substrate_of`/`product_of` pair for every reversible reaction, both at sign
+1, and the encoder derives `depletes_substrate` (-1) from every `substrate_of`. A metabolite that is both
substrate and product of one reaction therefore receives +1 and -1 **from that same reaction**, and the mean is
zero.

**Cost, measured here for the first time.** **15,867 metabolite-reaction pairs are both substrate and product of
the same reaction, over 5,729 distinct metabolites and 5,725 reactions.** Every one contributes an exactly
cancelling pair under the in-degree mean. This is a large and previously unnamed share of the 7,811 exact
cancellations in section 1, and its own consequence is sharper: **no perturbation can drive a reversible reaction
one way rather than the other.** Roughly two fifths of the reaction set is reversible, so this is not a corner.

**How much of the graph is mute because of it.** The **median affected metabolite has 67 percent of its incoming
traffic in self-cancelling mirror pairs**, 4,221 have more than half, and **1,809 have nothing else at all**: their
field is identically zero under every perturbation, whatever is done upstream.

**What the physiology asks for** (the owner, 13 October: "The reaction should only be able to follow the reversal
path under certain conditions: i.e. a competitor ligand with higher binding affinity for the target comes in
contact"). Exactly so: reversal is conditional on the local state, and an unconditional bidirectional pair is not
a weaker version of that but a different thing, which here cancels to nothing.

**An earlier draft of this section said the faithful version was unavailable to the encoder. That was wrong, in
two ways.**

- **The nominal direction needs no new data and no new machinery.** `build_physiology_graph` already writes it -
  `is_substrate = coefficient < 0` sets each edge from the stoichiometric sign - and then adds the mirror edge on
  top for every reversible reaction. Dropping the mirror leaves the direction Human-GEM writes the reaction in,
  and recovers the 1,809 silent metabolites, at the rebuild, for nothing.
- **The resting direction is not state-dependent, so it can be set before the propagation.** It is a property of
  the cell at baseline, and the propagation computes deviations *around* baseline; a thermodynamic estimate or a
  flux solution could set it for the 5,713 reversible reactions. That is the route competition already takes -
  "it can enter only as a model input computed before the propagation and not as a relation inside it"
  (`docs/label_source_review.md` 13.4) - and it is the template here rather than an obstacle.

**And a perturbation can reverse one - this review was wrong a second time.** It said flipping the direction
under perturbation was a state-dependent adjacency and therefore a nonlinearity. It is not: in a deviation model
reversal is the **sign of the rate**, not a switch on the edges. Linearised mass action, `dv = k_f dS - k_r dP`,
puts the product into the rate at a negative sign, so a perturbation that accumulates the product makes the rate
deviation negative and the reaction runs backwards, regenerating the substrate. One derived edge per product of a
reversible reaction, product -> reaction at -1, carries it, with no gate on the state. Built as amendment item 27
(`product_opposes`); the build now writes the nominal direction only, and the encoder refuses a graph that still
has the mirror. What remains is the limit of linearity everywhere: nothing saturates, so a reaction can run
backwards faster than its reverse capacity would allow.

*(The pralidoxime example is a near neighbour rather than this case: reactivating a covalently inhibited enzyme is
a displacement at a binding site, which belongs with the competition work of section 12 and is equally
unrepresented.)*

## 3. The limiting-step rule is registered with no measurement at all — open

**Biology.** The throughput of a chain at steady state is set by its weakest step, not by the product of its steps.
The rule is right, and it is the owner's (13 October: "min along the path makes WAY more sense to me").

**Model.** `b(t) = max over routes of min over nodes on the route`, applied **once** to the settled state.

**Three costs the encoder states and nothing measures.** The bottleneck **overstates**, because the settled state
already sums contributions from every route including thin ones and `b` then scales the whole sum by the widest
route's ceiling. The product rule **understates** in the other direction. And `b` is computed on the **unsigned,
unweighted** structure, so a route through an inhibitory edge raises the ceiling exactly as an activating route
does. The code's justification is "a bottleneck cares only whether a route exists, so signs, weights and the
in-degree divisor play no part in it". For **existence** that holds, and excluding the in-degree divisor is right
because it is a modelling artifact rather than physiology. Excluding the **sign** does not follow from it: `b`
scales a signed sum, so the ceiling can be set by a well-expressed inhibitory route while the response travels an
activating one, and those are different routes with different expression. The fix sits where `b` already does -
the bottleneck is applied after the fixed point, so the settled sign at the target is known and the max-min can
be taken over routes whose sign product matches it, with no nonlinearity entering the loop. Decision 4.12.

**No measurement decides whether the rule is used, and this review was wrong to imply one should** (the owner, 13
October: "that doesn't need measured for confirmation to be in the final models. It should just be in there to
mirror the physiology"). The rule is chosen because the throughput of a chain is set by its weakest step, which is
a fact about pathways and not a hypothesis this study tests. A comparison against the product rule would say how
much the choice moves predictions - a different question, and not a reason to hold the rule back.

**What does need checking is whether the implementation is the rule.** The first of the three costs above is not a
property of a limiting step at all: a bottleneck computed on unsigned structure lets an inhibitory route raise the
ceiling. That is a correctness question answerable by a unit test, and
`experiments/measure_limiting_step_rule.py` is repurposed to it. The second and third are properties of the rule
as registered and are disclosed rather than fixed.

## 4. Cell-class expression no longer conditions the route — open, and the decision log is stale

**Biology.** Whether a perturbation reaches a symptom depends on which cells express the intervening steps, so the
*route* is cell-type-specific, not only the gain.

**Model.** Under `limiting_step` the weights are taken out of the propagation entirely; the class's expression
re-enters only as the post-hoc ceiling. The encoder says so plainly: "The rule therefore conditions the magnitude of
a class's response on its expression and no longer conditions the route."

**Cost.** Any prediction that turns on a route existing in one cell type and not another is unidentifiable.
**Unquantified** - and `measure_limiting_step_rule.py` item 4 was written for exactly this ("the loss of
class-specific routing is a number and not an assertion").

**Compounding it: 34.0 percent of nodes are flat at 1 in every class** - every metabolite, every protein entity,
`VM_c`, and every gene or reaction with unknown expression (median weight over the 32,815 varying nodes is 0.372).
Under a max-min rule a route through a flat node costs nothing, so a long metabolic chain between two expressed
enzymes is free. The flat default is deliberate, "so a class channel does not block what nothing is known about",
but nothing sizes what it does to the bottleneck.

**Three documents describe the superseded product rule as current**: `docs/open_decisions.md` 4.3,
`docs/label_source_review.md` 28.4 (whose correction was itself reversed by item 20) and the module docstring of
`brain_expression_weights.py`.

## 5. The compartment pH ceiling limits no throughput — undisclosed

**Claim.** `compartment_ph.py`: "the enzymes of an acidified compartment simply pass less". Item 22a: "a multiplier
on the throughput of every node in that compartment".

**What the code does.** `state = state * compartment_ph_weight(state)`, applied **once, after the fixed point and
after the bottleneck**. Every node in an acidified compartment passed its full message to every neighbour
throughout the propagation; only its own final value is scaled. Unlike the cell-class bottleneck, which genuinely
carries its ceiling along routes through the max-min recurrence, **the pH weight carries nowhere: a collapse of
lysosomal pH cannot change the field at any cytosolic node.**

**Reach.** 18,928 of 49,796 nodes take a weight. **25,741 nodes have no compartment at all** - every gene, protein
and drug node - and **5,127 transport reactions carry a multi-compartment label and so take none**, which is
precisely the class of reaction whose rate a proton gradient most directly sets.

**And the node the ceiling is read from is the worst-conditioned in the graph.** Cytosolic H+ (`MAM02039c`) has
2,048 incoming `product_of` and 1,281 outgoing `substrate_of` edges, so one perturbed reaction moves the reading by
1/4,096 or 1/2,562 of its own change - a dilution about ten times worse than the 405 at `VM_c` that the review
called a defect - and it is one of the exactly-cancelling metabolites of section 1. This is why the settled field
there is of order 1e-11 and why the sensitivity has to be calibrated to a numerical target rather than a
physiological one.

## 6. The pH weight in code is not the form the preregistration registers — undisclosed, and the registration is wrong

**Item 22a registers** `w_pH = exp(-(sensitivity * deviation) ** 2)` and registers as a property that "the form is
even in the deviation, so a gain-of-function stays the exact negative of a loss-of-function".

**The code implements the log form**: `ΔpH = -log10(1 + sensitivity * deviation)`, `w = exp(-(ΔpH / width) ** 2)`.
The encoder's own docstring says what follows: "A pH shift is not even ... So a gain and a loss now take different
weights at a node whose compartment the perturbation moves, and the two stop being exact negatives there."

**Cost.** The mirror property is not cosmetic - items 20 and 22 lean on it ("the response must mirror: a gain of
function has to produce the negative of the field a loss of function produces"), and the sign-and-typing ablations
read it. The code is right and the physiology is on its side: acidifying by a fraction and alkalinising by it do
move pH by different amounts. **The registration is what is wrong**, and it asserts a property the model does not
have. Corrected in the preregistration today as item 22b.

## 7. The membrane potential averages the currents — measured, disclosed

**Biology.** A membrane potential is the integral of the currents crossing the membrane, a sum, and each channel
then reads the one potential.

**Model.** `VM_c` has in-degree 409 and its output to each of 7 channels passes undivided.

**Cost, measured.** Deleting the 405 reaction edges **lowers** the field at `VM_c` by a factor of 26 (3.72e-07 to
1.44e-08), because the binding divisor is the in-degree at the ion pool - cytosolic Na+ at 364 - and not the 405 at
`VM_c`. Both normalisation alternatives the encoder carries are worse. `VM_c` also carries none of the 138
descriptor columns. Disclosed and open (decision 4.1).

**Two things about the electrical layer that are not disclosed.** First, the linear-response encoder **never reads
`edge_weight`** except on drug entry edges, so item 13's six fixed conductance weights - potassium at 0.318 per
side against 0.091 for sodium and chloride, "the 3.5-fold conductance share that a weight of 1 on every edge threw
away" - are **inert in that encoder**: K+, Na+ and Cl- reach `VM_c` with equal magnitude. The message-passing
encoder does apply them, so the two registered encoders disagree about the physics. Second, the graph on disk
carries the superseded `learned_from_reactions` build - 4 unit-weight `driving_force` edges, and neither of the two
neuronal extracellular pools item 13 adds. Item 13 reads as describing the graph; it describes the rebuild.

## 8. One species, two disconnected pools — measured, with one case undisclosed

**Protons.** Human-GEM `MAM*` against Reactome `CHEBI*`, by compartment: cytosol 3,329 against 44, extracellular
719 against 8, mitochondrion 686 against 14, ER 559 against 1. The Human-GEM copies are currency and send nothing;
the Reactome copies propagate. **The only proton signal that travels does so through one or two percent of the
reactions that touch a proton.**

**Ions.** **Zero edges join the five `NEURON_` pools to the eight metabolic ion pools.** A perturbation that changes
sodium through metabolism cannot move the membrane potential, and the reverse.

**Chloride, which no document names.** `MAM01442c` and `MAM01442e` are **not** currency and do propagate, in
parallel with `NEURON_MAM01442c` and with no edge between them. So chloride runs **two live, disconnected pools**
rather than one live and one silent. That matters more than it sounds: the GABA-A reversal potential is set by
intracellular chloride through NKCC1 and KCC2, and it is the mechanism behind the developmental switch and several
drug effects this study has labels for.

## 9. Currency metabolites receive but never send — measured

**Model.** Edges whose **source** is a currency node are dropped, so the node computes a deviation and passes
nothing on. **222 currency nodes; 13,013 of 280,966 edges (4.63 percent) leave one and are dropped**, against
14,178 that arrive and are kept.

**What is computed and discarded**, by reaction degree: pH 6,181, redox about 4,600, acyl carriers about 2,900,
energy charge about 2,700, oxygen 877, methylation 368. The rule is right - a proton hub with 3,329 edges would
make everything adjacent to everything - but the cell-state quantities it defines are read nowhere.

**Four tagging defects were fixed in code on 13 October and are in no graph yet**: a charge suffix defeated the
name list, so cytosolic NAD+ sent nothing while peroxisomal NAD+ sent freely; nodes added after tagging were never
tagged; cytosolic glycine was tagged currency by degree although it is an inhibitory transmitter and the NMDA
co-agonist; and the fix would have tagged `NEURON_MAM02519c` and silenced the electrical layer, which until now was
protected only by the second defect. The set goes from 222 to 290 at the rebuild.

## 10. Protein complexes: a flat ceiling, and AND is indistinguishable from OR — measured

**The cell-class weight reaches a reaction through its catalysing gene**, and a Reactome reaction's catalyst is
often a complex, which carries no weight. Of 7,769 reactions regulated by a gene or protein node, 98.9 percent have
a class-varying weight; of the **413 regulated only by a complex, none does** - every one is flat at 1.0 in every
class, adenylyl cyclase among them. 1,527 of the 1,681 complexes have a member that already carries a weight, so
the fix is available and is open decision 4.2.

**An obligate complex is a conjunction and a DefinedSet is a disjunction, and the graph gives both one relation.**
**542 entities need all their members and 238 need any one, and nothing in the graph tells them apart.** The
message-passing encoder has a per-channel minimum for `member_of`; the linear-response encoder has no conjunction at
all, so losing one member of a two-member complex leaves **0.50** of the input instead of abolishing it. 147
two-member conjunctions, 1,202 member protein nodes, 201 of the perturbed genes. The owner accepted this on 9
October as imposed by the encoder, but it is a live difference between the two registered encoders and it is not in
the preregistration's open list.

## 11. A knockout and an inhibitor are nearly the same input — open

**Biology.** A germline knockout removes the protein from development onward in every cell; an inhibitor blocks
activity, partially, reversibly, in the cells the drug reaches, in an adult.

**Model.** Both arrive at the protein with sign -1 and magnitude 1, differing only in the entry hop. No dose,
affinity, occupancy or brain exposure enters the registered input; 154 of 222 drugs carry a total input of 1.

**Cost.** The three fixes the owner directed - occupancy, exposure, amount against activity - are unbuilt. The
nearest measured number sizes the consequence rather than the conflation: the development base rate is 0.159 for
drugs against 0.063 for genes, and a predictor that knows only whether a perturbation is a drug gains +0.059 micro
[+0.018, +0.091] over the better baseline.

**And nothing anywhere represents adaptation.** No receptor desensitisation, internalisation, downregulation or
tolerance appears in the code or the documents, while the labels are largely **chronic**-dosing adverse events.
Transcriptional feedback is partly representable through the 41,638 `regulates_transcription_of` edges; receptor
and transporter adaptation is not. This is the largest unrepresented piece of pharmacology in the study and it has
no entry in the open list.

## 12. Competition, partial agonism and cooperativity — open, with counts

**Competition is not modelled.** No source labels it for these drugs: none of the 647 human Guide to PHARMACOLOGY
rows for the study's drugs is marked competitive or non-competitive. Computed from a shared binding site, of 9,908
(drug, drug, gene) triples 1,184 share a site, 932 are at different sites and 7,792 have at least one site unknown.
Three registered limits: competition is saturable and so cannot be a relation inside a linear propagation; the
drug-against-endogenous half needs a synaptic ligand concentration no source constrains to within orders of
magnitude; and **of 80 positive two-drug triples only 2 have opposing actions at one target**, so there is almost
nothing to test it on.

**Partial agonism is a constant.** Partial agonist, PAM, positive modulator and stabiliser map to +0.5 and NAM to
-0.5: **357 entry rows at +0.5 and 19 at -0.5, over 41 drugs, 33 of which have no other kind of target.** Two
errors are registered with it: a partial agonist's direction depends on ambient tone (aripiprazole, varenicline,
buspirone, memantine), and a PAM does nothing without the endogenous ligand (every benzodiazepine). The fix -
scaling the endogenous ligand's state rather than adding to it - is not linear, so it cannot sit inside the
linear-response encoder without ending the linearity items 20 and 22 rest on.

**Cooperativity is not carried.** ChEMBL's 350 cooperativity activities run 0.2 to 3000, and 8 of 535
cooperativity-style records touch any of the study's molecules, none of them a cooperativity.

**And 49 of the 1,058 mechanism edges carry sign 0.0** - SUBSTRATE, MODULATOR, BINDING AGENT, CROSS-LINKING AGENT,
OTHER. Those drug-target mechanisms are in the graph and deliver exactly nothing. Undisclosed.

## 13. The head assumes equifinality with no interaction between modules — measured

**Claim.** "Any single module whose link times activation is near 1 explains the symptom on its own, so the
functional form encodes equifinality without a cross-module interaction term."

**Cost, measured.** The sufficiency test finds the own-ablation interval excludes zero in **0 of 62 symptom-module
rows** for the mechanistic model and 2 of 44 for the linear response. **No module is sufficient for any symptom**,
which is the premise of the functional form failing on its own terms.

**A second, structural cost.** Every symptom is a monotone function of the same K module activations, so symptoms
cannot differ in how modules combine - only in how strongly. Measured as an identity: with links at their common
initial value every untrained symptom ranks the perturbations identically, correlation 1.000000 and zero
perturbations ranked differently.

## 14. Smaller departures, each with its cost

- **One elasticity per relation, shared by every node.** The encoder cites metabolic control analysis, where an
  elasticity is specific to a reaction and a metabolite; the model has 26 gains per channel for 13,793 reactions
  and 280,966 edges. This is the mechanism behind the "one global ratio" problem of section 1. Undisclosed as a
  departure.
- **The linear response is applied outside the regime that justifies it.** The law is derived for "a small
  enzyme-level change"; every gene perturbation is a complete knockout at magnitude 1 and drug blockade is full.
  No saturation term, no stated validity range. Priced only jointly with section 1.
- **A propagation step is not a unit of time.** "The graph has no rate constants, so a step is not a unit of time."
  Transcription (hours), trafficking (minutes) and binding (milliseconds) collapse onto one fixed point, separated
  only by learned gains. Tissues fill at median step 2 and reach 0.631 of their settled level by step 8, so the
  step index separates by graph distance and not by rate; a peak-concentration reading cannot be taken from it.
- **The subthreshold claim suppresses real regenerative feedback.** The seven `voltage_gates` loops through `VM_c`
  are genuine positive feedback and are suppressed by the contraction bound. Disclosed in item 13 and correctly
  flagged as not to be read as excitability.
- **The resting pH table is curated, cited and never read.** Ten values with provenance; `compartment_ph_weight`
  computes a deviation from the compartment's own rest, so the resting value cancels and no absolute-pH prediction
  - ionisation, lysosomal trapping of bases - is possible. `docs/data_sources.md` describes them as values the
  weight reads. Undisclosed.

---

## What this adds up to

**The finding that should change something.** Reversible reactions cancelling against themselves: 15,867 pairs,
and 1,809 metabolites whose field is identically zero under every perturbation. It is undisclosed, it sits
underneath the study's worst measured result, and the nominal-direction build recovers it at the rebuild.

**Two things this review got wrong, corrected by the owner on 13 October.** It treated the limiting-step rule's
missing measurement as a gap; no measurement decides whether a rule that mirrors the physiology is used, and what
needs checking is narrower - whether the implementation is the rule. And it reported item 22a's gain-versus-loss
mirror - the claim that a gain of function produces the exact negative of the field a loss of function produces -
as a divergence between registration and code, when the property should never have been registered at all: the
owner never accepted it, and it is withdrawn rather than corrected.

**The remedy for most of it is already written and is not the default.** `total_in_degree` divides each entry by
the destination's in-degree summed over every relation, so the unsigned row sums are one exactly as under
`in_degree` - the same contraction, no global shrink - while the counts survive: five producers against one
consumer give (5-1)/6 rather than 0. The per-node mean cannot be justified as physiology; it is there so the
iteration contracts, which is a numerical requirement. Decision 4.13.

**The pattern worth naming.** Most of the costed deviations are the same deviation: the in-degree mean. It erases
reaction counts (section 1), cancels reversible reactions (2), dilutes the proton node to 1/4,096 (5) and the
membrane potential to the pool's in-degree (7). The study's one external check of the propagation - direction
agreement 0.615 against a 0.806 majority baseline - is a check of that choice more than of anything else, and it
fails it. That is the deviation to measure next, and it is more consequential than any individual item below it.

**What is well handled.** The currency rule, the two-pool diagnosis, the complex conjunction semantics, the
competition survey and the partial-agonism sign are all documented with counts and with the objection stated
before the decision. The limiting-step rule is better physiology than what it replaced. The subthreshold caveat is
disclosed rather than glossed.
