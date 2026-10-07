# What limits the electrical layer (generated from the measurement in this file's history, 7 October 2026)

The neuronal variant adds a membrane potential node, `VM_c`, so that a perturbation of a channel, a pump or a
receptor can reach a symptom through excitability (design section 5.2, component 2). This note records what actually
limits that route, because the first number reported for it was the wrong one and led to the wrong priority.

## The correction

An earlier report said "only 43 of 451 slice perturbation genes have a path to VM". That was the count of genes whose
**own catalysed reaction** carries a `changes_membrane_potential` edge, which is a one-reaction route, not a path. By
directed path the picture is the opposite:

| graph | perturbation genes | reach VM by their own catalysed reaction | reach VM by any directed path |
|---|---|---|---|
| graph_neuronal (slice) | 451 | 40 (8.9%) | 445 (98.7%) |
| graph_full_neuronal | 451 | 40 (8.9%) | 450 (99.8%) |

Nearly every perturbation gene reaches the membrane potential. Reach is not the limitation.

## What the limitation is

Shortest directed distance from a perturbation gene to `VM_c`:

| edges | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 10 | 12 | unreachable |
|---|---|---|---|---|---|---|---|---|---|---|
| slice | 40 | 3 | 190 | 4 | 185 | 1 | 19 | 1 | 2 | 6 |
| full | 40 | 13 | 236 | 41 | 105 | 0 | 14 | 0 | 1 | 1 |

The mode is four edges, and a large group sits at six. Against the measured decay of the linear-response field under
in-degree averaging (the change at L-dopa is about 1e-5 of the response at the perturbed gene and at dopamine about
1e-8, two and three edges out; design section 5.2), a four to six edge path puts the arriving signal at roughly 1e-8
to 1e-12 of the perturbation. These genes reach the membrane potential in the graph and do not reach it in the field.
**The binding constraint is magnitude decay over path length, not connectivity.**

The second constraint is on the output side. `VM_c` has 409 inputs and **7** `voltage_gates` edges, to 7 reactions
catalysed by 18 genes in total, on both graphs. Even a membrane potential that moves can only act on those 7
reactions, so the feedback arm is narrow independently of how well the input arm carries signal.

## What the full graph does and does not do

The full graph was expected to be the lever here. It is not. It leaves the one-reaction route unchanged at 40 genes
and the gating fan-out unchanged at 7 edges and 18 genes, because gating is set by the curated channel-family patterns
rather than by graph size. What it does do is shorten paths: about 80 genes move from six edges to four or five, so
the median distance falls from five to four. That is a real gain of roughly one order of magnitude in arriving
signal, and it is worth having, but it does not change what limits the layer.

## Consequences for the order of work

1. The magnitude fix comes first, not the full graph: replacing per-node in-degree averaging with a global spectral
   scaling in the linear-response encoder, which keeps the map linear and so keeps the difference-field property, and
   the per-relation softmax mixture over mean, max, min and standard deviation in the message-passing encoder, where
   nonlinearity costs nothing. A four-edge path is the case these have to carry.
2. Widening the gating fan-out moves up. Splitting the lumped channel reactions by subfamily adds gated reactions and
   gated genes; at 7 reactions and 18 genes the feedback arm is the smaller of the two constraints only because the
   input arm currently delivers nothing.
3. The full graph stays worth building for its shorter paths, and it is built (data/processed/graph_full_neuronal,
   36,865 nodes), but it is no longer the first thing to try.

## Every slice perturbation is a loss of function, so the mirror constraint is untested (7 October 2026)

All 451 perturbations in the monogenic slice carry perturbation sign -1.0; the evidence table holds 861
records, every one of perturbation_type `gene`, relation `induces`, evidence_class `monogenic`, grade
`A`. The model never sees a perturbation of the other direction.

The design's requirement that a gain of function mirror a loss of function (section 5.2) is therefore
untestable rather than satisfied: with the input sign constant, nothing in the data distinguishes an
odd encoder from one that is not. The constraint earns its place when a direction of the other sign
enters, which is the drug arm, where an agonist and an antagonist are opposite-signed perturbations of
the same node.

The objection to a plain maximum has to be restated on grounds that survive this. It was first stated
as a mirror violation, which a constant input sign makes moot. The reason it still fails is that it
clips the falling half of the field: the messages at a node carry both signs whatever the input sign
is, so a maximum over negative messages returns the least negative, which at the first step is the
zero of an untouched neighbour. That is about message signs rather than input signs.

## Why the midpoint between the sum and the mean is dominated (7 October 2026)

Taking both scalings and splitting the difference is the natural compromise, and the geometric midpoint
between a sum S and a mean S / d is sign(S) |S| / sqrt(d), which is in-degree normalisation at exponent
0.5: the arm removed earlier the same day. Two things came out of measuring it properly.

The removal rested on a confounded measurement. The implemented arm divided by in-degree to the
exponent and then rescaled globally by the largest row sum. Contraction does not need that; it needs
the spectral radius of the unsigned aggregate below one, which is smaller. On the metabolic graph:

| scheme | largest row sum | spectral radius |
|---|---|---|
| mean, exponent 1, with the relations-feeding divisor | 1.000 | 1.000 |
| mean, exponent 1, without it | 1.000 | 1.423 |
| midpoint, exponent 0.5 | 45.465 | 12.377 |
| sum, exponent 0 | 2067.000 | 41.166 |

So the arm was dividing by 45.5 where 12.4 suffices, paying 3.7 times more than contraction requires.
The recorded field of 2.9e-15 at VM_c measured that waste rather than the exponent.

The midpoint is still dominated, and provably. One perturbed producing reaction at in-degree d
contributes 1 / d under the mean, 1 / (12.4 sqrt(d)) at the midpoint and 1 / 41.2 under the sum. The
midpoint beats the mean when sqrt(d) > 12.4, so when d > 154, and beats the sum when
sqrt(d) < 41.2 / 12.4, so when d < 11. Those conditions cannot both hold, so no in-degree makes the
midpoint the best of the three, and d > 154 is out of reach for any propagating metabolite: the
heaviest non-currency one is arginine at 114 producing reactions.

The compromise is also aimed at the wrong quantity. The global constant is uniform over nodes, so it
sets the overall scale and leaves the relative weighting between nodes untouched, and the overall scale
is what the learned per-relation gain and the signed-logarithm response scale already absorb. A
per-fold macro AUPRC is a ranking, so a uniform factor of 41 costs nothing that a gain cannot undo.
What the sum buys is the relative stoichiometric weighting between a metabolite with five producers and
one with one, which no rescaling recovers. The crossover at d of about 41 also says the choice is
properly per node rather than global, and in-degree is a node property, so a learned exponent per
relation or per node type is a weaker instrument than it looks: metabolites alone span in-degree 1 to
114.

## In-degree normalisation dilutes, and it dilutes the readouts hardest (7 October 2026)

Within a relation the encoder takes a mean, so the message a node receives is independent of how many
edges feed it and each edge contributes one over that count. The intuition that a well-connected
metabolite should feel a change more strongly, directly rather than inversely with its edge count, is
the sum, which is what dropping the per-node divisor gives.

The dilution is modest at a typical metabolite and severe in a tail. Producing reactions per metabolite
run at a median of 2 and a 90th percentile of 4, so one perturbed producer contributes a half or a
quarter of its own signal. Excluding the 229 currency metabolites, the heaviest are arginine at 114
producing reactions, tryptophan at 111, phenylalanine at 108, lysine at 103, alanine at 100 and
tyrosine at 99, from protein turnover; 3.4% of non-currency metabolites are diluted at least tenfold
and 0.4% at least fiftyfold.

The tail is not where it would do least harm. Taking the 157 base metabolites the laboratory labels
read out, which map to 515 graph nodes:

| diluted by at least | laboratory readout nodes | all metabolite nodes |
|---|---|---|
| 10x | 135 of 501 (26.9%) | 4.6% |
| 50x | 36 of 501 (7.2%) | 1.0% |

A readout node is diluted tenfold about six times as often as a metabolite node in general, and the
reason is selection rather than anything about the model: a clinical laboratory measures the abundant,
central metabolites because those are the measurable ones, and central is the same thing as
well-connected here. In-degree normalisation therefore attenuates the signal most at the nodes the
laboratory sign agreement and AUROC are computed on, and phenylalanine and tyrosine, the readouts for
the one mechanism the design validates by hand, sit at 108 and 99 producing reactions. This is a
further argument for scoring the spectral arm, independent of the cancellation argument above.

**Not a defect, checked and withdrawn.** Of the 515 readout nodes, 73 are marked `is_currency`,
covering 16 of the 157 base metabolites and 869 of 4,930 label rows. That does not make those labels
unpredictable: `signed_stacked_adjacency` drops edges whose source is a non-propagating node, so a
currency metabolite receives the response and only declines to pass it on, which is what assumption A9
intends. The field at those nodes is intact.

## What can break an exact cancellation, and what cannot (7 October 2026)

An exact cancellation cannot be resolved from the message values: every statistic in the mixture is
odd, and a permutation-invariant odd function must vanish on a multiset symmetric under negation,
since negating it returns the same multiset. Breaking the tie therefore needs information from outside
the multiset. Four candidates were measured.

**Message sign, preferring the negative on the grounds that negatives are rare.** The premise holds
globally and fails where the rule would fire. Negative relation-messages are 21.9% of all messages,
but at a mixed-sign destination there are exactly 7,811 negative against 7,811 positive on the
metabolic graph, precisely one each. The global rarity comes from the single-sign destinations, which
are reactions receiving substrate_of and catalyzed_by. There is no rarity to exploit at the decision
point, and with the input sign constant the rule reduces to a fixed preference for whichever relation
carries the negative message, which the learned per-relation gain already parameterises.

**Relation index, which is what argmax does today.** It is odd, because it ignores the message signs,
but the index is the relation's position in the stack, an arbitrary convention. It resolves every tie
and means nothing.

**The stoichiometric counts, which is to say the sum rather than the mean.** This is the one that
works, and it is not a tie-break bolted on but a consequence of not dividing by in-degree. Five
producing reactions against one consuming reaction is a different node from one against five, and
in-degree normalisation erases exactly that difference: the five contribute their average, so the
count divides out and the row sum is mean(+1, -1) = 0 whatever the stoichiometry. Measured on the
7,811 both-produced-and-consumed metabolites of the metabolic graph, the producing and consuming
counts differ at 3,722 of them (47.7%), and 3,806 of 8,123 (46.9%) on the neuronal graph. So the sum
resolves about half the cancellations for free, and it is the spectral normalisation arm rather than
new machinery.

**The stoichiometric coefficients, for the half the counts leave tied.** Largely a dead end, and the
builder is discarding them: build_physiology_graph.py reads each coefficient only for its sign, as
`is_substrate = coefficient < 0`, and writes every edge with sign 1.0. Reading the raw Human-GEM
model, 6,172 metabolites are both produced and consumed, 3,271 of them with equal reaction counts, and
the production and consumption coefficient totals differ at only 92 of those (2.8%); just 14.5% of the
6,172 touch any non-unit coefficient. Carrying the coefficients is a correctness improvement worth
making on its own terms and it is not the tie-breaker.

That leaves roughly half the cancellations with nothing to break them, and for those zero is the right
answer rather than a failure. A metabolite made by one reaction and consumed by one reaction at equal
stoichiometry genuinely has no net response to a change that moves both equally. The encoder returning
zero there is correct; the defect was returning zero at the other half as well.

## Sign cancellation: the first measurement was wrong, and a third of destinations annihilate (7 October 2026)

**Retracted.** An earlier version of this section reported that 65 of 23,709 nodes (0.3%) received messages of
both signs and concluded that there was almost nothing for a richer aggregator to resolve. That count omitted the
encoder's derived depletes_substrate relation, which is where nearly every negative sign in the model comes from:
signed_stacked_adjacency appends the reverse of every substrate_of edge with sign -1, and the count was taken from
edges.parquet, where the only stored negatives are the 59 inhibits edges and the handful of mixed curated relations.
Recomputed on the relation stack the encoder actually builds:

| | graph | graph_neuronal |
|---|---|---|
| destinations receiving any message | 21,337 | 24,121 |
| destinations receiving both signs | 7,811 (36.6%) | 8,712 (36.1%) |
| destinations seeing both signs within one relation | 0 | 18 of 32,793 pairs (0.1%) |

The within-relation figure of 18 reproduces, so only the across-relation count was wrong. Every relation carries a
single sign; no edge differs in sign from its relation. The sign is therefore a property of the relation rather than
of the edge, and all the mixing happens across relations.

The cancellation is exact, and it is a statement about the operator rather than about the field. Read the row sums
of sum_r gain_r S_r at equal gains, which is the node's response to a change shared by all of its inputs: 7,437 of
the 21,337 nodes with any input have a row sum of zero to within 1e-9, and every one of them is a metabolite, 7,437
of 8,460 (88%). No reaction has one. The reason is arithmetic. A metabolite that is both produced and consumed
receives product_of with mean +1 and depletes_substrate with mean -1, and their mean is zero.

What that does and does not say. It does not say the field is zero at those metabolites: in the dynamics the
producing and the consuming reactions carry different upstream states, so what arrives is their difference, not
zero. A test written on the stronger reading failed and was rewritten. What it does say is that 88% of metabolites
are blind to the part of their input that production and consumption share, and pass on only the difference between
them. Under attenuation that is where precision goes: the two inputs are small and similar, and their difference is
smaller still.

Learned gains do not remove this, they only move it. relation_gain() maps signed relations through a sigmoid, so
their gains are positive, and they are shared across every node. The net direction at a both-produced-and-consumed
metabolite is set by one global ratio, the product_of gain against the depletes_substrate gain, so the model cannot
raise one metabolite and lower another from the same pair of relations. A difference field whose sign at a third of
its destinations is fixed by one scalar carries little that a readout can use, which is a candidate cause of both the
attenuation recorded above and the null result in docs/graph_content_null_results.md.

Two consequences follow, and they reverse earlier rankings.

1. A mixture over several aggregators per node type, deprioritised on the retracted measurement, is reinstated and
   implemented (CrossRelationAggregator). Building it ruled out three of the four statistics first proposed. The
   response has to mirror, so that a gain of function gives the negative of what a loss of function gives, which
   requires every statistic to be odd. A plain maximum is not merely even but sign-biased: over negative messages it
   returns the least negative, which at the first step is the zero of an untouched neighbour, so it suppresses the
   whole falling branch and a metabolite that should drop stays at exactly zero. Measured on a four-node test graph,
   a maximum-dominated mixture gave 3.8e-15 where the mean gave 0.048. The minimum fails in mirror image and a
   standard deviation is even. The odd replacements keep what each was wanted for: the signed maximum magnitude
   takes the largest message with its sign, the median is odd already, and the standard deviation multiplied by the
   sign of the mean is the asymmetry-attached-to-dispersion idea in the only form that respects the mirror.
   The mixture fixes blindness to the shared input, not direction: its weights are per node type, so all 8,460
   metabolites still share one choice, exactly as they share one gain ratio.
2. Averaging within a relation before combining relations is the questionable step, not the choice of statistic.
   Mass balance at a metabolite is a signed sum of production and consumption, not a mean of two means, so dividing
   by in-degree per relation discards the stoichiometric weighting that decides the direction. The spectral
   normalisation arm, which applies one global constant and no per-node divisor, preserves that weighting; it was
   ranked last on field magnitude at VM_c (0.70x the in-degree arm) and that ranking measured magnitude, not sign
   preservation, so it does not settle this. The arm needs a run scored on symptoms before either normalisation is
   preferred.

The decay is plain attenuation, and a median divisor of two does not produce 1e-10 at four edges, so the paths that
reach the membrane potential run through the hub destinations of the same distribution whose maximum is 2048.

## The normalisation arms

`NORMALISATIONS` in models/linear_response_encoder.py selects among three:

- `in_degree`, the default: divide each entry by the destination's in-degree under its relation times the number of
  relations feeding it, so every row of sum_r |S_r| sums to at most one.
- `spectral`: one global divisor, the spectral radius of the unsigned aggregate, and no per-node divisor. Kept as a
  measured negative result. That radius is 36.3 on graph_neuronal while the in-degree divisor is below it at 98.8% of
  destinations, so it divides harder than what it replaces nearly everywhere; the median field arriving at the
  membrane potential fell to 0.7 of its in-degree value, 0.6 at distance four and lower still beyond.
- `in_degree_power` (removed 7 October 2026): divided by in-degree raised to an exponent, then rescaled the whole matrix once so the
  row sums return below one. Note what this does and does not do. Where in-degrees are uniform the rescale exactly
  undoes the softening and the exponent changes nothing, which a test asserts. It acts only where in-degrees differ,
  moving division off the hub destinations and onto the sparse ones. Since the paths to the membrane potential are the
  ones through hubs, that is the direction wanted, but the mechanism is redistribution rather than reduction and the
  gain is bounded by how heterogeneous the in-degrees along a path are.
