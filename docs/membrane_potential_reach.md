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

The cancellation is not partial but bimodal. Taking all gains equal, so that the net message at a destination is the
mean over relations of the mean sign within each relation, the ratio |net| / mean |sign| is 1 at every destination
that receives one sign and below 0.01 at every destination that receives both: 7,811 of 21,337 on the metabolic graph
and 8,389 of 24,121 on the neuronal graph. The reason is arithmetic. A metabolite that is both produced and consumed
receives product_of with mean +1 and depletes_substrate with mean -1, and the cross-relation mean of those is zero.

Learned gains do not remove this, they only move it. relation_gain() maps signed relations through a sigmoid, so
their gains are positive, and they are shared across every node. The net direction at a both-produced-and-consumed
metabolite is set by one global ratio, the product_of gain against the depletes_substrate gain, so the model cannot
raise one metabolite and lower another from the same pair of relations. A difference field whose sign at a third of
its destinations is fixed by one scalar carries little that a readout can use, which is a candidate cause of both the
attenuation recorded above and the null result in docs/graph_content_null_results.md.

Two consequences follow, and they reverse earlier rankings.

1. A mixture over several aggregators per node or relation type, deprioritised on the retracted measurement, is
   reinstated. It is the stated fix for exactly this: a maximum, a minimum or a signed asymmetry term survives where
   the mean annihilates, and a third moment is available at the destinations that matter, since a destination
   receiving both signs receives at least two messages by construction.
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
- `in_degree_power`: divide by in-degree raised to `normalisation_exponent`, then rescale the whole matrix once so the
  row sums return below one. Note what this does and does not do. Where in-degrees are uniform the rescale exactly
  undoes the softening and the exponent changes nothing, which a test asserts. It acts only where in-degrees differ,
  moving division off the hub destinations and onto the sparse ones. Since the paths to the membrane potential are the
  ones through hubs, that is the direction wanted, but the mechanism is redistribution rather than reduction and the
  gain is bounded by how heterogeneous the in-degrees along a path are.
