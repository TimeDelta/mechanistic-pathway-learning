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
