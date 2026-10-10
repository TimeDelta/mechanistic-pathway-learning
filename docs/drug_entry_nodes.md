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
