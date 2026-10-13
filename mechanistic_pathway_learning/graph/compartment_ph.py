"""Compartment pH as a throughput ceiling, not as a message (the owner, 13 October 2026: "i don't think we've
incorporated compartment pH yet, which is important", and of this design: "perfect solution. implement it").

What the graph already has. Protons are nodes, one per compartment, with 6,271 reaction edges between them and the
reactions that produce or consume them: on graph_full_neuronal_split_binders_drugs_v4_enzymes, 2,048 producing and
1,281 consuming reactions at the cytosolic proton alone. So a perturbation's effect on a compartment's proton balance
is already computed. Two things stop it being used.

- **It goes nowhere.** The nine Human-GEM proton nodes carry `is_currency`, and a currency node's outgoing edges are
  dropped (linear_response_encoder.signed_stacked_adjacency), so the node still *receives* and its deviation is
  computed, but it sends nothing. That is the right call for a message: a node touching 3,329 reactions that sent
  freely would couple the graph to itself through H+ and make every pair of reactions neighbours.
- **There are two pools for one species.** Human-GEM contributes `MAM02039*` and Reactome `CHEBI15378*`, with
  different base metabolite ids, so they never merged. The Human-GEM copy holds 50 to 100 times the reactions of the
  Reactome copy in every compartment both describe, and only the Reactome copy propagates. Whatever pH signal exists
  today therefore travels through about one to two percent of the reactions that touch a proton.

What this module does instead. A proton deviation is read **after** the propagation has settled, at whichever pool
node governs each compartment, and becomes a multiplier on the throughput of every node in that compartment:

    w_pH(compartment) = exp(-(sensitivity * proton deviation) ** 2)

A weight does not propagate, so the hub problem does not arise; the enzymes of an acidified compartment simply pass
less. The form is registered and the magnitude is learned, which is the honest division, and three things follow from
it being even in the deviation and applied after the fixed point: a gain-of-function stays the exact negative of a
loss-of-function, since the deviation changes sign and the square does not; the fixed point is reached before the
weight is applied, so the contraction the stopping rule rests on is untouched; and the weight lies in (0, 1], so it
cannot raise a norm. What it does break is positive homogeneity of degree one - doubling the input no longer doubles
the response, because the weight reads the deviation's size - and that is the registered cost.

**Why no pH constants are needed, which is not obvious.** Converting a deviation to pH units would need the resting
proton concentration per compartment and a scale from the model's arbitrary units to molar, and a pH-activity curve
would need each enzyme's pH optimum. Enzyme-specific optima are not available at the scale of 13,793 reactions. The
defensible simplification is that each reaction's optimum is its own compartment's resting state, which is where the
compartment sits when nothing is perturbed. The weight is then a function of the deviation from rest alone, the
resting value cancels out of the form, and no compartment pH constant enters. This is also why the module claims
nothing about absolute pH: it says that moving a compartment away from its resting proton balance costs throughput,
and how much is learned.

The constants would be needed for a different question, the ionisation of a drug against a compartment's pH and so
the trapping of bases in the lysosome, which is not this module.
"""
from __future__ import annotations

import numpy as np

# the pool to prefer where a compartment has more than one proton node: the one that the reactions actually use
PROTON_DISPLAY_NAME = "H+"
PREFERRED_PROTON_PREFIX = "MAM"  # Human-GEM, 50 to 100 times the reaction count of the Reactome copy per compartment


def proton_node_of_compartment(node_ids, node_display_names, node_compartments, node_degrees=None) -> dict[str, int]:
    """compartment -> the index of the proton node that governs it.

    Where both pools describe a compartment the Human-GEM node wins, because the reactions are there; where only one
    does, that one governs. A compartment named by a transport reaction ("c;m") is not a compartment of its own and
    never appears here."""
    candidates: dict[str, list[int]] = {}
    for index, (node_id, display_name, compartment) in enumerate(zip(node_ids, node_display_names, node_compartments)):
        if str(display_name).strip() != PROTON_DISPLAY_NAME or ";" in str(compartment) or not str(compartment):
            continue
        candidates.setdefault(str(compartment), []).append(index)
    governing: dict[str, int] = {}
    for compartment, indices in candidates.items():
        preferred = [index for index in indices if str(node_ids[index]).startswith(PREFERRED_PROTON_PREFIX)]
        pool = preferred or indices
        if len(pool) > 1 and node_degrees is not None:  # same source twice should not happen; take the one the reactions use
            pool = [max(pool, key=lambda index: node_degrees[index])]
        governing[compartment] = pool[0]
    return governing


def compartment_index_per_node(node_compartments, compartments_in_order) -> np.ndarray:
    """[N] index into compartments_in_order, -1 for a node that is in none of them (no compartment, or several, as a
    transport reaction has). A node with no governing compartment takes no pH weight."""
    position = {compartment: index for index, compartment in enumerate(compartments_in_order)}
    return np.asarray([position.get(str(compartment), -1) for compartment in node_compartments], dtype=np.int64)


def compartment_ph_inputs(data) -> dict:
    """The encoder keywords for the pH weight, or {} when the graph names no proton node in any single compartment."""
    if getattr(data, "node_compartment", None) is None:
        return {}
    degrees = getattr(data, "node_degree", None)
    governing = proton_node_of_compartment(data.node_ids, data.node_display_name, data.node_compartment, degrees)
    if not governing:
        return {}
    compartments = sorted(governing)
    return {"ph_compartment_names": compartments,
            "ph_proton_node_index": np.asarray([governing[compartment] for compartment in compartments], dtype=np.int64),
            "ph_compartment_of_node": compartment_index_per_node(data.node_compartment, compartments)}
