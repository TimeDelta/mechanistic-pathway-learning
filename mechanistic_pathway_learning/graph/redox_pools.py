"""Reducing capacity as its own nodes: the glutathione and NADPH pools (graph variant of design section 5.7, with
graph/reactome_import.py and graph/oxidative_regulation.py; experiments/build_neuronal_graph_variant.py).

Glutathione and NADPH are currency metabolites in the graph (cytosolic GSH has 209 edges, NADPH 637): they receive a
signal but do not pass it on, because a metabolite that touches hundreds of reactions would let every reaction reach
every other in two steps (assumption A9). That rule costs the oxidative layer its central mechanism, since a fall in
reducing capacity is how an oxidant load reaches the enzymes that defend against it, and how a defect in glutathione
synthesis or in the pentose phosphate pathway reaches oxidant handling.

This module adds one node per compartment for each pool and connects it to the reactions that set it, from Human-GEM's
own stoichiometry rather than from a gene list:

  - REDOX_GSH_<compartment>, the reduced glutathione pool. Every reaction that consumes glutathione in that
    compartment without producing it there (the peroxidases, the transferases, the conjugations, and a transport
    reaction that carries it out) lowers the pool, changes_redox_pool -1; every reaction that produces it there without
    consuming it (glutathione reductase, the synthesis step, the hydrolases that release it, and a transport reaction
    that carries it in) raises it, +1. A reaction with it on both sides of the same compartment changes neither.
  - REDOX_NADPH_<compartment>, the NADPH pool. Every reaction that reduces NADP+ to NADPH raises it (+1) and the
    reactions of the antioxidant demand, those that consume NADPH while reducing glutathione or thioredoxin, lower it
    (-1). NADPH's other consumers (fatty acid and sterol synthesis, the cytochromes P450: around 500 reactions in the
    cytosol alone) are left out: they are real competing demands, but at that number the in-degree average would make
    each one negligible and the pool would barely move.

Each pool then feeds back, redox_capacity +1, to the reactions that depend on it for oxidant handling only: for
glutathione the reactions that reduce a peroxide or a hydroperoxide with it, for NADPH the antioxidant demand above.
Restricting the feedback keeps the pool's out-degree in the tens, so it carries a mechanism without becoming the hub
the currency rule exists to prevent; the conjugation reactions still deplete it, so a xenobiotic or quinone load
reaches the defence.
"""
from __future__ import annotations

from collections import defaultdict

import pandas as pd

POOL_CHANGE_RELATION, POOL_CAPACITY_RELATION = "changes_redox_pool", "redox_capacity"
GLUTATHIONE, GLUTATHIONE_DISULPHIDE = "MAM02026", "MAM02027"
NADPH, NADP = "MAM02555", "MAM02554"
THIOREDOXINS = frozenset({"MAM02990", "MAM02487"})  # reduced thioredoxin, cytosolic and mitochondrial forms
# the oxidants whose reduction the glutathione pool supports: hydrogen peroxide, the lipid hydroperoxides Human-GEM
# names (HPETE and HPODE species, 4-hydroperoxy-2-nonenal), dehydroascorbate and the oxidised thioredoxins
PEROXIDE_SUBSTRATES = frozenset({"MAM02041", "MAM00983", "MAM01655", "MAM02666", "MAM02486"})
PEROXIDE_NAME_WORDS = ("hydroperoxy", "peroxide", "hpete", "hpode", "hydroperoxide")
POOL_COMPARTMENTS = ("c", "m", "r", "n", "x")


def pool_node(pool: str, compartment: str) -> str:
    return f"REDOX_{pool}_{compartment}"


def reaction_participants(edges: pd.DataFrame, base_of: dict, compartment_of: dict) -> tuple[dict, dict]:
    """({reaction: {(base metabolite, compartment)}} consumed, the same produced)."""
    consumed, produced = defaultdict(set), defaultdict(set)
    for source, target, relation in zip(edges.source_id, edges.target_id, edges.relation_type):
        if relation == "substrate_of" and base_of.get(source):
            consumed[target].add((base_of[source], compartment_of.get(source)))
        elif relation == "product_of" and base_of.get(target):
            produced[source].add((base_of[target], compartment_of.get(target)))
    return consumed, produced


def add_redox_pools(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The graph with the glutathione and NADPH pools of each compartment and their edges."""
    base_of = dict(zip(nodes.node_id, nodes.base_metabolite_id))
    compartment_of = dict(zip(nodes.node_id, nodes.compartment.astype(str)))
    name_of = dict(zip(nodes.node_id, nodes.display_name.astype(str)))
    consumed, produced = reaction_participants(edges, base_of, compartment_of)
    blank = {column: None for column in nodes.columns}
    new_rows, edge_rows = [], []
    counts: dict[str, int] = defaultdict(int)

    def handles_an_oxidant(reaction: str, compartment: str) -> bool:
        if any(base in PEROXIDE_SUBSTRATES for base, letter in consumed[reaction] if letter == compartment):
            return True
        return any(word in name_of.get(reaction, "").lower() for word in PEROXIDE_NAME_WORDS)

    for compartment in POOL_COMPARTMENTS:
        glutathione_pool, nadph_pool = pool_node("GSH", compartment), pool_node("NADPH", compartment)
        glutathione_rows, nadph_rows = [], []
        for reaction in set(consumed) | set(produced):
            takes = {base for base, letter in consumed[reaction] if letter == compartment}
            makes = {base for base, letter in produced[reaction] if letter == compartment}
            if GLUTATHIONE in takes and GLUTATHIONE not in makes:
                glutathione_rows.append((reaction, -1.0, "glutathione consumed"))
                if handles_an_oxidant(reaction, compartment):
                    glutathione_rows.append((reaction, 0.0, "oxidant reduced by glutathione"))
            elif GLUTATHIONE in makes and GLUTATHIONE not in takes:
                glutathione_rows.append((reaction, 1.0, "glutathione produced"))
            antioxidant_demand = NADPH in takes and (GLUTATHIONE in makes or makes & THIOREDOXINS)
            if NADP in takes and NADPH in makes:
                nadph_rows.append((reaction, 1.0, "NADP+ reduced"))
            elif antioxidant_demand:
                nadph_rows.append((reaction, -1.0, "NADPH spent on reducing glutathione or thioredoxin"))
                nadph_rows.append((reaction, 0.0, "antioxidant demand for NADPH"))
        for pool, rows, name in ((glutathione_pool, glutathione_rows, f"reduced glutathione pool ({compartment})"),
                                 (nadph_pool, nadph_rows, f"NADPH pool ({compartment})")):
            if not rows:
                continue
            new_rows.append({**blank, "node_id": pool, "node_type": "metabolite", "display_name": name, "compartment": compartment,
                             "base_metabolite_id": GLUTATHIONE if "glutathione" in name else NADPH, "is_currency": False, "in_metabolic_layer": True, "degree": 0})
            for reaction, sign, evidence in rows:
                if sign:
                    edge_rows.append({"source_id": reaction, "target_id": pool, "relation_type": POOL_CHANGE_RELATION, "sign": sign, "evidence_source": evidence})
                    counts[evidence] += 1
                else:  # the pool's own feedback on the reactions that draw on it
                    edge_rows.append({"source_id": pool, "target_id": reaction, "relation_type": POOL_CAPACITY_RELATION, "sign": 1.0, "evidence_source": evidence})
                    counts[evidence] += 1

    new_nodes = pd.concat([nodes, pd.DataFrame(new_rows, columns=nodes.columns)], ignore_index=True)
    new_edges = pd.concat([edges, pd.DataFrame(edge_rows, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict())
    new_edges = new_edges.drop_duplicates(["source_id", "target_id", "relation_type"], ignore_index=True)
    degrees = new_edges.source_id.value_counts().add(new_edges.target_id.value_counts(), fill_value=0)
    new_nodes["degree"] = degrees.reindex(new_nodes.node_id).fillna(0).astype("int64").to_numpy()
    new_relations = list(relation_types) + [relation for relation in (POOL_CHANGE_RELATION, POOL_CAPACITY_RELATION) if relation not in relation_types]
    summary = {"pools": [row["node_id"] for row in new_rows], "edges_by_role": dict(sorted(counts.items())),
               "largest_pool_out_degree": int(max((new_edges.source_id == row["node_id"]).sum() for row in new_rows)) if new_rows else 0}
    return new_nodes, new_edges, new_relations, summary
