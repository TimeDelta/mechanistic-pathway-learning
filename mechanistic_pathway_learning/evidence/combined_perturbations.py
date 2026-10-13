"""The combined perturbation type: one labelled example that carries two or more perturbations at once.

What it is for. Rows 4 to 6 of section 8 of docs/label_source_review.md name three label kinds that no single
perturbation can express: a drug given with the enzyme that clears it less active, a drug given with its plasma
carrier lower, and two drugs given together. Section 8 also states the property that makes them safe to add after a
held-out set is drawn: "a combination inherits the groups of its parts". This module is that type, and the loader and
the trainer read it through this module alone.

How it is written. A combination is an ordinary row of an evidence table with

  perturbation_type   "combination"
  perturbation_id     the part ids joined by "+", in the order the parts are given (combined_perturbation_id)
  perturbation_nodes  the union of the parts' triples, so a reader that knows nothing of parts still seeds it correctly
  perturbation_parts  this module's JSON: one entry per part, with the part's id, its kind and its own triples

The parts column is what the union of triples cannot replace. Three things are read per part and not per row:

  - on a gene/protein split graph a drug acts on proteins and a gene perturbation on its gene node
    (evaluation/experiment_data.py, drug_seeds_on_proteins), so the routing follows each part's kind;
  - on a graph with drug entry nodes a drug part is seeded on the drug's own node and its mechanism edges carry its
    sign and magnitude on, while a gene part keeps its own nodes;
  - a leakage group is formed from the parts, which is what makes the inherited group exact rather than inferred from
    which nodes happen to be seeded.

What this module does not decide. Whether a combination enters the loss at all is the owner's decision and is a
command-line choice in experiments/run_main_model.py (--combined-perturbations), not a default.
"""
from __future__ import annotations

import json
from dataclasses import dataclass

COMBINATION_PERTURBATION_TYPE = "combination"
PART_KINDS: tuple[str, ...] = ("drug", "gene")
PART_SEPARATOR = "+"
PARTS_COLUMN = "perturbation_parts"
COMBINED_RECORD_FILE = "combined_perturbation_records.parquet"  # written by experiments/build_combined_perturbations.py
LOWERED_GENE_SIGN = -1.0  # a less active enzyme or a lower carrier: the same sign a loss of function carries
# The row kinds that hold a drug, for the rules written for drugs: the equal loss shares and the drug strata. A
# combination of a gene and a gene would not be one of these, and the builder writes none.
DRUG_ROW_KINDS: tuple[str, ...] = ("drug", COMBINATION_PERTURBATION_TYPE)

# A receptor activation row (preregistration item 12c) is a curated statement that agonism at one receptor causes a
# symptom, seeded +1 at that gene node. It is a perturbation the study has no other way to express, and it is unlike
# the rest in two ways that decide how it is used: there are 23 of them carrying 35 labels, and they carry no
# pharmacology at all - no exposure, occupancy, route or plasma binding - so the row is identifiable from its inputs
# alone while 43 percent of its labels are one symptom (cognitive impairment, against 14.4 percent of all kept
# positives). A model could rank that symptom from the row's type rather than from any mechanism. So the type is
# never scored: see never_scored_rows below for what that means precisely, which is not the same as "always trained".
RECEPTOR_ACTIVATION_PERTURBATION_TYPE = "receptor_activation"
PERTURBATION_TYPES_NEVER_SCORED: tuple[str, ...] = (RECEPTOR_ACTIVATION_PERTURBATION_TYPE,)


@dataclass(frozen=True)
class PerturbationPart:
    """One part of a combination: a drug or a gene, with the node triples it seeds.

    part_id is the id the part would carry as a perturbation of its own (a drug's perturbation id, a gene's symbol), so
    a part that is also a perturbation of the table is recognised by its id and shares that perturbation's leakage
    group. nodes are [node_id, sign, magnitude] triples in the graph, the same form as perturbation_nodes.
    """

    part_id: str
    part_kind: str
    nodes: tuple[tuple[str, float, float], ...]

    def __post_init__(self) -> None:
        if self.part_kind not in PART_KINDS:
            raise ValueError(f"part {self.part_id!r}: kind {self.part_kind!r} is not one of {PART_KINDS}")
        if not self.part_id:
            raise ValueError(f"a part of kind {self.part_kind!r} has no id")
        if PART_SEPARATOR in self.part_id:
            raise ValueError(f"part id {self.part_id!r} holds {PART_SEPARATOR!r}, which joins part ids into a combination id")

    def as_json_entry(self) -> dict:
        return {"part_id": self.part_id, "part_kind": self.part_kind,
                "nodes": [[node_id, float(sign), float(magnitude)] for node_id, sign, magnitude in self.nodes]}


def combined_perturbation_id(part_ids: list[str] | tuple[str, ...]) -> str:
    """The perturbation id of a combination: its part ids joined by PART_SEPARATOR, in the given order.

    The order is the one the builder chose (the drug first for an enzyme or a carrier statement, the first drug first
    for a pair), so the id reads as the statement does and two statements over the same parts in the same order give
    the same id."""
    if len(part_ids) < 2:
        raise ValueError(f"a combination has at least two parts, got {list(part_ids)}")
    return PART_SEPARATOR.join(part_ids)


def parts_json(parts: list[PerturbationPart]) -> str:
    return json.dumps([part.as_json_entry() for part in parts])


def parts_with_disjoint_nodes(parts: list[PerturbationPart]) -> tuple[list[PerturbationPart], list[str]]:
    """(the parts with each node held by one part alone, the node ids two parts had seeded).

    Two parts of a combination can name the same node: two drugs that both inhibit one transporter, or a drug given with
    a lowered gene that is also one of its targets. The node then carries the sum of the signed magnitudes, held by the
    first part that names it, and sign and magnitude are read back off that sum. Summing is the composition the
    linear-response encoder gives a node anyway, since its input is sign x magnitude at the node; the
    message-passing encoder injects a learned function of (sign, magnitude) per seed and sums those, which is a
    different composition, so the sum is taken here, in the data, where it is the same for both encoders and visible.

    A node whose signed magnitudes cancel exactly keeps magnitude 0 and the first part's sign, which is an input of
    nothing at that node rather than an input of unknown direction."""
    signed_total: dict[str, float] = {}
    first_sign: dict[str, float] = {}
    order: list[tuple[str, int]] = []  # node id and the position of the part that holds it
    shared: list[str] = []
    for position, part in enumerate(parts):
        for node_id, sign, magnitude in part.nodes:
            if node_id in signed_total:
                shared.append(node_id)
            else:
                order.append((node_id, position))
                first_sign[node_id] = float(sign)
            signed_total[node_id] = signed_total.get(node_id, 0.0) + float(sign) * float(magnitude)
    held_by_part: dict[int, list[tuple[str, float, float]]] = {position: [] for position in range(len(parts))}
    for node_id, position in order:
        total = signed_total[node_id]
        sign = first_sign[node_id] if total == 0 else (1.0 if total > 0 else -1.0)
        held_by_part[position].append((node_id, sign, abs(total)))
    disjoint = [PerturbationPart(part.part_id, part.part_kind, tuple(held_by_part[position])) for position, part in enumerate(parts)]
    return disjoint, sorted(set(shared))


def nodes_json_of_parts(parts: list[PerturbationPart]) -> str:
    """The parts' triples in part order, each node once (parts_with_disjoint_nodes)."""
    disjoint, _ = parts_with_disjoint_nodes(parts)
    return json.dumps([[node_id, sign, magnitude] for part in disjoint for node_id, sign, magnitude in part.nodes])


def read_parts(perturbation_id: str, perturbation_type: str, perturbation_nodes: str, perturbation_parts: str | None = None) -> list[PerturbationPart]:
    """The parts of a row, for a single perturbation as well as for a combination.

    A row that is not a combination has one part: itself, with its own kind and its own triples. So every caller can
    read parts and none needs to branch on the type. A combination without a parts column is an error rather than a
    guess, because the parts cannot be recovered from the union of the triples: which triples belong to which part, and
    each part's kind, are exactly what the union loses.
    """
    triples = tuple((str(node_id), float(sign), float(magnitude)) for node_id, sign, magnitude in json.loads(perturbation_nodes or "[]"))
    if perturbation_type != COMBINATION_PERTURBATION_TYPE:
        return [PerturbationPart(perturbation_id, perturbation_type, triples)]
    if not perturbation_parts:
        raise ValueError(f"{perturbation_id} is a {COMBINATION_PERTURBATION_TYPE} with no {PARTS_COLUMN}; the parts cannot be read from the union of the triples")
    parts = [PerturbationPart(entry["part_id"], entry["part_kind"],
                              tuple((str(node_id), float(sign), float(magnitude)) for node_id, sign, magnitude in entry["nodes"]))
             for entry in json.loads(perturbation_parts)]
    expected = combined_perturbation_id([part.part_id for part in parts])
    if expected != perturbation_id:
        raise ValueError(f"{perturbation_id}: its parts name {expected}; a combination's id is its part ids joined by {PART_SEPARATOR!r}")
    return parts
