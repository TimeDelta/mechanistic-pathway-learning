"""Comparative checks: a drug given with a second, varied factor, compared with the drug alone.

What a check tests, and what it cannot. Only the two-drug kind reads parameters that a label trains; the three kinds
that go through a lowered carrier or a lowered enzyme do not, because no label of the study varies either. The report
keeps them apart for that reason (KINDS_THAT_READ_A_TRAINED_RELATION below, and the user's objection of 12 October
2026 that it records).

What a check is. A statement from a source that a drug's effects are at least as large when a second factor is present:
its plasma carrier is lower, a second drug displaces it from the carrier, the enzyme that clears it is less active, or
a second drug is given with it. The model is asked the same question: do the symptom probabilities of the combined
input lie at or above those of the reference input? Nothing is trained on a check. It is reported apart from every
score (the user's decisions of 10 and 11 October 2026: "Build, report as checks"; "Add the valproate, etc displacements
as checks"; docs/label_source_review.md, sections 8 and 13).

The four kinds, and the two inputs each compares:

  carrier_level     reference: the drug.            combined: the drug, with its carrier's nodes lowered.
  enzyme_activity   reference: the drug.            combined: the drug, with the enzyme's nodes lowered.
  displacement      reference: the two drugs.       combined: the two drugs, with the carrier's nodes lowered.
  second_drug       reference: the first drug.      combined: the two drugs.

Displacement is written as a lowered carrier because the graph has no edge from a drug to the carrier it occupies. Such
an edge would make every drug with a known carrier move the carrier when taken alone, which is the extra input the
change-carrying edges were built to avoid (docs/drug_entry_nodes_measured.md, section 6). So a displacement check tests
the same edge as a carrier check, with two drugs present; it says so in its kind.

A check needs the drug's own node and, for the first three kinds, the edge from the carrier or the enzyme to that node
(relations sequesters and clears, graph/drug_entry_nodes.py). On a graph without drug nodes, or under the loader's
"targets" mode, lowering a carrier does not reach the drug, and those checks are reported as not expressible.

A combined input belongs to the leakage groups of all its parts. A check is "held out" when one of its drugs, or its
carrier or enzyme where that gene is itself a perturbation of the study, is held out, so it never counts for a model
that saw part of it only on the other side of the split.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.drug_entry_nodes import CARRIER_GENES, CLEARS_RELATION, SEQUESTERS_RELATION

CHECK_KINDS: tuple[str, ...] = ("carrier_level", "enzyme_activity", "displacement", "second_drug")
# Which kinds read parameters that a label trains, and so can be read as a test of a model (the user, 12 October 2026:
# the checks "are almost certain to fail without having trained any enzymes so I don't see the point in checking them").
# A second drug enters through its mechanism edges, the relation every drug perturbation trains. A lowered carrier or
# enzyme reaches the drug through sequesters or clears, and no label of the study varies either, so:
#   - in the message-passing encoder those relation weights get no gradient, and the check reads untrained weights;
#   - in the linear-response encoder the gain of a signed relation is positive by construction and the edge carries the
#     sign, so the direction of the field at the drug's node is built in and the check reads the head's monotonicity.
# Either way the result says nothing about carriage or clearance, so those kinds are reported as a description of the
# model's behaviour and not as a test it passes or fails.
KINDS_THAT_READ_A_TRAINED_RELATION: tuple[str, ...] = ("second_drug",)
EXPECTATIONS: tuple[str, ...] = ("at_least_as_high", "none")  # "none": a negative control, reported and not judged
CHECK_COLUMNS: tuple[str, ...] = ("statement_id", "kind", "drug_perturbation_id", "drug", "lowered_genes", "second_drug_perturbation_id", "second_drug", "symptom",
                                  "expectation", "source", "source_id", "source_version", "quote", "second_mechanism_or_limit")
LOWERED_SIGN = -1.0
DEFAULT_TOLERANCE = 1e-6


@dataclass(frozen=True)
class ComparativeCheck:
    statement_id: str
    kind: str
    drug_perturbation_id: str
    lowered_genes: tuple[str, ...] = ()  # the carrier's or the enzyme's gene symbols; empty for second_drug
    second_drug_perturbation_id: str = ""  # the displacing drug, or the second drug of a pair
    symptom: str = ""  # one study symptom, or "" for the drug's own kept positive symptoms
    expectation: str = "at_least_as_high"

    def __post_init__(self) -> None:
        if self.kind not in CHECK_KINDS:
            raise ValueError(f"{self.statement_id}: kind {self.kind!r} is not one of {CHECK_KINDS}")
        if self.expectation not in EXPECTATIONS:
            raise ValueError(f"{self.statement_id}: expectation {self.expectation!r} is not one of {EXPECTATIONS}")
        if self.kind in ("carrier_level", "enzyme_activity", "displacement") and not self.lowered_genes:
            raise ValueError(f"{self.statement_id}: a {self.kind} check names the gene that is lowered")
        if self.kind in ("displacement", "second_drug") and not self.second_drug_perturbation_id:
            raise ValueError(f"{self.statement_id}: a {self.kind} check names its second drug")

    @property
    def drug_perturbation_ids(self) -> tuple[str, ...]:
        return (self.drug_perturbation_id, self.second_drug_perturbation_id) if self.second_drug_perturbation_id else (self.drug_perturbation_id,)

    @property
    def reads_a_trained_relation(self) -> bool:
        """True when the check runs through a relation that some label of the study trains (see KINDS_THAT_READ_A_TRAINED_RELATION)."""
        return self.kind in KINDS_THAT_READ_A_TRAINED_RELATION

    @property
    def edge_relation(self) -> str | None:
        """The relation whose edge the check passes through, None for a pair of drugs."""
        return {"carrier_level": SEQUESTERS_RELATION, "displacement": SEQUESTERS_RELATION, "enzyme_activity": CLEARS_RELATION}.get(self.kind)


def carrier_genes(carrier: str) -> tuple[str, ...]:
    """The gene symbols of a carrier as the carrier table names it: albumin, orosomucoid or both."""
    if carrier not in CARRIER_GENES:
        raise ValueError(f"carrier {carrier!r} is not one of {sorted(CARRIER_GENES)}")
    return CARRIER_GENES[carrier]


def load_checks(path: Path) -> list[ComparativeCheck]:
    """The checks of a table written by experiments/build_comparative_checks.py; the provenance columns stay in the table."""
    table = pd.read_csv(path, dtype=str).fillna("")
    return [ComparativeCheck(row.statement_id, row.kind, row.drug_perturbation_id, tuple(gene for gene in row.lowered_genes.split(";") if gene),
                             row.second_drug_perturbation_id, row.symptom, row.expectation) for row in table.itertuples(index=False)]


def nodes_of_gene(data, gene_symbol: str, protein_of_gene: dict[str, list[str]] | None = None) -> list[int]:
    """The node indices that stand for a gene's product: its protein nodes on a split graph, its gene node otherwise."""
    node_ids = (protein_of_gene or {}).get(f"GENE:{gene_symbol}") or [f"GENE:{gene_symbol}"]
    return [data.node_index[node_id] for node_id in node_ids if node_id in data.node_index]


def edge_exists(data, relation: str, source_nodes: list[int], target_node: int) -> bool:
    if relation not in data.relation_types:
        return False
    of_relation = np.asarray(data.edge_relation) == data.relation_types.index(relation)
    return bool((of_relation & np.isin(np.asarray(data.edge_source), source_nodes) & (np.asarray(data.edge_target) == target_node)).any())


def why_not_expressible(check: ComparativeCheck, data, protein_of_gene: dict[str, list[str]] | None = None) -> str:
    """Empty when the check can be put to a model on this data; otherwise the reason it cannot."""
    position_of = {perturbation_id: position for position, perturbation_id in enumerate(data.perturbation_ids)}
    for perturbation_id in check.drug_perturbation_ids:
        if perturbation_id not in position_of:
            return f"{perturbation_id} is not a perturbation of this table"
    if check.symptom and check.symptom not in data.symptoms:
        return f"{check.symptom} is not a symptom of this table"
    if check.edge_relation is None:
        return ""
    seeds = np.asarray(data.perturbation_seeds[position_of[check.drug_perturbation_id]])
    entry_mask = getattr(data, "entry_node_mask", None)
    if entry_mask is None or len(seeds) != 1 or not bool(np.asarray(entry_mask)[seeds[0]]):
        return "the drug has no node of its own here, so a lowered carrier or enzyme does not reach it"
    lowered_nodes = [node for gene in check.lowered_genes for node in nodes_of_gene(data, gene, protein_of_gene)]
    if not lowered_nodes:
        return f"no node for {', '.join(check.lowered_genes)}"
    if not edge_exists(data, check.edge_relation, lowered_nodes, int(seeds[0])):
        return f"no {check.edge_relation} edge from {', '.join(check.lowered_genes)} to the drug's node"
    return ""


def seeds_of(data, perturbation_id: str) -> list[tuple[int, float, float]]:
    position = data.perturbation_ids.index(perturbation_id)
    return [(int(node), float(sign), float(magnitude))
            for node, sign, magnitude in zip(data.perturbation_seeds[position], data.perturbation_signs[position], data.perturbation_magnitudes[position])]


def check_inputs(check: ComparativeCheck, data, protein_of_gene: dict[str, list[str]] | None = None) -> tuple[list[tuple[int, float, float]], list[tuple[int, float, float]]]:
    """(reference seeds, combined seeds) of a check, each a list of (node index, sign, magnitude).

    A lowered gene is seeded with sign -1 and a total magnitude of 1 shared among its nodes, the way a drug's magnitude
    is shared among the components of one target."""
    drug = seeds_of(data, check.drug_perturbation_id)
    second_drug = seeds_of(data, check.second_drug_perturbation_id) if check.second_drug_perturbation_id else []
    lowered = []
    for gene in check.lowered_genes:
        gene_nodes = nodes_of_gene(data, gene, protein_of_gene)
        lowered += [(node, LOWERED_SIGN, 1.0 / len(gene_nodes)) for node in gene_nodes]
    if check.kind == "second_drug":
        return drug, drug + second_drug
    if check.kind == "displacement":
        return drug + second_drug, drug + second_drug + lowered
    return drug, drug + lowered


def symptoms_read(check: ComparativeCheck, data) -> list[int]:
    """The symptom columns a check reads: the one it names, or the drug's kept positive symptoms."""
    if check.symptom:
        return [data.symptoms.index(check.symptom)]
    position = data.perturbation_ids.index(check.drug_perturbation_id)
    kept = np.asarray(data.outcomes[position]) > 0
    if data.label_mask is not None:
        kept &= np.asarray(data.label_mask[position])
    return [int(column) for column in np.flatnonzero(kept)]


def side_of(check: ComparativeCheck, held_out_perturbation_ids: set[str]) -> str:
    """"held out" when one of the check's parts is held out, "development" otherwise (the union rule for a combination).

    The parts are its drugs and its lowered genes: a gene perturbation of the study carries its gene symbol as its id,
    so a carrier or an enzyme that is itself a held-out perturbation puts the check on the held-out side."""
    parts = (*check.drug_perturbation_ids, *check.lowered_genes)
    return "held out" if any(part in held_out_perturbation_ids for part in parts) else "development"


def evaluate_checks(checks: list[ComparativeCheck], data, predict: Callable[[list[list[tuple[int, float, float]]], list[int]], np.ndarray],
                    protein_of_gene: dict[str, list[str]] | None = None, tolerance: float = DEFAULT_TOLERANCE) -> pd.DataFrame:
    """One row per check with the difference between the combined and the reference input on the symptoms it reads.

    predict(seed lists, perturbation positions) returns the symptom probabilities of each seed list, [number of lists,
    number of symptoms]; the position is that of the check's drug in data, for a head that reads a covariate of the
    perturbation, so that the two inputs of a check differ in their seeds alone.

    holds is True when no symptom read falls by more than the tolerance, and None for a negative control, a check that
    cannot be expressed, a check whose drug has no symptom to read, and every kind that does not read a trained
    relation: there the difference is reported and nothing is judged."""
    rows, seed_lists, positions, pending = [], [], [], []
    for check in checks:
        reason = why_not_expressible(check, data, protein_of_gene)
        columns = [] if reason else symptoms_read(check, data)
        if not reason and not columns:
            reason = "the drug has no kept positive symptom to read"
        row = {"statement_id": check.statement_id, "kind": check.kind, "expectation": check.expectation, "expressible": not reason, "reason": reason,
               "reads_a_trained_relation": check.reads_a_trained_relation,
               "symptoms_read": len(columns), "mean_difference": np.nan, "smallest_difference": np.nan, "share_not_lower": np.nan, "holds": None}
        rows.append(row)
        if reason:
            continue
        reference, combined = check_inputs(check, data, protein_of_gene)
        position = data.perturbation_ids.index(check.drug_perturbation_id)
        pending.append((row, columns, len(seed_lists)))
        seed_lists += [reference, combined]
        positions += [position, position]
    if pending:
        probabilities = np.asarray(predict(seed_lists, positions))
        for row, columns, start in pending:
            difference = probabilities[start + 1, columns] - probabilities[start, columns]
            row.update(mean_difference=float(difference.mean()), smallest_difference=float(difference.min()), share_not_lower=float((difference >= -tolerance).mean()))
            # holds is read only where the comparison tests something a label trained; elsewhere the difference is reported alone
            if row["expectation"] == "at_least_as_high" and row["reads_a_trained_relation"]:
                row["holds"] = bool((difference >= -tolerance).all())
    return pd.DataFrame(rows)
