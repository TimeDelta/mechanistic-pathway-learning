"""Minimum description length selection of the module count (design section 5.5).

Total code length for a fitted model = model cost + data cost, where the model
cost is node_cost * active support nodes + link_cost * nonzero links (from
NoisyOrPathwayModuleHead.description_length_penalty) and the data cost is the
weighted negative log-likelihood on the training observations. Candidate values
of K are trained separately and the smallest total wins. Costs are in nats when
the likelihood is in nats; a node cost of log2(num_graph_nodes) bits is the
natural choice (the cost of naming a node), converted to nats.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class CandidateDescriptionLength:
    num_pathway_modules: int
    model_cost_nats: float
    data_cost_nats: float

    @property
    def total_nats(self) -> float:
        return self.model_cost_nats + self.data_cost_nats


def node_naming_cost_nats(num_graph_nodes: int) -> float:
    return math.log(num_graph_nodes)


def select_minimum_description_length(candidates: list[CandidateDescriptionLength]) -> CandidateDescriptionLength:
    if not candidates:
        raise ValueError("no candidates")
    return min(candidates, key=lambda candidate: candidate.total_nats)
