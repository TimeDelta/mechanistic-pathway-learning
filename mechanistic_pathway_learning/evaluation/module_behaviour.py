"""Name a learned module from how it behaves, not only from what its support holds.

A support-side name says which part of the reconstruction a module's gates ended up on. It is silent about the two
things a reader of a mechanism card wants: which perturbations switch the module on, and which symptoms it carries.
Both are already in a finished run, so reading them costs no training (the user, 9 October 2026: "is there an
efficient way to derive learned module names from the nodes, edges and training examples themselves ... i.e. this
module seems to encode for movement-related symptoms or encodes for effects on protein complexes when a knockout of
x,y or z happens"):

- the **drivers** are the held-out perturbations in the top fifth of the module's activation
  (`test_module_activations.npy`), named by their own labels and by any group they share;
- the **reach** is the share of those drivers whose seed nodes have a directed path into the support no longer than
  the propagation the encoder actually runs, which is what makes a driver name mechanistic rather than a correlation;
- the **symptom side** is the module's link vector (`module_symptom_links`) read against the other modules' and then
  grouped by the RDoC domain the crosswalk already assigns each symptom.

Every reading refuses to name rather than name weakly, and each refusal carries the number behind it. Two refusals are
expected rather than exceptional, because both are already recorded failures of the module layer: a module whose
activation does not vary across perturbations has no drivers (docs/module_health.md: "Modules that add the same amount
to every perturbation cannot be read as pathways"), and a module whose link vector matches the other modules' carries
no symptom family of its own (the same document's `median_link_correlation`, whose target is at most 0.7).
"""
from __future__ import annotations

from collections import Counter, deque

import numpy as np
import pandas as pd

# the reading of experiments/module_health.py, kept at the same value so a module that is "responding" there has
# drivers here
RESPONDING_STANDARD_DEVIATION = 0.05
DRIVER_ACTIVATION_QUANTILE = 0.8
DRIVERS_NAMED_IN_A_LABEL = 3
DRIVER_GROUP_SHARE_FOR_A_NAME = 0.5
# a group that holds half the drivers but half the scored perturbations too says nothing about the module
DRIVER_GROUP_ENRICHMENT_FOR_A_NAME = 2.0
CURATED_GENES_FOR_A_DRIVER_NAME = 3
# docs/module_health.md gives 0.7 as the median link correlation a module set read as pathways has to be under
LINK_CORRELATION_FOR_A_DISTINCT_PROFILE = 0.7
# and the module's own link for a symptom has to stand this far above what the other modules give that symptom, so a
# name does not rest on a difference that standardising magnified
SYMPTOM_LINK_MARGIN_FOR_A_NAME = 0.2
DOMAIN_SHARE_FOR_A_NAME = 0.5
DOMAIN_ENRICHMENT_FOR_A_NAME = 1.5
# and it has to hold two of the module's symptoms, so a domain name never rests on a single symptom
SYMPTOMS_IN_A_DOMAIN_NAME = 2
# a driver name is mechanistic only if the drivers reach the support more often than the perturbations that are not
# drivers; this is the ratio of the two shares a name needs
REACH_ENRICHMENT_FOR_A_MECHANISTIC_NAME = 1.5
DRIVERS_RECORDED = 50
# the backward walk from a support is bounded, because a support holding a hub reaction reaches much of the graph
REACH_VISITED_NODE_CAP = 20000


def reverse_adjacency(edge_source: np.ndarray, edge_target: np.ndarray, excluded_nodes: set[str]) -> dict[str, list[str]]:
    """Target -> sources, for walking backwards from a support to the perturbations that can reach it."""
    reverse: dict[str, list[str]] = {}
    for source, target in zip(edge_source, edge_target):
        if source in excluded_nodes or target in excluded_nodes:
            continue
        reverse.setdefault(target, []).append(source)
    return reverse


def hops_into_support(support_node_ids: set[str], reverse_edges: dict[str, list[str]], max_hops: int) -> tuple[dict[str, int], bool]:
    """For each node within `max_hops` of the support along a directed path, the length of the shortest such path.

    Returns the distances and whether the walk stopped at REACH_VISITED_NODE_CAP, in which case a node absent from the
    distances may still reach the support and the reported reach is a lower bound.
    """
    distance = {node_id: 0 for node_id in support_node_ids}
    frontier = deque((node_id, 0) for node_id in support_node_ids)
    truncated = False
    while frontier:
        node_id, hops = frontier.popleft()
        if hops >= max_hops:
            continue
        for source in reverse_edges.get(node_id, ()):
            if source in distance:
                continue
            if len(distance) >= REACH_VISITED_NODE_CAP:
                truncated = True
                frontier.clear()
                break
            distance[source] = hops + 1
            frontier.append((source, hops + 1))
    return distance, truncated


def propagation_hops(recorded_arguments: dict) -> int:
    """How far the encoder the run used actually propagates, which bounds a driver's reach into a support."""
    if recorded_arguments.get("encoder") == "linear_response":
        return int(recorded_arguments.get("propagation_steps", 0) or 0)
    return int(recorded_arguments.get("num_layers", 0) or 0)


def symptom_rdoc_domains(crosswalk: pd.DataFrame) -> dict[str, list[str]]:
    """The crosswalk's own domain per symptom, which is where a symptom family comes from rather than a new grouping."""
    domains = {}
    for row in crosswalk.itertuples(index=False):
        listed = [part.strip() for part in str(getattr(row, "rdoc_domain", "") or "").split(";") if part.strip()]
        domains[row.target_symptom] = listed
    return domains


def domain_shares(symptoms: list[str], domains: dict[str, list[str]]) -> dict[str, float]:
    """Each domain's share of a set of symptoms, a symptom with two domains counting half to each."""
    shares: dict[str, float] = {}
    total = 0.0
    for symptom in symptoms:
        listed = domains.get(symptom) or []
        if not listed:
            continue
        for domain in listed:
            shares[domain] = shares.get(domain, 0.0) + 1.0 / len(listed)
        total += 1.0
    return {domain: count / total for domain, count in shares.items()} if total else {}


def module_drivers(module_activations: np.ndarray, perturbation_labels: list[str], perturbation_types: list[str],
                   curated_module_genes: dict[str, set[str]] | None = None,
                   hops_of_perturbation: dict[int, int] | None = None, hop_budget: int | None = None) -> dict:
    """The held-out perturbations that switch one module on, any group they share and how far they are from its support."""
    standard_deviation = float(np.std(module_activations))
    reading: dict = {"activation_standard_deviation": round(standard_deviation, 4),
                     "responds": standard_deviation > RESPONDING_STANDARD_DEVIATION,
                     "num_perturbations": int(len(module_activations))}
    if not reading["responds"]:
        return reading
    threshold = float(np.quantile(module_activations, DRIVER_ACTIVATION_QUANTILE))
    positions = [int(position) for position in np.argsort(-module_activations) if module_activations[position] >= threshold]
    reading["num_drivers"] = len(positions)
    reading["top_drivers"] = [perturbation_labels[position] for position in positions[:DRIVERS_NAMED_IN_A_LABEL]]
    driver_types = Counter(perturbation_types[position] for position in positions)
    background_types = Counter(perturbation_types)
    reading["by_perturbation_type"] = dict(driver_types)
    reading["group"] = None
    driver_genes = {perturbation_labels[position] for position in positions
                    if perturbation_types[position] == "gene"}
    scored_genes = {label for label, kind in zip(perturbation_labels, perturbation_types) if kind == "gene"}
    for module, genes in sorted((curated_module_genes or {}).items()):
        shared = driver_genes & genes
        # the same enrichment test the support side needs: a curated module holding three drivers says nothing if it
        # holds as large a share of the perturbations that are not drivers
        background = len(scored_genes & genes) / len(perturbation_labels) if perturbation_labels else 0.0
        share = len(shared) / len(positions) if positions else 0.0
        if len(shared) >= CURATED_GENES_FOR_A_DRIVER_NAME and share >= DRIVER_GROUP_ENRICHMENT_FOR_A_NAME * background:
            reading["group"] = {"kind": "curated_module", "name": module, "drivers": len(shared), "of": len(positions),
                                "share": round(share, 2), "background": round(background, 3)}
            break
    if reading["group"] is None and positions:
        perturbation_type, count = driver_types.most_common(1)[0]
        share = count / len(positions)
        background = background_types[perturbation_type] / len(perturbation_types)
        if share >= DRIVER_GROUP_SHARE_FOR_A_NAME and share >= DRIVER_GROUP_ENRICHMENT_FOR_A_NAME * background:
            reading["group"] = {"kind": "perturbation_type", "name": f"{perturbation_type} perturbations",
                                "drivers": count, "of": len(positions)}
    reading["drivers"] = [perturbation_labels[position] for position in positions[:DRIVERS_RECORDED]]
    if hops_of_perturbation is not None:
        others = [position for position in range(len(module_activations)) if position not in set(positions)]
        reached = [hops_of_perturbation[position] for position in positions if position in hops_of_perturbation]
        reached_by_others = [position for position in others if position in hops_of_perturbation]
        driver_share = len(reached) / len(positions) if positions else 0.0
        other_share = len(reached_by_others) / len(others) if others else 0.0
        reading.update({
            "drivers_reaching_the_support": len(reached),
            "others_reaching_the_support": len(reached_by_others), "num_others": len(others),
            "median_hops_into_the_support": float(np.median(reached)) if reached else None,
            "hop_budget": hop_budget,
            "reach_is_enriched": bool(other_share == 0.0 and driver_share > 0.0
                                      or other_share > 0.0 and driver_share >= REACH_ENRICHMENT_FOR_A_MECHANISTIC_NAME * other_share),
        })
    return reading


def module_symptom_side(links: np.ndarray, module_index: int, symptoms: list[str], domains: dict[str, list[str]]) -> dict:
    """The symptoms one module feeds more than the other modules do, and the domain they fall in."""
    profile = links[module_index]
    others = np.delete(links, module_index, axis=0)
    reading: dict = {"distinct_profile": True, "link_correlation_with_the_other_modules": None,
                     "symptoms": [], "domain": None}
    if len(others):
        other_mean = others.mean(axis=0)
        if profile.std() > 0 and other_mean.std() > 0:
            correlation = float(np.corrcoef(profile, other_mean)[0, 1])
            reading["link_correlation_with_the_other_modules"] = round(correlation, 3)
            reading["distinct_profile"] = correlation < LINK_CORRELATION_FOR_A_DISTINCT_PROFILE
        above = [symptoms[s] for s in range(len(symptoms))
                 if other_mean[s] > 0 and profile[s] >= (1.0 + SYMPTOM_LINK_MARGIN_FOR_A_NAME) * other_mean[s]]
    else:
        above = [symptoms[s] for s in np.argsort(-profile)[:1]]
    reading["symptoms"] = above
    if not reading["distinct_profile"] or not above:
        return reading
    shares = domain_shares(above, domains)
    background = domain_shares(symptoms, domains)
    if shares:
        domain, share = max(shares.items(), key=lambda item: item[1])
        # the share is weighted, a symptom with two domains counting half to each; the count is not, because "two of
        # its symptoms are sensorimotor" is the claim a reader checks
        in_the_domain = sum(1 for symptom in above if domain in (domains.get(symptom) or []))
        if (in_the_domain >= SYMPTOMS_IN_A_DOMAIN_NAME and share >= DOMAIN_SHARE_FOR_A_NAME
                and share >= DOMAIN_ENRICHMENT_FOR_A_NAME * background.get(domain, 0.0)):
            reading["domain"] = {"name": domain, "share": round(share, 2), "of_symptoms": len(above),
                                 "symptoms_in_the_domain": in_the_domain}
    return reading


def behaviour_label(behaviour: dict) -> str:
    """One readable clause per reading, each carrying the number it rests on, or the reason there is no name."""
    if not behaviour.get("available", True):
        return behaviour.get("reason", "no behaviour recorded")
    drivers, symptom_side = behaviour.get("drivers", {}), behaviour.get("symptom_side", {})
    parts = []
    if not drivers.get("responds", False):
        parts.append("no driver: its activation varies by "
                     f"{drivers.get('activation_standard_deviation', float('nan')):.3f} over the held-out perturbations")
    else:
        named, number = ", ".join(drivers.get("top_drivers", [])), drivers.get("num_drivers", 0)
        group = drivers.get("group")
        head = (f"driven by {group['name']} ({group['drivers']} of {group['of']} drivers, {named} highest)" if group
                else f"driven by {named} and {max(number - len(drivers.get('top_drivers', [])), 0)} others with no group in common")
        reached, budget = drivers.get("drivers_reaching_the_support"), drivers.get("hop_budget")
        if reached is not None and budget:
            head += (f"; {reached} of {number} reach the support within {budget} hops against "
                     f"{drivers.get('others_reaching_the_support', 0)} of {drivers.get('num_others', 0)} others"
                     + ("" if drivers.get("reach_is_enriched") else ", so the activation does not follow the path"))
        parts.append(head)
    domain = symptom_side.get("domain")
    if domain:
        parts.append(f"feeds {domain['name']} symptoms ({domain['symptoms_in_the_domain']:g} of its "
                     f"{domain['of_symptoms']} symptoms)")
    elif not symptom_side.get("distinct_profile", True):
        correlation = symptom_side.get("link_correlation_with_the_other_modules")
        parts.append(f"no symptom family of its own: its link profile matches the other modules (correlation {correlation})")
    elif not symptom_side.get("symptoms"):
        parts.append(f"no symptom above {1 + SYMPTOM_LINK_MARGIN_FOR_A_NAME:.1f} times the other modules' link")
    else:
        parts.append("symptoms " + ", ".join(symptom_side["symptoms"][:3]) + ", in no one domain")
    return "; ".join(parts)
