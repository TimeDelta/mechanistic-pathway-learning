"""Currency metabolite tagging (experiment design, assumption A9).

Currency metabolites (ATP, NADH, water, protons and similar) connect to a large
share of all reactions and make any path- or propagation-based feature a
measure of hub proximity rather than mechanism. They are tagged by a name list
plus a degree rule, removed from path-based features and message passing, and
kept for flux computation where stoichiometry needs them.

Identifiers differ by reconstruction: Recon3D uses names with a compartment
suffix such as "atp[c]" or "atp_c"; Human-GEM uses opaque identifiers such as
"MAM01371c" with names in a separate table. The functions here take a mapping
from identifier to a display name so either style works.

Three defects found on 13 October 2026, when the owner asked what else the currency mask holds back, and fixed here.
A graph built before this date carries them; the v4 rebuild is what clears them.

1. **A charge suffix defeated the name list.** It held `nad`, `nadp`, `hco3` and `glutathione` while the display names
   are `NAD+`, `NADP+`, `HCO3-` and `GSH`, so no copy of those matched by name and only the high-degree copies were
   caught by the degree rule. Cytosolic NAD+ sent nothing while peroxisomal NAD+, with 135 edges, sent freely. A name
   is now matched with its charge suffix stripped as well as whole, so `NAD+`, `Ca2+` and `HCO3-` all reach the list.
2. **The degree rule tagged things that are not cofactors.** Glycine's cytosolic copy was tagged at 171 edges although
   glycine is an inhibitory neurotransmitter and the NMDA co-agonist, and every other transmitter in the graph
   propagates: glutamate at 296 edges, glutamine at 258, aspartate at 107, dopamine, serotonin, acetylcholine,
   noradrenaline, histamine, adenosine, taurine and D-serine. The owner, 13 October 2026: "glycine can be removed as
   currency as long as there's still a mechanism for the amount of glycine reaching a receptor to vary ... it should
   also be consistent with other NTs". NEVER_CURRENCY_NAMES is that consistency: a transmitter is never tagged, by
   either rule, so the amount of it reaching a receptor varies the way every other transmitter's does.
3. **Nodes added after the tagging ran were never tagged.** `tag_currency` ran on the base graph and
   `reactome_import.py` then added nodes with `is_currency=False` written in, so thirteen species had copies that
   disagreed: H2O in the vacuole propagated while its other nine copies did not, and likewise O2, Na+ in the
   intermembrane space and mitochondrion, Ca2+ in five compartments, Mg2+ and FMN. `retag_currency_by_name` re-applies
   the name rule to every metabolite once the graph is complete, which a build calls last. The degree rule is not
   re-applied, because a degree is only meaningful against the whole graph's distribution and the base graph's is the
   one the threshold was chosen on.

A fourth, found on 17 October 2026 by experiments/measure_perturbation_reach.py: the fix of defect 3 matched the name
list against every node, not only metabolites, so the gene of carbonic anhydrase 2, displayed "CA2", matched the list's
calcium ("ca2") and GENE:CA2 was tagged on the v4 and v5 graphs. A currency node does not propagate in the
linear-response encoder, so the CA2 perturbation reached no node past its gene and its transcriptional regulators did
not reach its protein. The name rule now tags a node only when node_types names it a metabolite, and removes a tag any
other node carries; load_experiment_data reads a graph built before this the same way (currency_flags).
"""
from __future__ import annotations

import re
from collections.abc import Mapping

DEFAULT_CURRENCY_METABOLITE_NAMES: frozenset[str] = frozenset(
    {
        "atp", "adp", "amp", "gtp", "gdp", "gmp", "utp", "udp", "ump", "ctp", "cdp", "cmp",
        "nad", "nadh", "nadp", "nadph", "fad", "fadh2", "fmn", "fmnh2",
        "coa", "acetyl-coa", "accoa",
        "h2o", "h", "h+", "proton", "o2", "co2", "hco3", "h2o2",
        "pi", "ppi", "phosphate", "diphosphate",
        "nh3", "nh4", "ammonia", "ammonium",
        "na1", "na+", "k", "k+", "cl", "cl-", "ca2", "ca2+", "mg2", "mg2+", "fe2", "fe3", "zn2",
        "so4", "sulfate", "h2s",
        "sam", "amet", "ahcys", "sah",
        # "gsh" and "gssg" are the display names Human-GEM uses; they are not charge variants of "glutathione", so
        # stripping a charge does not reach them and they are named (defect 1)
        "thf", "glutathione", "gthrd", "gthox", "gsh", "gssg",
    }
)

# Never tagged currency, by either rule: a neurotransmitter or signalling molecule whose amount reaching a receptor is
# the mechanism, not a background cofactor. Glycine is the one the degree rule caught; the rest are here so the rule
# cannot catch them later, and because every one of them already propagates (defect 2 above).
NEVER_CURRENCY_NAMES: frozenset[str] = frozenset(
    {
        "glycine", "gly", "d-serine", "l-glutamate", "glutamate", "glu", "l-aspartate", "aspartate",
        "gaba", "4-aminobutanoate", "4-aminobutyrate", "dopamine", "serotonin", "5-hydroxytryptamine",
        "acetylcholine", "noradrenaline", "norepinephrine", "adrenaline", "epinephrine", "histamine",
        "adenosine", "taurine", "glutamine", "gln", "l-glutamine", "melatonin", "anandamide",
    }
)
# Node identifiers a build created in order to propagate, which no rule may tag. The neuronal ion pools
# (reactome_import.py) are the case: the membrane potential reads them, the cytosolic sodium pool carries 369 edges
# and so sits above the degree rule's 0.995 quantile, and tagging it would stop it reaching VM_c and silence the
# electrical layer. Until defect 3 was fixed these nodes were protected only by accident, because they are added after
# the tagging ran and so were written is_currency=False whatever the rules said.
NEVER_CURRENCY_IDENTIFIER_PREFIXES: tuple[str, ...] = ("NEURON_",)
# a trailing charge, so "NAD+" also matches the list's "nad" and "HCO3-" its "hco3" (defect 1 above)
CHARGE_SUFFIX_PATTERN = re.compile(r"[+-]+$")
COMPARTMENT_SUFFIX_PATTERN = re.compile(r"(\[[a-z]\]|_[a-z])$")


def strip_compartment_suffix(metabolite_identifier: str) -> str:
    """'atp[c]' -> 'atp', 'atp_c' -> 'atp'; identifiers without a suffix are returned unchanged."""
    return COMPARTMENT_SUFFIX_PATTERN.sub("", metabolite_identifier)


def tag_currency_metabolites(
    metabolite_display_name_by_identifier: Mapping[str, str],
    metabolite_degree_by_identifier: Mapping[str, int],
    currency_metabolite_names: frozenset[str] = DEFAULT_CURRENCY_METABOLITE_NAMES,
    degree_quantile_threshold: float = 0.995,
) -> set[str]:
    """Return the identifiers to exclude from path-based features.

    A metabolite is tagged when its normalized name is on the currency list, or when
    its degree is at or above the given quantile of all metabolite degrees. The
    degree rule catches reconstruction-specific currency species the list misses.
    """
    never = {identifier for identifier, display_name in metabolite_display_name_by_identifier.items()
             if never_currency(identifier, display_name)}
    tagged: set[str] = set()
    for metabolite_identifier, display_name in metabolite_display_name_by_identifier.items():
        if metabolite_identifier in never:
            continue
        if (normalised_names(display_name) | normalised_names(metabolite_identifier)) & currency_metabolite_names:
            tagged.add(metabolite_identifier)
    if metabolite_degree_by_identifier:
        sorted_degrees = sorted(metabolite_degree_by_identifier.values())
        quantile_position = min(len(sorted_degrees) - 1, int(degree_quantile_threshold * len(sorted_degrees)))
        degree_cutoff = sorted_degrees[quantile_position]
        for metabolite_identifier, degree in metabolite_degree_by_identifier.items():
            if degree >= degree_cutoff and degree > 0 and metabolite_identifier not in never:
                tagged.add(metabolite_identifier)
    return tagged


def never_currency(identifier: str, display_name: str) -> bool:
    """Whether no rule may tag this node: a neurotransmitter by name, or a pool a build made in order to propagate."""
    if str(identifier).startswith(NEVER_CURRENCY_IDENTIFIER_PREFIXES):
        return True
    return bool((normalised_names(display_name) | normalised_names(identifier)) & NEVER_CURRENCY_NAMES)


def normalised_names(value: str) -> set[str]:
    """The forms of a name the lists are matched against: as written, without a compartment suffix, and without a
    trailing charge, each lowercased. "NAD+" gives {"nad+", "nad"}, so the list's "nad" matches it (defect 1)."""
    without_compartment = strip_compartment_suffix(str(value)).strip().lower()
    return {without_compartment, CHARGE_SUFFIX_PATTERN.sub("", without_compartment)}


CURRENCY_NODE_TYPE = "metabolite"  # only a metabolite can be currency (the fourth defect above)


def currency_flags(nodes) -> "np.ndarray":
    """A node table's is_currency, read for metabolites only, as a boolean array in the table's order.

    A graph built before 17 October 2026 tags GENE:CA2 (the fourth defect above); reading it through this drops that
    tag, so a model and an analysis of an older graph see what a rebuilt one carries. An absent column or value is
    False."""
    import numpy as np

    if "is_currency" not in nodes.columns:
        return np.zeros(len(nodes), dtype=bool)
    tagged = nodes["is_currency"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    return tagged & (nodes["node_type"].astype(str).to_numpy() == CURRENCY_NODE_TYPE)


def currency_tags_on_other_node_types(nodes) -> list[str]:
    """The nodes a table tags as currency that are not metabolites, which currency_flags reads as untagged."""
    if "is_currency" not in nodes.columns:
        return []
    tagged = nodes["is_currency"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    other_type = nodes["node_type"].astype(str).to_numpy() != CURRENCY_NODE_TYPE
    return [str(node_id) for node_id in nodes["node_id"].to_numpy()[tagged & other_type]]


def retag_currency_by_name(node_identifiers, display_names, already_tagged, *, node_types,
                           currency_metabolite_names: frozenset[str] = DEFAULT_CURRENCY_METABOLITE_NAMES) -> set[str]:
    """The currency set with the name rule applied to every metabolite given, not only those present when the degree
    rule ran.

    A build calls this once its graph is complete (defect 3). already_tagged is what the degree rule produced; a node
    whose name is on NEVER_CURRENCY_NAMES is removed from it even if the degree rule tagged it, and so is a node that is
    not a metabolite: node_types gives each node's type, and a gene or protein whose symbol spells a listed name (CA2,
    the fourth defect) is never tagged."""
    node_types = list(node_types)
    node_identifiers, display_names = list(node_identifiers), list(display_names)
    if not len(node_types) == len(node_identifiers) == len(display_names):
        raise ValueError("node_identifiers, display_names and node_types must have one entry per node")
    tagged = set(already_tagged)
    for identifier, display_name, node_type in zip(node_identifiers, display_names, node_types):
        forms = normalised_names(display_name) | normalised_names(identifier)
        if node_type != CURRENCY_NODE_TYPE or never_currency(identifier, display_name):
            tagged.discard(identifier)
        elif forms & currency_metabolite_names:
            tagged.add(identifier)
    return tagged
