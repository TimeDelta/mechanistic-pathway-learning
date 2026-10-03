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
        "thf", "glutathione", "gthrd", "gthox",
    }
)

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
    tagged: set[str] = set()
    for metabolite_identifier, display_name in metabolite_display_name_by_identifier.items():
        normalized_name = strip_compartment_suffix(display_name).strip().lower()
        if normalized_name in currency_metabolite_names or strip_compartment_suffix(metabolite_identifier).lower() in currency_metabolite_names:
            tagged.add(metabolite_identifier)
    if metabolite_degree_by_identifier:
        sorted_degrees = sorted(metabolite_degree_by_identifier.values())
        quantile_position = min(len(sorted_degrees) - 1, int(degree_quantile_threshold * len(sorted_degrees)))
        degree_cutoff = sorted_degrees[quantile_position]
        for metabolite_identifier, degree in metabolite_degree_by_identifier.items():
            if degree >= degree_cutoff and degree > 0:
                tagged.add(metabolite_identifier)
    return tagged
