"""The three currency-tag defects of 13 October 2026, and the interaction fixing one of them exposed.

The owner asked what else the currency mask holds back ("are there any other currency nodes in the graph that are
getting masked out?"), which found them, and then directed the glycine half ("glycine can be removed as currency as
long as there's still a mechanism for the amount of glycine reaching a receptor to vary ... it should also be
consistent with other NTs")."""
import pandas as pd
import pytest

from mechanistic_pathway_learning.graph.tag_currency_metabolites import (
    DEFAULT_CURRENCY_METABOLITE_NAMES,
    NEVER_CURRENCY_IDENTIFIER_PREFIXES,
    NEVER_CURRENCY_NAMES,
    never_currency,
    normalised_names,
    retag_currency_by_name,
    tag_currency_metabolites,
)


def test_a_charge_suffix_no_longer_defeats_the_name_list():
    """Defect 1: the list holds `nad`, `nadp` and `hco3` while the display names are NAD+, NADP+ and HCO3-, so no copy
    matched by name and only high-degree copies were caught by the degree rule."""
    for display_name in ("NAD+", "NADP+", "HCO3-", "Ca2+", "Mg2+", "Cl-", "H+", "Na+", "GSH", "GSSG"):
        assert normalised_names(display_name) & DEFAULT_CURRENCY_METABOLITE_NAMES, display_name
    assert normalised_names("NAD+") == {"nad+", "nad"}
    # a name that is not a charge variant of a listed one still does not match
    assert not normalised_names("dopamine") & DEFAULT_CURRENCY_METABOLITE_NAMES


def test_a_neurotransmitter_is_never_currency_by_either_rule():
    """Defect 2: glycine's cytosolic copy was tagged at 171 edges although every other transmitter in the graph
    propagates - glutamate at 296 edges, glutamine at 258, aspartate at 107."""
    names = {"glycine_c": "glycine", "glutamate_c": "glutamate", "atp_c": "ATP", "dopamine_c": "dopamine"}
    degrees = {"glycine_c": 10_000, "glutamate_c": 10_000, "atp_c": 10_000, "dopamine_c": 10_000}
    tagged = tag_currency_metabolites(names, degrees)
    assert tagged == {"atp_c"}  # the degree rule cannot reach a transmitter, however many edges it has
    for transmitter in ("glycine", "glutamate", "GABA", "dopamine", "serotonin", "acetylcholine", "histamine"):
        assert normalised_names(transmitter) & NEVER_CURRENCY_NAMES, transmitter


def test_the_name_rule_is_reapplied_to_nodes_added_after_the_degree_rule_ran():
    """Defect 3: tag_currency ran on the base graph and later nodes were written is_currency=False whatever their name,
    so thirteen species had copies that disagreed."""
    base = tag_currency_metabolites({"h2o_c": "H2O"}, {"h2o_c": 10_000})
    assert base == {"h2o_c"}
    complete = retag_currency_by_name(["h2o_c", "h2o_v", "glycine_v"], ["H2O", "H2O", "glycine"], base)
    assert complete == {"h2o_c", "h2o_v"}  # the late vacuolar copy is tagged, the late transmitter is not


def test_a_pool_built_to_propagate_is_never_tagged_however_high_its_degree():
    """The interaction fixing defect 3 exposed. The neuronal sodium pool carries 369 edges, above the degree rule's
    quantile, and the membrane potential reads it: tagging it would silence the electrical layer. Before defect 3 was
    fixed it was protected only by accident, because it is added after the tagging ran."""
    names = {"NEURON_MAM02519c": "Na+ (neuronal pool, cytosol)", "MAM02519c": "Na+"}
    degrees = {"NEURON_MAM02519c": 369, "MAM02519c": 601}
    tagged = tag_currency_metabolites(names, degrees)
    assert "NEURON_MAM02519c" not in tagged and "MAM02519c" in tagged
    assert retag_currency_by_name(list(names), list(names.values()), {"NEURON_MAM02519c"}) == {"MAM02519c"}
    assert never_currency("NEURON_MAM01413c", "Ca2+ (neuronal pool, cytosol)")
    assert NEVER_CURRENCY_IDENTIFIER_PREFIXES == ("NEURON_",)


def test_the_fix_on_the_live_node_table_releases_glycine_and_keeps_every_neuronal_pool():
    """A guard on the graph the rebuild will start from, since none of this reaches it until then."""
    nodes_path = pd.io.common.stringify_path("data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes/nodes.parquet")
    try:
        nodes = pd.read_parquet(nodes_path)
    except (FileNotFoundError, OSError):
        pytest.skip("the v4 graph is not on disk")
    metabolites = nodes[nodes.node_type == "metabolite"]
    fixed = retag_currency_by_name(
        list(nodes.node_id), list(nodes.display_name),
        tag_currency_metabolites(dict(zip(metabolites.node_id, metabolites.display_name)),
                                 dict(zip(metabolites.node_id, metabolites.degree))))
    neuronal = [node_id for node_id in nodes.node_id if str(node_id).startswith("NEURON_")]
    assert neuronal and not any(node_id in fixed for node_id in neuronal)
    glycine = nodes[(nodes.display_name.astype(str).str.strip() == "glycine")]
    assert not any(node_id in fixed for node_id in glycine.node_id)
    # every copy of a cofactor species now agrees, which is what defects 1 and 3 broke
    for species in ("NAD+", "NADP+", "HCO3-", "GSH", "H2O", "O2", "Ca2+", "Mg2+", "FMN"):
        copies = nodes[nodes.display_name.astype(str).str.strip() == species]
        if len(copies) > 1:
            assert all(node_id in fixed for node_id in copies.node_id), species
