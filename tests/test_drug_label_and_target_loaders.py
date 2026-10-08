"""Offline tests for the SIDER loader and the ChEMBL drug-target mapping using synthetic files."""
import gzip
import json
from pathlib import Path

from mechanistic_pathway_learning.evidence.load_drug_label_events import (
    is_nervous_system_atc,
    load_sider_events,
    pubchem_cid_from_stitch_flat,
)
import pandas as pd

from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import (
    ChemblCaches,
    DrugTarget,
    check_graph_compounds_have_no_mechanism,
    GraphNodeLookup,
    drug_targets_for_pubchem_cid,
    graph_compound_targets,
    has_dominant_target,
    load_chembl_caches,
    load_drugs_acting_as_graph_compounds,
    mechanism_targets,
    mechanisms_by_parent_and_molecule,
)


def write_synthetic_sider(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    with gzip.open(directory / "meddra_all_se.tsv.gz", "wt", encoding="utf-8") as side_effects:
        side_effects.write("CID100003016\tCID000003016\tC0033975\tLLT\tC0033975\tPsychotic disorder\n")
        side_effects.write("CID100003016\tCID000003016\tC0033975\tPT\tC0033975\tPsychotic disorder\n")
        side_effects.write("CID100003016\tCID000003016\tC0917801\tPT\tC0917801\tInsomnia\n")
        side_effects.write("CID100003016\tCID000003016\tC0015230\tPT\tC0015230\tRash\n")
    with gzip.open(directory / "meddra_freq.tsv.gz", "wt", encoding="utf-8") as frequencies:
        frequencies.write("CID100003016\tCID000003016\tC0917801\t\t1-10%\t0.01\t0.10\tPT\tC0917801\tInsomnia\n")
        frequencies.write("CID100003016\tCID000003016\tC0917801\tplacebo\t2%\t0.02\t0.02\tPT\tC0917801\tInsomnia\n")
    with gzip.open(directory / "meddra_all_indications.tsv.gz", "wt", encoding="utf-8") as indications:
        indications.write("CID100003016\tC0003467\tNLP_indication\tAnxiety\tPT\tC0003467\tAnxiety\n")
        indications.write("CID100003016\tC0033975\ttext_mention\tPsychosis\tPT\tC0033975\tPsychotic disorder\n")
    (directory / "drug_names.tsv").write_text("CID100003016\tdiazepam\n", encoding="utf-8")
    (directory / "drug_atc.tsv").write_text("CID100003016\tN05BA01\n", encoding="utf-8")


def write_synthetic_chembl(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "pubchem_to_chembl_with_parents.json").write_text(json.dumps({"3016": ["CHEMBL12"]}))
    (directory / "mechanisms.json").write_text(json.dumps([
        {"molecule_chembl_id": "CHEMBL12", "target_chembl_id": "CHEMBL2093872", "action_type": "POSITIVE ALLOSTERIC MODULATOR", "mechanism_of_action": "GABA-A receptor; anion channel positive allosteric modulator"}
    ]))
    (directory / "targets.json").write_text(json.dumps({
        "CHEMBL2093872": {"pref_name": "GABA-A receptor; anion channel", "target_type": "PROTEIN COMPLEX GROUP", "organism": "Homo sapiens", "gene_symbols": ["GABRA1", "GABRA2"], "accessions": ["P14867", "P47869"]}
    }))


def test_sider_loader_maps_terms_filters_indication_methods_and_reads_frequency(tmp_path: Path) -> None:
    sider_directory = tmp_path / "sider"
    write_synthetic_sider(sider_directory)
    events = load_sider_events(sider_directory, Path("docs/symptom_crosswalk.csv"))
    by_relation_and_term = {(event.relation, event.meddra_preferred_term): event for event in events}
    assert ("induces", "Psychotic disorder") in by_relation_and_term
    assert by_relation_and_term[("induces", "Psychotic disorder")].target_symptom == "psychosis"
    assert ("induces", "Rash") not in by_relation_and_term  # unmapped term dropped
    insomnia = by_relation_and_term[("induces", "Insomnia")]
    assert abs(insomnia.label_frequency - 0.055) < 1e-9  # placebo row excluded, range midpoint kept
    assert ("relieves", "Anxiety") in by_relation_and_term  # NLP_indication kept
    assert ("relieves", "Psychotic disorder") not in by_relation_and_term  # text_mention dropped
    assert insomnia.drug_name == "diazepam" and insomnia.pubchem_cid == 3016
    assert is_nervous_system_atc(insomnia.atc_codes)
    assert pubchem_cid_from_stitch_flat("CID100003016") == 3016


def test_drug_target_mapping_sign_and_dominance(tmp_path: Path) -> None:
    chembl_directory = tmp_path / "chembl"
    write_synthetic_chembl(chembl_directory)
    caches = load_chembl_caches(chembl_directory)
    drug_targets = drug_targets_for_pubchem_cid(3016, caches)
    assert len(drug_targets) == 1 and has_dominant_target(drug_targets)
    assert drug_targets[0].sign == 0.5
    triples = GraphNodeLookup({"GABRA1": "GENE:GABRA1", "GABRA2": "GENE:GABRA2", "HMBS": "GENE:HMBS"}, {}, {}).perturbation_nodes(drug_targets)
    assert sorted(map(tuple, triples)) == [("GENE:GABRA1", 0.5, 0.5), ("GENE:GABRA2", 0.5, 0.5)]
    assert drug_targets_for_pubchem_cid(999999, caches) == []


def synthetic_graph_nodes() -> pd.DataFrame:
    """Two genes and the compartment copies of two Human-GEM metabolites, in the columns of nodes.parquet."""
    return pd.DataFrame([
        {"node_id": "GENE:TTR", "node_type": "gene", "gene_symbol": "TTR", "ensembl_gene_id": "ENSG00000118271", "base_metabolite_id": None},
        {"node_id": "GENE:GABRA1", "node_type": "gene", "gene_symbol": "GABRA1", "ensembl_gene_id": "ENSG00000022355", "base_metabolite_id": None},
        {"node_id": "MAM01821c", "node_type": "metabolite", "gene_symbol": None, "ensembl_gene_id": None, "base_metabolite_id": "MAM01821"},
        {"node_id": "MAM01821m", "node_type": "metabolite", "gene_symbol": None, "ensembl_gene_id": None, "base_metabolite_id": "MAM01821"},
        {"node_id": "MAM01821s", "node_type": "metabolite", "gene_symbol": None, "ensembl_gene_id": None, "base_metabolite_id": "MAM01821"},
        {"node_id": "MAM02382c", "node_type": "metabolite", "gene_symbol": None, "ensembl_gene_id": None, "base_metabolite_id": "MAM02382"},
        {"node_id": "MAM02382e", "node_type": "metabolite", "gene_symbol": None, "ensembl_gene_id": None, "base_metabolite_id": "MAM02382"},
    ])


def test_mechanisms_recorded_on_a_salt_form_are_found_through_the_parent() -> None:
    mechanisms = [{"molecule_chembl_id": "CHEMBL1200964", "parent_molecule_chembl_id": "CHEMBL629", "target_chembl_id": "CHEMBL222", "action_type": "INHIBITOR"}]
    targets = {"CHEMBL222": {"pref_name": "Norepinephrine transporter", "target_type": "SINGLE PROTEIN", "organism": "Homo sapiens", "gene_symbols": ["SLC6A2"], "accessions": ["P23975"]}}
    by_molecule = mechanisms_by_parent_and_molecule(mechanisms)
    assert [target.target_chembl_id for target in mechanism_targets(["CHEMBL629"], by_molecule, targets)] == ["CHEMBL222"]  # the parent
    assert [target.target_chembl_id for target in mechanism_targets(["CHEMBL1200964"], by_molecule, targets)] == ["CHEMBL222"]  # the salt itself
    caches = ChemblCaches({"2160": ["CHEMBL629"]}, by_molecule, targets)
    assert [target.gene_symbols for target in drug_targets_for_pubchem_cid(2160, caches)] == [["SLC6A2"]]


def test_extra_chembl_ids_supply_targets_when_unichem_gives_none() -> None:
    mechanisms = [{"molecule_chembl_id": "CHEMBL1201", "parent_molecule_chembl_id": "CHEMBL1201", "target_chembl_id": "CHEMBL222", "action_type": "INHIBITOR"}]
    targets = {"CHEMBL222": {"pref_name": "Norepinephrine transporter", "target_type": "SINGLE PROTEIN", "organism": "Homo sapiens", "gene_symbols": ["SLC6A2"], "accessions": []}}
    caches = ChemblCaches({}, mechanisms_by_parent_and_molecule(mechanisms), targets)
    assert drug_targets_for_pubchem_cid(5000, caches) == []
    assert [target.target_chembl_id for target in drug_targets_for_pubchem_cid(5000, caches, extra_chembl_ids=["CHEMBL1201"])] == ["CHEMBL222"]


def test_nucleic_acid_target_maps_to_its_gene_by_ensembl_id() -> None:
    lookup = GraphNodeLookup.from_nodes(synthetic_graph_nodes(), None)
    transthyretin_mrna = DrugTarget("CHEMBL3885585", "RNAI INHIBITOR", [], "NUCLEIC-ACID", "Transthyretin mRNA", ["ENSG00000118271"])
    assert lookup.perturbation_nodes([transthyretin_mrna]) == [["GENE:TTR", -1.0, 1.0]]
    generic_dna = DrugTarget("CHEMBL2311221", "CROSS-LINKING AGENT", [], "NUCLEIC-ACID", "DNA", [])
    assert lookup.perturbation_nodes([generic_dna]) == []


def test_metal_target_seeds_every_compartment_copy_from_the_curated_table(tmp_path: Path) -> None:
    table_path = tmp_path / "non_protein_drug_targets.csv"
    table_path.write_text("target_chembl_id,target_pref_name,target_type,base_metabolite_ids,note\nCHEMBL2363058,Iron,METAL,MAM01821,\nCHEMBL2366381,Aluminium,METAL,,not in Human-GEM\n")
    lookup = GraphNodeLookup.from_nodes(synthetic_graph_nodes(), table_path)
    iron = DrugTarget("CHEMBL2363058", "CHELATING AGENT", [], "METAL", "Iron", [])
    triples = lookup.perturbation_nodes([iron])
    assert sorted(node_id for node_id, _, _ in triples) == ["MAM01821c", "MAM01821m", "MAM01821s"]
    assert all(sign == -1.0 and abs(magnitude - 1 / 3) < 1e-12 for _, sign, magnitude in triples)
    assert lookup.perturbation_nodes([DrugTarget("CHEMBL2366381", "CHELATING AGENT", [], "METAL", "Aluminium", [])]) == []


def test_drug_that_is_a_graph_compound_raises_that_metabolite_only_without_a_mechanism(tmp_path: Path) -> None:
    table_path = tmp_path / "drugs_acting_as_graph_compounds.csv"
    table_path.write_text("drug_name,pubchem_cid,chembl_ids,base_metabolite_ids,note\nlithium,28486,CHEMBL1200826,MAM02382,Li+\n")
    by_pubchem_cid, by_chembl_id = load_drugs_acting_as_graph_compounds(table_path)
    lookup = GraphNodeLookup.from_nodes(synthetic_graph_nodes(), None)
    by_cid = graph_compound_targets([], 28486, [], by_pubchem_cid, by_chembl_id)
    by_chembl = graph_compound_targets([], None, ["CHEMBL1200826"], by_pubchem_cid, by_chembl_id)
    assert len(by_cid) == 1 and has_dominant_target(by_cid) and by_cid[0].sign == 1.0
    assert sorted(map(tuple, lookup.perturbation_nodes(by_chembl))) == [("MAM02382c", 1.0, 0.5), ("MAM02382e", 1.0, 0.5)]
    receptor = DrugTarget("CHEMBL2093872", "AGONIST", ["GABRA1"], "PROTEIN COMPLEX GROUP")
    assert graph_compound_targets([receptor], 28486, [], by_pubchem_cid, by_chembl_id) == [receptor]  # a mechanism target wins
    assert graph_compound_targets([], 1, ["CHEMBL1"], by_pubchem_cid, by_chembl_id) == []


def test_sider_route_falls_back_to_the_bridge_parent_then_to_the_graph_compound(tmp_path: Path) -> None:
    from mechanistic_pathway_learning.evidence.assemble_evidence_table import bridge_chembl_parents_by_sider_cid, sider_drug_targets

    bridge_path = tmp_path / "ingredient_identifier_bridge.json"
    bridge_path.write_text(json.dumps({"ingredients": {
        "100": {"ingredient_id": "100", "ingredient_name": "drug without unichem", "identifier_system": "RXCUI", "bridge_method": "unii", "chembl_parent": "CHEMBL1201",
                "unified_with_sider": True, "sider_pubchem_cid": 5000, "perturbation_id": "CID000005000"},
        "200": {"ingredient_id": "200", "ingredient_name": "not unified", "identifier_system": "RXCUI", "bridge_method": "unii", "chembl_parent": "CHEMBL9999",
                "unified_with_sider": False, "sider_pubchem_cid": None, "perturbation_id": "ONSIDES:200"},
    }}))
    parents_by_cid = bridge_chembl_parents_by_sider_cid(bridge_path)
    assert parents_by_cid == {5000: ["CHEMBL1201"]}
    mechanisms = [{"molecule_chembl_id": "CHEMBL1201", "parent_molecule_chembl_id": "CHEMBL1201", "target_chembl_id": "CHEMBL222", "action_type": "INHIBITOR"}]
    targets = {"CHEMBL222": {"pref_name": "Norepinephrine transporter", "target_type": "SINGLE PROTEIN", "organism": "Homo sapiens", "gene_symbols": ["SLC6A2"], "accessions": []}}
    caches = ChemblCaches({"28486": ["CHEMBL1200826"]}, mechanisms_by_parent_and_molecule(mechanisms), targets)
    assert [target.target_chembl_id for target in sider_drug_targets(5000, caches, parents_by_cid, {}, {})] == ["CHEMBL222"]
    lithium = sider_drug_targets(28486, caches, parents_by_cid, {}, {"CHEMBL1200826": ["MAM02382"]})
    assert [(target.target_chembl_id, target.action_type) for target in lithium] == [("COMPOUND:MAM02382", "EXOGENOUS SUPPLY")]
    assert bridge_chembl_parents_by_sider_cid(None) == {}


def test_a_listed_graph_compound_with_a_chembl_mechanism_is_refused() -> None:
    mechanisms = [{"molecule_chembl_id": "CHEMBL1200826", "parent_molecule_chembl_id": "CHEMBL1200826", "target_chembl_id": "CHEMBL1786", "action_type": "INHIBITOR"}]
    by_molecule = mechanisms_by_parent_and_molecule(mechanisms)
    check_graph_compounds_have_no_mechanism({"CHEMBL96": ["MAM00970"]}, by_molecule)
    try:
        check_graph_compounds_have_no_mechanism({"CHEMBL1200826": ["MAM02382"]}, by_molecule)
    except ValueError as error:
        assert "CHEMBL1200826" in str(error)
    else:
        raise AssertionError("a graph compound with a mechanism was accepted")


def test_the_single_target_rule_can_be_lifted() -> None:
    from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import passes_single_target_rule

    two_targets = [DrugTarget("CHEMBL1786", "INHIBITOR", ["IMPA1"]), DrugTarget("CHEMBL2095188", "INHIBITOR", ["GSK3A", "GSK3B"])]
    assert not has_dominant_target(two_targets) and not passes_single_target_rule(two_targets)
    assert has_dominant_target(two_targets, max_targets=0) and has_dominant_target(two_targets, max_targets=None)
    assert passes_single_target_rule(two_targets, max_targets=None) and has_dominant_target(two_targets, max_targets=2)
    assert not has_dominant_target([], max_targets=None)  # a drug with no mechanism target still has no entry point
    lookup = GraphNodeLookup({"IMPA1": "GENE:IMPA1", "GSK3A": "GENE:GSK3A", "GSK3B": "GENE:GSK3B"}, {}, {})
    assert sorted(map(tuple, lookup.perturbation_nodes(two_targets))) == [("GENE:GSK3A", -1.0, 0.5), ("GENE:GSK3B", -1.0, 0.5), ("GENE:IMPA1", -1.0, 1.0)]
