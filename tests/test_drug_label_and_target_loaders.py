"""Offline tests for the SIDER loader and the ChEMBL drug-target mapping using synthetic files."""
import gzip
import json
from pathlib import Path

from mechanistic_pathway_learning.evidence.load_drug_label_events import (
    is_nervous_system_atc,
    load_sider_events,
    pubchem_cid_from_stitch_flat,
)
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import (
    drug_targets_for_pubchem_cid,
    has_dominant_target,
    load_chembl_caches,
    perturbation_from_drug_targets,
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
    triples = perturbation_from_drug_targets(drug_targets, {"GABRA1": 0, "GABRA2": 1, "HMBS": 2})
    assert sorted(triples) == [(0, 0.5, 0.5), (1, 0.5, 0.5)]
    assert drug_targets_for_pubchem_cid(999999, caches) == []
