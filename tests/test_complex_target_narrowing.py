"""A drug's ChEMBL complex targets narrowed to the subunits the Guide to PHARMACOLOGY names, wired into the seeding of
both label routes (decision 8.11; the owner, 13 October 2026: "agreed on both 1 and 2" and "ya absolutely add back
obligate partners"). The rule itself, replace_complex_rows, is tested in tests/test_gtopdb_drug_targets.py; this file
tests what the evidence build does with it."""
import json

import pandas as pd
import pytest

from mechanistic_pathway_learning.perturbation.gtopdb_drug_targets import ComplexTargetNarrowing
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import DrugTarget, GraphNodeLookup
from tests.test_gtopdb_drug_targets import write_gtopdb

GABA_A_SUBUNITS = ["GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA6", "GABRB2", "GABRD", "GABRG2"]


def gaba_a_complex(action_type="POSITIVE ALLOSTERIC MODULATOR"):
    return DrugTarget("CHEMBL2093872", action_type, list(GABA_A_SUBUNITS), "PROTEIN COMPLEX GROUP", "GABA-A receptor; anion channel")


def interaction(ligand, gene, action, species="Human"):
    return {"Ligand ID": ligand, "Target Species": species, "Target Gene Symbol": gene, "Action": action, "Type": "Allosteric modulator", "Primary Target": "t"}


def write_release(tmp_path):
    rows = [interaction(11, gene, "Positive") for gene in ("GABRA1", "GABRA2", "GABRA3")]  # a benzodiazepine: three alpha subunits
    rows += [interaction(11, "GPR68", "Agonist")]  # a GtoPdb target outside every ChEMBL row of the drug
    rows += [interaction(22, "GABRA1", "Positive")]  # one named subunit: an assay, not a selectivity claim
    rows += [interaction(33, gene, "Positive") for gene in ("GABRB2", "GABRD")]  # a drug matched by its name
    mapping = [{"Ligand id": 11, "Name": "diazepam", "PubChem CID": "3016"}, {"Ligand id": 22, "Name": "a one-row drug", "PubChem CID": "5000"},
               {"Ligand id": 33, "Name": "ganaxolone", "PubChem CID": ""}]
    directory = write_gtopdb(tmp_path / "gtopdb", rows, mapping)
    sites = tmp_path / "binding_sites.csv"
    pd.DataFrame([{"drug": "diazepam", "gene": "GABRA1", "site": "benzodiazepine site", "source": "ChEMBL", "evidence": "", "site_id": ""}]).to_csv(sites, index=False)
    return ComplexTargetNarrowing.from_gtopdb(directory, sites)


def test_a_complex_narrows_to_the_named_subunits_and_its_site_s_obligate_subunit(tmp_path):
    narrowing = write_release(tmp_path)
    serotonin = DrugTarget("CHEMBL214", "AGONIST", ["HTR1A"], "SINGLE PROTEIN")
    targets, record = narrowing.narrow([gaba_a_complex(), serotonin], pubchem_cid=3016, drug_name="Diazepam")
    assert record["ligand"] == 11 and record["matched_by"] == "pubchem_cid"
    narrowed = targets[0]
    assert narrowed.target_chembl_id == "CHEMBL2093872" and narrowed.gene_symbols == ["GABRA1", "GABRA2", "GABRA3", "GABRG2"]
    assert narrowed.gene_signs == {"GABRA1": 0.5, "GABRA2": 0.5, "GABRA3": 0.5, "GABRG2": 0.5}, "GtoPdb's Positive is 0.5; the obligate GABRG2 takes ChEMBL's sign"
    assert record["narrowed"]["CHEMBL2093872"]["obligate_partners"] == ["GABRG2"]
    assert targets[1] is serotonin, "a single-gene ChEMBL target is kept whatever GtoPdb says"
    assert all("GPR68" not in target.gene_symbols for target in targets) and len(targets) == 2, "GtoPdb's additions are not added"


def test_one_named_subunit_keeps_the_complex_whole_and_says_why(tmp_path):
    narrowing = write_release(tmp_path)
    targets, record = narrowing.narrow([gaba_a_complex()], pubchem_cid=5000, drug_name="a one-row drug")
    assert targets[0].gene_symbols == GABA_A_SUBUNITS and not targets[0].gene_signs
    assert record["kept_whole"] == {"CHEMBL2093872": "chembl complex, gtopdb names only 1 of its subunits"} and record["narrowed"] == {}


def test_a_drug_is_matched_by_its_name_when_gtopdb_has_no_identifier_for_it(tmp_path):
    narrowing = write_release(tmp_path)
    targets, record = narrowing.narrow([gaba_a_complex()], pubchem_cid=999999, drug_name=" Ganaxolone ")
    assert record["matched_by"] == "name" and targets[0].gene_symbols == ["GABRB2", "GABRD"]


def test_a_drug_gtopdb_does_not_hold_keeps_every_target(tmp_path):
    narrowing = write_release(tmp_path)
    original = [gaba_a_complex()]
    targets, record = narrowing.narrow(original, pubchem_cid=1, drug_name="unknown")
    assert targets == original and record["ligand"] is None and record["kept_whole"] == {"CHEMBL2093872": "no GtoPdb ligand"}


def test_the_narrowing_refuses_a_missing_release_rather_than_seeding_every_subunit(tmp_path):
    with pytest.raises(FileNotFoundError, match="no-narrow-complex-targets"):
        ComplexTargetNarrowing.from_gtopdb(tmp_path / "absent")
    narrowing = write_release(tmp_path)
    assert len(narrowing.provenance["interactions.csv_sha256"]) == 64 and narrowing.provenance["gtopdb_additions"] == "not added"


def test_the_seeds_spread_one_unit_over_the_kept_subunits_with_their_own_signs(tmp_path):
    nodes = pd.DataFrame({"node_id": [f"GENE:{gene}" for gene in GABA_A_SUBUNITS] + ["GENE:HTR1A"], "node_type": "gene",
                          "gene_symbol": GABA_A_SUBUNITS + ["HTR1A"]})
    lookup = GraphNodeLookup.from_nodes(nodes, non_protein_targets_path=None)
    whole = lookup.perturbation_nodes([gaba_a_complex()])
    assert len(whole) == 8 and {sign for _, sign, _ in whole} == {0.5} and whole[0][2] == pytest.approx(1 / 8)
    mixed = DrugTarget("CHEMBL1", "ANTAGONIST", ["GABRA1", "GABRB2", "NOT_A_NODE"], gene_signs={"GABRA1": 0.5})
    seeds = lookup.perturbation_nodes([mixed])
    assert seeds == [["GENE:GABRA1", 0.5, 0.5], ["GENE:GABRB2", -1.0, 0.5]], "a subunit without its own sign takes the target's; magnitude 1/n over the mapped nodes"
    narrowing = write_release(tmp_path)
    narrowed, _ = narrowing.narrow([gaba_a_complex()], pubchem_cid=3016, drug_name="diazepam")
    assert [node for node, _, _ in lookup.perturbation_nodes(narrowed)] == ["GENE:GABRA1", "GENE:GABRA2", "GENE:GABRA3", "GENE:GABRG2"]
    assert all(magnitude == pytest.approx(0.25) for _, _, magnitude in lookup.perturbation_nodes(narrowed))


class RecordingNarrowing:
    """A stand-in with the narrowing's interface, recording which drug each call named."""

    provenance = {"gtopdb_directory": "stand-in"}

    def __init__(self):
        self.calls = []

    def narrow(self, drug_targets, pubchem_cid=None, drug_name=None):
        self.calls.append((pubchem_cid, drug_name))
        narrowed = [DrugTarget(target.target_chembl_id, target.action_type, target.gene_symbols[:1]) for target in drug_targets]
        return narrowed, {"ligand": 1, "matched_by": "pubchem_cid", "narrowed": {target.target_chembl_id: {} for target in drug_targets}, "kept_whole": {}}


def test_the_onsides_route_narrows_once_per_ingredient_by_the_unified_sider_drug():
    from mechanistic_pathway_learning.evidence.load_onsides_label_events import onsides_reports
    from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import IngredientBridge

    bridge = IngredientBridge("83367", "diazepam", "RXCUI", "unii", chembl_parent="CHEMBL12", atc_codes_rxnav=["N05BA01"], perturbation_id="CID100003016",
                              unified_with_sider=True, sider_pubchem_cid=3016, sider_drug_name="diazepam")
    statements = pd.DataFrame([{"ingredient_id": "83367", "ingredient_name": "diazepam", "meddra_id": 1, "meddra_name": "Somnolence", "meddra_term_type": "PT",
                                "target_symptom": "somnolence_or_hypersomnia", "label_section": section, "source_country": "US", "label_count": 3,
                                "single_ingredient_label_count": 3, "labels_of_ingredient_in_country": 3, "max_pred1": float("nan"),
                                "mean_pred1": float("nan"), "match_methods": "SM"} for section in ("AR", "WP")])
    mechanisms = {"CHEMBL12": [{"molecule_chembl_id": "CHEMBL12", "target_chembl_id": "CHEMBL2093872", "action_type": "POSITIVE ALLOSTERIC MODULATOR"}]}
    targets = {"CHEMBL2093872": {"organism": "Homo sapiens", "gene_symbols": list(GABA_A_SUBUNITS), "target_type": "PROTEIN COMPLEX GROUP"}}
    node_lookup = {gene: f"GENE:{gene}" for gene in GABA_A_SUBUNITS}
    plain, plain_counts = onsides_reports(statements, {"83367": bridge}, mechanisms, targets, node_lookup)
    assert len(json.loads(plain.perturbation_nodes.iloc[0])) == 8 and "complex_target_narrowing" not in plain_counts
    recorder = RecordingNarrowing()
    narrowed, counts = onsides_reports(statements, {"83367": bridge}, mechanisms, targets, node_lookup, complex_target_narrowing=recorder)
    assert recorder.calls == [(3016, "diazepam")], "one call for the ingredient's two statements, by the unified SIDER drug"
    assert json.loads(narrowed.perturbation_nodes.iloc[0]) == [["GENE:GABRA1", 0.5, 1.0]]
    assert counts["complex_target_narrowing"]["ingredients_narrowed"] == 1 and "83367" in counts["complex_target_narrowing"]["by_ingredient"]


def test_both_label_routes_offer_the_narrowing_and_take_it_by_default():
    import os
    import subprocess
    import sys
    from pathlib import Path

    repository_root = Path(__file__).resolve().parents[1]
    environment = {**os.environ, "PYTHONPATH": str(repository_root)}  # this checkout's package, not an installed one
    for command in ([sys.executable, "-m", "mechanistic_pathway_learning.evidence.assemble_evidence_table", "--help"],
                    [sys.executable, "experiments/build_onsides_reports.py", "--help"]):
        help_text = subprocess.run(command, capture_output=True, text=True, check=True, cwd=repository_root, env=environment).stdout
        assert "--narrow-complex-targets" in help_text and "--no-narrow-complex-targets" in help_text and "--gtopdb-dir" in help_text, command[-2]
    from mechanistic_pathway_learning.evidence import assemble_evidence_table
    import inspect

    for source in (inspect.getsource(assemble_evidence_table.main), (repository_root / "experiments/build_onsides_reports.py").read_text(encoding="utf-8")):
        assert '"--narrow-complex-targets", action=argparse.BooleanOptionalAction, default=True' in source, "the two routes must make the same choice"
