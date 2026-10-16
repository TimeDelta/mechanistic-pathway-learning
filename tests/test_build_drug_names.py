"""drug_names.csv is derived from the evidence table, the OnSIDES bridge and PubChem (experiments/label_review/build_drug_names.py)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments" / "label_review"))
import build_drug_names  # noqa: E402


def test_rxnorm_names_come_from_the_unified_ingredients_only():
    document = {"ingredients": {
        "1": {"ingredient_name": "fluoxetine", "unified_with_sider": True, "perturbation_id": "CID100003386"},
        "2": {"ingredient_name": "fluoxetine hydrochloride", "unified_with_sider": True, "perturbation_id": "CID100003386"},
        "3": {"ingredient_name": "atenolol", "unified_with_sider": False, "perturbation_id": "CID100002249"},
        "4": {"ingredient_name": "unbridged", "unified_with_sider": True, "perturbation_id": ""},
    }}
    assert build_drug_names.rxnorm_names_by_perturbation(document) == {"CID100003386": "fluoxetine;fluoxetine hydrochloride"}
    assert build_drug_names.rxnorm_names_by_perturbation({}) == {}


def test_titles_are_read_from_the_cache_and_fetched_only_when_allowed(tmp_path, monkeypatch):
    cache = tmp_path / "titles.json"
    cache.write_text(json.dumps({"3386": "Fluoxetine"}))
    assert build_drug_names.pubchem_titles([3386, 2249], cache, fetch=False) == {3386: "Fluoxetine"}
    fetched = []

    def fake_fetch(cids):
        fetched.append(list(cids))
        return {2249: "Atenolol"}

    monkeypatch.setattr(build_drug_names, "fetch_titles", fake_fetch)
    assert build_drug_names.pubchem_titles([3386, 2249], cache, fetch=True) == {3386: "Fluoxetine", 2249: "Atenolol"}
    assert fetched == [[2249]], "only the uncached CID is requested"
    assert json.loads(cache.read_text()) == {"2249": "Atenolol", "3386": "Fluoxetine"}, "the cache is extended in place"


def test_the_table_carries_names_titles_and_selection_membership(tmp_path):
    drugs = pd.DataFrame({"perturbation_id": ["CID100002249", "CID100003386"], "perturbation_label": ["atenolol", "fluoxetine"]})
    selection = tmp_path / "selection.parquet"
    pd.DataFrame({"perturbation_id": ["CID100003386", "CID100002249"], "symptom": ["anxiety", "fatigue"], "keep": [True, False]}).to_parquet(selection)
    selections = {"in_current": build_drug_names.kept_drugs(selection), "in_approved": build_drug_names.kept_drugs(tmp_path / "absent.parquet")}
    table = build_drug_names.build_table(drugs, {"CID100003386": "fluoxetine"}, {3386: "Fluoxetine", 2249: "Atenolol"}, selections)
    assert list(table.columns) == ["perturbation_id", "perturbation_label", "rxnorm_names", "pubchem_title", "in_current", "in_approved"]
    assert table.set_index("perturbation_id").loc["CID100003386"].tolist() == ["fluoxetine", "fluoxetine", "Fluoxetine", True, False]
    assert table.set_index("perturbation_id").loc["CID100002249"].tolist() == ["atenolol", "", "Atenolol", False, False], "a drug with no kept row, and an absent selection, are False"
