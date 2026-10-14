"""Decision 9.26's two extra routes. The invariant that matters is the one the EC route exists to protect.

UniProt's EC column names every human protein of an EC class, and one class often names many isoforms - EC
1.14.14.1 names 25. A BRENDA statement about the class does not say which isoform, so joining it onto all of them
would manufacture 25 edges from one statement. The route is taken only where the class names exactly one accession.
"""
import json
from pathlib import Path

import pandas as pd
import pytest

from experiments.map_brenda_regulation_to_graph import human_accessions_by_ec_number

UNIPROT_PATH = Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz")
MAPPED_PAIRS_PATH = Path("configs/brenda_metabolite_enzyme_regulation.csv")
MAPPING_SUMMARY_PATH = Path("docs/brenda_regulation_mapping.json")


@pytest.mark.skipif(not UNIPROT_PATH.exists(), reason="the pinned UniProt dump is not present")
def test_an_ec_class_can_name_many_human_accessions_which_is_why_the_route_is_restricted():
    accessions_of_ec = human_accessions_by_ec_number(UNIPROT_PATH)
    assert len(accessions_of_ec["1.14.14.1"]) > 10, "the cytochrome P450 class names many isoforms"
    classes_naming_one = sum(1 for accessions in accessions_of_ec.values() if len(accessions) == 1)
    assert 0 < classes_naming_one < len(accessions_of_ec), "some classes name one accession and some name several"


@pytest.mark.skipif(not MAPPED_PAIRS_PATH.exists(), reason="the mapped table has not been built")
def test_the_mapped_table_carries_only_the_three_regulation_relations_with_matching_signs():
    pairs = pd.read_csv(MAPPED_PAIRS_PATH)
    assert set(pairs.relation_type) <= {"inhibits", "activates", "cofactor_of"}
    assert set(pairs[pairs.relation_type == "inhibits"].sign) == {-1.0}
    assert set(pairs[pairs.relation_type == "activates"].sign) == {1.0}
    assert set(pairs[pairs.relation_type == "cofactor_of"].sign) == {0.0}
    assert not pairs.duplicated(subset=["base_metabolite_id", "uniprot_accession", "relation_type"]).any()


@pytest.mark.skipif(not MAPPING_SUMMARY_PATH.exists(), reason="the mapping summary has not been built")
def test_the_summary_records_the_records_the_ec_route_refused():
    summary = json.loads(MAPPING_SUMMARY_PATH.read_text())
    refused = summary["enzyme_resolution"]["EC class names none or several human accessions"]
    assert refused > 0, "the refusal count is what shows the cross-product was not created"
    assert summary["pairs_written"] > 0
