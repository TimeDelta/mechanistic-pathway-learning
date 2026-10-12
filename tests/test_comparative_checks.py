"""Tests: comparative checks (evaluation/comparative_checks.py; experiments/build_comparative_checks.py and
run_comparative_checks.py). A check compares a drug with a second, varied factor against a reference input. The tests
cover what each kind compares, when a graph can express a check, the union rule that places a combination on a side of
the held-out set, the reading of a result, and one pass through an encoder and head built the way the trainer builds
them (untrained: the pass checks the plumbing, not the pharmacology)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch

from experiments.build_comparative_checks import carrier_and_displacement_rows, enzyme_rows
from experiments.run_comparative_checks import predict_seed_lists, structure_lines, structure_table
from experiments.run_main_model import build_argument_parser, build_models
from mechanistic_pathway_learning.evaluation.comparative_checks import (
    CHECK_COLUMNS,
    ComparativeCheck,
    check_inputs,
    evaluate_checks,
    load_checks,
    side_of,
    symptoms_read,
    why_not_expressible,
)
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder
from tests.test_clearance_edges import variant_with_enzymes
from tests.test_drug_entry_nodes import DRUG_SEEDS, evidence_row, write_graph

CARRIER = ComparativeCheck("carrier", "carrier_level", "drug_one", ("ALB",))
ENZYME = ComparativeCheck("enzyme", "enzyme_activity", "drug_two", ("C",), symptom="anxiety")
DISPLACEMENT = ComparativeCheck("displacement", "displacement", "drug_one", ("ALB",), "drug_two")
PAIR = ComparativeCheck("pair", "second_drug", "drug_one", (), "drug_two", "anxiety")
NEGATIVE_PAIR = ComparativeCheck("negative_pair", "second_drug", "drug_two", (), "drug_one", "anxiety", "none")
ALL_CHECKS = [CARRIER, ENZYME, DISPLACEMENT, PAIR, NEGATIVE_PAIR]


def toy_data(tmp_path: Path, drug_entry: str = "nodes", with_enzymes: bool = True):
    _, _, new_nodes, new_edges, relations, entry_table, _ = variant_with_enzymes() if with_enzymes else variant_with_enzymes(None)
    write_graph(tmp_path / "graph", new_nodes, new_edges, relations, entry_table)
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir(exist_ok=True)
    rows = [evidence_row("ALB", "gene", "cluster:ALB", [["GENE:ALB", -1.0, 1.0]])]
    rows += [evidence_row(perturbation, "drug", perturbation, [list(triple) for triple in triples]) for perturbation, (_, triples) in DRUG_SEEDS.items()]
    rows.append({**evidence_row("drug_one", "drug", "drug_one", [list(triple) for triple in DRUG_SEEDS["drug_one"][1]]), "symptom": "fatigue"})
    pd.DataFrame(rows).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    return load_experiment_data(tmp_path / "graph", evidence_directory, drug_entry=drug_entry)


def test_each_kind_compares_the_inputs_it_says(tmp_path: Path) -> None:
    data = toy_data(tmp_path)
    node = data.node_index
    drug_one, drug_two = (node["DRUG:drug_one"], 1.0, 1.0), (node["DRUG:drug_two"], 1.0, 1.0)
    assert check_inputs(CARRIER, data) == ([drug_one], [drug_one, (node["GENE:ALB"], -1.0, 1.0)])
    assert check_inputs(ENZYME, data) == ([drug_two], [drug_two, (node["GENE:C"], -1.0, 1.0)])
    # displacement: both drugs are in the reference, so the difference is the carrier's edge alone
    assert check_inputs(DISPLACEMENT, data) == ([drug_one, drug_two], [drug_one, drug_two, (node["GENE:ALB"], -1.0, 1.0)])
    assert check_inputs(PAIR, data) == ([drug_one], [drug_one, drug_two])
    # a gene with two protein nodes shares a total magnitude of 1 between them
    shared = check_inputs(CARRIER, data, {"GENE:ALB": ["GENE:ALB", "GENE:A"]})[1][1:]
    assert shared == [(node["GENE:ALB"], -1.0, 0.5), (node["GENE:A"], -1.0, 0.5)]
    assert [data.symptoms[column] for column in symptoms_read(CARRIER, data)] == ["anxiety", "fatigue"]  # the drug's own kept positives
    assert [data.symptoms[column] for column in symptoms_read(ENZYME, data)] == ["anxiety"]  # the one the statement names
    with pytest.raises(ValueError, match="names the gene that is lowered"):
        ComparativeCheck("broken", "carrier_level", "drug_one")
    with pytest.raises(ValueError, match="names its second drug"):
        ComparativeCheck("broken", "second_drug", "drug_one")
    with pytest.raises(ValueError, match="not one of"):
        ComparativeCheck("broken", "dose", "drug_one")


def test_a_check_is_expressible_only_where_the_drug_has_a_node_and_the_edge_exists(tmp_path: Path) -> None:
    with_enzymes = toy_data(tmp_path / "enzymes")
    assert [why_not_expressible(check, with_enzymes) for check in ALL_CHECKS] == [""] * 5
    without_enzymes = toy_data(tmp_path / "plain", with_enzymes=False)
    assert why_not_expressible(ENZYME, without_enzymes) == "no clears edge from C to the drug's node"
    assert why_not_expressible(CARRIER, without_enzymes) == ""
    # drug_two has no known carrier
    assert why_not_expressible(ComparativeCheck("no_carrier", "carrier_level", "drug_two", ("ALB",)), with_enzymes) == "no sequesters edge from ALB to the drug's node"
    assert why_not_expressible(ComparativeCheck("no_node", "carrier_level", "drug_one", ("NOT_A_GENE",)), with_enzymes) == "no node for NOT_A_GENE"
    assert why_not_expressible(ComparativeCheck("absent", "carrier_level", "another_drug", ("ALB",)), with_enzymes) == "another_drug is not a perturbation of this table"
    assert why_not_expressible(ComparativeCheck("symptom", "second_drug", "drug_one", (), "drug_two", "catatonia"), with_enzymes) == "catatonia is not a symptom of this table"
    # seeded on its targets a drug has no node, so a lowered carrier does not reach it; a pair of drugs needs no node
    as_targets = toy_data(tmp_path / "targets", drug_entry="targets")
    assert why_not_expressible(CARRIER, as_targets).startswith("the drug has no node of its own here")
    assert why_not_expressible(PAIR, as_targets) == ""


def test_a_combination_is_held_out_when_one_of_its_parts_is() -> None:
    assert side_of(PAIR, {"drug_two"}) == "held out" and side_of(PAIR, {"drug_three"}) == "development"
    assert side_of(DISPLACEMENT, {"drug_two"}) == "held out"
    assert side_of(CARRIER, {"ALB"}) == "held out"  # the carrier is itself a gene perturbation of the study
    assert side_of(CARRIER, set()) == "development"


def test_the_reading_of_a_check_and_what_is_left_unjudged(tmp_path: Path) -> None:
    data = toy_data(tmp_path)
    anxiety, fatigue = data.symptoms.index("anxiety"), data.symptoms.index("fatigue")

    def predict(seed_lists, positions):
        """A stand-in model: every lowered node adds 0.1 to anxiety and takes 0.05 from fatigue; a second drug adds 0.2 to anxiety."""
        probabilities = np.full((len(seed_lists), len(data.symptoms)), 0.3)
        for row, seeds in enumerate(seed_lists):
            lowered = sum(1 for _, sign, _ in seeds if sign < 0)
            drugs = sum(1 for node, _, _ in seeds if bool(data.entry_node_mask[node]))
            probabilities[row, anxiety] += 0.1 * lowered + 0.2 * (drugs - 1)
            probabilities[row, fatigue] -= 0.05 * lowered
        assert positions == [data.perturbation_ids.index(name) for name in ("drug_one", "drug_one", "drug_two", "drug_two", "drug_one", "drug_one", "drug_one", "drug_one", "drug_two", "drug_two")]
        return probabilities

    absent = ComparativeCheck("absent", "carrier_level", "another_drug", ("ALB",))
    table = evaluate_checks([*ALL_CHECKS, absent], data, predict).set_index("statement_id")
    assert table.loc["carrier", ["symptoms_read", "holds"]].tolist() == [2, False]  # anxiety rises and fatigue falls: one symptom lower
    assert table.loc["carrier", "share_not_lower"] == 0.5 and table.loc["carrier", "mean_difference"] == pytest.approx(0.025) and table.loc["carrier", "smallest_difference"] == pytest.approx(-0.05)
    assert bool(table.loc["enzyme", "holds"]) and table.loc["enzyme", "mean_difference"] == pytest.approx(0.1)
    assert bool(table.loc["pair", "holds"]) and table.loc["pair", "mean_difference"] == pytest.approx(0.2)
    assert table.loc["negative_pair", "holds"] is None and table.loc["negative_pair", "mean_difference"] == pytest.approx(0.2)  # reported, not judged
    assert not table.loc["absent", "expressible"] and table.loc["absent", "holds"] is None and np.isnan(table.loc["absent", "mean_difference"])
    assert evaluate_checks([absent], data, lambda seed_lists, positions: pytest.fail("no input to predict")).expressible.tolist() == [False]


def test_the_checks_table_keeps_its_provenance_and_the_report_counts_by_kind_and_side(tmp_path: Path) -> None:
    data = toy_data(tmp_path)
    perturbation_id_of = {"drug one": "drug_one", "drug two": "drug_two"}
    statements = pd.DataFrame([
        {"statement_id": "carrier_statement", "kind": "carrier_level", "drug": "drug one", "carrier": "albumin", "displacing_drug": "", "symptom": "anxiety", "source": "a label",
         "source_id": "set_id 1", "source_version": "2026", "quote": "lower albumin", "second_mechanism_or_limit": "a trend"},
        {"statement_id": "displacement_statement", "kind": "displacement", "drug": "drug one", "carrier": "albumin", "displacing_drug": "drug two", "symptom": "", "source": "a label",
         "source_id": "set_id 1", "source_version": "2026", "quote": "displaces", "second_mechanism_or_limit": "also inhibits metabolism"},
        {"statement_id": "outside", "kind": "displacement", "drug": "drug one", "carrier": "albumin", "displacing_drug": "warfarin", "symptom": "", "source": "", "source_id": "",
         "source_version": "", "quote": "", "second_mechanism_or_limit": ""}])
    rows, left_out = carrier_and_displacement_rows(statements, perturbation_id_of)
    assert left_out == ["outside: warfarin not a study drug"]
    enzymes = pd.DataFrame([{"drug": "drug two", "enzyme_gene": "C", "gives_an_edge": "True", "fda_table_section": "1", "fda_table_drug": "Drug Two", "fda_table_gene": "C",
                             "fda_table_subgroup": "poor metabolizers", "fda_table_description": "Results in higher systemic concentrations.", "source": "https://www.fda.gov/table", "page_saved_on": "2026-10-10"},
                            {"drug": "drug one", "enzyme_gene": "", "gives_an_edge": "False", "fda_table_section": "3", "fda_table_drug": "Drug One", "fda_table_gene": "C",
                             "fda_table_subgroup": "poor metabolizers", "fda_table_description": "May alter systemic concentrations.", "source": "https://www.fda.gov/table", "page_saved_on": "2026-10-10"}])
    rows += enzyme_rows(enzymes, perturbation_id_of)
    checks_path = tmp_path / "checks.csv"
    pd.DataFrame(rows, columns=list(CHECK_COLUMNS)).to_csv(checks_path, index=False)
    checks = load_checks(checks_path)
    assert [(check.statement_id, check.kind, check.lowered_genes, check.second_drug_perturbation_id, check.symptom) for check in checks] == [
        ("carrier_statement", "carrier_level", ("ALB",), "", "anxiety"), ("displacement_statement", "displacement", ("ALB",), "drug_two", ""), ("enzyme_drug_two_C", "enzyme_activity", ("C",), "", "")]
    written = pd.read_csv(checks_path).fillna("")
    assert written.quote.tolist() == ["lower albumin", "displaces", "Drug Two, C, poor metabolizers: Results in higher systemic concentrations."]

    structure = structure_table(checks, data, {}, {"drug_two"})
    assert structure.expressible.tolist() == [True, True, True] and structure.side.tolist() == ["development", "held out", "held out"]
    assert structure_table(checks, data, {}, None).side.tolist() == ["not assigned"] * 3
    report = "\n".join(structure_lines(structure, Path("graph"), Path("evidence"), checks_path, "nodes", Path("held_out.csv")))
    assert "| carrier_level | carrier lowered with the drug | 1 | 1 | 1 | 0 |" in report and "| enzyme_activity | clearing enzyme lowered with the drug | 1 | 1 | 0 | 1 |" in report


def test_the_checks_pass_through_an_encoder_and_head_built_as_the_trainer_builds_them(tmp_path: Path) -> None:
    data = toy_data(tmp_path)
    arguments = build_argument_parser().parse_args(["--graph-dir", str(tmp_path / "graph"), "--evidence-dir", str(tmp_path / "evidence"), "--encoder", "message_passing",
                                                    "--node-state-dim", "6", "--num-layers", "2", "--num-modules", "2", "--batch-size", "3"])
    torch.manual_seed(0)
    device = torch.device("cpu")
    encoder, head = build_models(data, arguments, device)
    assert encoder.change_only_relation_indices == (data.relation_types.index("sequesters"), data.relation_types.index("clears"))
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(torch.as_tensor(np.stack([data.edge_source, data.edge_target]), dtype=torch.long),
                                                                             torch.as_tensor(data.edge_relation, dtype=torch.long), len(data.node_ids), len(data.relation_types))

    def predict(seed_lists, positions):
        return predict_seed_lists(encoder, head, data, seed_lists, positions, adjacencies, arguments, device)

    table = evaluate_checks(ALL_CHECKS, data, predict)
    assert table.expressible.all() and table.mean_difference.notna().all() and set(table.holds[table.expectation == "at_least_as_high"]) <= {True, False}
    # the reference of a check is the trainer's own prediction for that drug: same seeds, same call
    from experiments.run_main_model import predict as trainer_predict

    position = data.perturbation_ids.index("drug_one")
    reference, _ = check_inputs(CARRIER, data)
    np.testing.assert_allclose(predict([reference], [position])[0], trainer_predict(encoder, head, data, np.array([position]), adjacencies, arguments, device)[0], atol=1e-6)
    assert json.loads(table.to_json(orient="records"))[0]["statement_id"] == "carrier"
