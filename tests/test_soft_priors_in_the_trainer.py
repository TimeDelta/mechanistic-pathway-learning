"""The soft priors in the loader and the trainer (experiments/run_main_model.py): off by default and then invisible,
on only for the rows a fit trains on, recorded in the fingerprint and results.json, and the permutation control."""
from __future__ import annotations

import json
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from experiments.run_main_model import (
    SOFT_PRIOR_SETTINGS,
    VALUES_BEFORE_THE_V4_DEFAULTS,
    build_argument_parser,
    build_soft_priors,
    configuration_fingerprint,
    fit_targets_with_soft_priors,
)
from mechanistic_pathway_learning.evaluation.experiment_data import PER_PERTURBATION_ARRAYS, load_experiment_data, restrict_to_perturbations

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SYMPTOMS = ("anxiety", "psychosis", "insomnia")


def write_toy_dataset(root: Path, num_genes: int = 60) -> dict[str, Path]:
    """A chain of reactions catalysed by genes, grade-A labels, grade-C rows on some unlabelled pairs, and a report
    table with IMPC calls and absences on anxiety and PubTator relations on psychosis."""
    generator = np.random.default_rng(0)
    genes = [f"G{index}" for index in range(num_genes)]
    nodes = [{"node_id": f"GENE:{gene}", "node_type": "gene", "gene_symbol": gene, "is_currency": False, "degree": 1.0} for gene in genes]
    nodes += [{"node_id": f"R{index}", "node_type": "reaction", "gene_symbol": None, "is_currency": False, "degree": 4.0} for index in range(10)]
    edges = [{"source_id": f"GENE:{gene}", "target_id": f"R{index % 10}", "relation_type": "catalyzed_by", "sign": 1.0} for index, gene in enumerate(genes)]
    edges += [{"source_id": f"R{index}", "target_id": f"R{(index + 1) % 10}", "relation_type": "product_of", "sign": 1.0} for index in range(10)]
    (root / "graph").mkdir(parents=True)
    pd.DataFrame(nodes).to_parquet(root / "graph/nodes.parquet", index=False)
    pd.DataFrame(edges).to_parquet(root / "graph/edges.parquet", index=False)
    (root / "graph/relation_types.json").write_text(json.dumps(["catalyzed_by", "product_of"]))

    def record(gene, index, symptom, grade="A", weight=1.0, references=1):
        return {"perturbation_id": gene, "perturbation_type": "gene", "perturbation_label": gene, "group_id": f"group{index // 3}", "disease_cluster_id": f"cluster{index // 3}",
                "symptom": symptom, "relation": "induces", "evidence_class": "monogenic", "grade": grade, "weight": weight, "label_frequency": None, "in_metabolic_layer": True,
                "evidence_date": None, "perturbation_nodes": json.dumps([[f"GENE:{gene}", -1.0, 1.0]]), "positive_report_count": 1, "negative_report_count": 0,
                "report_count": 1, "distinct_pubmed_reference_count": references, "source": "HPO"}

    records = []
    for index, gene in enumerate(genes):
        labelled = [symptom for symptom in SYMPTOMS if generator.random() < 0.4] or [SYMPTOMS[index % 3]]
        records += [record(gene, index, symptom) for symptom in labelled]
        records += [record(gene, index, symptom, grade="C", weight=0.0, references=2) for symptom in SYMPTOMS if symptom not in labelled and generator.random() < 0.3]
    (root / "evidence").mkdir()
    pd.DataFrame(records).to_parquet(root / "evidence/evidence_records.parquet", index=False)
    reports = []
    for index, gene in enumerate(genes):
        reports.append({"report_id": f"IMPC|{gene}|anxiety", "perturbation_id": gene, "symptom": "anxiety", "relation": "induces", "evidence_class": "animal", "source": "IMPC",
                        "report_value": int(index % 2), "source_record_id": f"MGI:{index}", "references": "", "pmid": None})
        for offset in range(1 + index % 3):
            pmid = str(1000 + 10 * index + offset)
            reports.append({"report_id": f"PubTator3-cause|{gene}|psychosis|{pmid}", "perturbation_id": gene, "symptom": "psychosis", "relation": "induces",
                            "evidence_class": "literature", "source": "PubTator3-cause", "report_value": 1, "source_record_id": f"PMID:{pmid}", "references": f"PMID:{pmid}", "pmid": pmid})
    pd.DataFrame(reports).to_parquet(root / "reports.parquet", index=False)
    table = pd.DataFrame(records)
    is_grade_c = table.grade == "C"
    pd.DataFrame({"perturbation_id": table.perturbation_id, "symptom": table.symptom, "perturbation_type": "gene", "label_frequency": None, "keep": ~is_grade_c,
                  "reason": np.where(is_grade_c, "grade C only", "kept"), "masks_a_negative": is_grade_c}).to_parquet(root / "selection.parquet", index=False)
    return {"graph": root / "graph", "evidence": root / "evidence", "reports": root / "reports.parquet", "selection": root / "selection.parquet"}


def test_the_prior_is_off_by_default_and_registered_as_absent_before_it_existed():
    arguments = build_argument_parser().parse_args(["--graph-dir", "g", "--evidence-dir", "e"])
    assert arguments.soft_prior_strength == 0.0 and arguments.soft_prior_reports == [] and arguments.permute_soft_priors is False
    assert arguments.soft_prior_tiers == Path("configs/soft_prior_tiers.csv") and arguments.soft_prior_tempering == 0.5 and arguments.soft_prior_count_cap == 2
    assert VALUES_BEFORE_THE_V4_DEFAULTS["soft_prior_strength"] == 0.0 and VALUES_BEFORE_THE_V4_DEFAULTS["soft_prior_reports"] == []
    assert {"soft_prior_tiers", "soft_prior_reports", "permute_soft_priors"} <= SOFT_PRIOR_SETTINGS


def test_the_fingerprint_names_the_prior_only_while_it_is_on(tmp_path):
    data = Namespace(outcomes=np.zeros((2, 2)), weights=np.zeros((2, 2)), edge_source=np.array([0]), edge_target=np.array([1]), edge_relation=np.array([0]),
                     perturbation_seeds=[np.array([0]), np.array([1])], perturbation_signs=[np.array([-1.0]), np.array([-1.0])], perturbation_magnitudes=[np.array([1.0]), np.array([1.0])])
    tiers = tmp_path / "tiers.csv"
    tiers.write_text("evidence_class,direction,prior_probability,unit_strength,note\nanimal,positive,0.6,0.3,x\n")
    reports = tmp_path / "reports.parquet"
    pd.DataFrame({"a": [1]}).to_parquet(reports)
    base = dict(graph_dir="g", soft_prior_tiers=tiers, soft_prior_tempering=0.5, soft_prior_count_cap=2, soft_prior_record_grades=["C"], soft_prior_reports=[reports],
                permute_soft_priors=True)
    off = configuration_fingerprint(Namespace(**base, soft_prior_strength=0.0), data, None)
    assert not [key for key in off if "soft_prior" in key], "at strength 0 the settings mean nothing and an old run's fingerprint stays as it was"
    on = configuration_fingerprint(Namespace(**base, soft_prior_strength=1.0), data, None)
    assert on["argument:soft_prior_strength"] == "1.0" and on["argument:permute_soft_priors"] == "True" and on["argument:soft_prior_tempering"] == "0.5"
    assert len(on["file:soft_prior_tiers"]) == 64 and len(on["file:soft_prior_reports:reports.parquet"]) == 64


def test_fit_targets_leave_the_tensors_alone_without_priors_and_touch_only_fit_rows_with_them():
    outcomes = torch.tensor([[1.0, 0.0], [0.0, 0.0], [0.0, 1.0]])
    weights = torch.tensor([[1.0, 0.2], [0.2, 0.2], [0.2, 1.0]])
    plain, applied = fit_targets_with_soft_priors(outcomes, weights, None, np.array([0, 1]), None, None, None, 1.0)
    assert plain["outcomes"] is outcomes and plain["weights"] is weights and applied is None
    probability = np.array([[np.nan, 0.7], [0.0, 0.7], [0.7, np.nan]])
    strength = np.array([[0.0, 0.5], [0.3, 0.5], [0.5, 0.0]])
    eligible = (outcomes.numpy() == 0)
    fitted, applied = fit_targets_with_soft_priors(outcomes, weights, None, np.array([0, 1]), probability, strength, eligible, 0.5)
    assert np.allclose(fitted["outcomes"].numpy(), [[1.0, 0.7], [0.0, 0.7], [0.0, 1.0]]) and np.allclose(fitted["weights"].numpy(), [[1.0, 0.25], [0.2, 0.25], [0.2, 1.0]])
    assert outcomes.tolist() == [[1.0, 0.0], [0.0, 0.0], [0.0, 1.0]], "the plain tensors, which the validation loss reads, are unchanged"
    assert applied["positive_cells"] == 2 and applied["absence_cells"] == 1, "row 2 is not fitted"
    frequency_targets = torch.tensor([[0.3, 0.0], [0.0, 0.0], [0.0, 0.8]])
    fitted, _ = fit_targets_with_soft_priors(outcomes, weights, frequency_targets, np.array([0, 1, 2]), probability, strength, eligible, 1.0)
    assert np.allclose(fitted["positive_targets"].numpy(), [[0.3, 0.7], [0.0, 0.7], [0.7, 0.8]]) and fitted["outcomes"] is outcomes, "under frequency targets the prior writes into the frequency tensor"


def test_the_loader_marks_the_cells_set_aside_for_lacking_a_grade_and_slices_them_with_the_rows(tmp_path):
    paths = write_toy_dataset(tmp_path, num_genes=12)
    data = load_experiment_data(paths["graph"], paths["evidence"], label_selection=paths["selection"])
    assert data.masked_for_soft_evidence is not None and data.masked_for_soft_evidence.shape == data.outcomes.shape
    assert np.array_equal(data.masked_for_soft_evidence, ~data.label_mask), "every masked cell of this selection is a grade-C set-aside"
    assert (data.outcomes[data.masked_for_soft_evidence] == 0).all()
    assert "masked_for_soft_evidence" in PER_PERTURBATION_ARRAYS
    kept = np.zeros(len(data.perturbation_ids), dtype=bool)
    kept[:5] = True
    assert np.array_equal(restrict_to_perturbations(data, kept).masked_for_soft_evidence, data.masked_for_soft_evidence[:5])
    assert load_experiment_data(paths["graph"], paths["evidence"]).masked_for_soft_evidence is None


def test_build_soft_priors_reads_the_records_and_the_reports_into_aligned_matrices(tmp_path):
    paths = write_toy_dataset(tmp_path, num_genes=12)
    data = load_experiment_data(paths["graph"], paths["evidence"], label_selection=paths["selection"])
    arguments = Namespace(soft_prior_strength=1.0, soft_prior_tiers=REPOSITORY_ROOT / "configs/soft_prior_tiers.csv", evidence_dir=paths["evidence"],
                          soft_prior_record_grades=["C"], soft_prior_reports=[paths["reports"]], soft_prior_tempering=0.5, soft_prior_count_cap=2, soft_prior_probability="tier")
    probability, strength, summary, literature, cells = build_soft_priors(data, arguments)
    assert cells is None and summary["probability_source"] == "tier"
    assert probability.shape == data.outcomes.shape and literature is not None and (literature.evidence_class == "literature").any()
    psychosis = data.symptoms.index("psychosis")
    unlabelled_psychosis = data.outcomes[:, psychosis] == 0
    assert np.isin(probability[unlabelled_psychosis, psychosis], [0.55, 0.75]).all() and (strength[unlabelled_psychosis, psychosis] > 0).all(), \
        "every gene has a literature relation on psychosis; a grade-C record outranks it where there is one"
    anxiety = data.symptoms.index("anxiety")
    absences = [index for index, gene in enumerate(data.perturbation_ids) if int(gene[1:]) % 2 == 0 and data.outcomes[index, anxiety] == 0 and not data.masked_for_soft_evidence[index, anxiety]]
    assert absences and (probability[absences, anxiety] == 0.0).all() and np.allclose(strength[absences, anxiety], 0.3)
    literature_pairs = summary["by_class_and_direction"]["literature:positive"]["pairs"]
    grade_c_on_psychosis = int((pd.read_parquet(paths["evidence"] / "evidence_records.parquet").query("grade == 'C' and symptom == 'psychosis'")).shape[0])
    assert literature_pairs + grade_c_on_psychosis == 12 and summary["classes_without_a_tier"] == {}, "one prior per gene on psychosis: the record where there is one, else the literature"
    assert summary["strength"] == 1.0 and summary["reports"] == [str(paths["reports"])]
    assert build_soft_priors(data, Namespace(soft_prior_strength=0.0)) == (None, None, None, None, None)
    calibrated = Namespace(**{**vars(arguments), "soft_prior_probability": "calibrated"})
    probability_again, strength_again, _, _, cells = build_soft_priors(data, calibrated)
    assert np.array_equal(np.nan_to_num(probability_again), np.nan_to_num(probability)) and np.array_equal(strength_again, strength), "the tiers' matrices; a fit calibrates them on its rows"
    assert cells is not None and len(cells) == int((strength > 0).sum()) and {"row", "column", "log_records", "log_popularity", "perturbation:drug"} <= set(cells.columns)


def run_trainer(paths: dict[str, Path], run_directory: Path, *extra: str) -> dict:
    command = [sys.executable, "experiments/run_main_model.py", "--graph-dir", str(paths["graph"]), "--evidence-dir", str(paths["evidence"]), "--run-dir", str(run_directory),
               "--encoder", "message_passing", "--max-epochs", "2", "--num-bootstrap", "3", "--refit-on-validation", "--validation-fraction", "0.3",
               "--no-keep-large-groups-in-training", "--fold", "0", "--seed", "0", "--batch-size", "8", "--placebo-negative-weights-file", "none",
               "--group-by", "disease_cluster", "--label-selection", str(paths["selection"]), *extra]
    completed = subprocess.run(command, cwd=REPOSITORY_ROOT, capture_output=True, text=True, env={"PYTHONPATH": ".", "PATH": "/usr/bin:/bin:/usr/local/bin", "OMP_NUM_THREADS": "1"})
    assert completed.returncode == 0, completed.stderr[-3000:]
    split_directory = next(path for path in run_directory.iterdir() if path.is_dir())
    return {"results": json.loads((split_directory / "results.json").read_text()), "predictions": np.load(split_directory / "test_predictions.npy")}


def test_end_to_end_the_prior_off_is_the_plain_run_and_on_it_trains_scores_and_records(tmp_path):
    paths = write_toy_dataset(tmp_path)
    plain = run_trainer(paths, tmp_path / "plain")
    off = run_trainer(paths, tmp_path / "off", "--soft-prior-strength", "0", "--soft-prior-reports", str(paths["reports"]), "--permute-soft-priors")
    assert np.array_equal(plain["predictions"], off["predictions"]), "strength 0 with every other flag set is the registered run, bit for bit"
    assert plain["results"]["configuration_fingerprint"] == off["results"]["configuration_fingerprint"] and off["results"]["soft_priors"] is None
    on = run_trainer(paths, tmp_path / "on", "--soft-prior-strength", "1.0", "--soft-prior-reports", str(paths["reports"]))
    assert not np.array_equal(plain["predictions"], on["predictions"])
    soft = on["results"]["soft_priors"]
    assert soft["applied_in_fit"]["positive_cells"] > 0 and soft["applied_in_fit"]["absence_cells"] > 0 and soft["permutation"] is None
    assert soft["applied_in_refit"]["fit_rows"] == on["results"]["num_train"] + on["results"]["num_validation"] > soft["applied_in_fit"]["fit_rows"] == on["results"]["num_train"]
    assert set(soft["by_class_and_direction"]) == {"animal:negative", "animal:positive", "literature:positive", "monogenic:positive"}
    assert "macro_auprc_by_literature_popularity_bin" in on["results"] and len(on["results"]["macro_auprc_by_literature_popularity_bin"]) == 5
    assert on["results"]["arguments"]["soft_prior_reports"] == [str(paths["reports"])]
    permuted = run_trainer(paths, tmp_path / "permuted", "--soft-prior-strength", "1.0", "--soft-prior-reports", str(paths["reports"]), "--permute-soft-priors", "--permute-labels")
    record = permuted["results"]["soft_priors"]["permutation"]
    assert record["seed"] == 0 and len(record["layout_sha256"]) == 64 and record["pairs"] > 0 and permuted["results"]["labels_permuted"]
    calibrated = run_trainer(paths, tmp_path / "calibrated", "--soft-prior-strength", "1.0", "--soft-prior-reports", str(paths["reports"]), "--soft-prior-probability", "calibrated",
                             "--soft-prior-calibration-minimum-cases", "2")
    calibration = calibrated["results"]["soft_priors"]["calibration"]
    assert set(calibration) == {"fit", "refit"} and calibration["refit"]["training_cells"] >= calibration["fit"]["training_cells"]
    assert calibrated["results"]["soft_priors"]["probability_source"] == "calibrated"
    if not calibration["fit"]["kept_tiers"]:
        assert "log_popularity" in calibration["fit"]["coefficients"] and calibration["fit"]["mean_prediction_at_median_popularity"]
        assert not np.array_equal(on["predictions"], calibrated["predictions"]), "a calibrated target moves the fit"
    tiers_kept = run_trainer(paths, tmp_path / "kept", "--soft-prior-strength", "1.0", "--soft-prior-reports", str(paths["reports"]), "--soft-prior-probability", "calibrated",
                             "--soft-prior-calibration-minimum-cases", "100000")
    kept = tiers_kept["results"]["soft_priors"]["calibration"]
    assert kept["fit"]["kept_tiers"] and "calibration cells" in kept["fit"]["reason"]
    assert np.array_equal(on["predictions"], tiers_kept["predictions"]), "a calibration that keeps the tiers is the tier run, bit for bit"


def test_a_drug_row_rescaled_by_equal_shares_scales_its_priors_the_same_way():
    from experiments.run_main_model import negative_share_scale

    outcomes = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    eligible = outcomes == 0
    weights = np.array([[1.0, 0.4, 0.4], [0.2, 0.2, 1.0], [0.0, 1.0, 0.0]])  # row 0's negatives were doubled; row 2 has no weighted negative
    assert negative_share_scale(weights, eligible, 0.2).tolist() == [2.0, 1.0, 1.0]
    probability = np.array([[np.nan, 0.7, np.nan], [0.7, np.nan, np.nan], [0.7, np.nan, np.nan]])
    strength = np.array([[0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])
    fitted, _ = fit_targets_with_soft_priors(torch.as_tensor(outcomes, dtype=torch.float32), torch.as_tensor(weights, dtype=torch.float32), None, np.arange(3),
                                             probability, strength, eligible, 1.0, row_scale=negative_share_scale(weights, eligible, 0.2))
    assert np.allclose(fitted["weights"].numpy()[[0, 1, 2], [1, 0, 0]], [1.0, 0.5, 0.5]), "row 0's prior is doubled with its negatives"
