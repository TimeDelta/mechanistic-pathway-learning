"""Put the comparative checks to a graph and, when finished runs are given, to their trained models.

The checks are statements that a drug's effects are at least as large with a second factor present: its carrier lowered,
a second drug displacing it from the carrier, its clearing enzyme less active, or a second drug given with it
(configs/comparative_checks.csv, experiments/build_comparative_checks.py; what each kind compares is written in
mechanistic_pathway_learning/evaluation/comparative_checks.py). Nothing is trained on them and no score of the study
reads them. They are reported here, apart from every score.

Without --run-dirs the script reports structure only: how many checks of each kind the graph can express (the drug has
its own node, and the carrier's or the enzyme's edge to it exists) and on which side of the held-out set each falls. It
loads no model and trains nothing.

With --run-dirs each finished fold is rebuilt as experiments/evaluate_on_rewired_graph.py rebuilds it (the arguments of
its results.json, its best weights) and asked every check that its data can express. Its data are the ones it was
trained on, with the lockbox removed, so a check on a held-out drug is not asked. Runs on the lockbox, on permuted
labels and on rewired graphs are refused there, and so here. A fold writes comparative_checks.json beside its results.

Two things to keep in mind when reading a result.

  - In the linear-response encoder a lowered carrier or enzyme reaches the drug through one signed edge with one gain,
    so "at least as high" holds or fails for every drug at once: it tests the sign the model learned for the relation,
    once, and not each statement.
  - A pair of drugs that share a target with the same sign adds inputs at one node. Any encoder that is monotone in its
    input passes that comparison; the pairs that share no target are the informative ones.

Usage:
  PYTHONPATH=. python experiments/run_comparative_checks.py
  PYTHONPATH=. python experiments/run_comparative_checks.py --run-dirs runs/<a finished run> --markdown-output none
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from mechanistic_pathway_learning.evaluation.comparative_checks import (
    CHECK_KINDS,
    KINDS_THAT_READ_A_TRAINED_RELATION,
    ComparativeCheck,
    evaluate_checks,
    load_checks,
    side_of,
    why_not_expressible,
)
from mechanistic_pathway_learning.evaluation.experiment_data import (
    load_experiment_data,
    read_lockbox,
    read_protein_of_gene,
    restrict_to_perturbations,
)

OUTPUT_NAME = "comparative_checks.json"
KIND_TEXT = {"carrier_level": "carrier lowered with the drug", "enzyme_activity": "clearing enzyme lowered with the drug",
             "displacement": "two drugs, carrier lowered", "second_drug": "two drugs against one"}


def predict_seed_lists(encoder, head, data, seed_lists: list[list[tuple[int, float, float]]], positions: list[int], adjacencies, arguments, device) -> np.ndarray:
    """Symptom probabilities of arbitrary seed lists, with the encoder and the head called as experiments/run_main_model.py
    calls them. positions names, for each list, the perturbation whose degree covariate the head reads under --degree-offset."""
    from experiments.run_main_model import encode, perturbation_covariate

    width = max(1, max(len(seeds) for seeds in seed_lists))
    node_index = torch.full((len(seed_lists), width), -1, dtype=torch.long)
    sign_and_magnitude = torch.zeros((len(seed_lists), width, 2))
    for row, seeds in enumerate(seed_lists):
        for column, (node, sign, magnitude) in enumerate(seeds):
            node_index[row, column] = node
            sign_and_magnitude[row, column, 0], sign_and_magnitude[row, column, 1] = sign, magnitude
    covariate = torch.as_tensor(perturbation_covariate(data)[np.asarray(positions)], device=device) if arguments.degree_offset else None
    encoder.eval()
    head.eval()
    predictions = []
    with torch.no_grad():
        for start in range(0, len(seed_lists), arguments.batch_size):
            rows = slice(start, start + arguments.batch_size)
            field = encode(encoder, node_index[rows].to(device), sign_and_magnitude[rows].to(device), adjacencies, arguments.field)
            predictions.append(head(field, relation_index=0, perturbation_covariate=None if covariate is None else covariate[rows]).symptom_probability.cpu().numpy())
    return np.concatenate(predictions, axis=0)


def structure_table(checks: list[ComparativeCheck], data, protein_of_gene: dict[str, list[str]], held_out_perturbation_ids: set[str] | None) -> pd.DataFrame:
    """One row per check: whether this graph can express it, the reason when it cannot, and its side of the held-out set."""
    return pd.DataFrame([{"statement_id": check.statement_id, "kind": check.kind, "expectation": check.expectation,
                          "reason": (reason := why_not_expressible(check, data, protein_of_gene)), "expressible": not reason,
                          "reads_a_trained_relation": check.reads_a_trained_relation,
                          "side": "not assigned" if held_out_perturbation_ids is None else side_of(check, held_out_perturbation_ids)} for check in checks])


def structure_lines(structure: pd.DataFrame, graph_directory: Path, evidence_directory: Path, checks_path: Path, drug_entry: str, held_out_path: Path | None = None) -> list[str]:
    lines = ["# Comparative checks: what the graph can express", "",
             "Generated by `experiments/run_comparative_checks.py`. Do not edit by hand.", "",
             f"Checks: `{checks_path}` ({len(structure)}). Graph: `{graph_directory}`, drugs seeded on their {'own nodes' if drug_entry == 'nodes' else 'targets'}. Table: `{evidence_directory}`.",
             "No model is loaded and nothing is trained. A check is reported apart from every score of the study.", "",
             "A check is a test of a model only where it reads a relation that some label trains, which is the two-drug kind alone: a second drug enters through its",
             "mechanism edges, the relation every drug perturbation trains. A lowered carrier or enzyme reaches the drug through `sequesters` or `clears`, and no label",
             "of the study varies either, so in the message-passing encoder the check reads untrained weights and in the linear-response encoder it reads a direction",
             "that is built into the edge's sign. Those three kinds are reported as a description of what the model does, with no verdict (the user's objection of",
             "12 October 2026; `docs/label_source_review.md`, 13.7).", "",
             "| kind | what is compared | tests a trained relation | checks | expressible | development | held out |", "| --- | --- | --- | --- | --- | --- | --- |"]
    for kind in CHECK_KINDS:
        for expectation, label in (("at_least_as_high", ""), ("none", ", negative controls")):
            of_kind = structure[(structure.kind == kind) & (structure.expectation == expectation)]
            if not len(of_kind):
                continue
            expressible = of_kind[of_kind.expressible]
            lines.append(f"| {kind}{label} | {KIND_TEXT[kind]} | {'yes' if kind in KINDS_THAT_READ_A_TRAINED_RELATION else 'no'} | {len(of_kind)} | {len(expressible)}"
                         f" | {int((expressible.side == 'development').sum())} | {int((expressible.side == 'held out').sum())} |")
    reasons = Counter(structure[~structure.expressible].reason.str.replace(r"^\S+ is not a perturbation of this table$", "a drug of the statement is not in this table", regex=True))
    if reasons:
        lines += ["", "Not expressible, by reason:", ""] + [f"- {count}: {reason}" for reason, count in reasons.most_common()]
    if held_out_path is None:
        lines += ["", "No list of held-out perturbations was given, so no check is assigned a side."]
    else:
        lines += ["", f"Sides are read from `{held_out_path}` (membership only). A check is held out when one of its drugs is, or its carrier or enzyme where that gene is itself a perturbation of the study.",
                  "A finished run is asked the development checks only: its data have the lockbox removed."]
    return lines + [""]


def evaluate_fold(fold_directory: Path, checks: list[ComparativeCheck], data_cache: dict) -> dict:
    """The checks put to one finished fold; writes comparative_checks.json beside its results and returns its content."""
    from experiments.evaluate_on_rewired_graph import refuse_unsupported, stored_arguments, trained_models

    results = json.loads((fold_directory / "results.json").read_text())
    arguments = stored_arguments(results)
    refuse_unsupported(results, arguments, fold_directory)
    data_key = (str(arguments.graph_dir), str(arguments.evidence_dir), arguments.group_by, tuple(arguments.label_grades or ()), str(arguments.label_selection), str(arguments.lockbox),
                arguments.metabolic_layer_only, getattr(arguments, "drug_entry", None))
    if data_key not in data_cache:
        # the data the run was trained on, as the trainer loads it: the lockbox removed. A check whose drug is in the
        # lockbox is then "not a perturbation of this table" and is not asked, which keeps the lockbox blind
        data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only, group_by=arguments.group_by,
                                    label_grades=tuple(arguments.label_grades) if arguments.label_grades else None, label_selection=arguments.label_selection,
                                    drug_entry=getattr(arguments, "drug_entry", None) or "nodes")
        if arguments.lockbox is not None:
            data = restrict_to_perturbations(data, ~read_lockbox(arguments.lockbox, data, arguments.group_by, arguments.evidence_dir))
        data_cache[data_key] = data
    data = data_cache[data_key]
    device = torch.device("cpu")
    refit = (fold_directory / "refit_checkpoint.pt").exists() and bool(getattr(arguments, "refit_on_validation", False))
    if refit:
        refit_checkpoint = torch.load(fold_directory / "refit_checkpoint.pt", map_location=device, weights_only=False)
        checkpoint = {"state": {}, "encoder": refit_checkpoint["encoder"], "head": refit_checkpoint["head"]}
    else:
        checkpoint = torch.load(fold_directory / "checkpoint.pt", map_location=device, weights_only=False)
    torch.manual_seed(arguments.seed)
    encoder, head, adjacencies = trained_models(data, arguments, checkpoint, device)
    protein_of_gene = read_protein_of_gene(arguments.graph_dir)
    table = evaluate_checks(checks, data, lambda seed_lists, positions: predict_seed_lists(encoder, head, data, seed_lists, positions, adjacencies, arguments, device), protein_of_gene)
    judged = table[table.expressible & (table.expectation == "at_least_as_high") & table.holds.notna()]
    summary = {"fold": str(fold_directory), "weights": "refit" if refit else "best epoch", "checks": len(table), "expressible": int(table.expressible.sum()),
               "kinds_that_read_a_trained_relation": list(KINDS_THAT_READ_A_TRAINED_RELATION),
               "by_kind": {kind: {"judged": int((judged.kind == kind).sum()), "hold": int(judged[judged.kind == kind].holds.astype(bool).sum()),
                                  "mean_difference": float(judged[judged.kind == kind].mean_difference.mean()) if (judged.kind == kind).any() else None} for kind in CHECK_KINDS},
               "rows": json.loads(table.to_json(orient="records"))}
    (fold_directory / OUTPUT_NAME).write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--checks", type=Path, default=Path("configs/comparative_checks.csv"))
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes"), help="for the structure report")
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
    parser.add_argument("--drug-entry", choices=("nodes", "targets"), default="nodes")
    parser.add_argument("--held-out-perturbations", type=Path, default=None,
                        help="a table with a perturbation_id column naming the held-out perturbations (membership only; never a prediction). Without it no check is assigned a side")
    parser.add_argument("--run-dirs", type=Path, nargs="*", default=[], help="finished runs of experiments/run_main_model.py whose folds are asked the checks")
    parser.add_argument("--markdown-output", default="docs/comparative_checks.md", help="the structure report; none writes no file")
    arguments = parser.parse_args()

    checks = load_checks(arguments.checks)
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster", label_selection=arguments.label_selection, drug_entry=arguments.drug_entry)
    held_out = None
    if arguments.held_out_perturbations is not None:
        held_out_table = pd.read_csv(arguments.held_out_perturbations)
        held_out = set(held_out_table.perturbation_id[held_out_table.held_out.astype(bool)] if "held_out" in held_out_table.columns else held_out_table.perturbation_id)
    structure = structure_table(checks, data, read_protein_of_gene(arguments.graph_dir), held_out)
    lines = structure_lines(structure, arguments.graph_dir, arguments.evidence_dir, arguments.checks, arguments.drug_entry, arguments.held_out_perturbations)
    print("\n".join(lines[6:]))
    if str(arguments.markdown_output).lower() != "none":
        Path(arguments.markdown_output).write_text("\n".join(lines))
        print(f"wrote {arguments.markdown_output}")

    data_cache: dict = {}
    for run_directory in arguments.run_dirs:
        for results_path in sorted(run_directory.rglob("results.json")):
            summary = evaluate_fold(results_path.parent, checks, data_cache)
            print(f"{results_path.parent}: {summary['expressible']} of {summary['checks']} checks expressible; " + "; ".join(
                f"{kind} {by_kind['hold']} of {by_kind['judged']} hold" for kind, by_kind in summary["by_kind"].items() if by_kind["judged"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
