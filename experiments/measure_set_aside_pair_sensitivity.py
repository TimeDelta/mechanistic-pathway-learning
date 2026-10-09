"""Metric-side sensitivity of the better-label selection: how the model does on the pairs the rule set aside.

The selection rule of experiments/build_label_selection.py keeps a gene pair at a reported frequency of 0.30 and a
drug pair at 0.01 or when both SIDER and OnSIDES list it, and masks the rest: a set-aside pair is neither a positive
nor a negative, so it leaves both the loss and the metric. The user asked what that costs (9 October 2026: "i was
thinking both the metric and training sensitivities, each stratified by reason for dropping and only report the
model-minus-baseline gap within each set ... the whole point is to determine the effects on the model").

This is the metric arm, and it needs no new training: the trainer saves a prediction for every pair of its held-out
perturbations, masked or not, so a set-aside pair already has a score. Per stratum the pairs of that stratum are
scored as positives against the pairs no source reports, and only the model-minus-baseline gap is reported, because
the level of an AUPRC inside a stratum is set by that stratum's prevalence and is not comparable across strata.

Reading it. A gap on a set-aside stratum close to the gap on the kept positives says the rule is not flattering the
model: the pairs it removed are ones the model orders as well as the ones it kept. A gap near zero says the model has
nothing to say about those pairs, which is a reason to leave them out. A gap above the kept one would say the rule is
removing the pairs the model is best at, which would make the headline number pessimistic rather than optimistic.

The pairs no source reports are used as the negatives, as the study's main metric does; they are unlabelled rather
than known negatives, which is why the gap rather than the level is the reading.

The training arm (train with a stratum admitted as positives, score on the kept metric) is a separate run and is
specified in docs/set_aside_pair_sensitivity.md rather than started here.

Writes docs/set_aside_pair_sensitivity.md. Reads no lockbox outcome: the runs it reads are development runs with the
lockbox perturbations removed. Idempotent.

Usage:
  python experiments/measure_set_aside_pair_sensitivity.py --run-dir runs/full/<run> --fold fold0_seed0
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

FREQUENCY_IN_A_REASON = re.compile(r"frequency \d+(\.\d+)?")
KEPT_STRATUM = "kept positives"
# the comparator the preregistration names for the confirmatory readings: "Against the better of popularity and
# degree_popularity in each draw"
ENDPOINT_BASELINES = ("popularity", "degree_popularity")
BEST_ENDPOINT_BASELINE = "the better of popularity and degree_popularity"
MINIMUM_POSITIVES_IN_A_STRATUM = 10
MINIMUM_POSITIVES_FOR_A_SYMPTOM = 5
OUTPUT_DOCUMENT = Path("docs/set_aside_pair_sensitivity.md")


def reason_family(reason: str) -> str:
    """The reason a pair was set aside, with the frequency value removed so the strata are families, not values."""
    return FREQUENCY_IN_A_REASON.sub("frequency", str(reason or ""))


def pair_strata(selection: pd.DataFrame, perturbation_ids: list[str], symptoms: list[str]) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """One boolean [perturbations, symptoms] mask per stratum, and the mask of pairs no source reports."""
    row_of = {perturbation_id: index for index, perturbation_id in enumerate(perturbation_ids)}
    column_of = {symptom: index for index, symptom in enumerate(symptoms)}
    shape = (len(perturbation_ids), len(symptoms))
    reported = np.zeros(shape, dtype=bool)
    strata: dict[str, np.ndarray] = {}
    for row in selection.itertuples(index=False):
        position, column = row_of.get(row.perturbation_id), column_of.get(row.symptom)
        if position is None or column is None:
            continue
        reported[position, column] = True
        stratum = KEPT_STRATUM if row.keep else reason_family(row.reason)
        strata.setdefault(stratum, np.zeros(shape, dtype=bool))[position, column] = True
    return strata, ~reported


def gaps_within_a_stratum(stratum: np.ndarray, unreported: np.ndarray, model_predictions: np.ndarray,
                          baseline_predictions: dict[str, np.ndarray], symptoms: list[str]) -> dict:
    """Micro and macro AUPRC of the model and each baseline on one stratum's pairs against the unreported pairs."""
    scored = stratum | unreported
    truth = stratum[scored].astype(int)
    reading = {"positives": int(stratum.sum()), "negatives": int(unreported.sum()),
               "prevalence": float(truth.mean()) if len(truth) else float("nan")}
    if reading["positives"] < MINIMUM_POSITIVES_IN_A_STRATUM or not truth.any():
        return reading | {"micro": {}, "macro": {}, "symptoms_in_the_macro": 0}

    def micro(predictions: np.ndarray) -> float:
        return float(average_precision_score(truth, predictions[scored]))

    symptom_columns = [index for index in range(len(symptoms)) if stratum[:, index].sum() >= MINIMUM_POSITIVES_FOR_A_SYMPTOM]

    def macro(predictions: np.ndarray) -> float:
        per_symptom = []
        for index in symptom_columns:
            rows = scored[:, index]
            column_truth = stratum[rows, index].astype(int)
            if column_truth.any() and not column_truth.all():
                per_symptom.append(average_precision_score(column_truth, predictions[rows, index]))
        return float(np.mean(per_symptom)) if per_symptom else float("nan")

    model_micro, model_macro = micro(model_predictions), macro(model_predictions)
    baseline_micro = {name: micro(predictions) for name, predictions in baseline_predictions.items()}
    baseline_macro = {name: macro(predictions) for name, predictions in baseline_predictions.items()}
    reading["micro"] = {name: model_micro - value for name, value in baseline_micro.items()}
    reading["macro"] = {name: model_macro - value for name, value in baseline_macro.items()}
    endpoint_micro = [baseline_micro[name] for name in ENDPOINT_BASELINES if name in baseline_micro]
    endpoint_macro = [baseline_macro[name] for name in ENDPOINT_BASELINES if name in baseline_macro]
    if endpoint_micro:
        reading["micro"][BEST_ENDPOINT_BASELINE] = model_micro - max(endpoint_micro)
        reading["macro"][BEST_ENDPOINT_BASELINE] = model_macro - max(endpoint_macro)
    reading["symptoms_in_the_macro"] = len(symptom_columns)
    return reading


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", type=Path, nargs="+", required=True, help="one or more run directories; each gets a section")
    parser.add_argument("--fold", default="fold0_seed0")
    parser.add_argument("--baselines-dir", type=Path, default=Path("runs/full/baselines_v2_development"))
    parser.add_argument("--selection", type=Path, default=None,
                        help="default: the label selection the run was trained with, from its results.json")
    parser.add_argument("--markdown-output", type=Path, default=OUTPUT_DOCUMENT)
    arguments = parser.parse_args()

    sections, baseline_names = [], []
    for run_directory in arguments.run_dir:
        section, baseline_names = section_for_a_run(run_directory, arguments)
        sections.append(section)

    lines = [
        "# The pairs the selection rule set aside: what the model does with them",
        "",
        "Generated by `experiments/measure_set_aside_pair_sensitivity.py`. Do not edit by hand.",
        "",
        f"Fold `{arguments.fold}`, baselines `{arguments.baselines_dir}`. Each stratum's pairs are scored as "
        "positives against the pairs no source reports; the other strata's pairs are left out of that comparison. "
        "Only the model-minus-baseline gap is reported, because the level of an AUPRC inside a stratum is set by its "
        "prevalence. The pairs no source reports are the negatives, as the study's main metric does; they are "
        "unlabelled rather than known negatives, which is the second reason to read the gap and not the level.",
        "",
    ]
    for section in sections:
        lines += section
    lines += [
        f"A stratum with fewer than {MINIMUM_POSITIVES_IN_A_STRATUM} pairs is left unscored, and the macro average "
        f"takes the symptoms with at least {MINIMUM_POSITIVES_FOR_A_SYMPTOM} of the stratum's pairs.",
        "",
        "## What the training arm needs, which this script does not do",
        "",
        "The metric arm asks whether the model already orders the set-aside pairs. The training arm asks whether "
        "admitting them would change what it learns, and that needs one run per stratum: the same configuration with "
        "that stratum's pairs admitted as positives in the loss, scored on the unchanged kept metric, reported as the "
        "model-minus-baseline gap on the same rows as the run that did not admit them. Those runs are training runs, "
        "and no new pilot is started while the user's instruction of 9 October 2026 (\"No new pilots yet\") stands.",
    ]
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(f"wrote {arguments.markdown_output}")
    return 0


def section_for_a_run(run_directory: Path, arguments) -> tuple[list[str], list[str]]:
    """One run's strata table, and the baselines it was compared with."""
    fold_directory = run_directory / arguments.fold
    results = json.loads((fold_directory / "results.json").read_text())
    symptoms = results["symptoms"]
    test_perturbations = results["test_perturbation_ids"]
    model_predictions = np.load(fold_directory / "test_predictions.npy")
    selection_path = arguments.selection or Path(results["arguments"]["label_selection"])
    selection = pd.read_parquet(selection_path)

    baseline_perturbations = json.loads((arguments.baselines_dir / "perturbation_ids.json").read_text())
    baseline_results = json.loads((arguments.baselines_dir / "results.json").read_text())
    baseline_rows = [baseline_perturbations.index(perturbation_id) for perturbation_id in test_perturbations
                     if perturbation_id in baseline_perturbations]
    if len(baseline_rows) != len(test_perturbations):
        raise SystemExit(f"{len(test_perturbations) - len(baseline_rows)} test perturbations are absent from {arguments.baselines_dir}")
    baseline_columns = [baseline_results["symptoms"].index(symptom) for symptom in symptoms]
    baseline_predictions = {}
    for path in sorted(arguments.baselines_dir.glob("predictions_grouped_*.npy")):
        name = path.stem[len("predictions_grouped_"):]
        if name.startswith("label_permutation_") or name.startswith("rewired_graph_"):
            continue
        baseline_predictions[name] = np.load(path)[np.ix_(baseline_rows, baseline_columns)]

    strata, unreported = pair_strata(selection, test_perturbations, symptoms)
    readings = {name: gaps_within_a_stratum(mask, unreported, model_predictions, baseline_predictions, symptoms)
                for name, mask in sorted(strata.items(), key=lambda item: (item[0] != KEPT_STRATUM, item[0]))}

    baseline_names = sorted(baseline_predictions) + [BEST_ENDPOINT_BASELINE]
    lines = [
        f"## `{run_directory.name}`",
        "",
        f"Selection `{selection_path}`; {len(test_perturbations)} held-out perturbations, {len(symptoms)} symptoms.",
        "",
        "| stratum | pairs | prevalence | symptoms in the macro | " + " | ".join(f"micro gap vs {name}" for name in baseline_names) + " |",
        "| --- | --- | --- |" + " --- |" * (1 + len(baseline_names)),
    ]
    for name, reading in readings.items():
        cells = [f"{reading['micro'].get(baseline, float('nan')):+.3f}" if reading["micro"] else "-"
                 for baseline in baseline_names]
        lines.append(f"| {name} | {reading['positives']} | {reading['prevalence']:.3f} | "
                     f"{reading['symptoms_in_the_macro']} | " + " | ".join(cells) + " |")
    lines += ["", "| stratum | " + " | ".join(f"macro gap vs {name}" for name in baseline_names) + " |",
              "| --- |" + " --- |" * len(baseline_names)]
    for name, reading in readings.items():
        cells = [f"{reading['macro'].get(baseline, float('nan')):+.3f}" if reading["macro"] else "-"
                 for baseline in baseline_names]
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    kept_gap = readings.get(KEPT_STRATUM, {}).get("micro", {}).get(BEST_ENDPOINT_BASELINE)
    if kept_gap is not None:
        leans = []
        for name, reading in readings.items():
            gap = reading.get("micro", {}).get(BEST_ENDPOINT_BASELINE)
            if name == KEPT_STRATUM or gap is None:
                continue
            leans.append(f"{name} {gap:+.3f} ({'above' if gap > kept_gap else 'below'} the kept {kept_gap:+.3f})")
        lines += [f"Against {BEST_ENDPOINT_BASELINE}, micro: " + "; ".join(leans) + ".", ""]
    lines.append("")
    for name, reading in readings.items():
        gap = ", ".join(f"{baseline} {reading['micro'][baseline]:+.3f}" for baseline in baseline_names) if reading["micro"] else "unscored"
        print(f"{run_directory.name} | {name}: {reading['positives']} pairs; micro gap {gap}")
    return lines, baseline_names


if __name__ == "__main__":
    raise SystemExit(main())
