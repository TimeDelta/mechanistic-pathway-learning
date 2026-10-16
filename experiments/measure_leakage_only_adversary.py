"""What can a predictor earn from the open leakage channel alone, under the item 25 split? (Open decision 9.35.)

The owner, 14 October 2026: "build the adversary. it's a good control", and then "you can remove the leakage ceiling
construction since that adversary is being used instead". So this script measures a control and sets no threshold: it
reports what a lookup-only predictor earns in the units of the registered reading, and the channel counts behind it.

What it does. Over --num-draws pseudo-draws of the item 25 rule on the development data, each draw gives a held-out set
and a training set. The adversary is fitted on the training labels and scores the held-out rows from the training
perturbations that share a seed node with them, nothing else. Its advantage over each registered baseline is read in
macro and micro AUPRC with the paired bootstrap of design section 7, the same reading the confirmatory test uses, so the
number is comparable with H1's 0.05 requirement and with the measured precision of the readings.

Why pseudo-draws and not the real one. The adversary needs the held-out labels to be *scored*, which the lockbox
forbids before `runs/confirmatory/SCORED` exists. Development pseudo-draws under the same rule need no lockbox label, so
this runs before the draw and tells the owner whether the channel is worth anything at all - which is what they would
otherwise learn only after the one scoring they get.

Usage:
  python experiments/measure_leakage_only_adversary.py --num-draws 20
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.evaluation.leakage_only_adversary import (
    LeakageOnlyAdversary, POOLING_RULES, channel_breakdown, leaked_neighbours, seed_nodes_of_perturbations,
    seen_target_share)
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    macro_statistic, micro_auprc, paired_bootstrap_macro_difference, paired_bootstrap_micro_difference)
from mechanistic_pathway_learning.evaluation.target_held_out_split import REGISTERED_HELD_OUT_FRACTION, held_out_perturbations, targets_of_drugs
from mechanistic_pathway_learning.models.baselines.popularity_baseline import PopularityBaseline

# a 95 percent bootstrap interval is about this many standard errors wide on each side, used only to turn the interval
# the registered machinery returns into the standard error the readings are compared against
STANDARD_ERRORS_IN_A_95_PERCENT_HALF_WIDTH = 1.959964


def registered_two_sided_draw(data, draw_seed: int) -> tuple[set[str], set[str]]:
    """The held-out set as the registration actually draws it: **each side at its own registered fraction** (item 25,
    "the genes by the registered disease-cluster grouping and the drugs by the rule below, one seventh of each rather
    than one seventh of the whole").

    Corrected 14 October 2026 after the owner asked "i thought gene perturbations were also being included in the
    lockbox draw ...". They are, and an earlier version of this script drew the drug side only: `held_out_perturbations`
    pulls a gene in only when its own node is a held-out drug **target**, and no gene node is a drug target on a
    gene/protein-split graph, so every pseudo-draw held out 39, 22 or 30 drugs and **zero** genes. Since 1,414 of the
    1,636 perturbations are genes, the adversary was being scored on 2 percent of the row mix it will meet.
    """
    held_out_ids, held_out_targets = held_out_perturbations(
        data.perturbation_ids, data.perturbation_types, data.perturbation_seeds, data.node_ids, draw_seed,
        REGISTERED_HELD_OUT_FRACTION)
    num_grouped_folds = int(round(1 / REGISTERED_HELD_OUT_FRACTION))
    fold_of_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, num_grouped_folds, draw_seed)
    for perturbation_id, kind in zip(data.perturbation_ids, data.perturbation_types):
        if kind != "drug" and fold_of_perturbation[perturbation_id] == 0:
            held_out_ids.add(str(perturbation_id))
    return held_out_ids, held_out_targets


def one_draw(data, draw_seed: int, pooling_rule: str, num_bootstrap: int) -> dict:
    """One pseudo-draw: the adversary's advantage over popularity and degree-scaled popularity, and the channel counts."""
    held_out_ids, held_out_targets = registered_two_sided_draw(data, draw_seed)
    held_out_rows = np.array([index for index, perturbation_id in enumerate(data.perturbation_ids) if str(perturbation_id) in held_out_ids], dtype=int)
    training_rows = np.array([index for index, perturbation_id in enumerate(data.perturbation_ids) if str(perturbation_id) not in held_out_ids], dtype=int)
    if held_out_rows.size == 0 or training_rows.size == 0:
        return {"draw_seed": draw_seed, "usable": False}

    held_out_identifiers = [str(data.perturbation_ids[row]) for row in held_out_rows]
    training_identifiers = [str(data.perturbation_ids[row]) for row in training_rows]
    label_mask = data.label_mask if data.label_mask is not None else np.ones(data.outcomes.shape, dtype=bool)

    seed_nodes_of_perturbation = seed_nodes_of_perturbations(data.perturbation_ids, data.perturbation_seeds, data.node_ids)
    neighbours = leaked_neighbours(held_out_identifiers, training_identifiers, seed_nodes_of_perturbation)
    perturbation_types_of_id = {str(perturbation_id): kind for perturbation_id, kind in zip(data.perturbation_ids, data.perturbation_types)}

    adversary = LeakageOnlyAdversary(pooling_rule).fit(training_identifiers, data.outcomes[training_rows], label_mask[training_rows])
    adversary_scores = adversary.predict(held_out_identifiers, neighbours)

    held_out_outcomes, held_out_mask = data.outcomes[held_out_rows], label_mask[held_out_rows]
    held_out_degrees = np.array([float(data.node_degree[data.perturbation_seeds[row]].sum()) for row in held_out_rows])
    training_degrees = np.array([float(data.node_degree[data.perturbation_seeds[row]].sum()) for row in training_rows])

    macro_average = macro_statistic()
    result = {"draw_seed": draw_seed, "usable": True, "pooling_rule": pooling_rule,
              "held_out_perturbations": int(held_out_rows.size), "held_out_targets": int(len(held_out_targets)),
              "held_out_drugs": int(sum(1 for row in held_out_rows if data.perturbation_types[row] == "drug")),
              "held_out_genes": int(sum(1 for row in held_out_rows if data.perturbation_types[row] != "drug")),
              "seen_target_share": seen_target_share(
                  targets_of_drugs(data.perturbation_ids, data.perturbation_types, data.perturbation_seeds, data.node_ids), held_out_identifiers),
              "adversary_macro_auprc": float(macro_average(adversary_scores, held_out_outcomes, held_out_mask)),
              "adversary_micro_auprc": float(micro_auprc(adversary_scores, held_out_outcomes, held_out_mask))}
    result.update(channel_breakdown(neighbours, perturbation_types_of_id))

    # The drug stratum, approved by the owner 14 October 2026 ("1. yes"). The whole-lockbox reading is diluted: a draw
    # holds out about 170 genes to 22-39 drugs, no gene shares a seed node with anything, and the adversary falls back
    # to the base rate on every gene row, so its whole-lockbox macro and micro sit almost on top of popularity whatever
    # the drug side is doing. The channel exists only among drugs, so the bound has force only where it is read.
    held_out_is_a_drug = np.array([data.perturbation_types[row] == "drug" for row in held_out_rows], dtype=bool)
    strata = {"all_held_out": np.ones(len(held_out_rows), dtype=bool), "drug_stratum": held_out_is_a_drug,
              "gene_stratum": ~held_out_is_a_drug}

    # type_popularity is the comparator the pooled reading needs, and reading against popularity alone was a defect.
    # Popularity gives every row the same score for a symptom, so it cannot order a drug above a gene; the adversary's
    # score differs from the base rate on drug rows only, so it orders drugs first, and with the drug base rate well
    # above the gene one that earns a large AUPRC gain for KNOWING WHICH PERTURBATIONS ARE DRUGS rather than for any
    # lookup. Against popularity the pooled macro advantage came out at +0.139, above the drug stratum's own +0.041,
    # which is the signature of that contamination. Within the drug stratum every row is a drug, so there
    # type_popularity and popularity coincide and the stratum numbers were never affected.
    training_kinds = np.asarray([data.perturbation_types[row] for row in training_rows])
    held_out_kinds = np.asarray([data.perturbation_types[row] for row in held_out_rows])

    for baseline_name, scale_by_degree in (("popularity", False), ("degree_popularity", True), ("type_popularity", False)):
        if baseline_name == "type_popularity":
            baseline_scores = np.zeros((len(held_out_rows), data.outcomes.shape[1]))
            for kind in np.unique(held_out_kinds):
                of_this_kind_in_training = training_kinds == kind
                if not of_this_kind_in_training.any():
                    of_this_kind_in_training = np.ones(len(training_rows), dtype=bool)
                fitted_on_kind = PopularityBaseline().fit(
                    data.outcomes[training_rows][of_this_kind_in_training], None,
                    label_mask[training_rows][of_this_kind_in_training])
                baseline_scores[held_out_kinds == kind] = fitted_on_kind.predict(int((held_out_kinds == kind).sum()))
        else:
            baseline = PopularityBaseline(scale_by_degree=scale_by_degree).fit(
                data.outcomes[training_rows], training_degrees, label_mask[training_rows])
            baseline_scores = baseline.predict(len(held_out_rows), held_out_degrees)
        for stratum_name, in_stratum in strata.items():
            if in_stratum.sum() < 2:
                continue
            prefix = "" if stratum_name == "all_held_out" else f"{stratum_name}_"
            macro = paired_bootstrap_macro_difference(adversary_scores[in_stratum], baseline_scores[in_stratum],
                                                      held_out_outcomes[in_stratum], num_bootstrap=num_bootstrap,
                                                      random_seed=draw_seed, mask=held_out_mask[in_stratum])
            micro = paired_bootstrap_micro_difference(adversary_scores[in_stratum], baseline_scores[in_stratum],
                                                      held_out_outcomes[in_stratum], num_bootstrap=num_bootstrap,
                                                      random_seed=draw_seed, mask=held_out_mask[in_stratum])
            for reading_name, reading in (("macro", macro), ("micro", micro)):
                half_width = (reading["upper"] - reading["lower"]) / 2.0
                result[f"{prefix}{reading_name}_advantage_over_{baseline_name}"] = reading["difference"]
                result[f"{prefix}{reading_name}_interval_over_{baseline_name}"] = [reading["lower"], reading["upper"]]
                result[f"{prefix}{reading_name}_standard_error_over_{baseline_name}"] = half_width / STANDARD_ERRORS_IN_A_95_PERCENT_HALF_WIDTH
    return result


def summarise(draws: list[dict]) -> dict:
    """Means and ranges over the usable draws, per reading."""
    usable = [draw for draw in draws if draw.get("usable")]
    summary = {"usable_draws": len(usable)}
    if not usable:
        return summary
    for key in sorted(usable[0]):
        values = [draw[key] for draw in usable if isinstance(draw.get(key), (int, float)) and not isinstance(draw.get(key), bool)]
        if len(values) == len(usable) and values:
            finite = [value for value in values if np.isfinite(value)]
            if finite:
                summary[key] = {"mean": float(np.mean(finite)), "minimum": float(np.min(finite)), "maximum": float(np.max(finite)),
                                "draws": len(finite)}
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
    parser.add_argument("--group-by", default="disease_cluster_and_targets")
    parser.add_argument("--num-draws", type=int, default=20)
    parser.add_argument("--pooling-rules", nargs="+", default=list(POOLING_RULES), choices=list(POOLING_RULES),
                        help="which pooling rules to measure; the default measures all of them, which is the cross "
                             "product this script is slowest on")
    parser.add_argument("--checkpoint-dir", type=Path, default=Path("runs/leakage_only_adversary_draws"),
                        help="one JSON per (pooling rule, draw). A draw already written here is skipped, so a run "
                             "stopped by the background time limit resumes where it left off instead of starting over")
    parser.add_argument("--num-bootstrap", type=int, default=1000)
    parser.add_argument("--lockbox", type=Path, default=None,
                        help="a drawn lockbox whose perturbations are excluded, so the calibration reads development rows only")
    parser.add_argument("--json-output", type=Path, default=Path("docs/leakage_only_adversary.json"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/leakage_only_adversary.md"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by,
                               label_selection=arguments.label_selection)
    excluded_lockbox_perturbations = set()
    if arguments.lockbox is not None:
        excluded_lockbox_perturbations = {str(perturbation_id) for perturbation_id in json.loads(arguments.lockbox.read_text())["perturbation_ids"]}
        keep = [index for index, perturbation_id in enumerate(data.perturbation_ids) if str(perturbation_id) not in excluded_lockbox_perturbations]
        data = _restrict(data, keep)

    results = {"pooling_rules": {}, "excluded_lockbox_perturbations": len(excluded_lockbox_perturbations),
               "evidence_dir": str(arguments.evidence_dir), "graph_dir": str(arguments.graph_dir),
               "label_selection": str(arguments.label_selection), "num_draws": arguments.num_draws}
    arguments.checkpoint_dir.mkdir(parents=True, exist_ok=True)
    completed_draws = 0
    for pooling_rule in arguments.pooling_rules:
        draws = []
        for draw_seed in range(arguments.num_draws):
            checkpoint_path = arguments.checkpoint_dir / f"{pooling_rule}_draw_{draw_seed:03d}.json"
            if checkpoint_path.exists():
                draws.append(json.loads(checkpoint_path.read_text()))
                continue
            draw = one_draw(data, draw_seed, pooling_rule, arguments.num_bootstrap)
            checkpoint_path.write_text(json.dumps(draw, indent=1) + "\n")
            draws.append(draw)
            completed_draws += 1
            print(f"  {pooling_rule} draw {draw_seed} done and checkpointed", flush=True)
        results["pooling_rules"][pooling_rule] = {"draws": draws, "summary": summarise(draws)}
        print(pooling_rule, json.dumps(results["pooling_rules"][pooling_rule]["summary"], indent=1), flush=True)
    print(f"draws computed this invocation: {completed_draws}; re-run to continue if the run was cut short")

    arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
    arguments.json_output.write_text(json.dumps(results, indent=1) + "\n")
    arguments.markdown_output.write_text(_markdown(results, arguments) + "\n")
    print(f"wrote {arguments.json_output} and {arguments.markdown_output}")


def _restrict(data, keep: list[int]):
    """A copy of the experiment data holding only the given perturbation rows."""
    import dataclasses
    per_perturbation_lists = ("perturbation_ids", "perturbation_labels", "perturbation_types", "group_ids",
                              "perturbation_seeds", "perturbation_signs", "perturbation_magnitudes")
    changes = {name: [getattr(data, name)[index] for index in keep] for name in per_perturbation_lists}
    for name in ("outcomes", "weights", "label_mask", "frequencies", "evidence_dates", "evidence_date_is_publication",
                 "evidenced_negative", "placebo_negative_weight_matrix"):
        value = getattr(data, name, None)
        if value is not None:
            changes[name] = np.asarray(value)[keep]
    return dataclasses.replace(data, **changes)


def _markdown(results: dict, arguments) -> str:
    lines = ["# The leakage-only adversary (generated by experiments/measure_leakage_only_adversary.py)", "",
             f"Evidence {arguments.evidence_dir}, label selection {arguments.label_selection}, graph {arguments.graph_dir}; "
             f"{arguments.num_draws} pseudo-draws of the item 25 rule at the registered held-out fraction, "
             f"{arguments.num_bootstrap} bootstrap resamples per reading.", "",
             "The adversary knows two things: which training perturbations share a seed node with the held-out "
             "perturbation, and the training labels. No edge, no propagation, no descriptor. It is handed the leaked "
             "labels directly, so its advantage is an upper bound on what any model could take from this channel by "
             "lookup. A held-out perturbation with no leaked neighbour takes the training base rate, which keeps the "
             "advantage a measure of the leak rather than a mixture of a leak and a handicap.", "",
             "Read against H1's registered requirement of a 0.05 macro and micro difference, and against the measured "
             "precision of the readings (paired-bootstrap 95 percent half-widths of 0.019 to 0.063 macro and 0.018 to "
             "0.041 micro, growing by about 1.4 at Holm's smallest level).", ""]
    for pooling_rule, block in results["pooling_rules"].items():
        summary = block["summary"]
        lines += [f"## {pooling_rule}", "", f"Usable draws: {summary.get('usable_draws', 0)}.", "",
                  "| reading | mean | minimum | maximum |", "|---|---|---|---|"]
        for key in sorted(summary):
            value = summary[key]
            if isinstance(value, dict):
                lines.append(f"| {key.replace('_', ' ')} | {value['mean']:.4f} | {value['minimum']:.4f} | {value['maximum']:.4f} |")
        lines.append("")
    lines += ["## The channel this does not cover", "",
              "Two perturbations sharing a seed node perturb the same node, so a label transfers by identity and a "
              "lookup substitutes for mechanism: that is the channel measured here. Two perturbations whose seeds sit "
              "in one pathway are the hypothesis under test rather than a leak, so transfer along a route is not "
              "measured and must not be: bounding it would bound the study's own subject. A channel outside both "
              "descriptions is a disclosed limitation."]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
