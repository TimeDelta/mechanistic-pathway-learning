"""Compare the neuronal and oxidative graph variant with the metabolic graph it extends, and read off what the variant's
own nodes do (design section 5.7; graph/reactome_import.py, neurotransmission_regulation.py, oxidative_regulation.py,
redox_pools.py).

Three readings, written to docs/neuronal_variant_comparison.md:
  - accuracy of the variant against the same configuration on the metabolic graph and against the variant without its
    curated receptor and oxidant layers: the per-fold macro AUPRC, the paired per-fold difference with a t interval over
    the folds both runs finished, the paired bootstrap over perturbations of the pooled out-of-fold difference (macro
    AUPRC and AUROC) on the rows both runs scored, and the same bootstrap inside degree strata, which gives no credit
    for ordering perturbations by degree (design section 5.2; aggregate_main_model_runs.within_degree_strata);
  - the predicted symptom probabilities of the genes the variant connects and the metabolic graph does not reach the
    same way: the voltage-gated channels, the pumps and cotransporters, the transmitter transporters and the redox
    enzymes, each shown against that symptom's base rate and marked with whether the pair is a labelled positive;
  - the rank of those genes among the held-out perturbations of their fold, per symptom, so a gene whose probability is
    high only because the symptom is common is visible as such.

Positives and base rates are the training labels themselves (load_experiment_data: relation "induces", grades A and B,
at least one positive report), not a re-filter of the evidence table. Degree strata use the metabolic graph's degrees
for both runs, so the variant's added edges do not move a perturbation between strata. A run with fewer than five folds
is reported with the folds it has. Idempotent; rewrites its output.

Usage:
  OMP_NUM_THREADS=1 python experiments/report_neuronal_variant_comparison.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from experiments.compare_twin_runs import cached_compare, describe_interval, fold_directories, format_interval, read_run
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import degree_strata

# the genes the variant connects through the electrical, vesicular and redox layers, with what each one tests
VARIANT_GENES = {
    "SCN1A": "voltage-gated sodium channel (curated split, depolarisation-activated)",
    "HCN1": "hyperpolarisation-activated channel (curated split, the opposite gating sign)",
    "ATP1A3": "sodium-potassium ATPase (hyperpolarising current, Human-GEM)",
    "SLC12A5": "KCC2, which sets the chloride pool and so the GABA-A response",
    "SLC12A2": "NKCC1, the opposing cotransporter",
    "CACNA1A": "voltage-gated calcium channel, the route from depolarisation to release",
    "SLC6A3": "dopamine transporter, uptake and the curated reverse transport",
    "SLC18A2": "VMAT2, now only at the vesicle",
    "SLC32A1": "VGAT, vesicular GABA and glycine loading",
    "SOD2": "mitochondrial superoxide dismutase (also the manganese cofactor route)",
    "GSR": "glutathione reductase, the NADPH demand of the glutathione pool",
    "G6PD": "pentose phosphate NADPH supply",
    "GPX4": "glutathione-dependent reduction of lipid hydroperoxides",
    "NFE2L2": "the NRF2 arm of the oxidative response",
    "MAOA": "monoamine oxidase, a hydrogen peroxide source",
}
COMPARISONS = (
    ("b3_linear_response_cofactors_laboratory_neuronal", "b3_linear_response_cofactors_laboratory", "the neuronal variant against the metabolic graph"),
    ("b3_linear_response_cofactors_laboratory_neuronal", "b3_linear_response_cofactors_laboratory_neuronal_no_curation", "the curated receptor and oxidant layers"),
    # b6_linear_response_gate_time_scales_cofactors is the same configuration on the metabolic graph, identical
    # argument for argument apart from --graph-dir, so it is the reference the noisy-OR variant is read against
    ("b6_linear_response_cofactors_neuronal", "b6_linear_response_gate_time_scales_cofactors", "the noisy-OR head on the variant against the metabolic graph"),
)
MINIMUM_POSITIVES_TO_SCORE = 5


def accuracy_lines(run: str, reference: str, description: str, data, strata: np.ndarray, num_bootstrap: int) -> list[str]:
    entry = cached_compare(run, reference, data, strata, num_bootstrap, require_all_folds=False)
    if "set_aside" in entry:
        return [f"- {description}: not yet comparable ({entry['set_aside']})."]
    num_shared_folds = len(entry["per_fold_differences"])
    per_fold, pooled_auprc, pooled_auroc, within = (entry["per_fold_difference"], entry["pooled_macro_auprc"], entry["pooled_macro_auroc"],
                                                    entry["within_degree_strata_macro_auprc"])
    return [(f"- {description}: {run} {entry['per_fold_a']:.3f} ± {entry['per_fold_a_sd']:.3f}, {reference} {entry['per_fold_b']:.3f} ± "
             f"{entry['per_fold_b_sd']:.3f}, over the {num_shared_folds} folds both finished."),
            (f"  Paired per fold (t interval, {num_shared_folds - 1} degrees of freedom): {format_interval(per_fold)}, the interval "
             f"{describe_interval(per_fold['lower'], per_fold['upper'])}. Per fold: {', '.join(f'{difference:+.3f}' for difference in entry['per_fold_differences'])}."),
            (f"  Pooled out-of-fold over the {entry['pooled_rows']} perturbations both scored (paired bootstrap, {pooled_auprc['num_resamples']} resamples): "
             f"macro AUPRC {format_interval(pooled_auprc)}, the interval {describe_interval(pooled_auprc['lower'], pooled_auprc['upper'])}; "
             f"macro AUROC {format_interval(pooled_auroc)}, the interval {describe_interval(pooled_auroc['lower'], pooled_auroc['upper'])}."),
            (f"  Within degree strata (scores ranked inside one stratum of one test fold): macro AUPRC {format_interval(within)}, the interval "
             f"{describe_interval(within['lower'], within['upper'])}.")]


def gene_prediction_table(run: str) -> pd.DataFrame:
    """Predicted probability and within-fold rank of the variant's genes, per symptom."""
    _, predictions = read_run(run)
    if predictions is None:
        return pd.DataFrame()
    symptoms = [column for column in predictions.columns if column not in ("perturbation_id", "fold")]
    ranks = predictions.copy()
    for symptom in symptoms:  # rank within the fold, 1 is the highest predicted probability
        ranks[symptom] = predictions.groupby("fold")[symptom].rank(ascending=False, method="average")
    rows = []
    for gene, role in VARIANT_GENES.items():
        match = predictions.perturbation_id == gene
        if not match.any():
            continue
        fold = int(predictions.loc[match, "fold"].iloc[0])
        held_out = int(match.sum())
        for symptom in symptoms:
            rows.append({"gene": gene, "role": role, "fold": fold, "symptom": symptom,
                         "probability": float(predictions.loc[match, symptom].iloc[0]),
                         "rank_in_fold": float(ranks.loc[match, symptom].iloc[0]),
                         "perturbations_in_fold": int((predictions.fold == fold).sum()), "rows": held_out})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"), help="the metabolic graph, whose degrees define the strata")
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/neuronal_variant_comparison.md"))
    parser.add_argument("--num-bootstrap", type=int, default=1000, help="the same count as compare_twin_runs.py, so its cached pairs are reused")
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster")
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, 5, 0)
    strata = np.array([fold_by_perturbation[p] for p in data.perturbation_ids]) * 1000 + degree_strata(data.perturbation_degrees)
    position_of = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    symptom_position = {symptom: index for index, symptom in enumerate(data.symptoms)}
    base_rates = dict(zip(data.symptoms, data.outcomes.mean(axis=0)))

    lines = ["# The neuronal and oxidative variant against the metabolic graph (generated by experiments/report_neuronal_variant_comparison.py)", ""]
    lines.append("Per-fold accuracy is each fold's macro AUPRC; the per-fold difference is paired by fold and given a t interval. The pooled "
                 "comparison scores the out-of-fold predictions of both runs on the perturbations both scored, with a paired bootstrap over "
                 "perturbations; the within-strata comparison ranks each score inside one degree stratum of one test fold first "
                 f"(metabolic-graph degrees, {len(np.unique(degree_strata(data.perturbation_degrees)))} strata). Symptoms with fewer than "
                 f"{MINIMUM_POSITIVES_TO_SCORE} positives among the rows are not scored.")
    lines.append("")
    lines.append("## Accuracy")
    lines.append("")
    for run, reference, description in COMPARISONS:
        lines += accuracy_lines(run, reference, description, data, strata, arguments.num_bootstrap)
    lines.append("")
    lines.append("## What the variant's own genes predict")
    lines.append("")
    lines.append("One row per gene, symptom and held-out fold: the predicted probability, its rank among that fold's held-out perturbations "
                 "(1 is the highest), the symptom's base rate (the fraction of the slice's perturbations labelled with it), and whether the "
                 "pair is a labelled positive. A high probability on a common symptom is not evidence; the rank is.")
    lines.append("")
    for run in dict.fromkeys(run for run, _, _ in COMPARISONS):
        table = gene_prediction_table(run)
        if table.empty:
            lines.append(f"### {run}: no out-of-fold predictions yet")
            lines.append("")
            continue
        table["base_rate"] = table.symptom.map(base_rates).astype(float)
        table["positive"] = [int(data.outcomes[position_of[gene], symptom_position[symptom]]) if gene in position_of else 0
                             for gene, symptom in zip(table.gene, table.symptom)]
        table = table.sort_values(["gene", "probability"], ascending=[True, False])
        positives_ranked = table[table.positive == 1]
        lines.append(f"### {run} ({table.gene.nunique()} of the variant's genes held out; {len(positives_ranked)} labelled positive pairs, "
                     f"median rank among them {positives_ranked.rank_in_fold.median():.0f})" if len(positives_ranked) else
                     f"### {run} ({table.gene.nunique()} of the variant's genes held out; no labelled positive pairs)")
        lines.append("")
        lines.append("| gene | role | fold | symptom | probability | rank in fold | perturbations in fold | base rate | positive |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for row in table.itertuples():
            lines.append(f"| {row.gene} | {row.role} | {row.fold} | {row.symptom} | {row.probability:.3f} | {row.rank_in_fold:.0f} | "
                         f"{row.perturbations_in_fold} | {row.base_rate:.3f} | {row.positive} |")
        lines.append("")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(f"wrote {arguments.markdown_output}")
    for run, reference, description in COMPARISONS:
        print(f"  {description}: {len(fold_directories(run))} and {len(fold_directories(reference))} folds")


if __name__ == "__main__":
    main()
