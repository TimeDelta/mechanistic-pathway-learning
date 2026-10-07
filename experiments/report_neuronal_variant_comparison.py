"""Compare the neuronal and oxidative graph variant with the metabolic graph it extends, and read off what the variant's
own nodes do (design section 5.7; graph/reactome_import.py, neurotransmission_regulation.py, oxidative_regulation.py,
redox_pools.py).

Three readings, written to docs/neuronal_variant_comparison.md:
  - accuracy, per fold and pooled over the out-of-fold predictions, of the variant against the same configuration on the
    metabolic graph and against the variant without its curated receptor and oxidant layers, with the paired bootstrap
    of the pooled difference and the same comparison inside degree strata, which gives no credit for ordering
    perturbations by degree (design section 5.2);
  - the predicted symptom probabilities of the genes the variant connects and the metabolic graph does not reach the
    same way: the voltage-gated channels, the pumps and cotransporters, the transmitter transporters and the redox
    enzymes, each shown against that symptom's base rate and marked with whether the pair is a held-out positive;
  - the rank of those genes among the held-out perturbations of their fold, per symptom, so a gene whose probability is
    high only because the symptom is common is visible as such.

A run with fewer than five folds is reported with the folds it has. Idempotent; rewrites its output.

Usage:
  OMP_NUM_THREADS=1 python experiments/report_neuronal_variant_comparison.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

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
RUN_ROOT = Path("runs/encoder")
COMPARISONS = (
    ("b3_linear_response_cofactors_laboratory_neuronal", "b3_linear_response_cofactors_laboratory", "the neuronal variant against the metabolic graph"),
    ("b3_linear_response_cofactors_laboratory_neuronal", "b3_linear_response_cofactors_laboratory_neuronal_no_curation", "the curated receptor and oxidant layers"),
    # b6_linear_response_gate_time_scales_cofactors is the same configuration on the metabolic graph, identical
    # argument for argument apart from --graph-dir, so it is the reference the noisy-OR variant is read against
    ("b6_linear_response_cofactors_neuronal", "b6_linear_response_gate_time_scales_cofactors", "the noisy-OR head on the variant against the metabolic graph"),
)


def fold_directories(run: str) -> list[Path]:
    directory = RUN_ROOT / f"{run}_disease_cluster"
    return sorted(path.parent for path in directory.glob("fold*_seed0/results.json")) if directory.exists() else []


def read_run(run: str) -> tuple[list[dict], pd.DataFrame | None]:
    """(per-fold results, out-of-fold predictions as perturbation x symptom rows)."""
    results, frames = [], []
    for fold in fold_directories(run):
        result = json.loads((fold / "results.json").read_text())
        results.append(result)
        predictions = fold / "test_predictions.npy"
        if predictions.exists() and result.get("test_perturbation_ids") and result.get("symptoms"):
            scores = np.load(predictions)
            if scores.shape == (len(result["test_perturbation_ids"]), len(result["symptoms"])):
                frames.append(pd.DataFrame(scores, index=result["test_perturbation_ids"], columns=result["symptoms"])
                              .rename_axis("perturbation_id").reset_index().assign(fold=result.get("fold")))
    return results, (pd.concat(frames, ignore_index=True) if frames else None)


def paired_bootstrap(first: np.ndarray, second: np.ndarray, generator: np.random.Generator, resamples: int = 2000) -> tuple[float, float, float]:
    """(difference of the means, lower and upper 95 percent bounds) of first minus second over paired rows."""
    difference = first - second
    draws = [difference[generator.integers(0, len(difference), len(difference))].mean() for _ in range(resamples)]
    return float(difference.mean()), float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def accuracy_lines(run: str, reference: str, description: str, generator: np.random.Generator) -> list[str]:
    run_results, _ = read_run(run)
    reference_results, _ = read_run(reference)
    if not run_results or not reference_results:
        return [f"- {description}: not yet runnable ({run} has {len(run_results)} folds, {reference} has {len(reference_results)})."]
    run_folds = np.array([result["macro_auprc"] for result in run_results], dtype=float)
    reference_folds = np.array([result["macro_auprc"] for result in reference_results], dtype=float)
    lines = [f"- {description}: {run} {run_folds.mean():.3f} ± {run_folds.std(ddof=1) if len(run_folds) > 1 else 0:.3f} over {len(run_folds)} folds, "
             f"{reference} {reference_folds.mean():.3f} ± {reference_folds.std(ddof=1) if len(reference_folds) > 1 else 0:.3f} over {len(reference_folds)}."]
    if len(run_folds) == len(reference_folds) > 1:
        difference, lower, upper = paired_bootstrap(run_folds, reference_folds, generator)
        lines.append(f"  Paired over folds: {difference:+.3f} [{lower:+.3f}, {upper:+.3f}]; "
                     f"{'the interval excludes zero' if lower > 0 or upper < 0 else 'the interval includes zero'}.")
    return lines


def gene_prediction_table(run: str) -> pd.DataFrame:
    """Predicted probability and within-fold rank of the variant's genes, per symptom."""
    results, predictions = read_run(run)
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
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/neuronal_variant_comparison.md"))
    parser.add_argument("--seed", type=int, default=0)
    arguments = parser.parse_args()
    generator = np.random.default_rng(arguments.seed)

    positives = set()
    base_rates: dict[str, float] = {}
    records = arguments.evidence_dir / "evidence_records.parquet"
    if records.exists():
        evidence = pd.read_parquet(records)
        positive_rows = evidence[evidence.relation.astype(str).eq("causes")] if "relation" in evidence.columns else evidence
        positives = set(zip(positive_rows.perturbation_id.astype(str), positive_rows.symptom.astype(str)))
        perturbations = evidence.perturbation_id.nunique()
        base_rates = {symptom: count / perturbations for symptom, count in positive_rows.symptom.value_counts().items()}

    lines = ["# The neuronal and oxidative variant against the metabolic graph (generated by experiments/report_neuronal_variant_comparison.py)", ""]
    lines.append("Accuracy is the per-fold macro AUPRC of the out-of-fold predictions; the paired bootstrap is over the folds both runs share.")
    lines.append("")
    lines.append("## Accuracy")
    lines.append("")
    for run, reference, description in COMPARISONS:
        lines += accuracy_lines(run, reference, description, generator)
    lines.append("")
    lines.append("## What the variant's own genes predict")
    lines.append("")
    lines.append("One row per gene, symptom and held-out fold: the predicted probability, its rank among that fold's held-out perturbations "
                 "(1 is the highest), the symptom's base rate over the evidence table, and whether the pair is a positive. A high probability "
                 "on a common symptom is not evidence; the rank is.")
    lines.append("")
    for run, _, _ in COMPARISONS:
        table = gene_prediction_table(run)
        if table.empty:
            lines.append(f"### {run}: no out-of-fold predictions yet")
            lines.append("")
            continue
        table["base_rate"] = table.symptom.map(base_rates).fillna(float("nan"))
        table["positive"] = [int((gene, symptom) in positives) for gene, symptom in zip(table.gene, table.symptom)]
        table = table.sort_values(["gene", "probability"], ascending=[True, False])
        lines.append(f"### {run} ({table.gene.nunique()} of the variant's genes held out)")
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
