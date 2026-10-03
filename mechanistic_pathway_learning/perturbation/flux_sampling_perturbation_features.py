"""Evidence-free mechanistic features from flux sampling (design section 5.2, route 2).

For each gene in the reconstruction, knock it out (its reactions become blocked
when the gene-protein-reaction rule evaluates false), sample the feasible flux
space and record, per reaction, the shift in sampled flux median and the change
in flux variability range relative to the unperturbed model. Sampling rather
than flux balance analysis is used because the brain has no defensible single
objective; flux variability is computed with fraction_of_optimum=0 for the same
reason.

Written for Slurm array jobs on a runtime-limited cluster: one gene per task,
idempotent (skips when the output file exists), with the unperturbed baseline
computed once by task 0 and read by the rest.

Usage:
  python -m mechanistic_pathway_learning.perturbation.flux_sampling_perturbation_features \
      --model-path data/raw/Human-GEM.xml --output-dir data/cache/flux_perturbation \
      --gene-index $SLURM_ARRAY_TASK_ID --num-samples 2000

Not yet executed against Human-GEM; the cobra API calls below follow cobrapy 0.29.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def sample_fluxes(model, num_samples: int, processes: int, random_seed: int):
    from cobra.sampling import sample  # imported lazily so the package imports without cobra

    return sample(model, num_samples, method="optgp", processes=processes, seed=random_seed)


def flux_variability(model):
    from cobra.flux_analysis import flux_variability_analysis

    return flux_variability_analysis(model, fraction_of_optimum=0.0)


def summarize_shift(baseline_samples, perturbed_samples, baseline_fva, perturbed_fva):
    import pandas as pd

    median_shift = perturbed_samples.median(axis=0) - baseline_samples.median(axis=0)
    baseline_range = baseline_fva["maximum"] - baseline_fva["minimum"]
    perturbed_range = perturbed_fva["maximum"] - perturbed_fva["minimum"]
    return pd.DataFrame(
        {
            "reaction_id": median_shift.index,
            "flux_median_shift": median_shift.values,
            "flux_range_change": (perturbed_range - baseline_range).reindex(median_shift.index).values,
        }
    )


def run_for_gene_index(model_path: Path, output_directory: Path, gene_index: int, num_samples: int, processes: int, random_seed: int) -> Path:
    from cobra.io import load_model, read_sbml_model

    output_directory.mkdir(parents=True, exist_ok=True)
    model = read_sbml_model(str(model_path)) if model_path.suffix.lower() == ".xml" else load_model(str(model_path))
    gene_list_path = output_directory / "gene_order.json"
    if not gene_list_path.exists():
        gene_list_path.write_text(json.dumps([gene.id for gene in model.genes]))
    gene_order = json.loads(gene_list_path.read_text())
    baseline_samples_path = output_directory / "baseline_samples.parquet"
    baseline_fva_path = output_directory / "baseline_fva.parquet"
    if not baseline_samples_path.exists():
        sample_fluxes(model, num_samples, processes, random_seed).to_parquet(baseline_samples_path)
        flux_variability(model).to_parquet(baseline_fva_path)
    import pandas as pd

    gene_id = gene_order[gene_index]
    output_path = output_directory / f"{gene_id}.parquet"
    if output_path.exists():
        return output_path
    baseline_samples = pd.read_parquet(baseline_samples_path)
    baseline_fva = pd.read_parquet(baseline_fva_path)
    with model:
        model.genes.get_by_id(gene_id).knock_out()
        perturbed_samples = sample_fluxes(model, num_samples, processes, random_seed)
        perturbed_fva = flux_variability(model)
    summarize_shift(baseline_samples, perturbed_samples, baseline_fva, perturbed_fva).to_parquet(output_path)
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gene-index", type=int, required=True)
    parser.add_argument("--num-samples", type=int, default=2000)
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--random-seed", type=int, default=0)
    arguments = parser.parse_args()
    output_path = run_for_gene_index(arguments.model_path, arguments.output_dir, arguments.gene_index, arguments.num_samples, arguments.processes, arguments.random_seed)
    print(output_path)


if __name__ == "__main__":
    main()
