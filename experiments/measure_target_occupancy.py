"""Occupancy of each study drug's mechanism targets at its therapeutic unbound concentration.

Why. A drug enters the model with magnitude 1 at each of its mechanism targets, so a drug with three targets carries
an input of 3 and the size of its input follows the number of targets curated for it. Pharmacology has a quantity for
"how strongly is this target engaged": the occupancy theta = C / (C + K), with C the unbound concentration at the
target and K the dissociation constant. This script computes it wherever both numbers exist, to show what an
occupancy-scaled input would look like beside the unit input, and to test the placeholder docs/off_target_scoping.md
uses ("the best-bound mechanism target is 90 percent occupied").

Inputs, none of them changed:
  data/processed/exposure/drug_exposure_<name>.csv     unbound plasma and brain concentrations (build_drug_exposure_table.py)
  DrugCentral and Guide to PHARMACOLOGY affinities     one value per drug and gene, the maximum (the user's rule of 9 October)

What the numbers leave out, so that they are read as a first measurement and not as occupancies to train on:
  - active metabolites are not added to the parent's concentration (nortriptyline to amitriptyline, norfluoxetine,
    paliperidone to risperidone), which understates engagement for those drugs;
  - an EC50 or IC50 is pooled with Ki and Kd, as in the scoping;
  - the concentration is at steady state within the therapeutic range, not a peak;
  - the brain concentration uses a Kp,uu that is from animals for most drugs.

Nothing here changes a label, a graph or a model input.

Usage:
  PYTHONPATH=. python experiments/measure_target_occupancy.py
"""
from __future__ import annotations

import argparse
import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data

SCOPING_SCRIPT = Path(__file__).resolve().parent / "scope_off_target_binding.py"
specification = importlib.util.spec_from_file_location("scope_off_target_binding", SCOPING_SCRIPT)
scoping = importlib.util.module_from_spec(specification)
specification.loader.exec_module(scoping)

EXPOSURE_DIRECTORY = Path("data/processed/exposure")
OUTPUT_DOCUMENT = Path("docs/target_occupancy_measured.md")
OCCUPANCY_BANDS = ((0.0, 0.1, "below 0.1"), (0.1, 0.5, "0.1 to 0.5"), (0.5, 0.9, "0.5 to 0.9"), (0.9, 1.01, "0.9 and above"))


def occupancy_at(concentration_micromolar: float, p_affinity: float) -> float:
    """One-site occupancy at a concentration in micromolar, for an affinity given as -log10 of a molar constant."""
    dissociation_constant_micromolar = 10.0 ** (6.0 - p_affinity)
    return concentration_micromolar / (concentration_micromolar + dissociation_constant_micromolar)


def affinities_by_drug(data, graph_directory: Path) -> dict[str, dict]:
    """perturbation id -> {"mechanism": {gene: pK}, "other": {gene: pK}} from both sources, the maximum per gene."""
    nodes = pd.read_parquet(graph_directory / "nodes.parquet", columns=["node_id", "gene_symbol"])
    symbol_of_node = dict(zip(nodes.node_id, nodes.gene_symbol))
    drugcentral = scoping.drugcentral_measurements(scoping.DRUGCENTRAL_INTERACTIONS)
    gtopdb = scoping.gtopdb_measurements(scoping.GTOPDB_INTERACTIONS, scoping.GTOPDB_LIGAND_MAPPING)
    drugcentral_id_of_name = drugcentral.groupby("drug_name").drugcentral_id.first().to_dict()
    result = {}
    for perturbation_id, kind, label, seeds in zip(data.perturbation_ids, data.perturbation_types, data.perturbation_labels, data.perturbation_seeds):
        if kind != "drug":
            continue
        pubchem_cid, name = scoping.drug_keys(perturbation_id, label)
        from_gtopdb = gtopdb[((gtopdb.pubchem_cid == pubchem_cid) if pubchem_cid else False) | (gtopdb.drug_name == name) | (gtopdb.inn == name)]
        drugcentral_ids = set(from_gtopdb.drugcentral_id.dropna().astype(int)) | ({drugcentral_id_of_name[name]} if name in drugcentral_id_of_name else set())
        measurements = pd.concat([drugcentral[drugcentral.drugcentral_id.isin(drugcentral_ids)], from_gtopdb], ignore_index=True)
        mechanism_genes = {symbol_of_node.get(data.node_ids[seed]) for seed in seeds} - {None}
        best = measurements.dropna(subset=["p_affinity"]).groupby("gene").p_affinity.max().to_dict() if len(measurements) else {}
        result[perturbation_id] = {"label": label, "mechanism_genes": sorted(mechanism_genes),
                                   "mechanism": {gene: value for gene, value in best.items() if gene in mechanism_genes},
                                   "other": {gene: value for gene, value in best.items() if gene not in mechanism_genes}}
    return result


def band_counts(values: np.ndarray) -> list[int]:
    return [int(((values >= low) & (values < high)).sum()) for low, high, _ in OCCUPANCY_BANDS]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", default="v4")
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"), help="a graph whose gene nodes carry the mechanism targets")
    parser.add_argument("--document", type=Path, default=OUTPUT_DOCUMENT)
    arguments = parser.parse_args()

    exposure = pd.read_csv(EXPOSURE_DIRECTORY / f"drug_exposure_{arguments.name}.csv").set_index("perturbation_id")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster", label_selection=arguments.label_selection)
    affinities = affinities_by_drug(data, arguments.graph_dir)

    rows = []
    for perturbation_id, entry in affinities.items():
        if perturbation_id not in exposure.index:
            continue
        drug = exposure.loc[perturbation_id]
        for site in ("plasma", "brain"):
            low, high = drug.get(f"unbound_{site}_low_micromolar"), drug.get(f"unbound_{site}_high_micromolar")
            if pd.isna(low) or pd.isna(high):
                continue
            for role in ("mechanism", "other"):
                for gene, p_affinity in entry[role].items():
                    rows.append({"perturbation_id": perturbation_id, "drug": entry["label"], "site": site, "role": role, "gene": gene, "p_affinity": p_affinity,
                                 "mechanism_targets_of_drug": len(entry["mechanism_genes"]),
                                 "occupancy_low": occupancy_at(low, p_affinity), "occupancy_high": occupancy_at(high, p_affinity)})
    occupancies = pd.DataFrame(rows)
    output_path = EXPOSURE_DIRECTORY / f"target_occupancy_{arguments.name}.csv"
    occupancies.to_csv(output_path, index=False)

    drugs_total = len(affinities)
    lines = ["# Occupancy of mechanism targets at therapeutic unbound concentrations", "",
             "Generated by `experiments/measure_target_occupancy.py`. Do not edit by hand.", "",
             f"Drugs: the {drugs_total} drugs of `{arguments.evidence_dir}`. Occupancy is C / (C + K) with C the unbound concentration within the therapeutic range",
             f"(`docs/drug_exposure_table.md`) and K from the larger of the DrugCentral and Guide to PHARMACOLOGY affinities of the drug at the gene. The table of pairs is `{output_path}`, not in git.",
             "Nothing here changes a label, a graph or a model input.", ""]
    for site, site_text in (("plasma", "unbound plasma concentration"), ("brain", "unbound brain concentration (plasma × Kp,uu,brain)")):
        at_site = occupancies[occupancies.site == site] if len(occupancies) else occupancies
        mechanism = at_site[at_site.role == "mechanism"] if len(at_site) else at_site
        if not len(mechanism):
            continue
        drugs_with_pairs = mechanism.perturbation_id.nunique()
        mechanism_targets_of_those = int(sum(len(affinities[perturbation_id]["mechanism_genes"]) for perturbation_id in mechanism.perturbation_id.unique()))
        lines += [f"## At the {site_text}", "",
                  f"{drugs_with_pairs} drugs have the concentration and an affinity at one or more of their mechanism genes: {len(mechanism)} drug-gene pairs, of the {mechanism_targets_of_those} mechanism genes those drugs have.", "",
                  "| occupancy | pairs at the low end of the therapeutic range | pairs at the high end |", "| --- | --- | --- |"]
        low_counts, high_counts = band_counts(mechanism.occupancy_low.to_numpy()), band_counts(mechanism.occupancy_high.to_numpy())
        lines += [f"| {label} | {low_count} | {high_count} |" for (_, _, label), low_count, high_count in zip(OCCUPANCY_BANDS, low_counts, high_counts)]
        best_per_drug = mechanism.groupby("perturbation_id").occupancy_high.max()
        summed = mechanism.groupby("perturbation_id").occupancy_high.sum()
        counted = mechanism.groupby("perturbation_id").gene.count()
        several = counted[counted >= 2].index
        spread = (mechanism[mechanism.perturbation_id.isin(several)].groupby("perturbation_id").occupancy_high.agg(lambda values: values.max() - values.min()))
        lines += ["",
                  f"- **The best-bound mechanism gene of a drug** has a median occupancy of {best_per_drug.median():.2f} at the high end of the range (quartiles {best_per_drug.quantile(0.25):.2f} and {best_per_drug.quantile(0.75):.2f}); "
                  f"it is at 0.9 or above for {int((best_per_drug >= 0.9).sum())} of the {len(best_per_drug)} drugs and below 0.5 for {int((best_per_drug < 0.5).sum())}. "
                  "The scoping's placeholder assumes 0.9 for every drug.",
                  f"- **Within a drug**, for the {len(several)} drugs with an affinity at two or more mechanism genes, the occupancies of its genes differ by a median of {spread.median():.2f} "
                  f"(0.5 or more for {int((spread >= 0.5).sum())} drugs). A unit input per target treats them as equal.",
                  f"- **Total input.** Summing the occupancies of a drug's mechanism genes gives a median total of {summed.median():.2f} per drug, against a median of {counted.median():.0f} such genes counted; "
                  f"the correlation between the two across drugs is {np.corrcoef(summed.to_numpy(), counted.to_numpy())[0, 1]:.2f}.", ""]
        other = at_site[at_site.role == "other"]
        if len(other):
            lines += [f"- **Genes outside the mechanism list**: {len(other)} drug-gene pairs have an affinity; {int((other.occupancy_high >= 0.25).sum())} of them, over {other[other.occupancy_high >= 0.25].perturbation_id.nunique()} drugs, "
                      f"reach an occupancy of 0.25 at the high end of the range and {int((other.occupancy_high >= 0.5).sum())} reach 0.5. These are the off-targets of `docs/off_target_scoping.md`, here with a measured concentration in place of its placeholder.", ""]
    lines += ["## What the numbers leave out", "",
              "- Active metabolites are not added to the parent's concentration, which understates engagement for drugs such as amitriptyline, fluoxetine and risperidone.",
              "- EC50 and IC50 values are pooled with Ki and Kd. An EC50 measures potency in an assay and not binding.",
              "- The concentration is a steady-state value inside the therapeutic range; it is not a peak after a dose.",
              "- The brain concentration rests on a Kp,uu that comes from animals for most drugs.",
              "- A drug that is not a reversible ligand (an irreversible enzyme inhibitor, a prodrug, a substrate that releases transmitter) is not described by one-site occupancy.", ""]
    arguments.document.write_text("\n".join(lines))
    print(f"wrote {output_path} ({len(occupancies)} rows) and {arguments.document}")
    print("\n".join(lines[8:40]))


if __name__ == "__main__":
    main()
