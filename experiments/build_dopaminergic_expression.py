"""Pseudobulk expression of the dopaminergic neuron cluster of the adult human brain atlas (Siletti et al. 2023), on the
Human Protein Atlas nCPM scale, for the class_dopaminergic_neuron column of the brain expression descriptors.

HPA's single-nucleus brain tables group the atlas's clusters into 34 cluster types, none of them dopaminergic, so the
cluster is summed here from the atlas's per-dissection count files on CZ CELLxGENE (collection
283d65eb-dd53-496d-adb7-7570c7caa443, CC BY 4.0; X holds raw UMI counts). Cluster 395 ("DA VGLUT2" in the atlas's
tables/cluster_annotation.xlsx, 998 nuclei from 3 donors, 93 percent of them in midbrain dissections) is the
dopaminergic cluster. Four dissections hold most of its nuclei: SN-RN, SN, PAG-DR and PAG.

Steps: sum the UMI counts of the cluster's nuclei per gene; keep the genes of HPA's cluster-type table (the same
Ensembl ids); counts per million over those genes; divide by trimmed_mean_scale_factor against HPA's splatter cluster
type (the supercluster cluster 395 belongs to), so the column sits on the scale of the other cell classes. Only the
cluster's rows are read from each file (h5py, row slices of the compressed sparse matrix), so memory stays small.

Usage:
  python experiments/build_dopaminergic_expression.py \
      --h5ad data/raw/siletti_cellxgene/dissection_sn_rn.h5ad data/raw/siletti_cellxgene/dissection_sn.h5ad \
             data/raw/siletti_cellxgene/dissection_pag_dr.h5ad data/raw/siletti_cellxgene/dissection_pag.h5ad \
      --output data/processed/brain_expression/dopaminergic_siletti_cluster395.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.brain_expression_descriptors import counts_per_million, read_zipped_table, trimmed_mean_scale_factor

DOPAMINERGIC_CLUSTER_ID = "395"
REFERENCE_CLUSTER_TYPE = "splatter"
ANCHOR_GENES = ["TH", "SLC6A3", "SLC18A2", "DDC", "KCNJ6", "NR4A2", "LMX1B", "EN1", "ALDH1A1", "SOX6", "CALB1", "SLC17A6",
                "GAD1", "GAD2", "AQP4", "MBP", "GCH1", "MAOA", "COMT", "DRD2"]


def categorical_column(group: h5py.Group, name: str) -> np.ndarray:
    """Values of an AnnData categorical (codes plus categories) or plain string column."""
    column = group[name]
    if isinstance(column, h5py.Group):
        categories = column["categories"][:].astype(str)
        return categories[column["codes"][:]]
    return column[:].astype(str)


def cluster_counts(h5ad_path: Path, cluster_id: str) -> tuple[pd.Series, dict]:
    """(UMI counts summed over the cluster's nuclei, indexed by Ensembl id; a description of what was read)."""
    with h5py.File(h5ad_path, "r") as handle:
        clusters = categorical_column(handle["obs"], "cluster_id")
        rows = np.flatnonzero(clusters == cluster_id)
        donors = categorical_column(handle["obs"], "donor_id")[rows]
        primary = categorical_column(handle["obs"], "is_primary_data")[rows] if "is_primary_data" in handle["obs"] else None
        matrix = handle["X"]
        if matrix.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{h5ad_path}: X is {matrix.attrs.get('encoding-type')}, expected csr_matrix")
        num_genes = int(matrix.attrs["shape"][1])
        indptr = matrix["indptr"][:]
        totals = np.zeros(num_genes)
        non_integer_values = 0
        for row in rows:
            start, end = int(indptr[row]), int(indptr[row + 1])
            values = matrix["data"][start:end].astype(float)
            non_integer_values += int((values != np.round(values)).sum())
            np.add.at(totals, matrix["indices"][start:end], values)
        if non_integer_values:
            raise ValueError(f"{h5ad_path}: {non_integer_values} non-integer values in X; expected raw UMI counts")
        gene_ids = handle["var"][handle["var"].attrs["_index"]][:].astype(str)
        symbols = categorical_column(handle["var"], "feature_name")
    counts = pd.Series(totals, index=pd.Index(gene_ids, name="ensembl_gene_id"))
    description = {"file": str(h5ad_path), "nuclei": int(len(rows)), "umis": float(totals.sum()),
                   "nuclei_by_donor": pd.Series(donors).value_counts().to_dict(),
                   "is_primary_data": None if primary is None else pd.Series(primary).value_counts().to_dict()}
    return counts, description | {"symbols": dict(zip(gene_ids, symbols))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5ad", type=Path, nargs="+", required=True)
    parser.add_argument("--cluster-id", default=DOPAMINERGIC_CLUSTER_ID)
    parser.add_argument("--hpa-cluster-type", type=Path, default=Path("data/raw/hpa/rna_single_nuclei_cluster_type/rna_single_nuclei_cluster_type.tsv.zip"))
    parser.add_argument("--reference-cluster-type", default=REFERENCE_CLUSTER_TYPE)
    parser.add_argument("--output", type=Path, default=Path("data/processed/brain_expression/dopaminergic_siletti_cluster395.parquet"))
    arguments = parser.parse_args()

    total_counts, descriptions, symbols = None, [], {}
    for path in arguments.h5ad:
        counts, description = cluster_counts(path, arguments.cluster_id)
        symbols.update(description.pop("symbols"))
        descriptions.append(description)
        total_counts = counts if total_counts is None else total_counts.add(counts, fill_value=0.0)
        print(json.dumps(description), flush=True)

    hpa = read_zipped_table(arguments.hpa_cluster_type)
    reference = hpa[hpa["Cluster type"] == arguments.reference_cluster_type].set_index("Gene")["nCPM"].astype(float)
    hpa_genes = pd.Index(hpa["Gene"].unique())
    counts_on_hpa_genes = total_counts.reindex(hpa_genes).fillna(0.0)
    cpm = counts_per_million(counts_on_hpa_genes)
    scale_factor = trimmed_mean_scale_factor(cpm, reference)
    table = pd.DataFrame({"gene_symbol": [symbols.get(gene, "") for gene in hpa_genes], "umi_count": counts_on_hpa_genes.to_numpy(),
                          "cpm": cpm.to_numpy(), "ncpm_aligned": (cpm / scale_factor).to_numpy()},
                         index=pd.Index(hpa_genes, name="ensembl_gene_id"))
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(arguments.output)

    hpa_wide = hpa.pivot_table(index="Gene", columns="Cluster type", values="nCPM", aggfunc="max")
    by_symbol = table.reset_index().set_index("gene_symbol")
    anchors = {}
    for symbol in ANCHOR_GENES:
        if symbol not in by_symbol.index:
            anchors[symbol] = None
            continue
        row = by_symbol.loc[[symbol]].iloc[0]
        hpa_row = hpa_wide.loc[row.ensembl_gene_id] if row.ensembl_gene_id in hpa_wide.index else None
        anchors[symbol] = {"dopaminergic_ncpm": round(float(row.ncpm_aligned), 1),
                           "hpa_max_ncpm": None if hpa_row is None else round(float(hpa_row.max()), 1),
                           "hpa_top_cluster_type": None if hpa_row is None else str(hpa_row.idxmax())}
    ranks = table.ncpm_aligned.rank(ascending=False)
    summary = {
        "output": str(arguments.output), "cluster_id": arguments.cluster_id, "files": descriptions,
        "nuclei": int(sum(description["nuclei"] for description in descriptions)),
        "umis_on_hpa_genes": float(counts_on_hpa_genes.sum()), "umis_all_genes": float(total_counts.sum()),
        "hpa_genes": int(len(hpa_genes)), "hpa_genes_with_a_count": int((counts_on_hpa_genes > 0).sum()),
        "reference_cluster_type": arguments.reference_cluster_type, "scale_factor_cpm_per_ncpm": scale_factor,
        "anchor_genes": anchors,
        "anchor_ranks": {symbol: int(ranks[by_symbol.loc[[symbol]].ensembl_gene_id.iloc[0]]) for symbol in ["TH", "SLC6A3", "SLC18A2"] if symbol in by_symbol.index},
    }
    arguments.output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "files"}, indent=2))


if __name__ == "__main__":
    main()
