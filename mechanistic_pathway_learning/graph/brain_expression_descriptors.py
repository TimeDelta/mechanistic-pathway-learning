"""Brain region and brain cell-class expression as node descriptors, for genes and for the reactions they catalyse.

The slice graphs carry no brain expression (only graph_full and graph_full_neuronal have the GTEx column), and no graph
carries regional or cell-type expression. This module gives both as two descriptor blocks in the convention of
node_descriptors.py (one block per node type, zero outside it, so the encoder's input layer acts as one linear map per
node type and nothing is learned per node). It adds no nodes and no edges; the cell-class layer copies are separate work.

Columns, per gene (log1p of the expression value, standardised over the nodes of the block's type that have a value):
  - gtex_brain_max and gtex_other_max: the largest GTEx v10 median TPM over the 13 brain tissues, and over every other
    tissue, so a peripheral anchor such as PAH (2.7 TPM in brain, 272 elsewhere) is told apart from a gene absent
    everywhere;
  - region_<name>: Human Protein Atlas consensus nTPM in each of 13 brain regions (rna_brain_region_hpa);
  - class_<name>: Human Protein Atlas single-nucleus nCPM (rna_single_nuclei_cluster_type, 34 cluster types from the
    adult human brain atlas of Siletti et al. 2023), the largest value over the cluster types of each of 10 cell
    classes in CELL_CLASS_OF_CLUSTER_TYPE. The 34 types are grouped to keep the column count low for 451 training
    perturbations; the grouping follows the atlas's own superclusters, with striatal medium spiny neurons kept apart
    from other inhibitory neurons because of their place in reward and motor circuits. The atlas has no dopaminergic
    cluster: SLC6A3 peaks at 0.6 nCPM, so midbrain dopamine neurons are not represented here;
  - has_expression: 1 when the gene (or a gene in the reaction's rule) is in the HPA region table.

Reactions take their value from the Human-GEM gene rule, minimum over "and" (every subunit is needed) and maximum over
"or" (any isozyme suffices). Machado and Herrgård 2014 (doi:10.1371/journal.pcbi.1003580), surveying methods that put
expression into metabolic models: "the expression level of reactions catalyzed by enzyme complexes ([AND] operator) is set
to the minimum expression level of the associated genes, and the expression level of reactions catalyzed by isozymes ([OR]
operator) is set to either the maximum or the sum" (operator names bracketed where the retrieved text dropped them). The maximum is used here so one high isozyme is not outweighed by
several silent ones. The same survey found that "for many conditions, the predictions obtained by simple flux balance
analysis using growth maximization and parsimony criteria are as good or better than those obtained using methods that
incorporate transcriptomic data", so expression is not assumed to help; the slice arm measures it. Genes absent from a
table are skipped inside the rule; a rule with no known gene gives a missing value, which becomes 0 (the mean) with
has_expression at 0. For genes expressed at a few nCPM, which cell class ranks first is noise: TH ranks medium spiny
neurons first at low counts.
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.node_descriptors import standardise

GTEX_BRAIN_TISSUE_PREFIX = "Brain_"
CELL_CLASS_OF_CLUSTER_TYPE = {
    "amygdala excitatory": "excitatory_neuron", "deep-layer corticothalamic and 6b": "excitatory_neuron",
    "deep-layer intratelencephalic": "excitatory_neuron", "deep-layer near-projecting": "excitatory_neuron",
    "hippocampal CA1-3": "excitatory_neuron", "hippocampal CA4": "excitatory_neuron",
    "hippocampal dentate gyrus": "excitatory_neuron", "thalamic excitatory": "excitatory_neuron",
    "upper-layer intratelencephalic": "excitatory_neuron", "lower rhombic lip": "excitatory_neuron",
    "upper rhombic lip": "excitatory_neuron", "mammillary body": "excitatory_neuron",
    "CGE interneuron": "cortical_interneuron", "MGE interneuron": "cortical_interneuron",
    "LAMP5-LHX6 and Chandelier": "cortical_interneuron",
    "medium spiny neuron": "medium_spiny_neuron", "eccentric medium spiny neuron": "medium_spiny_neuron",
    "cerebellar inhibitory": "other_inhibitory_neuron", "midbrain-derived inhibitory": "other_inhibitory_neuron",
    "splatter": "other_neuron", "miscellaneous": "other_neuron",
    "astrocyte": "astrocyte", "Bergmann glia": "astrocyte",
    "oligodendrocyte": "oligodendrocyte_lineage", "oligodendrocyte precursor cell": "oligodendrocyte_lineage",
    "committed oligodendrocyte precursor": "oligodendrocyte_lineage",
    "central nervous system macrophage": "immune", "leukocyte": "immune",
    "endothelial cell": "vascular", "pericyte": "vascular", "vascular associated smooth muscle cell": "vascular",
    "fibroblast": "vascular",
    "ependymal cell": "ependymal_choroid", "choroid plexus epithelial cell": "ependymal_choroid",
}
RULE_TOKEN = re.compile(r"\(|\)|\band\b|\bor\b|[^\s()]+")


def read_zipped_table(path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(path) as archive:
        with archive.open(archive.namelist()[0]) as handle:
            return pd.read_csv(handle, sep="\t")


def gtex_brain_and_other_maximum(gct_path: Path) -> pd.DataFrame:
    """Ensembl-indexed largest median TPM over the GTEx brain tissues and over every other tissue."""
    table = pd.read_csv(gct_path, sep="\t", skiprows=2)
    tissue_columns = [column for column in table.columns if column not in ("Name", "Description")]
    brain_columns = [column for column in tissue_columns if column.startswith(GTEX_BRAIN_TISSUE_PREFIX)]
    other_columns = [column for column in tissue_columns if column not in brain_columns]
    if not brain_columns or not other_columns:
        raise ValueError(f"{gct_path} needs both brain ({GTEX_BRAIN_TISSUE_PREFIX}) and other tissue columns")
    result = pd.DataFrame({
        "ensembl_gene_id": table["Name"].astype(str).str.split(".").str[0],
        "gtex_brain_max": table[brain_columns].max(axis=1).astype(float),
        "gtex_other_max": table[other_columns].max(axis=1).astype(float),
    })
    return result.groupby("ensembl_gene_id").max()


def hpa_region_table(region_table: pd.DataFrame) -> pd.DataFrame:
    """Ensembl-indexed nTPM per brain region, columns region_<name with underscores>."""
    wide = region_table.pivot_table(index="Gene", columns="Brain region", values="nTPM", aggfunc="max")
    wide.columns = ["region_" + str(name).replace(" ", "_") for name in wide.columns]
    return wide.rename_axis("ensembl_gene_id")


def hpa_cell_class_table(cluster_table: pd.DataFrame) -> pd.DataFrame:
    """Ensembl-indexed largest nCPM over the cluster types of each cell class, columns class_<class>."""
    unknown = set(cluster_table["Cluster type"].unique()) - set(CELL_CLASS_OF_CLUSTER_TYPE)
    if unknown:
        raise ValueError(f"cluster types with no cell class: {sorted(unknown)}")
    with_class = cluster_table.assign(cell_class=cluster_table["Cluster type"].map(CELL_CLASS_OF_CLUSTER_TYPE))
    wide = with_class.pivot_table(index="Gene", columns="cell_class", values="nCPM", aggfunc="max")
    wide.columns = ["class_" + str(name) for name in wide.columns]
    return wide.rename_axis("ensembl_gene_id")


def evaluate_gene_rule(rule: str, values: dict[str, float]) -> float:
    """Value of a gene rule: minimum over 'and', maximum over 'or', genes missing from values skipped (NaN if none known).

    'and' binds tighter than 'or', as in Boolean logic, so 'A or B and C' is A or (B and C).
    """
    tokens = RULE_TOKEN.findall(rule or "")
    if not tokens:
        return float("nan")
    position = 0

    def parse_or() -> float:
        nonlocal position
        operands = [parse_and()]
        while position < len(tokens) and tokens[position] == "or":
            position += 1
            operands.append(parse_and())
        known = [value for value in operands if not np.isnan(value)]
        return max(known) if known else float("nan")

    def parse_and() -> float:
        nonlocal position
        operands = [parse_atom()]
        while position < len(tokens) and tokens[position] == "and":
            position += 1
            operands.append(parse_atom())
        known = [value for value in operands if not np.isnan(value)]
        return min(known) if known else float("nan")

    def parse_atom() -> float:
        nonlocal position
        token = tokens[position]
        position += 1
        if token == "(":
            value = parse_or()
            if position >= len(tokens) or tokens[position] != ")":
                raise ValueError(f"unbalanced parentheses in gene rule {rule!r}")
            position += 1
            return value
        if token in ("and", "or", ")"):
            raise ValueError(f"unexpected {token!r} in gene rule {rule!r}")
        return float(values.get(token, float("nan")))

    value = parse_or()
    if position != len(tokens):
        raise ValueError(f"trailing tokens in gene rule {rule!r}")
    return value


def brain_expression_blocks(nodes: pd.DataFrame, gene_expression: pd.DataFrame) -> pd.DataFrame:
    """node_id-indexed gene_brain_* and reaction_brain_* blocks from an Ensembl-indexed table of raw expression values.

    gene_expression holds gtex_brain_max, gtex_other_max, region_* and class_* columns (raw TPM or nCPM). Values are
    log1p transformed, then each block is standardised over its own node type.
    """
    node_index = pd.Index(nodes.node_id, name="node_id")
    value_columns = list(gene_expression.columns)
    logged = np.log1p(gene_expression.astype(float).clip(lower=0.0))
    known_genes = set(gene_expression.index[gene_expression.filter(like="region_").notna().any(axis=1)])
    blocks = []

    is_gene = (nodes.node_type == "gene").to_numpy()
    per_gene = logged.reindex(nodes.ensembl_gene_id[is_gene].to_numpy())
    block = pd.DataFrame(0.0, index=node_index, columns=value_columns + ["has_expression"])
    block.loc[is_gene, value_columns] = standardise(per_gene).to_numpy()
    block.loc[is_gene, "has_expression"] = nodes.ensembl_gene_id[is_gene].isin(known_genes).astype(float).to_numpy()
    blocks.append(block.add_prefix("gene_brain_"))

    is_reaction = (nodes.node_type == "reaction").to_numpy()
    rules = nodes.gene_reaction_rule[is_reaction].fillna("").to_numpy()
    per_reaction = pd.DataFrame(np.nan, index=range(len(rules)), columns=value_columns)
    for column in value_columns:
        column_values = logged[column].dropna().to_dict()
        per_reaction[column] = [evaluate_gene_rule(rule, column_values) for rule in rules]
    has_known_gene = [any(token in known_genes for token in RULE_TOKEN.findall(rule)) for rule in rules]
    block = pd.DataFrame(0.0, index=node_index, columns=value_columns + ["has_expression"])
    block.loc[is_reaction, value_columns] = standardise(per_reaction).to_numpy()
    block.loc[is_reaction, "has_expression"] = np.asarray(has_known_gene, dtype=float)
    blocks.append(block.add_prefix("reaction_brain_"))
    return pd.concat(blocks, axis=1)
