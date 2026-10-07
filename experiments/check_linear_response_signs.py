"""Does the untrained linear response predict the direction of the metabolite changes measured in monogenic disease?

A parameter-free test of the stoichiometric part of the linear-response encoder (design section 5.2): for every slice
gene with laboratory-abnormality labels (mechanistic_pathway_learning/evidence/laboratory_abnormality_labels.py; e.g.
PAH: hyperphenylalaninemia, increased phenylalanine), the encoder at its initial gains (no training, no symptom label)
propagates the gene's perturbation (loss of function, sign -1) and the sign of the response at the labelled metabolite
is compared with the measured direction. The response at a metabolite is read in two ways: summed over its compartment
copies (the total pool) and at its extracellular copy alone (closest to a blood or urine measurement; pairs without
one are left out of that reading). Only labels with a direction count; a label the response does not reach (exactly
zero) is reported apart.

Agreement is compared with always predicting the majority direction ("increased") and summarised per pair and per gene
(the mean agreement of each gene's labels, so a gene with many labels does not dominate). The AUROC treats the response
as a score for "increased".

Usage:
  OMP_NUM_THREADS=1 python experiments/check_linear_response_signs.py --markdown-output docs/linear_response_sign_check.md
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evidence.laboratory_abnormality_labels import build_label_table, generic_metabolite_ids, read_chebi_relations
from mechanistic_pathway_learning.graph.cofactor_edges import CARRIER_RULES, cofactor_edge_mask
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder

STEP_COUNTS = (8, 32, 128)


def build_encoder(data, carrier_rule: str | None, steps: int) -> LinearResponseEncoder:
    """carrier_rule None gives the plain encoder, with carrier edges left as ordinary substrates and products."""
    torch.manual_seed(0)
    cofactor_edges = torch.as_tensor(cofactor_edge_mask(data.edge_source, data.edge_target, data.edge_relation, data.relation_types,
                                                         data.node_base_metabolite_id, data.node_display_name, data.is_currency,
                                                         carrier_rule=carrier_rule)) if carrier_rule is not None else None
    return LinearResponseEncoder(len(data.node_ids), data.relation_types, torch.as_tensor(data.edge_source), torch.as_tensor(data.edge_target),
                                 torch.as_tensor(data.edge_relation), torch.as_tensor(data.edge_sign), torch.as_tensor(data.structural_node_features()),
                                 32, non_propagating_nodes=torch.as_tensor(data.is_currency), cofactor_edges=cofactor_edges, num_propagation_steps=steps)


def responses_at_labels(encoder, data, labels: pd.DataFrame, batch_size: int = 32) -> pd.DataFrame:
    """Channel-0 response, divided by the sign of the channel's input weight, summed over compartment copies and at the
    extracellular copy, for every label row."""
    position_of_perturbation = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    copies_of_base: dict[str, list[int]] = {}
    extracellular_of_base: dict[str, int] = {}
    for index, (base_id, node_id) in enumerate(zip(data.node_base_metabolite_id, data.node_ids)):
        if isinstance(base_id, str) and base_id:
            copies_of_base.setdefault(base_id, []).append(index)
            if str(node_id).endswith("e"):
                extracellular_of_base[base_id] = index
    input_sign = float(torch.sign(encoder.input_weight[0]).item()) or 1.0
    totals, extracellular = np.full(len(labels), np.nan), np.full(len(labels), np.nan)
    genes = labels.gene_symbol.unique()
    for start in range(0, len(genes), batch_size):
        batch_genes = genes[start:start + batch_size]
        rows = [position_of_perturbation[gene] for gene in batch_genes]  # monogenic perturbations are identified by gene symbol
        max_nodes = max(len(data.perturbation_seeds[row]) for row in rows)
        node_index = torch.full((len(rows), max_nodes), -1, dtype=torch.long)
        sign_and_magnitude = torch.zeros((len(rows), max_nodes, 2))
        for slot, row in enumerate(rows):
            seeds = data.perturbation_seeds[row]
            node_index[slot, :len(seeds)] = torch.as_tensor(seeds)
            sign_and_magnitude[slot, :len(seeds), 0] = torch.as_tensor(data.perturbation_signs[row], dtype=torch.float32)
            sign_and_magnitude[slot, :len(seeds), 1] = torch.as_tensor(data.perturbation_magnitudes[row], dtype=torch.float32)
        with torch.no_grad():
            response = encoder.response(node_index, sign_and_magnitude)[:, :, 0].numpy() * input_sign
        for slot, gene in enumerate(batch_genes):
            for label_row in np.where(labels.gene_symbol.to_numpy() == gene)[0]:
                base_id = labels.base_metabolite_id.iloc[label_row]
                if base_id in copies_of_base:
                    totals[label_row] = response[slot, copies_of_base[base_id]].sum()
                if base_id in extracellular_of_base:
                    extracellular[label_row] = response[slot, extracellular_of_base[base_id]]
    return labels.assign(total_response=totals, extracellular_response=extracellular)


def agreement_summary(scored: pd.DataFrame, column: str) -> dict:
    signed = scored[(scored.direction != 0) & scored[column].notna()]
    reached = signed[signed[column] != 0]
    agree = np.sign(reached[column]) == reached.direction
    majority = (reached.direction == 1).mean() if len(reached) else float("nan")
    per_gene = agree.groupby(reached.gene_symbol).mean()
    per_gene_majority = (reached.direction == 1).groupby(reached.gene_symbol).mean()
    return {
        "signed pairs": len(signed),
        "reached (response not zero)": f"{len(reached)} ({len(reached) / max(len(signed), 1):.0%})",
        "sign agreement, pairs": f"{agree.mean():.3f}" if len(reached) else "n/a",
        "majority direction (increased), pairs": f"{majority:.3f}",
        "sign agreement, mean over genes": f"{per_gene.mean():.3f} ({len(per_gene)} genes)" if len(reached) else "n/a",
        "majority direction, mean over genes": f"{per_gene_majority.mean():.3f}" if len(reached) else "n/a",
        "AUROC of the response for 'increased'": f"{roc_auc_score(reached.direction == 1, reached[column]):.3f}" if reached.direction.nunique() == 2 else "n/a",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--hpo-owl", type=Path, default=Path("data/raw/hpo/hp-base.owl"))
    parser.add_argument("--hpo-annotations", type=Path, default=Path("data/raw/hpo/genes_to_phenotype.txt"))
    parser.add_argument("--chebi-obo", type=Path, default=Path("data/raw/chebi/chebi_core.obo.gz"))
    parser.add_argument("--human-gem-metabolites", type=Path, default=Path("data/raw/Human-GEM/model/metabolites.tsv"))
    parser.add_argument("--markdown-output", type=Path, default=None)
    arguments = parser.parse_args()
    torch.set_num_threads(1)
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster")
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    labels = build_label_table(arguments.hpo_owl.read_text(), pd.read_csv(arguments.hpo_annotations, sep="\t"), read_chebi_relations(arguments.chebi_obo),
                               pd.read_csv(arguments.human_gem_metabolites, sep="\t"), frozenset(generic_metabolite_ids(nodes)))
    slice_genes = set(data.perturbation_ids)
    labels = labels[labels.gene_symbol.isin(slice_genes) & labels.base_metabolite_id.isin(set(nodes.base_metabolite_id.dropna()))].reset_index(drop=True)
    signs = sorted({sign for row in data.perturbation_signs for sign in np.atleast_1d(row)})
    sections = []
    for carrier_rule in (None, *CARRIER_RULES):
        for steps in STEP_COUNTS:
            scored = responses_at_labels(build_encoder(data, carrier_rule, steps), data, labels)
            label = f"{'plain' if carrier_rule is None else f'carrier relations ({carrier_rule})'}, {steps} steps"
            sections.append((f"{label}, pool", agreement_summary(scored, "total_response")))
            sections.append((f"{label}, extracellular", agreement_summary(scored, "extracellular_response")))
            print(label, sections[-2][1]["sign agreement, pairs"], sections[-1][1]["sign agreement, pairs"], flush=True)
    keys = list(sections[0][1])
    lines = ["# Sign check of the untrained linear response against measured metabolite changes (generated by experiments/check_linear_response_signs.py)", "",
             f"{len(labels)} label rows for {labels.gene_symbol.nunique()} slice genes ({int((labels.direction != 0).sum())} with a direction); perturbation signs in the slice: {signs}. "
             "Initial gains, no training and no symptom label. 'pool': response summed over the metabolite's compartment copies; 'extracellular': its extracellular copy only.", "",
             "| reading | " + " | ".join(keys) + " |", "|---|" + "---|" * len(keys)]
    lines += [f"| {label} | " + " | ".join(str(values[key]) for key in keys) + " |" for label, values in sections]
    report = "\n".join(lines) + "\n"
    print(report)
    if arguments.markdown_output:
        arguments.markdown_output.write_text(report)


if __name__ == "__main__":
    main()
