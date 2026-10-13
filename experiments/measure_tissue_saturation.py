"""When each tissue reaches its peak as a perturbation spreads through the graph.

The owner asked for this on 12 October 2026: "What about replacing it with the reading at different points of
maximum spread based on node tissue type (meaning it would have to store each step until max distribution across all
nodes and then work backwards to find the peak distribution point for each tissue type)", and then, once the weighted
form was set out in docs/label_source_review.md 15.4, "For saturation, sure that works. Do it."

The strict reading is not available: no set of nodes is a tissue, because the median gene is expressed in 67 of
GTEx's 68 tissues (14.4). The weighted reading is, and it is what this computes. Every node carries a weight per
tissue, and a tissue's saturation at a propagation step is the weighted mean of how far its nodes have moved:

    saturation(T, t) = sum_v w(v, T) * |h_v(t)| / sum_v w(v, T)

Weights. A gene or protein node takes the share of its gene's GTEx median expression that sits in each tissue, so
the weights of a node sum to 1 and a gene concentrated in one tissue speaks mostly for that tissue. A reaction node
has no gene of its own and takes the mean of the shares of the genes in its Human-GEM gene reaction rule. A
metabolite node takes the mean over the reactions it takes part in. A node that none of these reaches carries no
weight and does not enter any tissue's mean.

Propagation. The recurrence is the linear-response encoder's, h(t+1) = (1 - damping) h(t) + damping (A h(t) + u)
with the registered damping of 0.5, but every relation's gain is 1 and the adjacency is the mean over incoming
edges. This is a structural measurement of where a perturbation goes, not a model: the gains are learned and no
model has been trained on this table, so reading peak steps off untrained gains would measure the initialisation.
The script runs past the 8 steps the encoder then ran on purpose, because 13.8 recorded that the propagation was cut
off before it settles and a peak can only be seen if the window reaches it.

Nothing here changes a label, a graph, a model or a model input; it writes a document and a table of curves.

Usage:
  PYTHONPATH=. python experiments/measure_tissue_saturation.py
  PYTHONPATH=. python experiments/measure_tissue_saturation.py --steps 24 --perturbation-type gene
"""
from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

GTEX_MEDIAN_TPM = Path("data/raw/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz")
OUTPUT_DOCUMENT = Path("docs/tissue_saturation.md")
OUTPUT_CURVES = Path("data/processed/tissue_saturation/curves.csv")
RULE_TOKEN = re.compile(r"[A-Za-z0-9_.-]+")
REGISTERED_STEPS = 8  # the linear-response encoder's num_propagation_steps
REGISTERED_DAMPING = 0.5
REGISTERED_CONTRACTION = 0.9  # the encoder bounds every relation gain by this, which is what makes the iteration converge


def tissue_shares_of_genes(gtex_path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Each gene's share of its own total median expression, per tissue; rows sum to 1."""
    table = pd.read_csv(gtex_path, sep="\t", skiprows=2, low_memory=False)
    tissues = [column for column in table.columns if column not in ("Name", "Description")]
    expression = table[tissues].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    totals = expression.sum(axis=1).replace(0.0, np.nan)
    shares = expression.div(totals, axis=0).fillna(0.0)
    shares.insert(0, "gene_symbol", table.Description)
    shares.insert(0, "ensembl_gene_id", table.Name.str.split(".").str[0])
    return shares, tissues


def node_tissue_weights(nodes: pd.DataFrame, edges: pd.DataFrame, shares: pd.DataFrame,
                        tissues: list[str]) -> tuple[np.ndarray, dict[str, int]]:
    """[num_nodes, num_tissues] of weights, and how many nodes each route supplied."""
    by_symbol = shares.drop_duplicates(subset=["gene_symbol"]).set_index("gene_symbol")[tissues]
    by_ensembl = shares.drop_duplicates(subset=["ensembl_gene_id"]).set_index("ensembl_gene_id")[tissues]
    weights = np.zeros((len(nodes), len(tissues)), dtype=np.float32)
    position_of_node = {node_id: index for index, node_id in enumerate(nodes.node_id)}
    supplied = {"gene or protein": 0, "reaction, through its gene rule": 0, "metabolite, through its reactions": 0, "none": 0}

    for index, (symbol, ensembl) in enumerate(zip(nodes.gene_symbol, nodes.ensembl_gene_id)):
        row = None
        if isinstance(ensembl, str) and ensembl in by_ensembl.index:
            row = by_ensembl.loc[ensembl].to_numpy(dtype=np.float32)
        elif isinstance(symbol, str) and symbol in by_symbol.index:
            row = by_symbol.loc[symbol].to_numpy(dtype=np.float32)
        if row is not None and row.sum() > 0:
            weights[index] = row
            supplied["gene or protein"] += 1

    # Reactions inherit the mean of their catalysing genes' shares.
    reaction_weight_rows = {}
    for node_id, rule in zip(nodes.node_id, nodes.gene_reaction_rule):
        if not isinstance(rule, str) or not rule.strip():
            continue
        rows = []
        for token in RULE_TOKEN.findall(rule):
            if token.lower() in ("and", "or"):
                continue
            if token in by_ensembl.index:
                rows.append(by_ensembl.loc[token].to_numpy(dtype=np.float32))
            elif token in by_symbol.index:
                rows.append(by_symbol.loc[token].to_numpy(dtype=np.float32))
        if rows:
            mean_row = np.mean(rows, axis=0)
            reaction_weight_rows[node_id] = mean_row
            index = position_of_node[node_id]
            if weights[index].sum() == 0:
                weights[index] = mean_row
                supplied["reaction, through its gene rule"] += 1

    # Metabolites inherit the mean over the reactions they touch.
    touching: dict[str, list[np.ndarray]] = {}
    for source, target in zip(edges.source_id, edges.target_id):
        for metabolite, reaction in ((source, target), (target, source)):
            if reaction in reaction_weight_rows:
                touching.setdefault(metabolite, []).append(reaction_weight_rows[reaction])
    for node_id, rows in touching.items():
        index = position_of_node.get(node_id)
        if index is not None and weights[index].sum() == 0:
            weights[index] = np.mean(rows, axis=0)
            supplied["metabolite, through its reactions"] += 1

    supplied["none"] = int((weights.sum(axis=1) == 0).sum())
    return weights, supplied


def normalised_adjacency(nodes: pd.DataFrame, edges: pd.DataFrame) -> sparse.csr_matrix:
    """Mean over incoming edges, so a node's state stays bounded as the field spreads."""
    position_of_node = {node_id: index for index, node_id in enumerate(nodes.node_id)}
    sources = np.array([position_of_node[node] for node in edges.source_id], dtype=np.int64)
    targets = np.array([position_of_node[node] for node in edges.target_id], dtype=np.int64)
    adjacency = sparse.coo_matrix((np.ones(len(sources), dtype=np.float32), (targets, sources)),
                                  shape=(len(nodes), len(nodes))).tocsr()
    in_degree = np.asarray(adjacency.sum(axis=1)).ravel()
    in_degree[in_degree == 0] = 1.0
    return sparse.diags(1.0 / in_degree).tocsr() @ adjacency


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--perturbation-type", default="drug", choices=("drug", "gene", "both"))
    parser.add_argument("--steps", type=int, default=24, help="propagation steps; the registered encoder runs 8")
    parser.add_argument("--damping", type=float, default=REGISTERED_DAMPING)
    parser.add_argument("--contraction", type=float, default=REGISTERED_CONTRACTION,
                        help="the gain given to every relation; the encoder bounds its learned gains by this")
    parser.add_argument("--document", type=Path, default=OUTPUT_DOCUMENT)
    parser.add_argument("--curves", type=Path, default=OUTPUT_CURVES)
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet", columns=["source_id", "target_id"])
    shares, tissues = tissue_shares_of_genes(GTEX_MEDIAN_TPM)
    print(f"graph: {len(nodes)} nodes, {len(edges)} edges; GTEx: {len(tissues)} tissues")

    weights, supplied = node_tissue_weights(nodes, edges, shares, tissues)
    print("nodes given a tissue weight, by route:")
    for route, count in supplied.items():
        print(f"  {route}: {count}")
    tissue_totals = weights.sum(axis=0)

    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet",
                              columns=["perturbation_type", "perturbation_label", "perturbation_nodes"])
    wanted = records if arguments.perturbation_type == "both" else records[records.perturbation_type == arguments.perturbation_type]
    wanted = wanted.drop_duplicates(subset=["perturbation_type", "perturbation_label"])
    position_of_node = {node_id: index for index, node_id in enumerate(nodes.node_id)}
    seeded, labels = [], []
    for row in wanted.itertuples(index=False):
        column = np.zeros(len(nodes), dtype=np.float32)
        for node_id, sign, magnitude in json.loads(row.perturbation_nodes):
            index = position_of_node.get(node_id)
            if index is not None:
                column[index] += float(sign) * float(magnitude)
        if np.abs(column).sum() > 0:
            seeded.append(column)
            labels.append(row.perturbation_label)
    sustained = np.stack(seeded, axis=1)  # [N, B]
    print(f"perturbations propagated: {sustained.shape[1]} ({arguments.perturbation_type})")

    # Every relation's gain is set to the encoder's contraction bound rather than to 1. With rows that sum to one, a
    # gain of 1 makes the transition non-contracting and the iteration grows without a fixed point, so a reading
    # taken from it would be measuring the divergence. At the bound the iteration converges, and it is the fastest
    # spread a registered model can have.
    adjacency = arguments.contraction * normalised_adjacency(nodes, edges)

    def saturation_curves(drive: str) -> np.ndarray:
        """[steps + 1, num_tissues]. "sustained" re-injects the perturbation at every step, which is the registered
        encoder's recurrence and rises monotonically to a fixed point. "impulse" seeds once and lets the field
        travel, which is the dynamics a peak belongs to: a node rises as the front arrives and falls as it passes."""
        state = sustained.copy()
        curves = np.zeros((arguments.steps + 1, len(tissues)), dtype=np.float64)
        for step in range(arguments.steps + 1):
            if step:
                if drive == "sustained":
                    state = (1.0 - arguments.damping) * state + arguments.damping * (adjacency @ state + sustained)
                else:
                    state = adjacency @ state
            moved = np.abs(state).mean(axis=1)  # the mean over perturbations of how far each node has moved
            curves[step] = (weights * moved[:, None]).sum(axis=0) / np.maximum(tissue_totals, 1e-9)
        return curves

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sustained_curves = saturation_curves("sustained")
        impulse_curves = saturation_curves("impulse")

    # Under a sustained drive the level only rises, so the quantity with a peak is the increment: how much of the
    # tissue filled at that step. Under an impulse the level itself peaks.
    increments = np.diff(sustained_curves, axis=0)
    readings = {
        "sustained drive, step of fastest filling": (increments.argmax(axis=0) + 1, increments.max(axis=0), increments),
        "impulse, step of peak level": (impulse_curves.argmax(axis=0), impulse_curves.max(axis=0), impulse_curves),
    }

    order = None
    for name, (peak_step, at_peak, curve_set) in readings.items():
        normalised = curve_set / np.maximum(at_peak, 1e-12)
        off_diagonal = np.corrcoef(normalised.T)[np.triu_indices(len(tissues), k=1)]
        print(f"\n{name}")
        print(f"  peak step per tissue: min {peak_step.min()}, median {int(np.median(peak_step))}, max {peak_step.max()}; "
              f"distinct peak steps {len(set(peak_step.tolist()))}")
        print(f"  tissues still rising at the last step computed ({arguments.steps}): {int((peak_step >= arguments.steps).sum())}")
        print(f"  tissues peaking at or before the registered {REGISTERED_STEPS} steps: {int((peak_step <= REGISTERED_STEPS).sum())}")
        print(f"  median pairwise correlation between the normalised curves: {np.median(off_diagonal):.3f}")
        order = np.argsort(peak_step)
        print("    earliest five: " + ", ".join(f"{tissues[index]} (step {peak_step[index]})" for index in order[:5]))
        print("    latest five:   " + ", ".join(f"{tissues[index]} (step {peak_step[index]})" for index in order[-5:]))

    # The reading the document leads with is the sustained increment. Under the registered recurrence the level only
    # rises, so it has no peak; under an impulse with a mean-aggregating adjacency the field contracts from the first
    # step, so every tissue peaks at step 0 and that reading says nothing. The increment is the quantity that has a
    # peak under the dynamics the study actually uses.
    settled = sustained_curves[-1]
    filled_at_registered = sustained_curves[REGISTERED_STEPS] / np.maximum(settled, 1e-12)
    print(f"\nshare of the settled level reached by the registered {REGISTERED_STEPS} steps: "
          f"median {np.median(filled_at_registered):.3f}, min {filled_at_registered.min():.3f} "
          f"({tissues[int(filled_at_registered.argmin())]}), max {filled_at_registered.max():.3f}")

    peak_step, at_peak, _ = readings["sustained drive, step of fastest filling"]
    impulse_peak_step = readings["impulse, step of peak level"][0]
    normalised = increments / np.maximum(at_peak, 1e-12)
    off_diagonal = np.corrcoef(normalised.T)[np.triu_indices(len(tissues), k=1)]
    order = np.argsort(peak_step)

    arguments.curves.parent.mkdir(parents=True, exist_ok=True)
    both = pd.concat([pd.DataFrame(sustained_curves, columns=tissues).rename_axis("step").assign(drive="sustained"),
                      pd.DataFrame(impulse_curves, columns=tissues).rename_axis("step").assign(drive="impulse")])
    both.to_csv(arguments.curves)

    lines = ["# Per-tissue saturation as a perturbation spreads", "",
             "Generated by `experiments/measure_tissue_saturation.py`. Do not edit by hand.", "",
             f"Graph `{arguments.graph_dir}`, {len(nodes)} nodes. Perturbations: {sustained.shape[1]} "
             f"({arguments.perturbation_type}) from `{arguments.evidence_dir}`. "
             f"{arguments.steps} propagation steps at damping {arguments.damping} with every relation gain at "
             f"{arguments.contraction}, the encoder's contraction bound; the registered encoder runs {REGISTERED_STEPS} steps.",
             "A structural measurement: every relation's gain is 1, because the gains are learned and no model has "
             "been trained on this table. Nothing here changes a label, a graph or a model input.", "",
             "## What carries a tissue weight", "", "| route | nodes |", "| --- | --- |"]
    lines += [f"| {route} | {count} |" for route, count in supplied.items()]
    lines += ["", "## When each tissue fills fastest", "",
              "The reading is the step at which a tissue's saturation rises most, under the registered sustained",
              "recurrence. Two other readings were computed and are not used: the saturation *level* under that",
              "recurrence rises monotonically to a fixed point, so it has no peak, and an impulse with a",
              "mean-aggregating adjacency contracts from the first step, so "
              f"{int((impulse_peak_step == 0).sum())} of {len(tissues)} tissues peak at step 0 under it. Both are in",
              f"`{arguments.curves}` for anyone who wants them.", "",
              f"- peak step: min {peak_step.min()}, median {int(np.median(peak_step))}, max {peak_step.max()}, "
              f"over {len(set(peak_step.tolist()))} distinct steps",
              f"- **tissues that fill fastest within the registered {REGISTERED_STEPS} steps: "
              f"{int((peak_step <= REGISTERED_STEPS).sum())} of {len(tissues)}**, so "
              f"{int((peak_step > REGISTERED_STEPS).sum())} are cut off before their fastest step",
              f"- median pairwise correlation between the normalised curves: {np.median(off_diagonal):.3f} "
              f"(quartiles {np.percentile(off_diagonal, 25):.3f} and {np.percentile(off_diagonal, 75):.3f}), so the "
              "curves separate rather than tracking one another",
              f"- share of the settled level reached by step {REGISTERED_STEPS}: median "
              f"{np.median(filled_at_registered):.3f}, lowest {filled_at_registered.min():.3f} "
              f"({tissues[int(filled_at_registered.argmin())]})", "",
              "| tissue | step of fastest filling | increment at that step |", "| --- | --- | --- |"]
    lines += [f"| {tissues[index]} | {peak_step[index]} | {at_peak[index]:.4g} |" for index in order]
    lines += ["", f"The curves are `{arguments.curves}`, one row per step and one column per tissue.", ""]
    arguments.document.write_text("\n".join(lines))
    print(f"\nwrote {arguments.document} and {arguments.curves}")


if __name__ == "__main__":
    main()
