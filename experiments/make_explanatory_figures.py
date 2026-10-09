"""Explanatory figures for docs/figures.md, drawn from the repository's own outputs so they can be regenerated.

  1. model_pipeline: perturbation, graph layers, linear-response propagation, output gate, pooling, head and labels;
  2. cell_class_channels: the cell-class channels of LinearResponseEncoder on a six-node chain, computed by the encoder;
  3. label_selection: the better_v2 label selection per symptom and the reasons pairs were set aside;
  4. dopaminergic_class: marker genes in the dopaminergic class against the highest Human Protein Atlas cluster type;
  5. slice_twin_comparisons: the one-change twin comparisons of the slice pilot (runs/twin_comparisons.json);
  6. slice_baselines_and_models: per-fold macro AUPRC of the slice baselines and model configurations;
  7. gene_protein_split: gene and protein nodes of graph_full_neuronal_split, with its counts;
  8. plasma_carriage: the binds edges that carry an extracellular cargo metabolite to its plasma binder protein.

Figures whose inputs are missing are skipped with a message. Writes PNG files to docs/figures/. Idempotent.

Usage:
  python experiments/make_explanatory_figures.py
  python experiments/make_explanatory_figures.py --only cell_class_channels label_selection
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder  # noqa: E402

OUTPUT_DIRECTORY = Path("docs/figures")
FULL_GRAPH_DIRECTORY = Path("data/processed/graph_full_neuronal")
SPLIT_GRAPH_DIRECTORY = Path("data/processed/graph_full_neuronal_split")
# the graph the two tested models read: the gene and protein split with the plasma carriage edges
# (docs/preregistration.md, amendments of 9 October 2026, first and fourth)
CONFIRMATORY_GRAPH_DIRECTORY = Path("data/processed/graph_full_neuronal_split_binders")
BINDER_SUMMARY = CONFIRMATORY_GRAPH_DIRECTORY / "plasma_binder_summary.json"
# the exact evidence_source values of a carriage edge. Matching the word "plasma" instead would also catch the
# electrical edges, whose sources name the plasma membrane ("net charge across the plasma membrane (Reactome)").
CARRIAGE_EVIDENCE_SOURCES = ("curated plasma carriage", "UniProt binding site (plasma binder carriage)")
LABEL_SELECTION_SUMMARY = Path("data/processed/label_selection/better_v2_full_v2.summary.json")
CELL_CLASS_WEIGHTS = Path("data/processed/cell_class_weights/full_neuronal_cell_class_weights.parquet")
MINIMUM_POSITIVES_FOR_MACRO = 5
DOPAMINERGIC_SUMMARY = Path("data/processed/brain_expression/dopaminergic_siletti_cluster395.summary.json")
TWIN_COMPARISONS = Path("runs/twin_comparisons.json")
SLICE_BASELINES = Path("runs/baselines_disease_cluster/results.json")
SLICE_AGGREGATE = Path("runs/phase3_aggregate.json")
BOX_COLOURS = {"input": "#dbe9f6", "graph": "#e3f1df", "propagation": "#fdebd3", "readout": "#efe3f4", "output": "#f4f4f4"}


def save(figure, name: str) -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIRECTORY / f"{name}.png"
    figure.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(figure)
    print(f"wrote {path}")


def graph_layer_counts(graph_directory: Path) -> dict:
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    edges = pd.read_parquet(graph_directory / "edges.parquet")
    source = edges.evidence_source.fillna("")
    layer_of_edge = np.select(
        [source.isin(["Human-GEM", "Human-GEM GPR"]), source.str.startswith("OmniPath"), source == "CollecTRI",
         source.isin(CARRIAGE_EVIDENCE_SOURCES), edges.relation_type.eq("encodes")],
        ["metabolic (Human-GEM)", "signaling (OmniPath)", "transcription (CollecTRI)",
         "plasma carriage (curated and UniProt)", "gene to protein (encodes)"],
        default="neuronal, electrical and redox (Reactome and curated)")
    return {"nodes": len(nodes), "edges": len(edges), "relations": edges.relation_type.nunique(),
            "node_types": nodes.node_type.value_counts().to_dict(), "edges_by_layer": pd.Series(layer_of_edge).value_counts().to_dict()}


def draw_box(axis, centre_x: float, centre_y: float, width: float, height: float, title: str, body: str, colour: str) -> None:
    axis.add_patch(FancyBboxPatch((centre_x - width / 2, centre_y - height / 2), width, height, boxstyle="round,pad=0.01,rounding_size=0.015",
                                  facecolor=colour, edgecolor="#555555", linewidth=1.0))
    axis.text(centre_x, centre_y + height / 2 - 0.035, title, ha="center", va="top", fontsize=9.5, fontweight="bold")
    axis.text(centre_x, centre_y + height / 2 - 0.085, body, ha="center", va="top", fontsize=7.4, linespacing=1.35)


def draw_arrow(axis, start: tuple[float, float], end: tuple[float, float]) -> None:
    axis.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=13, color="#444444", linewidth=1.2, connectionstyle="arc3,rad=0"))


def figure_model_pipeline() -> None:
    counts = graph_layer_counts(CONFIRMATORY_GRAPH_DIRECTORY)
    binder_summary = json.loads(BINDER_SUMMARY.read_text())
    label_summary = json.loads(LABEL_SELECTION_SUMMARY.read_text())
    cell_classes = [column for column in pd.read_parquet(CELL_CLASS_WEIGHTS).columns if column != "all_cells"]
    layer_lines = "\n".join(f"{layer}: {number:,} edges" for layer, number in sorted(counts["edges_by_layer"].items(), key=lambda item: -item[1]))
    node_lines = ", ".join(f"{number:,} {node_type.replace('_', ' ')}" for node_type, number in counts["node_types"].items() if number > 1)
    figure, axis = plt.subplots(figsize=(15, 8.4))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    top_row, bottom_row, height = 0.755, 0.245, 0.44
    draw_box(axis, 0.105, top_row, 0.19, height, "1. Perturbation",
             "monogenic: loss or gain of function\nat one gene (sign, magnitude)\n\ndrug: signed effect at its\nprotein targets (ChEMBL mechanisms)\n\n"
             "entered as a sustained input u\nat the perturbed nodes\n\non the split graph (figure 7):\nknockouts seed gene nodes,\ndrugs seed protein nodes",
             BOX_COLOURS["input"])
    draw_box(axis, 0.415, top_row, 0.38, height, "2. Graph (graph_full_neuronal_split_binders)",
             f"{counts['nodes']:,} nodes, {counts['edges']:,} edges, {counts['relations']} relations\n{node_lines}\n\n{layer_lines}\n\n"
             "signed edges: activation +1, inhibition and\nrepression -1, substrate depletion -1,\nplasma carriage 0\n\n"
             f"genes and proteins are separate nodes (figure 7)\nplasma carriage, {binder_summary['edges_added']} binds edges (figure 8)",
             BOX_COLOURS["graph"])
    draw_box(axis, 0.815, top_row, 0.35, height, "3. Encoder (one per tested model)",
             "message passing: node states from structural\nfeatures and descriptors, 3 layers of typed\nmessages (mean per relation), one more than on\na merged graph because a perturbed gene reaches\nits protein one hop later\n\n"
             "linear response, 8 steps:\nh(t+1) = (1 - d) h(t) + d w_k (sum_r g_r S_r h(t) + u)\n"
             "S_r: signed adjacency of relation r, normalised\ng_r: learned gain per relation and channel\n"
             f"w_k: node weight in cell class k ({len(cell_classes)} classes and all\ncells), extracellular pools shared across classes\n\n"
             "both read the difference field:\nperturbed minus unperturbed", BOX_COLOURS["propagation"])
    draw_box(axis, 0.2, bottom_row, 0.36, height, "4. Node properties",
             "fixed, per node type: protein (ESM-2 embeddings),\nmetabolite and reaction (EC) descriptors; brain\n"
             "region and cell-class expression, with the\ndopaminergic class\n\n"
             "message passing: in the node's starting state\nlinear response: in the output gate,\n"
             "field = (h W) * sigmoid(node properties W_g + b),\nwhich changes where the response is read,\nnot how it propagates\n\n"
             "treatments that keep the graph in charge (seed\nmasking, zero-initialised slow descriptor map):\nslice arms only, not applied", BOX_COLOURS["readout"])
    draw_box(axis, 0.52, bottom_row, 0.2, height, "5. Pooling and head",
             "sum over nodes\n\nnoisy-OR over pathway\nmodules (gated node\nsupports and links); each\nleak starts at its\nloss-optimal constant\n\n"
             "sigmoid head (one logit\nper symptom): slice\ncomparisons only, dropped\nfrom the tested models\non 8 October", BOX_COLOURS["readout"])
    draw_box(axis, 0.82, bottom_row, 0.32, height, "6. Symptoms and labels",
             f"{len(label_summary['all_by_symptom'])} psychiatric symptoms in the evidence\n\nlabels: HPO annotations (genes), SIDER and\n"
             "OnSIDES drug labels (drugs)\nbetter_v2: low-frequency positives and pairs\nwhose only evidence is grade C are masked,\n"
             "neither positive nor negative;\nunobserved pairs stay unlabelled\n\n"
             "folds grouped by disease cluster and drug\ntargets; lockbox_v2 held out; early-stopping\nvalidation rotated per seed; refit on\n"
             "training and validation together", BOX_COLOURS["output"])
    draw_arrow(axis, (0.2, top_row), (0.225, top_row))
    draw_arrow(axis, (0.605, top_row), (0.64, top_row))
    draw_arrow(axis, (0.81, top_row - height / 2), (0.2, bottom_row + height / 2))
    draw_arrow(axis, (0.38, bottom_row), (0.42, bottom_row))
    draw_arrow(axis, (0.62, bottom_row), (0.66, bottom_row))
    axis.set_title("The two tested models (confirmatory_v2, noisy-OR head on each encoder), from a perturbation to symptom probabilities",
                   fontsize=12, fontweight="bold")
    save(figure, "model_pipeline")


def toy_chain_encoder(cell_class_weights: torch.Tensor, extracellular_pool: torch.Tensor | None) -> LinearResponseEncoder:
    """The six-node chain of tests/test_linear_response_encoder.py with every gain and input weight equal, so the
    channels differ only by their cell-class weights."""
    relation_types = ["substrate_of", "product_of", "catalyzed_by"]
    edges = [(0, 1, "catalyzed_by", 1.0), (2, 1, "substrate_of", 1.0), (1, 3, "product_of", 1.0), (3, 4, "substrate_of", 1.0), (4, 5, "product_of", 1.0)]
    node_types = torch.tensor([0, 1, 2, 2, 1, 2])
    encoder = LinearResponseEncoder(
        6, relation_types, edge_source=torch.tensor([edge[0] for edge in edges]), edge_target=torch.tensor([edge[1] for edge in edges]),
        edge_relation=torch.tensor([relation_types.index(edge[2]) for edge in edges]), edge_sign=torch.tensor([edge[3] for edge in edges]),
        node_features=torch.nn.functional.one_hot(node_types, 3).float(), node_state_dim=4, num_propagation_steps=12,
        cell_class_weights=cell_class_weights, channels_per_cell_class=1, extracellular_pool_nodes=extracellular_pool)
    with torch.no_grad():
        encoder.gain_logit.zero_()
        encoder.input_weight.fill_(1.0)
    return encoder


def figure_cell_class_channels() -> None:
    node_names = ["gene G", "reaction R1", "substrate S", "product P\n(extracellular)", "reaction R2", "product Q"]
    weights = torch.ones(6, 2)
    weights[1, 1] = 0.0  # class B does not express reaction R1
    pool = torch.tensor([False, False, False, True, False, False])
    loss_of_function = (torch.tensor([[0]]), torch.tensor([[[-1.0, 1.0]]]))
    with torch.no_grad():
        isolated = toy_chain_encoder(weights, None).response(*loss_of_function)[0]
        pooled = toy_chain_encoder(weights, pool).response(*loss_of_function)[0]
    columns = {"class A\n(expresses R1)": isolated[:, 0], "class B\n(no R1)": isolated[:, 1],
               "class A,\nP shared": pooled[:, 0], "class B,\nP shared": pooled[:, 1]}
    values = torch.stack(list(columns.values()), dim=1).numpy()

    figure, (chain_axis, heat_axis) = plt.subplots(1, 2, figsize=(13.5, 4.8), gridspec_kw={"width_ratios": [1.1, 1]})
    positions = {0: (0.08, 0.82), 2: (0.08, 0.18), 1: (0.3, 0.5), 3: (0.6, 0.5), 4: (0.9, 0.5), 5: (0.9, 0.1)}
    edge_list = [(0, 1, "catalyzed by"), (2, 1, "substrate of"), (1, 3, "product of"), (3, 4, "substrate of"), (4, 5, "product of")]
    node_text = {}
    for node, (x, y) in positions.items():
        colour = "#fdebd3" if node in (1, 4) else ("#dbe9f6" if node == 0 else "#e3f1df")
        node_text[node] = chain_axis.text(x, y, node_names[node], ha="center", va="center", fontsize=8.5,
                                          bbox={"boxstyle": "round,pad=0.35", "facecolor": colour, "edgecolor": "#333333" if node != 3 else "#c0392b",
                                                "linewidth": 1.0 if node != 3 else 2.0})
    figure.canvas.draw()  # the text boxes need a size before arrows can be clipped to them
    for source, target, label in edge_list:
        chain_axis.annotate("", xy=(0.5, 0.5), xycoords=node_text[target], xytext=(0.5, 0.5), textcoords=node_text[source],
                            arrowprops={"arrowstyle": "-|>", "color": "#555555", "patchA": node_text[source].get_bbox_patch(),
                                        "patchB": node_text[target].get_bbox_patch(), "shrinkA": 2, "shrinkB": 2})
        middle = ((positions[source][0] + positions[target][0]) / 2, (positions[source][1] + positions[target][1]) / 2)
        offset = (0.07, 0.0) if source == 4 else (0.0, 0.07)
        chain_axis.text(middle[0] + offset[0], middle[1] + offset[1], label, fontsize=7.5, ha="center", color="#555555")
    chain_axis.text(0.6, 0.32, "shared pool\n(right two columns)", color="#c0392b", fontsize=8, ha="center")
    chain_axis.text(0.3, 0.6, "weight 0 in class B", color="#c0392b", fontsize=8, ha="center")
    chain_axis.set_xlim(-0.05, 1.05)
    chain_axis.set_ylim(0.0, 0.95)
    chain_axis.axis("off")
    chain_axis.set_title("Loss of function at G, two cell classes", fontsize=10)

    magnitude_floor = 1e-6
    signed_log = np.sign(values) * np.log10(1.0 + np.abs(values) / magnitude_floor)
    limit = float(np.abs(signed_log).max())
    image = heat_axis.imshow(signed_log, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0.0, vmin=-limit, vmax=limit), aspect="auto")
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            text = "0" if values[row, column] == 0 else f"{values[row, column]:+.2g}"
            heat_axis.text(column, row, text, ha="center", va="center", fontsize=8)
    heat_axis.set_xticks(range(len(columns)), list(columns), fontsize=8)
    heat_axis.set_yticks(range(len(node_names)), [name.replace("\n", " ") for name in node_names], fontsize=8)
    heat_axis.set_title("Response h after 12 steps (untrained, equal gains)", fontsize=10)
    colour_bar = figure.colorbar(image, ax=heat_axis, fraction=0.05)
    colour_bar.set_label("sign * log10(1 + |h| / 1e-6)", fontsize=8)
    figure.suptitle("Cell-class channels: one adjacency, each node's response scaled by its expression in the class", fontsize=11, fontweight="bold")
    figure.text(0.5, -0.04, "Class B does not express R1, so the loss at G stops at G in its channel. With P marked as a shared extracellular pool, "
                "P holds the mean of the two classes,\nso the change made in class A reaches R2 and Q in class B at half its size. Computed with LinearResponseEncoder (cell_class_weights, "
                "extracellular_pool_nodes) on the six-node chain of the unit tests.", ha="center", fontsize=8)
    save(figure, "cell_class_channels")


def figure_label_selection() -> None:
    summary = json.loads(LABEL_SELECTION_SUMMARY.read_text())
    symptoms = sorted(summary["all_by_symptom"], key=lambda symptom: summary["all_by_symptom"][symptom])
    kept = np.array([summary["kept_by_symptom"].get(symptom, 0) for symptom in symptoms])
    set_aside = np.array([summary["all_by_symptom"][symptom] for symptom in symptoms]) - kept
    figure, (symptom_axis, reason_axis) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.3, 1]})
    labels = [symptom.replace("_", " ") for symptom in symptoms]
    symptom_axis.barh(labels, kept, color="#2b6ca3", label="kept")
    symptom_axis.barh(labels, set_aside, left=kept, color="#bcd3e8", label="set aside (masked)")
    for row, (kept_count, set_aside_count) in enumerate(zip(kept, set_aside)):
        symptom_axis.text(kept_count + set_aside_count + 5, row, f"{kept_count} of {kept_count + set_aside_count}", va="center", fontsize=8)
    symptom_axis.set_xlabel("positive (perturbation, symptom) pairs, full evidence (evidence_full_v2)")
    symptom_axis.set_xlim(0, max(kept + set_aside) * 1.2)
    symptom_axis.legend(loc="lower right", fontsize=8)
    symptom_axis.set_title(f"{summary['kept_pairs']:,} of {summary['positive_pairs']:,} positive pairs kept", fontsize=10)
    reasons = summary["set_aside_reasons"]
    reason_labels = [reason.replace(" x ", " ").replace("frequency below", "frequency\nbelow").replace(", one label source", ",\none label source") for reason in reasons]
    reason_axis.barh(reason_labels, list(reasons.values()), color="#9aa5b1")
    for row, number in enumerate(reasons.values()):
        reason_axis.text(number + 8, row, str(number), va="center", fontsize=8)
    reason_axis.invert_yaxis()
    reason_axis.set_xlim(0, max(reasons.values()) * 1.2)
    reason_axis.tick_params(axis="y", labelsize=8)
    reason_axis.set_xlabel("positive pairs set aside")
    reason_axis.set_title("Why pairs were set aside", fontsize=10)
    too_few = [f"{symptom.replace('_', ' ')} ({summary['kept_by_symptom'].get(symptom, 0)})" for symptom in symptoms
               if summary["kept_by_symptom"].get(symptom, 0) < MINIMUM_POSITIVES_FOR_MACRO]
    figure.suptitle("Better (not more) training examples: selection better_v2 of the tested models (the user's decisions of 7 and 8 October 2026)",
                    fontsize=11, fontweight="bold")
    figure.text(0.5, -0.09, f"Kept: gene pairs at HPO frequency {summary['gene_minimum_frequency']:.2f} or more (Frequent and above); drug pairs at label "
                f"frequency {summary['drug_minimum_frequency']:.0%} or more, or listed by both SIDER and OnSIDES. A pair set aside stays out of the loss "
                f"and every metric;\nit is not a negative. Also masked: {summary['masked_negative_pairs']:,} pairs whose only evidence is grade C (a human "
                f"association), which were negatives before. Too few kept pairs for the macro average ({MINIMUM_POSITIVES_FOR_MACRO} needed): "
                f"{', '.join(too_few)}.", ha="center", fontsize=8)
    save(figure, "label_selection")


def figure_dopaminergic_class() -> None:
    summary = json.loads(DOPAMINERGIC_SUMMARY.read_text())
    genes = ["TH", "SLC6A3", "SLC18A2", "DDC", "KCNJ6", "LMX1B", "EN1", "NR4A2", "GAD1", "AQP4", "MBP"]
    anchors = summary["anchor_genes"]
    dopaminergic = np.array([anchors[gene]["dopaminergic_ncpm"] for gene in genes])
    hpa_highest = np.array([anchors[gene]["hpa_max_ncpm"] for gene in genes])
    positions = np.arange(len(genes))
    figure, axis = plt.subplots(figsize=(12, 4.6))
    axis.bar(positions - 0.2, dopaminergic + 0.1, width=0.4, color="#c0392b", label="dopaminergic class (atlas cluster 395)")
    axis.bar(positions + 0.2, hpa_highest + 0.1, width=0.4, color="#7f8c8d", label="highest of HPA's 34 cluster types (named below each gene)")
    for position, gene in zip(positions, genes):
        axis.annotate(anchors[gene]["hpa_top_cluster_type"].replace("vascular associated", "vascular"), xy=(position + 0.2, 0), xycoords=("data", "axes fraction"), xytext=(0, -22),
                      textcoords="offset points", ha="right", va="top", rotation=30, rotation_mode="anchor", fontsize=6.5, color="#555555")
    axis.set_yscale("log")
    axis.set_ylim(0.5, 2e5)
    axis.set_xticks(positions, genes)
    axis.set_ylabel("nCPM (log scale; 0.1 added)")
    axis.axvline(7.5, color="#bbbbbb", linestyle=":", linewidth=1)
    axis.text(5.2, 6e4, "dopamine neuron markers", ha="center", fontsize=8.5)
    axis.text(9, 6e4, "GABAergic and glial markers", ha="center", fontsize=8.5)
    axis.legend(fontsize=8, loc="upper left")
    nuclei = summary["nuclei"]
    donors = sorted({donor for entry in summary["files"] for donor in entry["nuclei_by_donor"]})
    axis.set_title(f"The dopaminergic cell class: {nuclei} nuclei from {len(donors)} donors, on HPA's nCPM scale "
                   f"(trimmed-mean factor {summary['scale_factor_cpm_per_ncpm']:.2f})", fontsize=10.5, fontweight="bold")
    figure.text(0.5, -0.27, "HPA's single-nucleus cluster types have no dopaminergic type (SLC6A3 peaks at 0.6 nCPM). Cluster 395 of the same atlas "
                "(Siletti et al. 2023, CC BY 4.0, via CZ CELLxGENE) is summed and scaled against HPA's splatter type.\nMBP at about 2 percent of its "
                "oligodendrocyte value is consistent with ambient RNA in single-nucleus data.", ha="center", fontsize=8)
    save(figure, "dopaminergic_class")


def short_change(change: str) -> str:
    return (change.replace("data/processed/node_descriptors/", "").replace("data/processed/", "").replace("--", "")
            .replace(".parquet", "").replace("_", " "))


def figure_slice_twin_comparisons() -> None:
    entries = [entry for entry in json.loads(TWIN_COMPARISONS.read_text()) if "set_aside" not in entry]
    entries.sort(key=lambda entry: entry["pooled_macro_auprc"]["difference"])
    figure, axis = plt.subplots(figsize=(11, 0.34 * len(entries) + 1.8))
    rows = np.arange(len(entries))
    for offset, key, colour, label in ((-0.17, "pooled_macro_auprc", "#2b6ca3", "pooled, paired bootstrap"),
                                       (0.17, "within_degree_strata_macro_auprc", "#e67e22", "within degree strata")):
        differences = np.array([entry[key]["difference"] for entry in entries])
        lower = np.array([entry[key]["lower"] for entry in entries])
        upper = np.array([entry[key]["upper"] for entry in entries])
        axis.errorbar(differences, rows + offset, xerr=[differences - lower, upper - differences], fmt="o", color=colour, markersize=3.5,
                      elinewidth=1.1, capsize=2, label=label)
    axis.axvline(0.0, color="#333333", linewidth=0.8)
    for threshold in (-0.05, 0.05):
        axis.axvline(threshold, color="#999999", linewidth=0.8, linestyle="--")
    axis.set_yticks(rows, [f"{entry['a']}\n  vs {entry['b']}: {short_change(entry['change'])}" for entry in entries], fontsize=6.3)
    axis.set_xlabel("A minus B, macro AUPRC (95 percent interval; dashed: the preregistered minimum difference of 0.05)")
    axis.legend(fontsize=8, loc="lower right")
    excluding_zero = sum(1 for entry in entries if entry["pooled_macro_auprc"]["lower"] > 0 or entry["pooled_macro_auprc"]["upper"] < 0)
    axis.set_title(f"Slice pilot: {len(entries)} one-change comparisons; {excluding_zero} pooled intervals exclude zero "
                   f"(about {0.05 * len(entries):.1f} expected by chance, uncorrected)", fontsize=10.5, fontweight="bold")
    save(figure, "slice_twin_comparisons")


def configuration_label(aggregate_name: str) -> str:
    """The aggregate names a run by its directory, with the path when two directories share a name: runs/encoder/ holds
    the reruns under the current code."""
    name = aggregate_name.removesuffix("_disease_cluster")
    if name.startswith("runs/encoder/"):
        return name.removeprefix("runs/encoder/") + " (current code)"
    if name.startswith("runs/"):
        return name.removeprefix("runs/") + " (earlier code)"
    return name


def figure_slice_baselines_and_models() -> None:
    baselines = json.loads(SLICE_BASELINES.read_text())["splits"]["grouped"]
    runs = json.loads(SLICE_AGGREGATE.read_text())["runs"]
    rows = [(name.replace("_", " ") + " (baseline)", entry["per_fold_macro_auprc_mean"], entry["per_fold_macro_auprc_sd"], True)
            for name, entry in baselines.items()]
    rows += [(configuration_label(name), entry["per_fold_macro_auprc_mean"], entry["per_fold_macro_auprc_sd"], False)
             for name, entry in runs.items() if entry.get("num_splits") == 5]
    rows.sort(key=lambda row: row[1])
    figure, axis = plt.subplots(figsize=(10, 0.27 * len(rows) + 1.6))
    positions = np.arange(len(rows))
    axis.barh(positions, [row[1] for row in rows], xerr=[row[2] for row in rows], color=["#7f8c8d" if row[3] else "#2b6ca3" for row in rows],
              error_kw={"elinewidth": 0.8, "capsize": 1.5})
    random_walk = baselines["random_walk_with_restart"]["per_fold_macro_auprc_mean"]
    axis.axvline(random_walk, color="#c0392b", linewidth=1, linestyle="--")
    axis.text(random_walk + 0.002, len(rows) - 0.5, "random walk", color="#c0392b", fontsize=8)
    axis.set_yticks(positions, [row[0] for row in rows], fontsize=6.5)
    axis.set_xlim(0.15, max(row[1] + row[2] for row in rows) + 0.02)
    axis.set_xlabel("macro AUPRC per disease-cluster fold, mean ± standard deviation over 5 folds (n divisor)")
    configurations = [row for row in rows if not row[3]]
    above_random_walk = sum(1 for row in configurations if row[1] > random_walk)
    axis.set_title(f"Slice pilot: {len(configurations)} five-fold configurations, {above_random_walk} with a mean above the random walk",
                   fontsize=10.5, fontweight="bold")
    save(figure, "slice_baselines_and_models")


def draw_node_link_diagram(axis, nodes: dict, edges: list) -> None:
    """nodes: name -> (x, y, label, colour); edges: (source, target, label). Arrows are clipped to the text boxes."""
    text_of_node = {}
    for name, (x, y, label, colour) in nodes.items():
        text_of_node[name] = axis.text(x, y, label, ha="center", va="center", fontsize=7.8,
                                       bbox={"boxstyle": "round,pad=0.35", "facecolor": colour, "edgecolor": "#333333"})
    axis.figure.canvas.draw()  # the text boxes need a size before arrows can be clipped to them
    for source, target, label in edges:
        axis.annotate("", xy=(0.5, 0.5), xycoords=text_of_node[target], xytext=(0.5, 0.5), textcoords=text_of_node[source],
                      arrowprops={"arrowstyle": "-|>", "color": "#555555", "patchA": text_of_node[source].get_bbox_patch(),
                                  "patchB": text_of_node[target].get_bbox_patch(), "shrinkA": 2, "shrinkB": 2})
        if label:
            source_x, source_y = nodes[source][:2]
            target_x, target_y = nodes[target][:2]
            axis.text((source_x + target_x) / 2 + 0.01, (source_y + target_y) / 2, label, fontsize=7, color="#555555", ha="left", va="center")


def figure_gene_protein_split() -> None:
    summary = json.loads((SPLIT_GRAPH_DIRECTORY / "split_summary.json").read_text())
    gene_to_protein = pd.read_parquet(SPLIT_GRAPH_DIRECTORY / "gene_to_protein.parquet")
    entries_of_protein = gene_to_protein.drop_duplicates("protein_node_id").set_index("protein_node_id").uniprot_entry
    calca_nodes = sorted(gene_to_protein.loc[gene_to_protein.gene_node_id == "GENE:CALCA", "protein_node_id"])
    gene_colour, protein_colour, seed_colour = "#dbe9f6", "#e3f1df", "#fdebd3"
    figure, axis = plt.subplots(figsize=(14, 5.6))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    nodes = {
        "knockout": (0.06, 0.88, "knockout", seed_colour),
        "transcription_factor": (0.06, 0.45, "transcription factor\n(protein node)", protein_colour),
        "gene_calca": (0.2, 0.67, "GENE:CALCA", gene_colour),
        "calcitonin": (0.36, 0.84, f"{calca_nodes[0]}\ncalcitonin", protein_colour),
        "cgrp": (0.36, 0.5, f"{calca_nodes[1]}\nCGRP", protein_colour),
        "drug": (0.52, 0.5, "drug whose ChEMBL\ntarget names CGRP", seed_colour),
        "gene_cdkn2a": (0.66, 0.84, "GENE:CDKN2A", gene_colour),
        "protein_cdkn2a": (0.66, 0.5, f"PROTEIN:CDKN2A\n{entries_of_protein['PROTEIN:CDKN2A']}", protein_colour),
        "gene_h33a": (0.84, 0.84, "GENE:H3-3A", gene_colour),
        "gene_h33b": (0.96, 0.84, "GENE:H3-3B", gene_colour),
        "protein_h33": (0.9, 0.5, f"PROTEIN:H3-3A\n{entries_of_protein['PROTEIN:H3-3A']} (H3.3)", protein_colour),
    }
    edges = [("knockout", "gene_calca", ""), ("transcription_factor", "gene_calca", "regulates_\ntranscription_of"),
             ("gene_calca", "calcitonin", "encodes"), ("gene_calca", "cgrp", "encodes"), ("drug", "cgrp", ""),
             ("gene_cdkn2a", "protein_cdkn2a", "encodes"), ("gene_h33a", "protein_h33", ""), ("gene_h33b", "protein_h33", "")]
    draw_node_link_diagram(axis, nodes, edges)
    captions = [(0.2, 0.3, "Entries the data tell apart: one protein\nnode each. The only such gene on the\nfull graph: a ChEMBL target names CGRP."),
                (0.66, 0.3, "Entries the data do not tell apart: one\nnode for all of them. OmniPath copies\nCDKN2A's rows to both entries."),
                (0.9, 0.3, "One entry from several genes: a shared\nnode, so a knockout of one gene\nleaves the other's route.")]
    for x, y, caption in captions:
        axis.text(x, y, caption, ha="center", va="top", fontsize=8)
    edge_ends = summary["edge_ends_at_genes_with_several_proteins"]
    axis.text(0.5, 0.02,
              f"{summary['gene_nodes']:,} gene nodes and {summary['protein_nodes']:,} protein nodes joined by {summary['encodes_edges']:,} encodes edges; "
              f"{summary['protein_nodes_of_several_entries']} protein nodes hold several entries, {summary['shared_protein_nodes']} are shared by "
              f"{summary['genes_on_a_shared_protein_node']} genes, {summary['genes_with_several_protein_nodes']} gene has several protein nodes and "
              f"{summary['gene_nodes_without_a_protein_node']} gene nodes have none.\n"
              f"Edges {summary['edges_before']:,} before, {summary['edges_after']:,} after. Edge ends at the split gene: {edge_ends['named_subset']} name one "
              f"entry, {edge_ends['named_all']} name both. Gene nodes carry the brain expression columns, protein nodes the protein descriptors.",
              ha="center", va="bottom", fontsize=8)
    axis.set_title("Gene and protein nodes (graph_full_neuronal_split): a gene's entries are split only where the data tell them apart",
                   fontsize=11, fontweight="bold")
    save(figure, "gene_protein_split")


def carriage_edges_with_names(graph_directory: Path) -> pd.DataFrame:
    """The plasma carriage edges of a graph, with the display name of the cargo metabolite and of the binder."""
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    display_name_of_node = nodes.set_index("node_id").display_name
    gene_symbol_of_node = nodes.set_index("node_id").gene_symbol
    edges = pd.read_parquet(graph_directory / "edges.parquet")
    carriage = edges[edges.evidence_source.isin(CARRIAGE_EVIDENCE_SOURCES)].copy()
    carriage["cargo_name"] = carriage.source_id.map(display_name_of_node)
    carriage["binder_name"] = carriage.target_id.map(gene_symbol_of_node).fillna(carriage.target_id.map(display_name_of_node))
    return carriage


def figure_plasma_carriage() -> None:
    summary = json.loads(BINDER_SUMMARY.read_text())
    carriage = carriage_edges_with_names(CONFIRMATORY_GRAPH_DIRECTORY)
    binder_of_most_carried_cargo = carriage.groupby("source_id").target_id.nunique().sort_values()
    worked_cargo_node = binder_of_most_carried_cargo.index[-1]
    worked_cargo_name = carriage.loc[carriage.source_id == worked_cargo_node, "cargo_name"].iloc[0]
    worked_binders = sorted(carriage.loc[carriage.source_id == worked_cargo_node, "target_id"])
    first_binder_symbol = carriage.loc[carriage.target_id == worked_binders[0], "binder_name"].iloc[0]

    cargo_colour, protein_colour, gene_colour, seed_colour = "#e8e2f4", "#e3f1df", "#dbe9f6", "#fdebd3"
    figure, (left_axis, right_axis) = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw={"width_ratios": [1.05, 1.0]})
    # the axes are placed before the node-link diagram is drawn: its arrows are anchored to the text artists, so a
    # later tight_layout would move the boxes out from under the arrows and bbox_inches="tight" would then save a
    # figure wide enough to hold them
    figure.subplots_adjust(left=0.04, right=0.97, top=0.86, bottom=0.08, wspace=0.18)
    left_axis.set_xlim(0, 1)
    left_axis.set_ylim(0, 1)
    left_axis.axis("off")

    worked_nodes = {"knockout": (0.17, 0.93, "knockout of the binder gene", seed_colour),
                    "binder_gene": (0.17, 0.74, f"GENE:{first_binder_symbol}\n(brain expression columns)", gene_colour),
                    # off the row of any binder, so that no arrow runs horizontally under its own label
                    "cargo": (0.17, 0.22, f"{worked_cargo_node}\n{worked_cargo_name}\n(extracellular)", cargo_colour)}
    worked_edges = [("knockout", "binder_gene", ""), ("binder_gene", worked_binders[0], "encodes")]
    for index, binder_node in enumerate(worked_binders):
        binder_symbol = carriage.loc[carriage.target_id == binder_node, "binder_name"].iloc[0]
        worked_nodes[binder_node] = (0.62, 0.86 - index * 0.28, f"{binder_node}\n({binder_symbol}, protein descriptors)", protein_colour)
        worked_edges.append(("cargo", binder_node, "binds, sign 0"))
    draw_node_link_diagram(left_axis, worked_nodes, worked_edges)
    left_axis.text(0.5, 0.06, f"One cargo with the most binders in the graph: {worked_cargo_name}. Carriage sits on the protein node, because binding is\n"
                              "the protein's property; a knockout seeds the gene node and reaches the carriage one encodes hop later, which is what moves\n"
                              "a binder gene's degree stratum (docs/plasma_binder_confirmatory_effect.md).", ha="center", va="bottom", fontsize=8)
    left_axis.set_title("How a carriage edge sits in the graph", fontsize=11, fontweight="bold")

    edges_per_binder = carriage.groupby("binder_name").size().sort_values()
    cargo_per_binder = carriage.groupby("binder_name").cargo_name.apply(lambda names: ", ".join(sorted(set(names))))
    positions = np.arange(len(edges_per_binder))
    right_axis.set_xlim(0, float(edges_per_binder.max()) * 3.4)
    right_axis.set_ylim(-0.8, len(edges_per_binder) - 0.2)
    right_axis.barh(positions, edges_per_binder.to_numpy(), color="#7a9cc6", height=0.62)
    right_axis.set_yticks(positions)
    right_axis.set_yticklabels(edges_per_binder.index, fontsize=9)
    right_axis.set_xlabel("binds edges")
    right_axis.set_xticks(np.arange(0, int(edges_per_binder.max()) + 1, 2))
    for position, binder_name in zip(positions, edges_per_binder.index):
        right_axis.text(edges_per_binder[binder_name] + 0.25, position, cargo_per_binder[binder_name], va="center", fontsize=6.8, color="#333333")
    sources = ", ".join(f"{number} {name.replace('_', ' ')}" for name, number in sorted(summary["edges_by_source"].items(), key=lambda item: -item[1]))
    right_axis.set_title(f"{summary['edges_added']} carriage edges over {len(edges_per_binder)} binders ({sources})", fontsize=11, fontweight="bold")
    right_axis.spines[["top", "right"]].set_visible(False)

    figure.suptitle("Plasma carriage in the confirmatory graph: an extracellular cargo metabolite binds its carrier protein (relation binds, sign 0)",
                    fontsize=12, fontweight="bold")
    save(figure, "plasma_carriage")


FIGURES = {
    "model_pipeline": (figure_model_pipeline, [CONFIRMATORY_GRAPH_DIRECTORY / "edges.parquet", BINDER_SUMMARY,
                                               LABEL_SELECTION_SUMMARY, CELL_CLASS_WEIGHTS]),
    "cell_class_channels": (figure_cell_class_channels, []),
    "label_selection": (figure_label_selection, [LABEL_SELECTION_SUMMARY]),
    "dopaminergic_class": (figure_dopaminergic_class, [DOPAMINERGIC_SUMMARY]),
    "slice_twin_comparisons": (figure_slice_twin_comparisons, [TWIN_COMPARISONS]),
    "slice_baselines_and_models": (figure_slice_baselines_and_models, [SLICE_BASELINES, SLICE_AGGREGATE]),
    "gene_protein_split": (figure_gene_protein_split, [SPLIT_GRAPH_DIRECTORY / "split_summary.json", SPLIT_GRAPH_DIRECTORY / "gene_to_protein.parquet"]),
    "plasma_carriage": (figure_plasma_carriage, [CONFIRMATORY_GRAPH_DIRECTORY / "edges.parquet", BINDER_SUMMARY]),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", nargs="*", choices=list(FIGURES), default=None)
    arguments = parser.parse_args()
    for name in arguments.only or FIGURES:
        draw, inputs = FIGURES[name]
        missing = [str(path) for path in inputs if not path.exists()]
        if missing:
            print(f"skipped {name}: missing {', '.join(missing)}")
            continue
        draw()


if __name__ == "__main__":
    main()
