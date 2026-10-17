"""Explanatory figures for docs/figures.md, drawn from the repository's own outputs so they can be regenerated.

  1. model_pipeline: perturbation, graph layers, linear-response propagation, output gate, pooling, head and labels;
  2. cell_class_channels: the cell-class channels of LinearResponseEncoder on a six-node chain, computed by the encoder;
  3. label_selection: the better_v2 label selection per symptom and the reasons pairs were set aside;
  4. dopaminergic_class: marker genes in the dopaminergic class against the highest Human Protein Atlas cluster type;
  5. slice_twin_comparisons: the one-change twin comparisons of the slice pilot (runs/twin_comparisons.json);
  6. slice_baselines_and_models: per-fold macro AUPRC of the slice baselines and model configurations;
  7. gene_protein_split: gene and protein nodes of graph_full_neuronal_split, with its counts;
  8. plasma_carriage: the binds edges that carry an extracellular cargo metabolite to its plasma binder protein;
  9. cell_state: what the graph computes and does not pass on (pH, energy charge, redox and the ion pools);
  10. graph_composition: the node types and the edge layers of the current graph, with each layer's signs;
  11. allele_perturbations: the allele as the perturbation (decision 9.54), HTT's two alleles and where each
      disease's mechanism came from;
  12. edge_priors: the prior that an edge is real by its evidence, the cost of closing it, and how a hard concrete
      gate reads at its prior (decisions 9.57 and 9.62);
  13. edge_gate_training: how the hard gates are fitted, stage by stage, the prior as the cost per edge in each pass,
      and the proxy gradient of the binary pass;
  14. gate_share_reading: a closed gate is a deleted edge in both encoders, computed with their adjacency functions,
      and the sole input whose cliff the slope cannot see;
  15. gate_simulation_recovery: in planted worlds, which closures to believe, the calibration of the posterior from
      the exact evidence, and the first-order flip estimate against the exact change;
  16. gate_simulation_conditions: the same when the favourable assumptions are dropped one at a time, and how much
      the labels say per edge in the worlds and on v5 (docs/edge_gate_simulation_results.json).

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
from matplotlib.path import Path as MatplotlibPath  # noqa: E402

from mechanistic_pathway_learning.graph.tag_currency_metabolites import currency_flags  # noqa: E402
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder  # noqa: E402

OUTPUT_DIRECTORY = Path("docs/figures")
FULL_GRAPH_DIRECTORY = Path("data/processed/graph_full_neuronal")
SPLIT_GRAPH_DIRECTORY = Path("data/processed/graph_full_neuronal_split")
# the graph the registered model reads with the v5 tables: the gene and protein split with the plasma carriage edges
# (docs/preregistration.md, amendments of 9 October 2026, first and fourth), the drug entry nodes (the user, 10 October
# 2026: "i definitely want drugs as nodes with an ablation on them") and the neomorphic allele node (decision 9.54)
CONFIRMATORY_GRAPH_DIRECTORY = Path("data/processed/graph_full_neuronal_split_binders_drugs_v5")
BINDER_SUMMARY = CONFIRMATORY_GRAPH_DIRECTORY / "plasma_binder_summary.json"
EVIDENCE_DIRECTORY = Path("data/processed/evidence_full_v5")
EDGE_PRIOR_TIERS = Path("configs/edge_prior_tiers.csv")
# the exact evidence_source values of a carriage edge. Matching the word "plasma" instead would also catch the
# electrical edges, whose sources name the plasma membrane ("net charge across the plasma membrane (Reactome)").
CARRIAGE_EVIDENCE_SOURCES = ("curated plasma carriage", "UniProt binding site (plasma binder carriage)")
LABEL_SELECTION_SUMMARY = Path("data/processed/label_selection/better_v2_full_v5.summary.json")
CELL_CLASS_WEIGHTS = Path("data/processed/cell_class_weights/full_neuronal_split_drugs_v5_cell_class_weights.parquet")
MINIMUM_POSITIVES_FOR_MACRO = 5
DOPAMINERGIC_SUMMARY = Path("data/processed/brain_expression/dopaminergic_siletti_cluster395.summary.json")
TWIN_COMPARISONS = Path("runs/twin_comparisons.json")
SLICE_BASELINES = Path("runs/baselines_disease_cluster/results.json")
SLICE_AGGREGATE = Path("runs/phase3_aggregate.json")
# 150 left the dense pipeline figure soft on a high-density display (the user, 9 October 2026)
# the cell-state figure is about what the graph holds, so it reads the graph the model reads
CELL_STATE_GRAPH_DIRECTORY = CONFIRMATORY_GRAPH_DIRECTORY
NEURONAL_POOL_PREFIX = "NEURON_"
# the species whose deviation is a piece of cell state rather than a metabolite of interest in itself; grouped as the
# quantity each one reports, since most of them are read as a ratio (graph/compartment_ph.py for the argument)
CELL_STATE_GROUPS = {
    "pH (H+)": ("H+",),
    "energy charge (ATP, ADP, AMP)": ("ATP", "ADP", "AMP"),
    "redox (NAD, NADP, FAD, GSH)": ("NAD+", "NADH", "NADP+", "NADPH", "FAD", "FADH2", "GSH", "GSSG", "FMN", "FMNH2"),
    "acyl carriers (CoA, carnitine)": ("CoA", "acetyl-CoA", "L-carnitine"),
    "acid-base (HCO3-, CO2)": ("HCO3-", "CO2"),
    "oxygen (O2)": ("O2",),
    "methylation (SAM, SAH, THF)": ("SAM", "SAH", "THF"),
    "metabolic ion pools": ("Na+", "K+", "Ca2+", "Mg2+", "Cl-"),
    "nitrogen and sulfur (NH3, H2S)": ("NH3", "H2S"),
}
SAVE_DPI = 240
BOX_COLOURS = {"input": "#dbe9f6", "graph": "#e3f1df", "propagation": "#fdebd3", "readout": "#efe3f4", "output": "#f4f4f4"}


def save(figure, name: str) -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIRECTORY / f"{name}.png"
    figure.savefig(path, dpi=SAVE_DPI, bbox_inches="tight")
    plt.close(figure)
    print(f"wrote {path}")


# Each edge's layer, read from its evidence_source, first rule that matches. The electrical rule comes before the Reactome
# one is reached because "net charge across the plasma membrane (Reactome)" names Reactome at its end, not its start.
OTHER_CURATED_LAYER = "neuronal, redox and other curated"


def edge_layers(edges: pd.DataFrame) -> np.ndarray:
    source = edges.evidence_source.fillna("")
    rules = [
        ("metabolic (Human-GEM)", source.isin(["Human-GEM", "Human-GEM GPR"])),
        ("signalling (OmniPath)", source.str.startswith("OmniPath")),
        ("transcription (CollecTRI)", source == "CollecTRI"),
        ("gene to protein (encodes)", edges.relation_type.eq("encodes")),
        ("Reactome complexes, sets and regulators", source.str.startswith("Reactome")),
        ("enzyme regulation (BRENDA)", source.str.startswith("BRENDA")),
        ("electrical (ions and charge)", source.str.contains("across the plasma membrane", regex=False)),
        ("drug mechanism (ChEMBL)", source.str.startswith("drug mechanism")),
        ("drug plasma sequestration", source.str.startswith("plasma carrier of the drug")),
        ("plasma carriage of metabolites", source.isin(CARRIAGE_EVIDENCE_SOURCES)),
        ("neomorphic allele mechanism", source.str.startswith("neomorphic allele mechanism")),
    ]
    return np.select([condition for _, condition in rules], [name for name, _ in rules], default=OTHER_CURATED_LAYER)


def graph_layer_counts(graph_directory: Path) -> dict:
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    edges = pd.read_parquet(graph_directory / "edges.parquet")
    layer_of_edge = edge_layers(edges)
    return {"nodes": len(nodes), "edges": len(edges), "relations": edges.relation_type.nunique(),
            "node_types": nodes.node_type.value_counts().to_dict(), "edges_by_layer": pd.Series(layer_of_edge).value_counts().to_dict()}


BOX_BODY_FONT_SIZE, BOX_TITLE_FONT_SIZE, BOX_LINE_SPACING = 8.8, 11.0, 1.4
# all in inches, because the pipeline figure's axis is drawn in inches (see figure_model_pipeline)
BOX_INNER_PAD, BOX_TITLE_GAP, SMALLEST_BOX_GAP, FIGURE_MARGIN, ROW_GUTTER = 0.22, 0.16, 0.5, 0.22, 1.0


def text_size_in_inches(figure, axis, text: str, font_size: float, weight: str = "normal") -> tuple[float, float]:
    """The width and height a text block occupies, in inches, measured by drawing it once and removing it.
    get_window_extent reports pixels at the canvas dpi, which is why the division is by figure.dpi and not by the
    dpi the file is finally saved at. The pipeline figure sizes its boxes from these measurements rather than from
    hand-set constants, so no box is left with dead space under its last line and no line can outgrow its box."""
    artist = axis.text(0.0, 0.0, text, fontsize=font_size, fontweight=weight, linespacing=BOX_LINE_SPACING)
    figure.canvas.draw()
    extent = artist.get_window_extent(figure.canvas.get_renderer())
    artist.remove()
    return extent.width / figure.dpi, extent.height / figure.dpi


def measured_box(figure, axis, title: str, body: str, colour: str) -> dict:
    title_width, title_height = text_size_in_inches(figure, axis, title, BOX_TITLE_FONT_SIZE, "bold")
    body_width, body_height = text_size_in_inches(figure, axis, body, BOX_BODY_FONT_SIZE)
    return {"title": title, "body": body, "colour": colour, "title_height": title_height,
            "width": max(title_width, body_width) + 2 * BOX_INNER_PAD,
            "height": title_height + BOX_TITLE_GAP + body_height + 2 * BOX_INNER_PAD}


def draw_box(axis, box: dict, left: float, top: float, height: float) -> None:
    axis.add_patch(FancyBboxPatch((left, top - height), box["width"], height, boxstyle="round,pad=0,rounding_size=0.1",
                                  facecolor=box["colour"], edgecolor="#555555", linewidth=1.1, zorder=1))
    centre_x = left + box["width"] / 2
    axis.text(centre_x, top - BOX_INNER_PAD, box["title"], ha="center", va="top", fontsize=BOX_TITLE_FONT_SIZE,
              fontweight="bold", zorder=2)
    axis.text(centre_x, top - BOX_INNER_PAD - box["title_height"] - BOX_TITLE_GAP, box["body"], ha="center", va="top",
              fontsize=BOX_BODY_FONT_SIZE, linespacing=BOX_LINE_SPACING, zorder=2)


def draw_arrow(axis, vertices: list[tuple[float, float]]) -> None:
    """A straight arrow between two points, an elbow through every point given. The step from the last box of one
    row to the first box of the next is routed through the gutter between the rows; drawn as a single diagonal it
    crossed the whole figure and read as an error."""
    path = MatplotlibPath(vertices, [MatplotlibPath.MOVETO] + [MatplotlibPath.LINETO] * (len(vertices) - 1))
    axis.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=18, color="#444444", linewidth=1.5,
                                   joinstyle="miter", zorder=3))


def edge_prior_counts(graph_directory: Path) -> dict:
    """How many edges the edge priors gate on this graph (decision 9.57), read from its edge evidence and the tiers."""
    from mechanistic_pathway_learning.graph.edge_evidence import edge_priors, read_edge_tiers

    evidence = pd.read_parquet(graph_directory / "edge_evidence.parquet")
    gated, _, _ = edge_priors(evidence, read_edge_tiers(EDGE_PRIOR_TIERS))
    return {"gated": int(gated.sum()), "edges": int(len(evidence))}


def figure_model_pipeline() -> None:
    counts = graph_layer_counts(CONFIRMATORY_GRAPH_DIRECTORY)
    binder_summary = json.loads(BINDER_SUMMARY.read_text())
    entry_summary = json.loads((CONFIRMATORY_GRAPH_DIRECTORY / "drug_entry_summary.json").read_text())
    label_summary = json.loads(LABEL_SELECTION_SUMMARY.read_text())
    gate_counts = edge_prior_counts(CONFIRMATORY_GRAPH_DIRECTORY)
    cell_classes = [column for column in pd.read_parquet(CELL_CLASS_WEIGHTS).columns if column != "all_cells"]
    layer_lines = "\n".join(f"{layer}: {number:,}" for layer, number in sorted(counts["edges_by_layer"].items(), key=lambda item: -item[1]))
    node_lines = ", ".join(f"{number:,} {node_type.replace('_', ' ')}" for node_type, number in counts["node_types"].items())
    # wrap the node counts at about 52 characters so the graph box keeps the width of its layer lines
    wrapped_node_lines, line = [], ""
    for part in node_lines.split(", "):
        candidate = f"{line}, {part}" if line else part
        if len(candidate) > 52 and line:
            wrapped_node_lines.append(line + ",")
            line = part
        else:
            line = candidate
    wrapped_node_lines.append(line)
    # the axis is drawn in inches, so a measured text size is directly a coordinate and the figure can be sized to
    # its contents; the size passed here is provisional and is replaced once the boxes have been measured
    figure = plt.figure(figsize=(16, 10))
    axis = figure.add_axes((0.0, 0.0, 1.0, 1.0))
    axis.axis("off")
    first_row = [
        measured_box(figure, axis, "1. Perturbation",
                     "a gene's allele (decision 9.54): a loss of\nfunction or an untyped allele seeds -1 at\nthe gene, a hypermorphic gain +1;\n"
                     f"{len(entry_summary.get('alleles', {}))} neomorph, HTT's CAG expansion, is its\nown node with "
                     f"{entry_summary.get('allele_mechanism_edges_added', 0)} cited edges\n\n"
                     f"a drug: its own node ({entry_summary['drug_nodes_added']}), joined to\nits targets by "
                     f"{entry_summary['mechanism_edges_added']:,} ChEMBL mechanism\nedges (ablation: seeded at the targets)\n\n"
                     "a sustained input u at the perturbed\nnodes; knockouts seed gene nodes,\nmechanism edges end at proteins (figure 7)",
                     BOX_COLOURS["input"]),
        measured_box(figure, axis, "2. Graph (the v5 drug and allele graph)",
                     f"{counts['nodes']:,} nodes, {counts['edges']:,} edges, {counts['relations']} relations\n" + "\n".join(wrapped_node_lines) + "\n\n"
                     f"edges by layer:\n{layer_lines}\n\n"
                     "signed edges: activation +1, inhibition and\nrepression -1, substrate depletion -1, binding 0\n\n"
                     f"edge priors (9.57): {gate_counts['gated']:,} edges carry a prior\nthat they are real; off in the registered\n"
                     "model; hard gates, a binary second pass and\nexact flip checks built (9.62; figures 13-16)",
                     BOX_COLOURS["graph"]),
        measured_box(figure, axis, "3. Encoder (one per tested model)",
                     "message passing: 3 layers of typed messages,\nrelations shared through a rank-2 factorisation,\n"
                     "edge signs from the graph; node states from\nstructural features and descriptors\n\n"
                     "linear response, to a fixed point (tolerance\n0.001, at most 60 rounds):\n"
                     "h(t+1) = (1 - d) h(t) + d (sum_r g_r S_r h(t) + u)\n"
                     "S_r: signed adjacency of relation r, normalised\nby total in-degree; g_r: learned gain per relation\n"
                     "with gates, both encoders read a closed gate as a\ndeleted edge (figure 14)\n"
                     "positive input weights, cofactor relations,\nreversible reactions, compartment pH\n"
                     f"class k reads b_k h*: b_k, the weakest node weight\non the widest route (the limiting step;\n"
                     f"{len(cell_classes)} classes and all cells); extracellular\npools shared across classes (figure 2)\n\n"
                     "both read the difference field:\nperturbed minus unperturbed", BOX_COLOURS["propagation"]),
    ]
    second_row = [
        measured_box(figure, axis, "4. Node properties",
                     "fixed, per node type: protein (ESM-2 embeddings),\nmetabolite and reaction (EC) descriptors; brain\n"
                     "region and cell-class expression, with the\ndopaminergic class (figure 4)\n\n"
                     "message passing: in the node's starting state\nlinear response: in the output gate,\n"
                     "field = (h W) * sigmoid(node properties W_g + b),\nwhich changes where the response is read,\n"
                     "not how it propagates\n\n"
                     "seed masked (registered): a perturbation's own\nnodes, their encodes partners and a drug node's\n"
                     "targets read only structural and expression\ncolumns, so no readout learns which protein it is",
                     BOX_COLOURS["readout"]),
        measured_box(figure, axis, "5. Pooling and head",
                     "sum over nodes\n\nnoisy-OR over pathway modules\n(gated node supports and links);\n"
                     "each leak starts at its\nloss-optimal constant; links,\nleaks and gates train at\ntheir own learning rates",
                     BOX_COLOURS["readout"]),
        measured_box(figure, axis, "6. Symptoms and labels",
                     f"{len(label_summary['all_by_symptom'])} fine symptom constructs (decision 9.58)\n\n"
                     "labels: HPO annotations (a gene's alleles),\nSIDER and OnSIDES drug labels (drugs); trial\n"
                     "arms, FAERS and IMPC knockouts enter the\nreliability fit as soft classes, never as labels\n"
                     "better_v2 (figure 3): low-frequency positives,\ngrade C pairs and drug pairs never above\n"
                     "placebo are masked, neither positive nor\nnegative; unobserved pairs stay unlabelled\n\n"
                     "folds grouped by disease cluster and drug\ntargets; the held-out set on v5 is drawn on\n"
                     "the owner's go-ahead (1.1); refit on\ntraining and validation together", BOX_COLOURS["output"]),
    ]
    title = "The registered model (the v4 flags, a noisy-OR head on each encoder) on the v5 tables, from a perturbation to symptom probabilities"
    title_width, title_height = text_size_in_inches(figure, axis, title, 14.0, "bold")
    rows = (first_row, second_row)
    row_heights = [max(box["height"] for box in row) for row in rows]
    figure_width = max([title_width] + [sum(box["width"] for box in row) + (len(row) - 1) * SMALLEST_BOX_GAP for row in rows]) + 2 * FIGURE_MARGIN
    figure_height = 2 * FIGURE_MARGIN + title_height + 0.35 + row_heights[0] + ROW_GUTTER + row_heights[1]
    figure.set_size_inches(figure_width, figure_height)
    axis.set_xlim(0, figure_width)
    axis.set_ylim(0, figure_height)
    axis.text(figure_width / 2, figure_height - FIGURE_MARGIN, title, ha="center", va="top", fontsize=14, fontweight="bold")

    row_tops = [figure_height - FIGURE_MARGIN - title_height - 0.35]
    row_tops.append(row_tops[0] - row_heights[0] - ROW_GUTTER)
    spans_by_row = []
    for row, row_top, row_height in zip(rows, row_tops, row_heights):
        # the slack goes into the gaps, so both rows share one left edge and one right edge
        gap = (figure_width - 2 * FIGURE_MARGIN - sum(box["width"] for box in row)) / (len(row) - 1)
        left, spans = FIGURE_MARGIN, []
        for box in row:
            draw_box(axis, box, left, row_top, row_height)
            spans.append((left, left + box["width"]))
            left += box["width"] + gap
        spans_by_row.append(spans)

    for spans, row_top, row_height in zip(spans_by_row, row_tops, row_heights):
        middle_y = row_top - row_height / 2
        for (_, right_edge), (next_left_edge, _) in zip(spans, spans[1:]):
            draw_arrow(axis, [(right_edge + 0.08, middle_y), (next_left_edge - 0.08, middle_y)])
    encoder_centre = sum(spans_by_row[0][-1]) / 2
    properties_centre = sum(spans_by_row[1][0]) / 2
    gutter_y = row_tops[1] + ROW_GUTTER / 2
    draw_arrow(axis, [(encoder_centre, row_tops[0] - row_heights[0] - 0.08), (encoder_centre, gutter_y),
                      (properties_centre, gutter_y), (properties_centre, row_tops[1] + 0.08)])
    save(figure, "model_pipeline")


def toy_chain_encoder(cell_class_weights: torch.Tensor, extracellular_pool: torch.Tensor | None, class_weight_rule: str) -> LinearResponseEncoder:
    """The six-node chain of tests/test_linear_response_encoder.py with every gain and input weight equal, so the
    channels differ only by their cell-class weights. Run to the registered fixed point (tolerance 0.001)."""
    relation_types = ["substrate_of", "product_of", "catalyzed_by"]
    edges = [(0, 1, "catalyzed_by", 1.0), (2, 1, "substrate_of", 1.0), (1, 3, "product_of", 1.0), (3, 4, "substrate_of", 1.0), (4, 5, "product_of", 1.0)]
    node_types = torch.tensor([0, 1, 2, 2, 1, 2])
    encoder = LinearResponseEncoder(
        6, relation_types, edge_source=torch.tensor([edge[0] for edge in edges]), edge_target=torch.tensor([edge[1] for edge in edges]),
        edge_relation=torch.tensor([relation_types.index(edge[2]) for edge in edges]), edge_sign=torch.tensor([edge[3] for edge in edges]),
        node_features=torch.nn.functional.one_hot(node_types, 3).float(), node_state_dim=4, class_weight_rule=class_weight_rule,
        cell_class_weights=cell_class_weights, channels_per_cell_class=1, extracellular_pool_nodes=extracellular_pool)
    with torch.no_grad():
        encoder.gain_logit.zero_()
        encoder.log_input_weight.fill_(0.0)  # an input weight of exp(0) = 1 in every channel
    return encoder


def figure_cell_class_channels() -> None:
    """The registered rule (limiting_step, the owner on 13 October 2026: "min along the path makes WAY more sense to me")
    beside the rule every run before it used (product), on the chain of the unit tests, computed by the encoder."""
    node_names = ["gene G", "reaction R1", "substrate S", "product P\n(extracellular)", "reaction R2", "product Q"]
    weights = torch.ones(6, 2)
    weights[1, 1] = 0.0  # class B does not express reaction R1
    weights[4, :] = 0.5  # R2 and Q at half their highest class, in both classes: two steps at 0.5 on one route
    weights[5, :] = 0.5
    pool = torch.tensor([False, False, False, True, False, False])
    loss_of_function = (torch.tensor([[0]]), torch.tensor([[[-1.0, 1.0]]]))
    with torch.no_grad():
        product = toy_chain_encoder(weights, None, "product").response(*loss_of_function)[0]
        limiting = toy_chain_encoder(weights, None, "limiting_step").response(*loss_of_function)[0]
        limiting_pooled = toy_chain_encoder(weights, pool, "limiting_step").response(*loss_of_function)[0]
    columns = {"product\nclass A": product[:, 0], "product\nclass B": product[:, 1],
               "limiting step\nclass A": limiting[:, 0], "limiting step\nclass B": limiting[:, 1],
               "limiting step,\nP shared, A": limiting_pooled[:, 0], "limiting step,\nP shared, B": limiting_pooled[:, 1]}
    values = torch.stack(list(columns.values()), dim=1).numpy()

    figure, (chain_axis, heat_axis) = plt.subplots(1, 2, figsize=(15.5, 5.0), gridspec_kw={"width_ratios": [1.0, 1.25]})
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
    chain_axis.text(0.9, 0.62, "R2 and Q: weight 0.5\nin both classes", color="#c0392b", fontsize=8, ha="center")
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
    heat_axis.set_xticks(range(len(columns)), list(columns), fontsize=7.5)
    heat_axis.set_yticks(range(len(node_names)), [name.replace("\n", " ") for name in node_names], fontsize=8)
    for boundary in (1.5, 3.5):
        heat_axis.axvline(boundary, color="#333333", linewidth=1.2)
    heat_axis.set_title("Response at the fixed point (untrained, equal gains); the registered rule is the limiting step", fontsize=9.5)
    colour_bar = figure.colorbar(image, ax=heat_axis, fraction=0.04)
    colour_bar.set_label("sign * log10(1 + |h| / 1e-6)", fontsize=8)
    figure.suptitle("Cell-class channels: one adjacency, each class's response limited by its expression", fontsize=11, fontweight="bold")
    figure.text(0.5, -0.07, "Product (every run before 13 October 2026): a node's weight multiplies what arrives at it on every round, so a route is charged "
                "the product of its weights, and Q, behind R2 at 0.5 and itself at 0.5,\nis cut by both. Limiting step (registered): the propagation runs unweighted "
                "to its fixed point and each node's state is then scaled once by b, the weakest weight on its widest route from the seed, so Q is limited\nby 0.5 "
                "once, and class B, which lacks R1, gets nothing past it. With P a shared extracellular pool, b at P is the largest over the pooled classes, so "
                "class B's R2 and Q receive what class A makes.\nComputed with LinearResponseEncoder (class_weight_rule, cell_class_weights, extracellular_pool_nodes) "
                "on the six-node chain of the unit tests.", ha="center", fontsize=8)
    save(figure, "cell_class_channels")


def figure_label_selection() -> None:
    summary = json.loads(LABEL_SELECTION_SUMMARY.read_text())
    symptoms = sorted(summary["all_by_symptom"], key=lambda symptom: summary["all_by_symptom"][symptom])
    kept = np.array([summary["kept_by_symptom"].get(symptom, 0) for symptom in symptoms])
    set_aside = np.array([summary["all_by_symptom"][symptom] for symptom in symptoms]) - kept
    figure, (symptom_axis, reason_axis) = plt.subplots(1, 2, figsize=(14, 0.2 * len(symptoms) + 2.6), gridspec_kw={"width_ratios": [1.3, 1]})
    labels = [symptom.replace("_", " ") for symptom in symptoms]
    symptom_axis.barh(labels, kept, color="#2b6ca3", label="kept")
    symptom_axis.barh(labels, set_aside, left=kept, color="#bcd3e8", label="set aside (masked)")
    for row, (kept_count, set_aside_count) in enumerate(zip(kept, set_aside)):
        symptom_axis.text(kept_count + set_aside_count + 5, row, f"{kept_count} of {kept_count + set_aside_count}", va="center", fontsize=8)
    symptom_axis.set_xlabel(f"positive (perturbation, symptom) pairs, the v5 tables ({EVIDENCE_DIRECTORY.name})")
    symptom_axis.tick_params(axis="y", labelsize=7.5)
    symptom_axis.set_xlim(0, max(kept + set_aside) * 1.2)
    symptom_axis.legend(loc="lower right", fontsize=8)
    symptom_axis.set_title(f"{summary['kept_pairs']:,} of {summary['positive_pairs']:,} positive pairs kept", fontsize=10)
    reasons = summary["set_aside_reasons"]
    reason_labels = [reason.replace(" x ", " ").replace("frequency below", "frequency\nbelow").replace(", one label source", ",\none label source")
                     .replace("never more frequent on the drug than on placebo", "never above\nplacebo").replace("(1 label entries with both arms)", "(1 label with both arms)").replace(" label entries with both arms", " labels with both arms")
                     for reason in reasons]
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
    figure.suptitle("Better (not more) training examples: selection better_v2 on the v5 tables (fine constructs, alleles as perturbations)",
                    fontsize=11, fontweight="bold")
    placebo = summary.get("mask_never_above_placebo") or {}
    figure.tight_layout(rect=(0, 0.07, 1, 0.97))
    figure.text(0.5, 0.055, f"Kept: gene pairs at HPO frequency {summary['gene_minimum_frequency']:.2f} or more (Frequent and above); drug pairs at label "
                f"frequency {summary['drug_minimum_frequency']:.0%} or more, or listed by both SIDER and OnSIDES. A pair set aside stays out of the loss "
                f"and every metric;\nit is not a negative. Also masked: {summary['masked_negative_pairs']:,} pairs whose only evidence is grade C (a human "
                f"association), and {placebo.get('kept_pairs_set_aside', 0)} kept drug pairs never more frequent on the drug than on placebo in a label "
                f"that gives both arms.\nToo few kept pairs for the macro average ({MINIMUM_POSITIVES_FOR_MACRO} needed): "
                f"{', '.join(too_few) if too_few else 'none'}.", ha="center", va="top", fontsize=8)
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
    """nodes: name -> (x, y, label, colour); edges: (source, target, label). Arrows are clipped to the text boxes.

    Returns the text object of each node, so a caller can draw further arrows clipped to the same boxes."""
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
    return text_of_node


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


def currency_state_counts(graph_directory: Path) -> dict:
    """Per group of currency species, the copies and the total degree: what the graph computes and does not pass on."""
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    nodes = nodes.assign(currency=currency_flags(nodes),
                         species=nodes.display_name.astype(str).str.strip())
    counts = {}
    for group, species_names in CELL_STATE_GROUPS.items():
        rows = nodes[nodes.species.isin(species_names)]
        counts[group] = {"copies": int(len(rows)),
                         "currency_copies": int(rows.currency.sum()),
                         "degree": int(rows.degree.sum()),
                         "currency_degree": int(rows[rows.currency].degree.sum())}
    neuronal = nodes[nodes.node_id.astype(str).str.startswith(NEURONAL_POOL_PREFIX)]
    counts["neuronal ion pools"] = {"copies": int(len(neuronal)), "currency_copies": int(neuronal.currency.sum()),
                                    "degree": int(neuronal.degree.sum()), "currency_degree": 0}
    return counts


def figure_cell_state() -> None:
    """What the model holds as the state of a cell, and which parts of it travel.

    The distinction the figure exists to make (the owner, 13 October 2026: "please include an explanatory figure for
    things like representation of the ph and other pieces of the cell state"): a node's state either propagates along
    edges or it caps what passes through a node, and the second is how everything that would otherwise be a hub
    enters. The counts are read from the graph rather than written here."""
    counts = currency_state_counts(CELL_STATE_GRAPH_DIRECTORY)
    figure, axes = plt.subplots(1, 2, figsize=(14.6, 6.6), gridspec_kw={"width_ratios": [1.1, 1.0]})
    axes[0].set_axis_off()  # only the diagram; the bar panel needs its axis for the group names
    for spine in ("top", "right"):
        axes[1].spines[spine].set_visible(False)

    travels, caps = "#e3f1df", "#fdebd3"
    nodes = {
        "seed": (0.5, 0.92, "perturbation\n(gene lowered, drug at its node)", BOX_COLOURS["input"]),
        "graph": (0.5, 0.66, "the response travels the graph\nto a fixed point: reactions, proteins,\nmetabolites, transcription", travels),
        "ion": (0.13, 0.40, "neuronal ion pools\n(Na+, K+, Cl-, Ca2+)", travels),
        "vm": (0.13, 0.17, "membrane potential\nV$_m$, computed", travels),
        "class": (0.84, 0.845, "cell-class expression\nw$_k$(node)", caps),
        "bottleneck": (0.84, 0.60, "limiting step: the weakest\nw$_k$ on the widest route", caps),
        "ph": (0.60, 0.30, "compartment pH\nfrom the proton pools", caps),
        "ratios": (0.86, 0.30, "energy charge, NAD+/NADH,\nNADPH/NADP+, GSH/GSSG, SAM/SAH\n(proposed, not built)", caps),
        "readout": (0.5, 0.045, "symptom probabilities", BOX_COLOURS["readout"]),
    }
    edges = [("seed", "graph", ""), ("graph", "ion", ""), ("ion", "vm", ""), ("graph", "readout", ""), ("vm", "readout", "")]
    text_of_node = draw_node_link_diagram(axes[0], nodes, edges)
    for capping, target in (("class", "graph"), ("bottleneck", "graph"), ("ph", "readout"), ("ratios", "readout")):
        axes[0].annotate("", xy=(0.5, 0.5), xycoords=text_of_node[target], xytext=(0.5, 0.5), textcoords=text_of_node[capping],
                         arrowprops={"arrowstyle": "-|>", "color": "#b36b00", "linestyle": (0, (3, 2)),
                                     "patchA": text_of_node[capping].get_bbox_patch(),
                                     "patchB": text_of_node[target].get_bbox_patch(), "shrinkA": 2, "shrinkB": 2})
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].set_title("What the state is, and what travels", fontsize=12, fontweight="bold")
    axes[0].text(0.5, 0.985, "solid: the state propagates along edges     dashed amber: the state caps throughput and does not propagate",
                 ha="center", va="top", fontsize=7.6, color="#555555")

    labels = list(counts)
    computed = [counts[label]["currency_degree"] for label in labels]
    passed = [counts[label]["degree"] - counts[label]["currency_degree"] for label in labels]
    positions = np.arange(len(labels))
    axes[1].barh(positions, computed, color="#fdebd3", edgecolor="#b36b00", label="reaction edges whose effect is computed and not passed on")
    axes[1].barh(positions, passed, left=computed, color="#e3f1df", edgecolor="#2f6b3a", label="reaction edges whose effect propagates")
    axes[1].set_yticks(positions)
    axes[1].set_yticklabels(labels, fontsize=8.5)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("edges at the nodes of that state, summed over compartmental copies", fontsize=9)
    axes[1].set_title("The state the graph computes and discards", fontsize=12, fontweight="bold")
    axes[1].legend(fontsize=7.6, loc="lower right", framealpha=0.95)
    for position, label in zip(positions, labels):
        total = counts[label]["degree"]
        axes[1].text(total + max(1, total * 0.015), position,
                     f"{counts[label]['currency_copies']}/{counts[label]['copies']} copies held back", va="center", fontsize=7.2, color="#555555")
    axes[1].set_xlim(0, max(counts[label]["degree"] for label in labels) * 1.46)
    figure.suptitle("The cell state the model carries", fontsize=14, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    save(figure, "cell_state")


def figure_graph_composition() -> None:
    """What the graph holds: its node types, its relations with their signs, and the source layers of its edges."""
    nodes = pd.read_parquet(CONFIRMATORY_GRAPH_DIRECTORY / "nodes.parquet")
    edges = pd.read_parquet(CONFIRMATORY_GRAPH_DIRECTORY / "edges.parquet")
    sign = edges.sign.astype(float)
    sign_classes = (("+1, activates or increases", "#2b6ca3", sign > 0), ("-1, inhibits, represses or depletes", "#c0392b", sign < 0),
                    ("0, binds or carries", "#9aa5b1", sign == 0))
    figure, (node_axis, relation_axis, layer_axis) = plt.subplots(1, 3, figsize=(18, 7.4), gridspec_kw={"width_ratios": [0.8, 1.15, 1.15]})
    node_counts = nodes.node_type.value_counts().sort_values()
    node_axis.barh([name.replace("_", " ") for name in node_counts.index], node_counts.to_numpy(), color="#7a9cc6")
    node_axis.set_xscale("log")
    for position, number in enumerate(node_counts.to_numpy()):
        node_axis.text(number * 1.15, position, f"{number:,}", va="center", fontsize=8)
    node_axis.set_xlim(0.8, node_counts.max() * 8)
    node_axis.set_xlabel("nodes (log scale)")
    node_axis.set_title(f"{len(nodes):,} nodes by type", fontsize=11, fontweight="bold")

    for axis, key, title in ((relation_axis, edges.relation_type, f"{len(edges):,} edges by relation ({edges.relation_type.nunique()} relations) and sign"),
                             (layer_axis, pd.Series(edge_layers(edges), index=edges.index), "the same edges by the source layer they come from")):
        order = key.value_counts().sort_values().index
        left = np.zeros(len(order))
        for label, colour, mask in sign_classes:
            counts = key[mask].value_counts().reindex(order, fill_value=0).to_numpy()
            axis.barh([name.replace("_", " ") for name in order], counts, left=left, color=colour, label=label)
            left += counts
        for position, total in enumerate(left):
            axis.text(total * 1.02 + 200, position, f"{int(total):,}", va="center", fontsize=7.5)
        axis.set_xlim(0, left.max() * 1.18)
        axis.tick_params(axis="y", labelsize=8)
        axis.set_xlabel("edges")
        axis.set_title(title, fontsize=11, fontweight="bold")
    relation_axis.legend(fontsize=8, loc="lower right")
    for axis in (node_axis, relation_axis, layer_axis):
        axis.spines[["top", "right"]].set_visible(False)
    figure.suptitle(f"The physiology graph the model reads ({CONFIRMATORY_GRAPH_DIRECTORY.name}), read from its tables when the figure is drawn",
                    fontsize=12.5, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    save(figure, "graph_composition")


MECHANISM_BASIS_LABELS = {
    "orphanet": "Orphanet's own type",
    "same": "the same disease's Orphanet record (exact)",
    "narrower": "a narrower Orphanet form (NTBT)",
    "gene2phenotype": "Gene2Phenotype, the second source",
    "adjudicated": "adjudicated from the primary paper",
    "conflict": "sources disagree: untyped",
    "untyped": "no source types it: untyped",
}
ALLELE_CLASS_COLOURS = {"loss_of_function": ("loss of function, -1 at the gene", "#2b6ca3"), "hypermorph": ("hypermorph, +1 at the gene", "#c0392b"),
                        "neomorph": ("neomorph, its own node", "#8e44ad"), "untyped": ("untyped, -1 at the gene by default", "#9aa5b1")}


def figure_allele_perturbations() -> None:
    reports = pd.read_parquet(EVIDENCE_DIRECTORY / "evidence_reports.parquet")
    monogenic = reports[reports.allele_class.fillna("") != ""]
    pairs = monogenic.assign(gene=monogenic.perturbation_id.str.split("|").str[0]).drop_duplicates(["gene", "source_record_id"])
    basis = pairs.allele_mechanism_basis.fillna("").map(lambda text: "conflict" if text.startswith("conflict") else text.split(" ")[0].rstrip(";"))
    summary = json.loads((EVIDENCE_DIRECTORY / "evidence_summary.json").read_text())["alleles"]
    entry_edges = pd.read_parquet(CONFIRMATORY_GRAPH_DIRECTORY / "edges.parquet")
    allele_edges = entry_edges[entry_edges.source_id.str.startswith("ALLELE:")]

    figure, (diagram_axis, source_axis) = plt.subplots(1, 2, figsize=(17, 7.2), gridspec_kw={"width_ratios": [1.15, 1.0]})
    figure.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.12, wspace=0.42)
    diagram_axis.set_xlim(0, 1)
    diagram_axis.set_ylim(0, 1)
    diagram_axis.axis("off")
    disease_colour, perturbation_colour, gene_colour, protein_colour, allele_colour = "#f4f4f4", "#fdebd3", "#dbe9f6", "#e3f1df", "#efe3f4"
    nodes = {
        "recessive": (0.14, 0.9, "Lopes-Maciel-Rodan syndrome\nOMIM:617435, recessive", disease_colour),
        "huntington": (0.62, 0.9, "Huntington disease\nORPHA:399, OMIM:143100,\njuvenile form ORPHA:248111", disease_colour),
        "loss": (0.14, 0.66, "perturbation HTT\n(loss of function or untyped)", perturbation_colour),
        "expansion": (0.62, 0.66, "perturbation\nHTT|CAG_repeat_expansion (neomorph)", perturbation_colour),
        "gene": (0.14, 0.4, "GENE:HTT\nseeded -1", gene_colour),
        "allele": (0.62, 0.42, "ALLELE:HTT_CAG_repeat_expansion\nseeded +1, present only when seeded", allele_colour),
    }
    edges = [("recessive", "loss", ""), ("huntington", "expansion", ""), ("loss", "gene", ""), ("expansion", "allele", "")]
    target_positions = [(0.32, 0.13), (0.55, 0.13), (0.78, 0.13), (0.98, 0.13)]
    for (row, position) in zip(allele_edges.sort_values("target_id").itertuples(index=False), target_positions):
        name = f"target_{row.target_id}"
        pmid = row.evidence_source.split("PMID ")[-1].rstrip(")")
        nodes[name] = (position[0] - 0.06, position[1], f"{row.target_id}\nPMID {pmid}", gene_colour if row.target_id.startswith("GENE:") else protein_colour)
        edges.append(("allele", name, f"{'+' if row.sign > 0 else '-'}1"))
    draw_node_link_diagram(diagram_axis, nodes, edges)
    diagram_axis.text(0.38, 0.535, "one leakage group, HTT: no split puts\none allele in training and the other in the test",
                      ha="center", va="center", fontsize=8, color="#555555", style="italic")
    diagram_axis.set_title("HTT's two alleles are two perturbations", fontsize=11, fontweight="bold")

    order = [key for key in MECHANISM_BASIS_LABELS if key in set(basis)]
    left = np.zeros(len(order))
    for allele_class, (label, colour) in ALLELE_CLASS_COLOURS.items():
        counts = np.array([int(((basis == key) & (pairs.allele_class == allele_class)).sum()) for key in order])
        source_axis.barh([MECHANISM_BASIS_LABELS[key] for key in order], counts, left=left, color=colour, label=label)
        left += counts
    for position, total in enumerate(left):
        source_axis.text(total + 25, position, f"{int(total):,}", va="center", fontsize=8)
    source_axis.invert_yaxis()
    source_axis.set_xlim(0, left.max() * 1.15)
    source_axis.set_xlabel("(gene, disease) pairs behind the monogenic reports, v5")
    source_axis.legend(fontsize=8, loc="upper right")
    source_axis.spines[["top", "right"]].set_visible(False)
    classes = summary["alleles_by_class"]
    source_axis.set_title(f"Where each disease's mechanism came from, in order: {classes.get('loss_of_function', 0):,} loss, "
                          f"{classes.get('hypermorph', 0)} hypermorph,\n{classes.get('neomorph', 0)} neomorph and {classes.get('untyped', 0):,} untyped alleles; "
                          f"{summary['genes_with_several_alleles']} genes with several", fontsize=10.5, fontweight="bold")
    figure.suptitle("The allele, not the gene, is the perturbation (decision 9.54)", fontsize=13, fontweight="bold")
    figure.text(0.5, 0.02, "Each step is consulted only where the ones before it typed nothing. A dominant-negative product counts as a loss of function, "
                "because the gene's net activity falls. Untyped pairs keep -1, which is expected for recessive diseases\nand wrong for an unknown share "
                "of dominant ones (docs/allele_perturbations.md, section 6). The allele's four edges each carry a PubMed id and a quote "
                "(configs/neomorphic_alleles.csv).", ha="center", fontsize=8)
    save(figure, "allele_perturbations")


def figure_edge_priors() -> None:
    """The prior that an edge is real, by its evidence (decision 9.57), in the binary form decided on 16 October 2026
    (decision 9.62: pi = sigmoid(n logit p), the cost of closing an edge n logit p nats), and how a hard concrete gate
    started at that prior reads, which is why the decided build ends with a binary pass."""
    from experiments.measure_edge_gate_readings import binary_priors, hard_concrete_masses, log_alpha_for_open_probability
    from mechanistic_pathway_learning.graph.edge_evidence import edge_priors, read_edge_tiers

    evidence = pd.read_parquet(CONFIRMATORY_GRAPH_DIRECTORY / "edge_evidence.parquet")
    gated, probability, pseudo_count = edge_priors(evidence, read_edge_tiers(EDGE_PRIOR_TIERS))
    gated_evidence = evidence[gated].assign(p=probability[gated], n=pseudo_count[gated])
    gated_evidence["prior"] = binary_priors(gated_evidence.p.to_numpy(), gated_evidence.n.to_numpy())["log_odds_times_n"]
    gated_evidence["tier"] = np.where(gated_evidence.contested, "lowered by an IntAct absence",
                                      gated_evidence.family.str.replace("_", " ") + ", " + gated_evidence.kind.str.replace("_", " "))
    tiers = gated_evidence.groupby("tier").agg(edges=("prior", "size"), low=("prior", "min"), high=("prior", "max")).sort_values("edges")

    figure, (tier_axis, gate_axis) = plt.subplots(1, 2, figsize=(16.5, 6.6), gridspec_kw={"width_ratios": [1.05, 1.0]})
    positions = np.arange(len(tiers))
    tier_axis.barh(positions, tiers.edges.to_numpy(), color="#7a9cc6")
    tier_axis.set_xscale("log")
    tier_axis.set_yticks(positions, tiers.index, fontsize=8.5)
    for position, (tier, row) in zip(positions, tiers.iterrows()):
        cost_low, cost_high = np.log(row.low / (1 - row.low)), np.log(row.high / (1 - row.high))
        prior_text = f"{row.low:.2f}" if abs(row.high - row.low) < 0.005 else f"{row.low:.2f}-{row.high:.2f}"
        cost_text = f"{cost_low:+.2f}" if abs(cost_high - cost_low) < 0.005 else f"{cost_low:+.2f} to {cost_high:+.2f}"
        tier_axis.text(row.edges * 1.12, position, f"{int(row.edges):,}; prior {prior_text}; closing {cost_text} nats", va="center", fontsize=7.6)
    tier_axis.set_xlim(1, tiers.edges.max() * 90)
    tier_axis.set_xlabel("gated edges (log scale)")
    tier_axis.spines[["top", "right"]].set_visible(False)
    tier_axis.set_title(f"{int(gated.sum()):,} of {len(evidence):,} edges carry a prior that they are real\n"
                        "prior = sigmoid(n logit p): p from the evidence tier, n from its tempered record count", fontsize=10.5, fontweight="bold")

    prior_grid = np.linspace(0.5, 0.995, 200)
    masses = hard_concrete_masses(log_alpha_for_open_probability(prior_grid))
    gate_axis.stackplot(prior_grid, masses["closed"], masses["fractional"], masses["open"], colors=["#d5dbe1", "#fdebd3", "#2b6ca3"],
                        labels=["gate exactly 0 (closed)", "gate between 0 and 1", "gate exactly 1 (open)"], alpha=0.95)
    gate_axis.plot(prior_grid, masses["test_time_gate"], color="#c0392b", linewidth=2, label="the paper's test-time gate (eq. 13)")
    gate_axis.plot(prior_grid, np.ones_like(prior_grid), color="#111111", linewidth=2, linestyle="--", label="the decided read-out: open, as P(open) is at least 1/2")
    gate_axis.set_xlim(0.5, 0.995)
    gate_axis.set_ylim(0, 1.05)
    gate_axis.set_xlabel("the edge's prior, where its hard concrete gate starts (P(open) = prior)")
    gate_axis.set_ylabel("share of sampled gates; gate value")
    gate_axis.legend(fontsize=7.8, loc="center left", bbox_to_anchor=(0.0, 0.32), framealpha=0.95)
    gate_axis.spines[["top", "right"]].set_visible(False)
    gate_axis.set_title("A hard concrete gate at its prior is not hard: most samples are fractional\n"
                        "(Louizos et al. 2018: gamma -0.1, zeta 1.1, beta 2/3), so the build ends with a binary pass", fontsize=10.5, fontweight="bold")
    figure.suptitle("Edge priors (decision 9.57) and the hard gates decided and built on 16 October 2026 (decision 9.62; figures 13 to 16)", fontsize=12.5, fontweight="bold")
    figure.text(0.5, 0.055, "Closing an edge costs logit(prior) = n logit p nats of fit; with the gate's entropy, the per-edge cost is the KL divergence to "
                "Bernoulli(prior). A positive closing cost keeps open an edge the data say nothing about;\nonly the 15 edges an IntAct absence lowered below "
                "one half are cheaper closed. The binary second pass trains the gates alone with the read-out in the forward pass and a surrogate gradient, a "
                "steep sigmoid at the threshold,\nin the backward pass. Numbers from experiments/measure_edge_gate_readings.py.", ha="center", va="top", fontsize=8)
    figure.tight_layout(rect=(0, 0.075, 1, 0.94))
    save(figure, "edge_priors")


def framed_row(axis, boxes: list[tuple[float, str, str, str]], top: float, title_size: float = 9.6, body_size: float = 7.9) -> list[tuple[float, float, float, float]]:
    """A row of framed boxes in axis coordinates, each a bold title over a body, centred at its x, sharing one top edge
    and one bottom edge (the lowest the texts need), so the row reads as one stage after another. Returns each frame's
    (left, bottom, right, top). The frames are measured from the drawn text, so no frame is narrower than its lines."""
    renderer = axis.figure.canvas.get_renderer()
    to_axis = axis.transAxes.inverted()
    pad_x, pad_y, title_gap = 0.007, 0.035, 0.025
    extents = []
    for centre_x, title, body, _ in boxes:
        title_text = axis.text(centre_x, top - pad_y, title, transform=axis.transAxes, ha="center", va="top", fontsize=title_size, fontweight="bold", zorder=3)
        title_extent = title_text.get_window_extent(renderer).transformed(to_axis)
        body_text = axis.text(centre_x, title_extent.y0 - title_gap, body, transform=axis.transAxes, ha="center", va="top", fontsize=body_size,
                              linespacing=1.38, zorder=3)
        body_extent = body_text.get_window_extent(renderer).transformed(to_axis)
        extents.append((min(title_extent.x0, body_extent.x0) - pad_x, body_extent.y0 - pad_y, max(title_extent.x1, body_extent.x1) + pad_x))
    bottom = min(extent[1] for extent in extents)
    frames = []
    for (left, _, right), (_, _, _, colour) in zip(extents, boxes):
        axis.add_patch(FancyBboxPatch((left, bottom), right - left, top - bottom, transform=axis.transAxes, boxstyle="round,pad=0,rounding_size=0.012",
                                      facecolor=colour, edgecolor="#555555", linewidth=1.1, zorder=1, clip_on=False))
        frames.append((left, bottom, right, top))
    for (_, low, right, high), (next_left, _, _, _) in zip(frames, frames[1:]):
        middle = (low + high) / 2
        axis.add_patch(FancyArrowPatch((right + 0.003, middle), (next_left - 0.003, middle), transform=axis.transAxes, arrowstyle="-|>", mutation_scale=16,
                                       color="#444444", linewidth=1.4, zorder=3))
    return frames


def figure_edge_gate_training() -> None:
    """How the hard edge gates of decision 9.62 (e) are fitted (models/edge_gates.py, models/edge_gate_passes.py): the
    stages in order; the per-edge cost the prior sets in each pass; and the proxy gradient of the binary pass, drawn
    from the functions the code uses (binary_edge_prior, HardConcreteEdgeGates.apply_flip_gradient's sigmoid)."""
    from mechanistic_pathway_learning.models.edge_gate_passes import DEFAULT_SURROGATE_WIDTH
    from mechanistic_pathway_learning.models.edge_gates import HARD_CONCRETE_BETA, HARD_CONCRETE_GAMMA, HARD_CONCRETE_ZETA, binary_edge_prior

    figure = plt.figure(figsize=(18.5, 8.8))
    grid = figure.add_gridspec(2, 2, height_ratios=[0.44, 1.0], hspace=0.24, wspace=0.18)
    flow_axis = figure.add_subplot(grid[0, :])
    flow_axis.axis("off")
    stages = [
        (0.095, "1. The prior that an edge is real",
         "p from the edge's evidence tier and n from\nthe tier's tempered record count (9.57):\n"
         "π = sigmoid(n · logit p), 9.62 sub-choice (iii)\n\nclosing the edge costs logit π nats of fit;\n"
         "reopening a closed edge pays them back", BOX_COLOURS["input"]),
        (0.297, "2. Stochastic pass, with the whole model",
         f"one hard concrete sample per gate and\nminibatch (γ {HARD_CONCRETE_GAMMA:g}, ζ {HARD_CONCRETE_ZETA:g}, β {HARD_CONCRETE_BETA:.2g};\n"
         "Louizos et al. 2018), read by both\nencoders as the sampled graph (figure 14)\n\n"
         "loss = data loss + Σ KL(Bern(q) ‖ Bern(π)) / W\nq = P(gate open), started at q = π", BOX_COLOURS["graph"]),
        (0.5, "3. Binary pass, the gates alone",
         "forward: the read-out, open where q ≥ 1/2,\nso the graph fitted is the graph read out\n\n"
         "backward: from D, the estimated change in\nthe objective if one edge alone flipped (its\n"
         "data part to first order, its prior part,\n± logit π / W, exact), passed to the gate\n"
         "through the sigmoid of the lower right", BOX_COLOURS["propagation"]),
        (0.703, "4. Exact flip checks",
         "the 30 edges whose D promises the largest\nfall are flipped one at a time; a flip is\n"
         "kept only where the exact objective falls\n\nin the simulations about one promised\n"
         "fall in three was real (figure 15)", BOX_COLOURS["readout"]),
        (0.905, "5. What a fit reports",
         "the read-out graph, and its fit;\nper closed edge, its exact evidence E: the\n"
         "rise in data loss × W, in nats, when it\nalone reopens, and P(false) =\nsigmoid(E − logit π);\n"
         "per edge, its closing frequency over\nrefits on halves (stability selection;\nsimulated, not yet in the trainer)", BOX_COLOURS["output"]),
    ]
    framed_row(flow_axis, stages, top=0.97)
    flow_axis.set_title("Fitting the hard edge gates (decision 9.62 (e), built 16 October 2026)", fontsize=12, fontweight="bold", loc="left")

    cost_axis = figure.add_subplot(grid[1, 0])
    open_share = np.linspace(1e-4, 1 - 1e-4, 400)
    tiers = (("curated", 0.9, 1.0, "#2b6ca3"), ("unreferenced", 0.7, 0.5, "#e08e2b"), ("text-mined only", 0.55, 0.5, "#8e44ad"))
    for offset, (name, p, n, colour) in enumerate(tiers):
        prior = float(binary_edge_prior(np.array([p]), np.array([n]))[0])
        kl = open_share * np.log(open_share / prior) + (1 - open_share) * np.log((1 - open_share) / (1 - prior))
        closed_cost, open_cost = -np.log(1 - prior), -np.log(prior)
        cost_axis.plot(open_share, kl, color=colour, linewidth=2,
                       label=f"{name}: p {p:g}, n {n:g}, π {prior:.3f}; closing costs logit π = {np.log(prior / (1 - prior)):.2f} nats")
        cost_axis.plot([0, 1], [closed_cost, open_cost], "o", color=colour, markersize=6)
        # the gap between the two end points, drawn at q = 1 where the open end sits
        cost_axis.annotate("", xy=(1.0 - 0.025 * (offset + 1), open_cost), xytext=(1.0 - 0.025 * (offset + 1), closed_cost),
                           arrowprops={"arrowstyle": "<->", "color": colour, "linewidth": 1.1})
        cost_axis.plot([1.0 - 0.025 * (offset + 1), 1.0], [closed_cost, closed_cost], color=colour, linewidth=0.8, linestyle=":")
    cost_axis.set_xlabel("q, the gate's probability of being open")
    cost_axis.set_ylabel("per-edge penalty, nats")
    cost_axis.set_xlim(0, 1)
    cost_axis.legend(fontsize=7.8, loc="upper center")
    cost_axis.spines[["top", "right"]].set_visible(False)
    cost_axis.set_title("The prior is the cost per edge. Stochastic pass: KL(Bern(q) ‖ Bern(π)), zero at q = π;\n"
                        "binary pass: the end points (dots). Closing costs the gap between them, logit π",
                        fontsize=10, fontweight="bold")

    proxy_axis = figure.add_subplot(grid[1, 1])
    log_odds = np.linspace(-4, 4, 801)
    proxy_axis.plot(log_odds, (log_odds >= 0).astype(float), color="#111111", linewidth=2.4, label="read-out: open when m ≥ 0 (q ≥ 1/2), the cliff")
    for width, style in ((0.25, ":"), (DEFAULT_SURROGATE_WIDTH, "-"), (1.0, "--")):
        proxy_axis.plot(log_odds, 1 / (1 + np.exp(-log_odds / width)), color="#2b6ca3", linestyle=style, linewidth=1.6 if width == DEFAULT_SURROGATE_WIDTH else 1.1,
                        label=f"sigmoid(m / T), T = {width:g}" + (" (default)" if width == DEFAULT_SURROGATE_WIDTH else ""))
    derivative_axis = proxy_axis.twinx()
    slope = 1 / (1 + np.exp(-log_odds / DEFAULT_SURROGATE_WIDTH))
    derivative_axis.fill_between(log_odds, slope * (1 - slope) / DEFAULT_SURROGATE_WIDTH, color="#e08e2b", alpha=0.3,
                                 label=f"σ'(m / T) / T at T = {DEFAULT_SURROGATE_WIDTH:g}: the share of D a gate receives")
    derivative_axis.set_ylim(0, 1.1 * 0.25 / DEFAULT_SURROGATE_WIDTH)
    derivative_axis.set_ylabel("proxy derivative, per unit of open log-odds", color="#9c5d12")
    proxy_axis.set_xlabel("m, the gate's open log-odds, logit q")
    proxy_axis.set_ylabel("gate as the forward pass reads it")
    proxy_axis.set_ylim(-0.03, 1.08)
    handles = proxy_axis.get_legend_handles_labels()[0] + derivative_axis.get_legend_handles_labels()[0]
    proxy_axis.legend(handles, [handle.get_label() for handle in handles], fontsize=7.8, loc="upper left")
    proxy_axis.spines[["top"]].set_visible(False)
    derivative_axis.spines[["top"]].set_visible(False)
    proxy_axis.set_title("The cliff and its proxy: dJ/dm = (J(open) − J(closed)) · σ'(m / T) / T,\n"
                         "so a gate near the threshold moves toward the side that lowers the objective", fontsize=10, fontweight="bold")
    figure.text(0.5, 0.015, "J(open) − J(closed) is −D for an open edge and +D for a closed one (HardConcreteEdgeGates.apply_flip_gradient). Adam scales each gate's step "
                "by its own gradient history, so the width T mostly decides which gates move at all;\nin the simulations T = 0.25, 0.5 and 1 gave the same read-out in world 0. "
                "W is the total loss weight of the fit (the number of training labels in the simulations), so one nat of prior weighs as one label's log-loss.",
                ha="center", fontsize=8)
    save(figure, "edge_gate_training")


def share_reading_values() -> dict:
    """What each encoder multiplies the inputs of a node by, computed with the encoders' own functions on a star of three
    inputs a, b and c into t (one relation, all signs +1, edge weights 1, 2 and 1): with every gate open, with b's gate
    closed, and on the graph without b; and the value of one input against its gate's value on a grid, among three
    inputs and as a sole input. Linear response: expanded_adjacency_entries and normalised_adjacency_values, the entries'
    weights times the gates (as the gated encoder passes them). Message passing: relation_adjacency_parts and
    gated_relation_adjacencies."""
    from mechanistic_pathway_learning.models.edge_gates import gated_relation_adjacencies, relation_adjacency_parts

    num_nodes, target_node = 4, 3
    source, target = torch.tensor([0, 1, 2]), torch.tensor([3, 3, 3])
    relation, sign, weight = torch.zeros(3, dtype=torch.long), torch.ones(3), torch.tensor([1.0, 2.0, 1.0])
    kept = torch.tensor([0, 2])

    def linear(sources, targets, relations, signs, weights, gates) -> list[float]:
        entries = LinearResponseEncoder.expanded_adjacency_entries(num_nodes, ["activates"], sources, targets, relations, signs, edge_weight=weights)
        values = LinearResponseEncoder.normalised_adjacency_values(entries, num_nodes, weight=entries["weight"] * gates[entries["origin"]])
        into_target = {int(entry): 0.0 for entry in range(3)}
        for value, entry_source, entry_target in zip(values.tolist(), entries["source"].tolist(), entries["target"].tolist()):
            if entry_target == target_node:
                into_target[entry_source] += value
        return [into_target[node] for node in range(3)]

    def message(sources, targets, relations, signs, weights, gates) -> list[float]:
        parts = relation_adjacency_parts(torch.stack([sources, targets]), relations, num_nodes, 1, edge_sign=signs, edge_weight=weights)
        adjacency = gated_relation_adjacencies(parts, gates, num_nodes)[0].to_dense()
        return [float(adjacency[target_node, node]) for node in range(3)]

    values = {}
    for name, function in (("linear_response", linear), ("message_passing", message)):
        values[name] = {"every_gate_open": function(source, target, relation, sign, weight, torch.ones(3)),
                        "b_closed": function(source, target, relation, sign, weight, torch.tensor([1.0, 0.0, 1.0])),
                        "graph_without_b": function(source[kept], target[kept], relation[kept], sign[kept], weight[kept], torch.ones(2))}
    gate_grid = np.linspace(0.0, 1.0, 201)
    one_of_three, sole = [], []
    for gate in gate_grid:
        one_of_three.append(linear(source, target, relation, sign, torch.ones(3), torch.tensor([float(gate), 1.0, 1.0]))[0])
        sole.append(linear(source[:1], target[:1], relation[:1], sign[:1], torch.ones(1), torch.tensor([float(gate)]))[0])
    values["gate_grid"], values["one_of_three_inputs"], values["sole_input"] = gate_grid, np.array(one_of_three), np.array(sole)
    return values


def figure_gate_share_reading() -> None:
    """Why a closed gate is a deleted edge in both encoders (decision 9.62), and why the binary pass estimates flips
    instead of following the slope: the numbers are computed by share_reading_values with the encoders' functions."""
    values = share_reading_values()
    figure = plt.figure(figsize=(17.5, 7.4))
    grid = figure.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.95], wspace=0.22, top=0.8, bottom=0.12)
    states = (("every_gate_open", "every gate open"), ("b_closed", "b's gate closed (0)"), ("graph_without_b", "the graph without b"))
    for column, (encoder, title, rule) in enumerate((
            ("linear_response", "Linear response: each input is a share", "value = sign · w · g / Σ (w · g) at t"),
            ("message_passing", "Message passing: a mean over open inputs", "value = sign · w · g / Σ g at t"))):
        axis = figure.add_subplot(grid[0, column])
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.axis("off")
        for row, (state, label) in enumerate(states):
            top = 0.84 - row * 0.3
            numbers = values[encoder][state]
            node_texts = {}
            for index, (name, x) in enumerate((("a  (w 1)", 0.1), ("b  (w 2)", 0.1), ("c  (w 1)", 0.1))):
                node_texts[index] = axis.text(x, top - index * 0.085, name, ha="center", va="center", fontsize=8,
                                              bbox={"boxstyle": "round,pad=0.25", "facecolor": "#dbe9f6" if not (state != "every_gate_open" and index == 1) else "#f4f4f4",
                                                    "edgecolor": "#333333", "linestyle": "-" if not (state == "graph_without_b" and index == 1) else ":"})
            target_text = axis.text(0.62, top - 0.085, "t", ha="center", va="center", fontsize=9, bbox={"boxstyle": "circle,pad=0.3", "facecolor": "#e3f1df", "edgecolor": "#333333"})
            figure.canvas.draw()
            for index in range(3):
                if state == "graph_without_b" and index == 1:
                    continue
                closed = state == "b_closed" and index == 1
                axis.annotate("", xy=(0.5, 0.5), xycoords=target_text, xytext=(0.5, 0.5), textcoords=node_texts[index],
                              arrowprops={"arrowstyle": "-|>", "color": "#c0392b" if closed else "#555555", "linestyle": "--" if closed else "-",
                                          "patchA": node_texts[index].get_bbox_patch(), "patchB": target_text.get_bbox_patch(), "shrinkA": 2, "shrinkB": 2})
                axis.text(0.36, top - index * 0.085 + 0.012, f"{numbers[index]:.3g}" + ("  (g = 0)" if closed else ""), fontsize=8.5, ha="center",
                          color="#c0392b" if closed else "#222222")
            axis.text(0.8, top - 0.085, label, fontsize=8.5, ha="left", va="center", fontweight="bold" if state != "every_gate_open" else "normal")
        same = np.allclose(values[encoder]["b_closed"], values[encoder]["graph_without_b"], atol=1e-6)
        axis.set_title(f"{title}\n{rule}\nb closed {'equals' if same else 'differs from'} the graph without b", fontsize=10, fontweight="bold")

    curve_axis = figure.add_subplot(grid[0, 2])
    grid_values = values["gate_grid"]
    curve_axis.plot(grid_values, values["one_of_three_inputs"], color="#2b6ca3", linewidth=2, label="one of three equal inputs: g / (g + 2)")
    curve_axis.plot(grid_values[1:], values["sole_input"][1:], color="#c0392b", linewidth=2, label="a sole input: 1 for any g > 0")
    curve_axis.plot([0], [values["sole_input"][0]], "o", color="#c0392b", markersize=7)
    curve_axis.plot([0], [1.0], "o", markerfacecolor="white", markeredgecolor="#c0392b", markersize=7)
    curve_axis.annotate("the sole input's cliff:\nits slope is 0 for every g > 0,\nso gradient descent cannot see\nthat closing it removes the input;\n"
                        "the flip estimate D does:\nthe value goes from 1 to 0", xy=(0.0, 0.0), xytext=(0.3, 0.5), fontsize=8,
                        arrowprops={"arrowstyle": "-|>", "color": "#555555"})
    curve_axis.set_xlabel("g, the input's gate")
    curve_axis.set_ylabel("the input's value at t (linear response)")
    curve_axis.set_ylim(-0.05, 1.1)
    curve_axis.legend(fontsize=8, loc="center right", bbox_to_anchor=(1.0, 0.76))
    curve_axis.spines[["top", "right"]].set_visible(False)
    curve_axis.set_title("What the slope sees and what it misses", fontsize=10, fontweight="bold")
    figure.suptitle("A closed gate is a deleted edge in both encoders (decision 9.62); computed with the encoders' adjacency functions", fontsize=12.5, fontweight="bold",
                    y=0.98)
    figure.text(0.5, 0.01, "Linear response divides by the gated sum of the input weights at the target (decision 4.11: weights are shares); message passing divides by "
                "the gated in-degree, the sum of the gates, and keeps the weights in the numerator.\nWith gates of exactly 0 and 1 both give the graph without the closed edges, "
                "which tests/test_hard_edge_gates.py checks end to end, the limiting-step bottleneck included. A soft gate that is its target's sole input "
                "therefore has no slope at all, which is why the binary pass estimates each flip instead.", ha="center", va="bottom", fontsize=8)
    save(figure, "gate_share_reading")


SIMULATION_DIGEST = Path("docs/edge_gate_simulation_results.json")
CALIBRATED_WORLDS = tuple(f"b3_s{seed}" for seed in range(5))
SELECTOR_COLOURS = {"read-out": "#111111", "exact evidence": "#2b6ca3", "posterior P(false)": "#8e44ad", "stability": "#e08e2b"}
SELECTOR_LEGEND = {"read-out": "the read-out: every closed edge", "exact evidence": "closures with exact evidence E ≥ the value shown, nats",
                   "posterior P(false)": "closures with P(false) ≥ the value shown", "stability": "closing frequency over refits ≥ the value shown"}
SELECTOR_LABEL_OFFSETS = {"read-out": (0, 8), "exact evidence": (0, -14), "posterior P(false)": (0, 8), "stability": (0, -14)}


def wilson_interval(successes: int, total: int, z: float = 1.96) -> tuple[float, float]:
    """The Wilson score interval of a binomial proportion; (nan, nan) for an empty total."""
    if total == 0:
        return float("nan"), float("nan")
    share = successes / total
    denominator = 1 + z * z / total
    centre = (share + z * z / (2 * total)) / denominator
    half_width = z * np.sqrt(share * (1 - share) / total + z * z / (4 * total * total)) / denominator
    return float(centre - half_width), float(centre + half_width)


def pooled_selector_points(runs: dict, worlds: tuple[str, ...]) -> dict[str, list[dict]]:
    """Per selector family, one point per threshold: closures per world (mean), the pooled precision over the worlds
    (false closures over closures) and its Wilson interval. The read-out is every closed edge."""
    families: dict[str, list[tuple[str, list[tuple[int, int]]]]] = {"read-out": [], "exact evidence": [], "posterior P(false)": [], "stability": []}

    def counts(entry: dict | None) -> tuple[int, int]:
        if entry is None or not entry.get("closed"):
            return 0, 0
        closed = int(entry["closed"])
        false = int(entry["false"]) if "false" in entry else int(round(entry["precision"] * closed))
        return closed, false

    families["read-out"].append(("every closed edge", [counts(runs[world]["fit"]) for world in worlds]))
    for threshold in ("0.5", "1.0", "2.0"):
        families["exact evidence"].append((f"{float(threshold):g}", [counts(runs[world]["closures_by_evidence"].get(f"at least {threshold} nats")) for world in worlds]))
    for level in ("0.7", "0.8", "0.9", "0.95"):
        families["posterior P(false)"].append((level, [counts(runs[world]["closures_by_posterior_false"].get(f"at least {level}")) for world in worlds]))
    for threshold in ("0.6", "0.7", "0.8", "0.9"):
        families["stability"].append((threshold, [counts(runs[world]["stability"]["thresholds"].get(threshold)) for world in worlds]))
    points: dict[str, list[dict]] = {}
    for family, thresholds in families.items():
        points[family] = []
        for label, per_world in thresholds:
            closed, false = sum(count[0] for count in per_world), sum(count[1] for count in per_world)
            if closed == 0:
                continue
            low, high = wilson_interval(false, closed)
            points[family].append({"label": label, "closures_per_world": closed / len(worlds), "precision": false / closed, "low": low, "high": high, "closed": closed})
    return points


def posterior_calibration(runs: dict, worlds: tuple[str, ...], edges: tuple[float, ...] = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) -> list[dict]:
    """The closed edges of the given worlds pooled and binned by their P(false); per bin, the mean P(false), the share
    that are false and its Wilson interval."""
    rows = [row for world in worlds for row in runs[world].get("closed_edges", []) if row[4] is not None]
    columns = runs[worlds[0]]["closed_edges_columns"]
    posterior_index, real_index = columns.index("posterior_false"), columns.index("real")
    bins = []
    for low, high in zip(edges, edges[1:]):
        inside = [row for row in rows if low <= row[posterior_index] < high or (high == edges[-1] and row[posterior_index] == high)]
        if not inside:
            continue
        false = sum(1 - row[real_index] for row in inside)
        interval = wilson_interval(false, len(inside))
        bins.append({"low": low, "high": high, "edges": len(inside), "mean_posterior_false": float(np.mean([row[posterior_index] for row in inside])),
                     "share_false": false / len(inside), "interval": interval})
    return bins


def figure_gate_simulation_recovery() -> None:
    """Five planted worlds, calibrated priors, the encoder's true parameters: which closures to believe (decision 9.62)."""
    digest = json.loads(SIMULATION_DIGEST.read_text())
    runs = digest["runs"]
    worlds = tuple(world for world in CALIBRATED_WORLDS if world in runs)
    figure, (selector_axis, calibration_axis, estimator_axis) = plt.subplots(1, 3, figsize=(19, 6.6), gridspec_kw={"width_ratios": [1.15, 0.95, 0.95]})
    base_rates = [runs[world]["world"]["false_edges"] / runs[world]["world"]["edges"] for world in worlds]
    for family, points in pooled_selector_points(runs, worlds).items():
        if not points:
            continue
        x = [point["closures_per_world"] for point in points]
        y = [point["precision"] for point in points]
        errors = np.clip(np.array([[point["precision"] - point["low"] for point in points], [point["high"] - point["precision"] for point in points]]), 0.0, None)
        selector_axis.errorbar(x, y, yerr=errors, color=SELECTOR_COLOURS[family], marker="o", markersize=5, linewidth=1.4, capsize=3, label=SELECTOR_LEGEND[family])
        offset = SELECTOR_LABEL_OFFSETS[family]
        for point in points:
            selector_axis.annotate(point["label"], (point["closures_per_world"], point["precision"]), textcoords="offset points", xytext=offset,
                                   fontsize=7.2, color=SELECTOR_COLOURS[family], ha="center")
    selector_axis.axhline(np.mean(base_rates), color="#999999", linestyle="--", linewidth=1)
    selector_axis.text(0.98, np.mean(base_rates) + 0.015, f"share of all edges that are false: {np.mean(base_rates):.2f}", fontsize=7.5, color="#666666",
                       ha="right", transform=selector_axis.get_yaxis_transform())
    selector_axis.set_xscale("log")
    selector_axis.set_xticks([3, 5, 10, 20, 50], ["3", "5", "10", "20", "50"])
    selector_axis.set_ylim(0, 1.05)
    selector_axis.set_xlabel("edges closed per world (log scale)")
    selector_axis.set_ylabel("share of the closed edges that are false (pooled, 95% Wilson interval)")
    selector_axis.legend(fontsize=8, loc="lower left")
    selector_axis.spines[["top", "right"]].set_visible(False)
    selector_axis.set_title(f"Which closures to believe, {len(worlds)} worlds pooled:\nthe read-out against three stricter rules", fontsize=10.5, fontweight="bold")

    bins = posterior_calibration(runs, worlds)
    calibration_axis.plot([0.5, 1.0], [0.5, 1.0], color="#999999", linestyle="--", linewidth=1, label="calibrated")
    for item in bins:
        calibration_axis.errorbar(item["mean_posterior_false"], item["share_false"], yerr=[[max(item["share_false"] - item["interval"][0], 0.0)], [max(item["interval"][1] - item["share_false"], 0.0)]],
                                  color="#8e44ad", marker="o", capsize=3)
        calibration_axis.annotate(f"{item['edges']}", (item["mean_posterior_false"], item["share_false"]), textcoords="offset points", xytext=(6, 4), fontsize=7.5, color="#8e44ad")
    calibration_axis.set_xlim(0.45, 1.02)
    calibration_axis.set_ylim(0, 1.05)
    calibration_axis.set_xlabel("P(false) = sigmoid(E − logit π), binned")
    calibration_axis.set_ylabel("share false in the bin (95% Wilson interval)")
    calibration_axis.legend(fontsize=8, loc="lower right")
    calibration_axis.spines[["top", "right"]].set_visible(False)
    calibration_axis.set_title("Is the posterior from the exact evidence calibrated?\n(numbers: closed edges in each bin)", fontsize=10.5, fontweight="bold")

    pairs = np.array([pair for world in worlds for pair in runs[world].get("estimator", {}).get("per_edge_data_part_nats", [])], dtype=float)
    if len(pairs):
        estimator_axis.scatter(pairs[:, 1], pairs[:, 0], s=7, alpha=0.45, color="#2b6ca3", edgecolors="none")
        limit = float(np.abs(pairs).max()) * 1.1
        estimator_axis.plot([-limit, limit], [-limit, limit], color="#999999", linestyle="--", linewidth=1)
        estimator_axis.set_xscale("symlog", linthresh=1e-3)
        estimator_axis.set_yscale("symlog", linthresh=1e-3)
        correlation = float(np.corrcoef(pairs[:, 0], pairs[:, 1])[0, 1])
        agree = float(np.mean(np.sign(pairs[:, 0]) == np.sign(pairs[:, 1])))
        estimator_axis.set_title(f"The first-order flip estimate against the exact change\n{len(pairs)} random edges; correlation {correlation:.2f}, "
                                 f"signs agree {agree:.0%}", fontsize=10.5, fontweight="bold")
    estimator_axis.set_xlabel("exact change in data loss × W if the edge alone flips, nats (symlog)")
    estimator_axis.set_ylabel("first-order estimate, nats (symlog)")
    estimator_axis.spines[["top", "right"]].set_visible(False)
    worlds_text = runs[worlds[0]]["world"]
    figure.suptitle(f"Edge-gate recovery in planted worlds ({worlds_text['nodes']} nodes, about {worlds_text['edges']:,} edges, {worlds_text['perturbations']} perturbations, "
                    f"{worlds_text['symptoms']} symptoms; docs/edge_gate_simulation.md)", fontsize=12, fontweight="bold")
    closed_rows = [row for world in worlds for row in runs[world].get("closed_edges", [])]
    below = sum(1 for row in closed_rows if row[4] is not None and row[4] < 0.5)
    figure.text(0.5, -0.04, "Every world is the most favourable case: priors calibrated by construction (each edge is real with its prior π) and the encoder's parameters "
                "the truth's, so only the gates learn. E is the exact rise in the training data loss, times W, when the edge alone is reopened with every other\nedge "
                f"as fitted. At a local optimum under single flips a closed edge has E > logit π, so P(false) > 1/2; {below} of the {len(closed_rows)} closures fell below one "
                "half and are not binned. Stability: closing frequency over 20 refits on complementary halves of the training perturbations (Shah and Samworth 2013). "
                "Drawn from docs/edge_gate_simulation_results.json.", ha="center", fontsize=8)
    figure.tight_layout()
    save(figure, "gate_simulation_recovery")


CONDITION_ROWS = (  # (label, runs pooled)
    ("calibrated priors, true parameters\n(5 worlds pooled)", CALIBRATED_WORLDS),
    ("calibrated priors, true parameters\n(world 0 alone)", ("b3_s0",)),
    ("without the binary pass and the\nflip checks (5 worlds pooled)", tuple(f"b3_s{seed}_nobinary" for seed in range(5))),
    ("32 binary epochs at a constant\nrate (5 worlds pooled; batch 4)", tuple(f"b4_s{seed}_epochs32" for seed in range(5))),
    ("32 binary epochs, rate x 0.8 each\n(5 worlds pooled; batch 4)", tuple(f"b4_s{seed}_epochs32_decay08" for seed in range(5))),
    ("curated edges taken as known, no\ngate (5 worlds pooled; 9.63)", tuple(f"b3_s{seed}_curated_fixed" for seed in range(5))),
    ("3,000 perturbations (world 0)", ("b3_s0_p3000",)),
    ("priors too optimistic: fewer edges\nreal than π says (world 0)", ("b3_s0_optimistic_priors",)),
    ("priors too pessimistic: more edges\nreal than π says (world 0)", ("b3_s0_pessimistic_priors",)),
    ("half the true positives observed\n(world 0)", ("b3_s0_pu",)),
    ("joint fit, fresh start, readout on\nthe true nodes (world 0)", ("b3_s0_joint_random_support",)),
    ("joint fit, fresh start, readout on\nevery node (world 0)", ("b2_s0_joint_random_dense",)),
    ("joint fit from the truth, readout\non every node (world 0)", ("b2_s0_joint_oracle_dense",)),
    ("v5's positives per gated edge\n(2,000 nodes, 20 symptoms)", ("b2_realratio",)),
)


def condition_summary(runs: dict, names: tuple[str, ...]) -> dict | None:
    """Pooled over the named runs that exist: closures, the read-out's precision and recall, the precision of the
    closures with at least 1 nat of exact evidence, the share of edges that are false, and the share of the attainable
    held-out gain the fitted graph realised, (open - fitted) / (open - true), averaged over the runs."""
    present = [runs[name] for name in names if name in runs]
    if not present:
        return None
    closed = sum(run["fit"]["closed"] for run in present)
    false_closed = sum(int(round((run["fit"]["precision"] or 0.0) * run["fit"]["closed"])) for run in present)
    false_edges = sum(run["fit"]["false_edges"] for run in present)
    edges = sum(run["world"]["edges"] for run in present)
    evidence = [run["closures_by_evidence"]["at least 1.0 nats"] for run in present if run.get("closures_by_evidence")]
    evidence_closed = sum(entry["closed"] for entry in evidence)
    evidence_false = sum(int(round((entry["precision"] or 0.0) * entry["closed"])) for entry in evidence)
    gains = []
    for run in present:
        loss = run["held_out_loss"]
        attainable = loss["every_edge_open"] - loss["true_graph"]
        if attainable > 0:
            gains.append((loss["every_edge_open"] - loss["fitted_graph"]) / attainable)
    return {"runs": len(present), "closed": closed, "precision": false_closed / closed if closed else None, "recall": false_closed / false_edges if false_edges else None,
            "evidence_closed": evidence_closed, "evidence_precision": evidence_false / evidence_closed if evidence_closed else None,
            "base_rate": false_edges / edges, "held_out_gain_share": float(np.mean(gains)) if gains else None}


def figure_gate_simulation_conditions() -> None:
    """How the read-out fares when the favourable assumptions are dropped one at a time, and how far v5 is from the worlds."""
    digest = json.loads(SIMULATION_DIGEST.read_text())
    runs = digest["runs"]
    rows = [(label, condition_summary(runs, names)) for label, names in CONDITION_ROWS]
    rows = [(label, summary) for label, summary in rows if summary is not None]
    figure, (precision_axis, gain_axis, density_axis) = plt.subplots(1, 3, figsize=(19, 7.6), gridspec_kw={"width_ratios": [1.25, 0.8, 0.95]})
    positions = np.arange(len(rows))[::-1]
    precision = [summary["precision"] or 0.0 for _, summary in rows]
    precision_axis.barh(positions, precision, color="#7a9cc6", height=0.62, label="read-out: every closed edge")
    for position, (_, summary) in zip(positions, rows):
        precision_axis.plot([summary["base_rate"]], [position], marker="|", color="#c0392b", markersize=16, markeredgewidth=2.2)
        if summary["evidence_precision"] is not None:
            precision_axis.plot([summary["evidence_precision"]], [position], marker="D", color="#2b6ca3", markersize=6)
        precision_axis.text(1.02, position, f"{summary['closed']} closed, recall {summary['recall']:.1%}; E ≥ 1 nat: {summary['evidence_closed']}",
                            va="center", fontsize=7.6)
    precision_axis.plot([], [], marker="|", color="#c0392b", linestyle="none", markersize=12, markeredgewidth=2.2, label="share of edges that are false")
    precision_axis.plot([], [], marker="D", color="#2b6ca3", linestyle="none", label="closures with E ≥ 1 nat")
    precision_axis.set_yticks(positions, [label for label, _ in rows], fontsize=8)
    precision_axis.set_xlim(0, 1.6)
    precision_axis.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    precision_axis.set_xlabel("share of the closed edges that are false")
    precision_axis.legend(fontsize=7.8, loc="upper center", bbox_to_anchor=(0.45, -0.1), ncol=3, frameon=False)
    precision_axis.spines[["top", "right"]].set_visible(False)
    precision_axis.set_title("Precision of the closures under each condition", fontsize=10.5, fontweight="bold")

    gains = [summary["held_out_gain_share"] for _, summary in rows]
    gain_axis.barh(positions, [gain if gain is not None else 0.0 for gain in gains], color=["#2b6ca3" if (gain or 0) >= 0 else "#c0392b" for gain in gains], height=0.62)
    for position, gain in zip(positions, gains):
        if gain is not None:
            gain_axis.text(gain + (0.03 if gain >= 0 else -0.03), position, f"{gain:+.2f}", va="center", ha="left" if gain >= 0 else "right", fontsize=7.6)
    gain_axis.axvline(0, color="#555555", linewidth=1)
    gain_axis.set_yticks(positions, [""] * len(rows))
    finite_gains = [gain for gain in gains if gain is not None]
    gain_axis.set_xlim(min(-0.6, min(finite_gains) * 1.3) if finite_gains else -1, 1.2)
    gain_axis.set_xlabel("(loss, every edge open − loss, fitted graph) /\n(loss, every edge open − loss, true graph), held out")
    gain_axis.spines[["top", "right"]].set_visible(False)
    gain_axis.set_title("Share of the attainable held-out gain\nthe fitted graph realised (1 = the true graph)", fontsize=10.5, fontweight="bold")

    density_rows = []
    for name, label in (("b3_s0", "world 0, 800 perturbations"), ("b3_s0_p3000", "world 0, 3,000 perturbations"), ("b2_realratio", "the v5-ratio world")):
        if name in runs and "training_positives_per_edge" in runs[name]["world"]:
            world = runs[name]["world"]
            density_rows.append((label, world["training_labels_per_edge"], world["training_positives_per_edge"]))
    if digest.get("v5"):
        density_rows.append(("v5, every perturbation\n(an upper bound for a fold)", digest["v5"]["cells_per_gated_edge"], digest["v5"]["positives_per_gated_edge"]))
    density_positions = np.arange(len(density_rows))[::-1]
    density_axis.barh(density_positions + 0.18, [row[1] for row in density_rows], height=0.34, color="#d5dbe1", label="label cells per gated edge")
    density_axis.barh(density_positions - 0.18, [row[2] for row in density_rows], height=0.34, color="#2b6ca3", label="positive labels per gated edge")
    for position, (_, cells, positives) in zip(density_positions, density_rows):
        density_axis.text(cells * 1.15, position + 0.18, f"{cells:.2g}", va="center", fontsize=7.6)
        density_axis.text(positives * 1.15, position - 0.18, f"{positives:.2g}", va="center", fontsize=7.6)
    density_axis.set_xscale("log")
    density_axis.set_yticks(density_positions, [row[0] for row in density_rows], fontsize=8)
    density_axis.set_xlim(0.01, 60)
    density_axis.legend(fontsize=7.8, loc="lower right")
    density_axis.spines[["top", "right"]].set_visible(False)
    density_axis.set_xlabel("per gated edge (log scale); the worlds count their training rows")
    density_axis.set_title("How much the labels say per edge", fontsize=10.5, fontweight="bold")
    figure.suptitle("Edge-gate recovery when the favourable assumptions are dropped one at a time (docs/edge_gate_simulation.md)", fontsize=12, fontweight="bold")
    figure.text(0.5, -0.08, "Each condition changes one thing from the calibrated world of figure 15; a world-0 condition is to be read against the world-0 row. Precision is "
                "the share of closed edges that are false; the red mark is that share among all edges, what closing at random would give. A held-out gain share\nbelow 0 "
                "means the fitted graph predicted the held-out perturbations worse than the database graph with every edge open; no bar where the true graph itself did not "
                "predict better than every edge open (with half the positives observed, the true graph is not the best\npredictor of the observed labels). The share is "
                "noisy: the attainable gain is 0.0013 to 0.0060 nats per held-out label, and over the five calibrated worlds the share ran from -1.20 (world 0) to +0.84. "
                "Joint fits learn the encoder's parameters and the readout\nwith the gates, with AdamW and early stopping on an inner validation share, as the trainer "
                "does. Drawn from docs/edge_gate_simulation_results.json.", ha="center", fontsize=8)
    figure.tight_layout()
    save(figure, "gate_simulation_conditions")


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
    "cell_state": (figure_cell_state, [CELL_STATE_GRAPH_DIRECTORY / "nodes.parquet"]),
    "graph_composition": (figure_graph_composition, [CONFIRMATORY_GRAPH_DIRECTORY / "nodes.parquet", CONFIRMATORY_GRAPH_DIRECTORY / "edges.parquet"]),
    "allele_perturbations": (figure_allele_perturbations, [EVIDENCE_DIRECTORY / "evidence_reports.parquet", EVIDENCE_DIRECTORY / "evidence_summary.json",
                                                           CONFIRMATORY_GRAPH_DIRECTORY / "edges.parquet"]),
    "edge_priors": (figure_edge_priors, [CONFIRMATORY_GRAPH_DIRECTORY / "edge_evidence.parquet", EDGE_PRIOR_TIERS]),
    "edge_gate_training": (figure_edge_gate_training, []),
    "gate_share_reading": (figure_gate_share_reading, []),
    "gate_simulation_recovery": (figure_gate_simulation_recovery, [SIMULATION_DIGEST]),
    "gate_simulation_conditions": (figure_gate_simulation_conditions, [SIMULATION_DIGEST]),
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
