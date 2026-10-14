"""Decision 4.2: what a protein complex taking its limiting member's cell-class weight would change.

The owner, 13 October 2026: "Taking the complex's limiting member weight per type makes sense i think but only if it
applies at the complex itself. That's what you mean, right?"

The wiring the question turns on. A complex is a `protein_entity` node; its members reach it by `member_of` edges
from gene nodes, and it reaches reactions by `catalyzed_by`, `activates`, `inhibits` and `substrate_of`. So a
response travels gene -> complex -> reaction and passes **through** the complex node. The encoder scales a
response by the destination node's class weight at each step, so a weight placed on the complex gates every path
through it exactly once, which is what "at the complex itself" means and is the only placement with the right
semantics:

  - **At the complex.** The gate is `min` over the members, taken within each class. An obligate subunit that a
    class does not express shuts the complex in that class however abundant the others are, which is what a complex
    is.
  - **On the member genes only**, which is close to today's behaviour, gates each path by the member it travelled
    through. A response entering through an abundant subunit passes even when another subunit is missing, so the
    complex behaves like a maximum over its members rather than a minimum. That is the wrong physiology.
  - **On the downstream reaction** would put the complex's weight somewhere the complex is not, double-gate any
    path that also crosses the complex node, and collide with the Human-GEM gene rule that already sets a
    reaction's weight.

So this script measures the first, and reports three things the decision needs: how many complexes can carry a
weight at all, how far the gate moves from the present flat 1, and how often the limiting member differs between
classes, because a member chosen once globally would not be the same rule as a minimum per class.

Counting only: no graph, label or model input is changed.

Usage:
  PYTHONPATH=. python experiments/label_review/complex_limiting_member_weight.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.brain_expression_descriptors import (
    gene_node_ensembl_ids,
    read_zipped_table,
    unambiguous_symbol_to_ensembl,
)
from mechanistic_pathway_learning.graph.brain_expression_weights import (
    COMPLEX_NODE_TYPE,
    MEMBER_RELATION,
    complex_limiting_member_expression,
    relative_expression_weights,
)

GRAPH_DIRECTORY = Path("data/processed/graph_full_neuronal")
CLASS_WEIGHTS = Path("data/processed/cell_class_weights/full_neuronal_cell_class_weights.parquet")
OUTPUT_DOCUMENT = Path("docs/complex_limiting_member_weight.md")
ALL_CELLS_CLASS = "all_cells"
REGULATOR_RELATION_SETS = {
    "catalyzed_by alone": ("catalyzed_by",),
    "catalysis and regulation": ("catalyzed_by", "activates", "inhibits", "regulates_release"),
    "any edge into the reaction": ("catalyzed_by", "activates", "inhibits", "regulates_release", "substrate_of"),
}


def the_other_order(member_edges: pd.DataFrame, node_weights: pd.DataFrame, classes: list[str]) -> pd.DataFrame:
    """The rejected aggregation: a minimum over the members' **relative weights** rather than their raw expression.

    Kept so the registered order can be scored against it. `node_weights` carries 1 for a node with no expression,
    which a minimum ignores, so this is computable wherever the registered order is."""
    members = member_edges.join(node_weights[classes], on="source_id", how="inner")
    return members.groupby("target_id")[classes].min()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=GRAPH_DIRECTORY)
    parser.add_argument("--class-weights", type=Path, default=CLASS_WEIGHTS)
    parser.add_argument("--class-expression", type=Path,
                        default=Path("data/processed/brain_expression/class_expression_for_weights.parquet"))
    parser.add_argument("--hpa-region", type=Path,
                        default=Path("data/raw/hpa/rna_brain_region_hpa/rna_brain_region_hpa.tsv.zip"))
    parser.add_argument("--document", type=Path, default=OUTPUT_DOCUMENT)
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet",
                            columns=["source_id", "target_id", "relation_type"])
    stored_weights = pd.read_parquet(arguments.class_weights)
    classes = [column for column in stored_weights.columns if column != ALL_CELLS_CLASS]
    class_expression = pd.read_parquet(arguments.class_expression)
    class_expression.columns = [column.removeprefix("class_") for column in class_expression.columns]
    class_expression = class_expression[[column for column in classes if column in class_expression.columns]]
    region_table = read_zipped_table(arguments.hpa_region)
    symbol_to_ensembl = unambiguous_symbol_to_ensembl(region_table["Gene name"], region_table["Gene"])

    complexes = set(nodes.node_id[nodes.node_type == COMPLEX_NODE_TYPE])
    member_edges = edges[(edges.relation_type == MEMBER_RELATION) & edges.target_id.isin(complexes)]
    limiting_expression = complex_limiting_member_expression(nodes, edges, class_expression, symbol_to_ensembl)
    registered = relative_expression_weights(limiting_expression)
    rejected = the_other_order(member_edges, stored_weights, classes).reindex(registered.index)

    reaction_nodes = set(nodes.node_id[nodes.node_type == "reaction"])
    reaction_counts = {}
    for label, relations in REGULATOR_RELATION_SETS.items():
        regulators = edges[edges.relation_type.isin(relations)]
        sources_of = regulators.groupby("target_id").source_id.apply(set)
        only_complexes = [reaction for reaction, sources in sources_of.items()
                          if reaction in reaction_nodes and sources and sources <= complexes]
        gated = [reaction for reaction in only_complexes
                 if any(source in registered.index for source in sources_of[reaction])]
        reaction_counts[label] = (len(only_complexes), len(gated))

    argmin_per_class = {}
    for class_name in class_expression.columns:
        joined = member_edges.join(
            class_expression[[class_name]].rename_axis("ensembl").reset_index()
            .merge(pd.DataFrame({"ensembl": gene_node_ensembl_ids(nodes, symbol_to_ensembl),
                                 "node_id": nodes.node_id}).dropna(), on="ensembl")
            .set_index("node_id")[[class_name]], on="source_id", how="inner").dropna(subset=[class_name])
        if len(joined):
            argmin_per_class[class_name] = joined.loc[joined.groupby("target_id")[class_name].idxmin()] \
                .set_index("target_id").source_id
    distinct_limiting_members = pd.DataFrame(argmin_per_class).nunique(axis=1)

    write_document(arguments.document, list(class_expression.columns), complexes, member_edges, registered, rejected,
                   reaction_counts, distinct_limiting_members)


def spread(table: pd.DataFrame, classes: list[str]) -> dict:
    values = table[classes].to_numpy().ravel()
    values = values[np.isfinite(values)]
    best_class = table[classes].max(axis=1).dropna()
    return {"median": float(np.median(values)), "low": float(np.quantile(values, .25)),
            "high": float(np.quantile(values, .75)), "below_a_tenth": float(np.mean(values < 0.1)),
            "half_or_more": float(np.mean(values >= 0.5)),
            "available_somewhere": float(np.mean(best_class >= 0.9))}


def write_document(document: Path, classes: list[str], complexes: set, member_edges: pd.DataFrame,
                   registered: pd.DataFrame, rejected: pd.DataFrame, reaction_counts: dict,
                   distinct_limiting_members: pd.Series) -> None:
    kept, dropped = spread(registered, classes), spread(rejected, classes)
    lines = [
        "# Decision 4.2: a complex taking its limiting member's cell-class weight", "",
        "Generated by `experiments/label_review/complex_limiting_member_weight.py`.", "",
        "## The question, and the answer", "",
        "The owner, 13 October 2026: \"Taking the complex's limiting member weight per type makes sense i think but "
        "only if it applies at the complex itself. That's what you mean, right?\"", "",
        "**Yes, and it is the only placement with the right semantics.** A complex is a `protein_entity` node that "
        "its members reach by `member_of` and that reaches reactions by `catalyzed_by`, `activates`, `inhibits` and "
        "`substrate_of`, so a response travels gene -> complex -> reaction and crosses the complex node. The "
        "encoder scales a response by the destination's class weight at each step, so a weight on the complex gates "
        "every path through it exactly once. The two alternatives are both wrong:", "",
        "- **On the member genes alone**, which is what happened before this decision, gates each path by whichever "
        "  member it arrived through. An abundant subunit then lets the response past while an obligate one is "
        "  missing: a maximum over the members where a complex is a minimum.",
        "- **On the downstream reaction** puts the weight where the complex is not, double-gates any path that also "
        "  crosses the complex node, and collides with the Human-GEM gene rule that already sets a reaction's "
        "  weight.", "",
        "## What the minimum is taken over, which is not what was first assumed", "",
        "A complex is the **\"and\"** of its subunits, and `evaluate_gene_rule` already reads \"and\" as a minimum over "
        "**raw** nCPM, with the relative weight formed afterwards. So the complex follows the same order: minimum "
        "on the raw expression, then `x_k / (max_j x_j + 1)`. The first version of this measurement took the "
        "minimum over the members' **relative weights** instead, which is a different quantity and a worse one: "
        "the relative weight is a specificity measure, so two subunits with different top classes pull each other "
        "down in every class and the complex ends up shut everywhere.", "",
        "| over the complexes that can carry a gate | registered: min of raw, then relative | rejected: min of relative weights |",
        "|---|---|---|",
        f"| median gate | **{kept['median']:.4f}** | {dropped['median']:.4f} |",
        f"| quartiles | {kept['low']:.4f} to {kept['high']:.4f} | {dropped['low']:.4f} to {dropped['high']:.4f} |",
        f"| below 0.1, largely shut in that class | {kept['below_a_tenth']:.1%} | {dropped['below_a_tenth']:.1%} |",
        f"| at or above 0.5, barely gated | {kept['half_or_more']:.1%} | {dropped['half_or_more']:.1%} |",
        f"| **complexes available in at least one class at 0.9 or above** | **{kept['available_somewhere']:.1%}** | {dropped['available_somewhere']:.1%} |",
        "",
        "The last row is the one that decides it. A complex that reaches 0.9 in no class is a complex nothing "
        f"propagates through anywhere, and the rejected order leaves {1 - dropped['available_somewhere']:.0%} of "
        f"them in that state against {1 - kept['available_somewhere']:.0%} under the registered order.", "",
        "## Coverage", "",
        "| reading | count |", "|---|---|",
        f"| complexes (`protein_entity` nodes) | {len(complexes)} |",
        f"| with at least one member | {member_edges.target_id.nunique()} |",
        f"| **with at least one member carrying class expression, so a gate is computable** | **{len(registered)}** |",
        f"| left flat at 1, as metabolites and ruleless reactions already are | {len(complexes) - len(registered)} |",
        "",
        "## The reactions it reaches, and a correction", "",
        "`docs/open_decisions.md` recorded **413** reactions regulated only by a complex. That number does not "
        "reproduce under any reading of \"regulated\", and is corrected here to the three that do:", "",
        "| what counts as regulating the reaction | reactions whose regulators are all complexes | of those, now gated |",
        "|---|---|---|",
    ]
    for label, (total, gated) in reaction_counts.items():
        lines.append(f"| {label} | {total} | **{gated}** |")
    lines += [
        "",
        "The middle row is the reading this decision is about, since `substrate_of` is a complex entering a "
        "reaction as a reactant rather than regulating it.", "",
        "## Why the minimum has to be taken per class", "",
        f"Over the {len(distinct_limiting_members)} complexes with a computable gate, the limiting member is the "
        f"same in every class for {int((distinct_limiting_members == 1).sum())} and differs across classes for "
        f"**{int((distinct_limiting_members > 1).sum())}**, up to {int(distinct_limiting_members.max())} distinct "
        "limiting subunits. Choosing one limiting subunit globally would therefore be a different rule and a "
        "weaker one, carrying the wrong subunit's profile in every class where another is scarcer. The minimum is "
        "taken within each class, which is what \"per type\" asks for.", "",
        "## Two consequences, one of them a problem", "",
        "**The compounding along a path is the encoder's standing convention, not a new double count.** A response "
        "on gene -> complex -> reaction is scaled by the complex's gate and also by the member gene's own weight, "
        "and that member's abundance is already inside the minimum. Every path in this encoder is scaled once per "
        "node it crosses, so a complex is no different from any other node on it; exempting complexes from the "
        "convention would need an argument that does not exist. Recorded so the choice is visible.", "",
        "**H2's rewired-expression control does not touch the complex gate, and should.** That control permutes the "
        "Human-GEM gene rules and recomputes the reactions' cell-class weights "
        "(`rewired_reaction_cell_class_weights`, which returns reaction rows only). A complex takes its gate from "
        "`member_of` edges and not from a rule, so the permutation leaves every complex gate intact and the "
        "\"rewired\" condition still carries real expression structure through "
        f"{len(registered)} complexes. The analogous permutation is of gene-to-complex membership, and extending "
        "the control that way is a change to a registered control: it is open decision 4.2b, not a build.", "",
    ]
    document.parent.mkdir(parents=True, exist_ok=True)
    document.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
