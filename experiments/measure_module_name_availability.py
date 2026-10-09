"""Name every analysed noisy-OR module from its support and from its behaviour, and count how often a name is earned.

The user, 9 October 2026: "can you please identify better names for the modules than a simple indexing?" The name is
built by `module_label` in experiments/analyze_pathway_modules.py, which takes the first name the support actually
carries: the curated module it overlaps by at least three genes, else the subsystem holding at least a quarter of the
support, else the node type holding at least three quarters, else "mixed support". Every name states the share it
rests on.

The user then asked whether a name could instead come from the nodes, edges and training examples: "this module seems
to encode for movement-related symptoms or encodes for effects on protein complexes when a knockout of x,y or z
happens". That second name is built by mechanistic_pathway_learning/evaluation/module_behaviour.py from the run's own
saved activations and links, and this script counts it beside the first.

A name is only as good as the reading behind it, so this script reads every module_analysis.json already written under
runs/ and reports which rule fired, per run and in total. A module that comes out "mixed support" with no driver is not
a naming failure: it is a module with no identity to name, which is a fact about the module and belongs beside the
module-health reading of docs/module_health.md.

Reads the analyses only. Trains nothing, changes no run and reads no lockbox outcome. Idempotent.

Usage:
  python experiments/measure_module_name_availability.py
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pandas as pd

from experiments.analyze_pathway_modules import (CURATED_GENES_FOR_A_NAME, NODE_TYPE_SHARE_FOR_A_NAME,
                                                 SUBSYSTEM_ENRICHMENT_FOR_A_NAME, SUBSYSTEM_SHARE_FOR_A_NAME,
                                                 module_label, subsystem_shares_of_graph)
from mechanistic_pathway_learning.evaluation.module_behaviour import behaviour_label

OUTPUT_DOCUMENT = Path("docs/module_names.md")
RULE_ORDER = ["curated module", "subsystem", "node type", "mixed", "empty"]
DRIVER_ORDER = ["group, path followed", "group, path not followed", "no group, path followed",
                "no group, path not followed", "no driver", "not recorded"]
SYMPTOM_ORDER = ["domain", "symptoms without a domain", "no symptom above the margin",
                 "profile matches the other modules", "not recorded"]
# runs/smoke trains a few steps to check the code path, so its gates mean nothing; counting its supports would pad
# the denominator of the headline count with a run no reading rests on
EXCLUDED_RUN_NAMES = ("smoke",)


def rule_that_named_it(label: str) -> str:
    if label == "empty support":
        return "empty"
    if label.startswith("mixed support"):
        return "mixed"
    if "across subsystems" in label:
        return "node type"
    if "of its genes" in label:
        return "curated module"
    return "subsystem"


def driver_rule(behaviour: dict) -> str:
    """Which driver reading one module earned: a shared group, and whether the drivers reach the support."""
    if not behaviour.get("available", False):
        return "not recorded"
    drivers = behaviour.get("drivers", {})
    if not drivers.get("responds", False):
        return "no driver"
    group = "group" if drivers.get("group") else "no group"
    if drivers.get("drivers_reaching_the_support") is None:
        return f"{group}, path not followed"
    return f"{group}, path {'followed' if drivers.get('reach_is_enriched') else 'not followed'}"


def symptom_rule(behaviour: dict) -> str:
    """Which symptom-side reading one module earned."""
    if not behaviour.get("available", False):
        return "not recorded"
    symptom_side = behaviour.get("symptom_side", {})
    if symptom_side.get("domain"):
        return "domain"
    if not symptom_side.get("distinct_profile", True):
        return "profile matches the other modules"
    if not symptom_side.get("symptoms"):
        return "no symptom above the margin"
    return "symptoms without a domain"


def names_of_run(analysis_path: Path, backgrounds: dict[str, dict[str, float]]) -> tuple[Counter, Counter, Counter, list[tuple[str, str, str]]]:
    """Which rule named each module of one run, and the named modules themselves.

    The subsystem rule needs the graph the run was trained on, because it asks whether a subsystem is over-represented
    in the support rather than merely large; each graph's background is read once and reused.
    """
    analysis = json.loads(analysis_path.read_text())
    graph_dir = analysis.get("graph_dir", "")
    if graph_dir and graph_dir not in backgrounds:
        nodes_path = Path(graph_dir) / "nodes.parquet"
        backgrounds[graph_dir] = (subsystem_shares_of_graph(pd.read_parquet(nodes_path, columns=["subsystem"]))
                                  if nodes_path.exists() else {})
    background = backgrounds.get(graph_dir, {})
    rules, drivers, symptom_sides, named = Counter(), Counter(), Counter(), []
    for split_name, split in analysis.get("splits", {}).items():
        for module_name, description in split.get("modules", {}).items():
            label = module_label(description, background)
            rule = rule_that_named_it(label)
            rules[rule] += 1
            behaviour = description.get("behaviour", {})
            drivers[driver_rule(behaviour)] += 1
            symptom_sides[symptom_rule(behaviour)] += 1
            if rule in ("curated module", "subsystem") or driver_rule(behaviour).startswith("group, path followed"):
                named.append((split_name, module_name, f"{label} | {behaviour_label(behaviour) if behaviour else 'no behaviour recorded'}"))
    return rules, drivers, symptom_sides, named


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--runs-dir", type=Path, default=Path("runs"))
    arguments = parser.parse_args()

    per_run, named_modules, totals, backgrounds = {}, {}, Counter(), {}
    per_run_drivers, per_run_symptoms, driver_totals, symptom_totals = {}, {}, Counter(), Counter()
    for analysis_path in sorted(arguments.runs_dir.rglob("module_analysis.json")):
        if analysis_path.parent.name in EXCLUDED_RUN_NAMES:
            continue
        rules, drivers, symptom_sides, named = names_of_run(analysis_path, backgrounds)
        run = str(analysis_path.parent)
        per_run[run], per_run_drivers[run], per_run_symptoms[run] = rules, drivers, symptom_sides
        totals.update(rules)
        driver_totals.update(drivers)
        symptom_totals.update(symptom_sides)
        if named:
            named_modules[run] = named

    modules = sum(totals.values())
    from_biology = totals["curated module"] + totals["subsystem"]
    lines = [
        "# Names for the learned modules, and how often a support carries one",
        "",
        "Generated by `experiments/measure_module_name_availability.py`. Do not edit by hand.",
        "",
        "A module is a learned gate, so its name is read off its support rather than assigned: the curated module it "
        f"overlaps by at least {CURATED_GENES_FOR_A_NAME} genes, else the subsystem holding at least "
        f"{SUBSYSTEM_SHARE_FOR_A_NAME:.0%} of the support and at least {SUBSYSTEM_ENRICHMENT_FOR_A_NAME:g} times its "
        "share of the graph (with the compartment when one holds half), else the node "
        f"type holding at least {NODE_TYPE_SHARE_FOR_A_NAME:.0%}, else \"mixed support\". The rule is `module_label` in "
        "`experiments/analyze_pathway_modules.py` and every name states the share it rests on. The index stays as the "
        "key, because it is what identifies a gate inside a split and what the other tables join on.",
        "",
        f"**{from_biology} of {modules} analysed module supports take a name from biology** (a curated module or a "
        f"subsystem); {totals['node type']} take only a node type and {totals['mixed']} are mixed. A mixed support is "
        "not a naming failure. It is a module with no identity to name, which is the same reading as the unused-module "
        "problem of `docs/module_health.md`, reached from the support rather than from the gate values.",
        "",
        "A second name comes from the run itself rather than from the support "
        "(`mechanistic_pathway_learning/evaluation/module_behaviour.py`): the held-out perturbations in the top fifth "
        "of the module's activation, whether they share a group, whether they reach the support through the graph more "
        "often than the perturbations that are not drivers, and the symptoms the module links to more strongly than "
        "the other modules do, grouped by the RDoC domain `docs/symptom_crosswalk.csv` already assigns each symptom. "
        f"**{driver_totals['group, path followed'] + driver_totals['no group, path followed']} of {modules} supports "
        f"have drivers that follow the graph** and {symptom_totals['domain']} carry a symptom domain of their own.",
        "",
        "| run | curated module | subsystem | node type | mixed | empty |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for run, rules in per_run.items():
        lines.append(f"| `{run}` | " + " | ".join(str(rules.get(rule, 0)) for rule in RULE_ORDER) + " |")
    lines.append("| **total** | " + " | ".join(f"**{totals.get(rule, 0)}**" for rule in RULE_ORDER) + " |")

    lines += ["", "## The drivers: which perturbations switch a module on", "",
              "`group` means the drivers share a curated module or a perturbation type that holds at least half of "
              "them and twice its share of the scored perturbations. `path followed` means the drivers have a "
              "directed path into the support, no longer than the propagation the encoder runs, at least 1.5 times as "
              "often as the perturbations that are not drivers; without that, the activation ordering is not coming "
              "through the graph and a driver name would be a correlation. `no driver` is a module whose activation "
              "varies by less than 0.05 across the held-out perturbations, the `responding_modules` reading of "
              "`docs/module_health.md`.", "",
              "| run | " + " | ".join(DRIVER_ORDER) + " |", "| --- |" + " --- |" * len(DRIVER_ORDER)]
    for run, drivers in per_run_drivers.items():
        lines.append(f"| `{run}` | " + " | ".join(str(drivers.get(rule, 0)) for rule in DRIVER_ORDER) + " |")
    lines.append("| **total** | " + " | ".join(f"**{driver_totals.get(rule, 0)}**" for rule in DRIVER_ORDER) + " |")

    lines += ["", "## The symptom side: what a module feeds that the others do not", "",
              "A module earns a symptom name only if its link vector is distinguishable from the other modules' "
              "(correlation with their mean below 0.7, the target `docs/module_health.md` sets) and at least one "
              "symptom's link stands 20 percent above what the other modules give that symptom. A domain name then "
              "needs two of those symptoms in one RDoC domain, half the weight and 1.5 times the domain's share of "
              "the run's symptoms.", "",
              "| run | " + " | ".join(SYMPTOM_ORDER) + " |", "| --- |" + " --- |" * len(SYMPTOM_ORDER)]
    for run, symptom_sides in per_run_symptoms.items():
        lines.append(f"| `{run}` | " + " | ".join(str(symptom_sides.get(rule, 0)) for rule in SYMPTOM_ORDER) + " |")
    lines.append("| **total** | " + " | ".join(f"**{symptom_totals.get(rule, 0)}**" for rule in SYMPTOM_ORDER) + " |")

    lines += ["", "## The modules a name fits", ""]
    if named_modules:
        lines += ["| run | split | module | support name | behaviour name |", "| --- | --- | --- | --- | --- |"]
        for run, named in named_modules.items():
            for split_name, module_name, label in named:
                lines.append(f"| `{run}` | {split_name} | {module_name} | {label} |")
    else:
        lines.append("None: no analysed support reaches a curated module or a subsystem threshold.")
    OUTPUT_DOCUMENT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUTPUT_DOCUMENT}")
    print(f"module supports: {modules}; named from biology: {from_biology}; node type only: {totals['node type']}; "
          f"mixed: {totals['mixed']}")
    print("drivers: " + ", ".join(f"{rule} {driver_totals.get(rule, 0)}" for rule in DRIVER_ORDER))
    print("symptom side: " + ", ".join(f"{rule} {symptom_totals.get(rule, 0)}" for rule in SYMPTOM_ORDER))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
