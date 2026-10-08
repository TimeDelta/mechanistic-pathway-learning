"""Build the OnSIDES later label slice as report rows and compare it with the SIDER drug pairs.

The comparison has two SIDER sides. The qualified side is the SIDER-label pairs of the assembled table
(data/processed/evidence_full, the 2015 pairs the time split trains on); a pair absent from it is
"absent from the qualified table". The raw side is every PT event of SIDER 4.1 (meddra_all_se.tsv.gz,
crosswalk symptoms) for any SIDER identifier linked to the perturbation (the unified CID and its
collisions), independent of SIDER's own ChEMBL mapping and qualification: a pair found there was on a
2015 label and is not a later-slice candidate even though the qualified table lacks it.

Usage:
  OMP_NUM_THREADS=1 python experiments/build_onsides_reports.py   # writes data/processed/onsides/onsides_reports.parquet, a summary and docs/onsides_label_slice.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.load_drug_label_events import load_sider_events, read_preferred_term_to_target_symptom
from mechanistic_pathway_learning.evidence.load_onsides_label_events import (
    ONSIDES_RELEASE,
    ONSIDES_RELEASE_DATE,
    aggregate_label_statements,
    load_onsides_tables,
    onsides_reports,
    sider_lower_level_to_preferred_terms,
)
from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import (
    IngredientBridge,
    load_ingredient_identifier_bridge,
    mechanisms_by_parent_and_molecule,
    stitch_flat_id_from_pubchem_cid,
)
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import (
    DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH,
    NON_PROTEIN_TARGETS_PATH,
    GraphNodeLookup,
    check_graph_compounds_have_no_mechanism,
    load_drugs_acting_as_graph_compounds,
)

SIDER_LABEL_SOURCE = "SIDER-label"


def qualified_sider_pairs(evidence_directory: Path) -> set[tuple[str, str]]:
    """(perturbation_id, symptom) of the SIDER-label reports in the assembled table; empty when the table is absent.

    The report table is read first since it names the source of every row; the record table of an older
    assembly is the fallback, with its drug rows that cite SIDER. OnSIDES rows appended to the same table are
    never counted as SIDER pairs.
    """
    reports_path, records_path = evidence_directory / "evidence_reports.parquet", evidence_directory / "evidence_records.parquet"
    if reports_path.exists():
        reports = pd.read_parquet(reports_path, columns=["perturbation_id", "symptom", "relation", "source"])
        sider_rows = reports[(reports.source == SIDER_LABEL_SOURCE) & (reports.relation == "induces")]
        return set(zip(sider_rows.perturbation_id, sider_rows.symptom))
    if records_path.exists():
        records = pd.read_parquet(records_path)
        drug_rows = records[(records.perturbation_type == "drug") & (records.relation == "induces") & records.source.astype(str).str.contains("SIDER")]
        return set(zip(drug_rows.perturbation_id, drug_rows.symptom))
    return set()


def sider_identifiers_by_perturbation(bridges: dict[str, IngredientBridge]) -> dict[str, set[str]]:
    """perturbation_id -> every STITCH flat identifier the bridge links to it (the unified CID and its collisions)."""
    identifiers: dict[str, set[str]] = {}
    for bridge in bridges.values():
        linked = {stitch_flat_id_from_pubchem_cid(cid) for cid in [bridge.sider_pubchem_cid, *bridge.sider_collision_pubchem_cids] if cid is not None}
        if linked:
            identifiers.setdefault(bridge.perturbation_id, set()).update(linked)
    return identifiers


def raw_sider_label_pairs(sider_directory: Path, crosswalk_path: Path) -> set[tuple[str, str]]:
    """(STITCH flat id, target symptom) of every SIDER 4.1 PT adverse-event row whose term is in the crosswalk, no qualification applied."""
    return {(event.stitch_flat_id, event.target_symptom) for event in load_sider_events(sider_directory, crosswalk_path) if event.relation == "induces" and event.target_symptom}


def pair_overlap(qualifying: pd.DataFrame, sider_pairs: set[tuple[str, str]], raw_sider_pairs: set[tuple[str, str]], bridges: dict[str, IngredientBridge]) -> dict:
    onsides_pairs = set(zip(qualifying.perturbation_id, qualifying.symptom))
    absent_from_qualified_table = onsides_pairs - sider_pairs
    linked_identifiers = sider_identifiers_by_perturbation(bridges)
    on_sider_labels = {(perturbation_id, symptom) for perturbation_id, symptom in absent_from_qualified_table if any((stitch_flat_id, symptom) in raw_sider_pairs for stitch_flat_id in linked_identifiers.get(perturbation_id, set()))}
    later_slice_candidates = absent_from_qualified_table - on_sider_labels
    unified_perturbations = {bridge.perturbation_id for bridge in bridges.values() if bridge.unified_with_sider}

    def by_symptom(pairs: set[tuple[str, str]]) -> dict[str, int]:
        return pd.Series([symptom for _, symptom in pairs]).value_counts().sort_index().to_dict() if pairs else {}

    return {
        "pairs_in_both": len(sider_pairs & onsides_pairs), "sider_only": len(sider_pairs - onsides_pairs), "sider_pairs": len(sider_pairs), "onsides_pairs": len(onsides_pairs),
        "onsides_absent_from_qualified_sider_table": len(absent_from_qualified_table),
        "onsides_absent_but_on_sider_labels": len(on_sider_labels),
        "onsides_later_slice_candidates": len(later_slice_candidates),
        "onsides_later_slice_candidates_on_unified_perturbations": sum(1 for perturbation_id, _ in later_slice_candidates if perturbation_id in unified_perturbations),
        "onsides_later_slice_candidates_on_onsides_only_perturbations": sum(1 for perturbation_id, _ in later_slice_candidates if perturbation_id not in unified_perturbations),
        "raw_sider_pairs": len(raw_sider_pairs),
        "onsides_only": len(later_slice_candidates),  # kept under the earlier key for readers of the summary: the later-slice candidates
        "onsides_absent_by_symptom": by_symptom(absent_from_qualified_table), "onsides_on_sider_labels_by_symptom": by_symptom(on_sider_labels), "onsides_only_by_symptom": by_symptom(later_slice_candidates),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--onsides-parquet-dir", type=Path, default=Path("data/raw/onsides/v3.1.1/parquet"))
    parser.add_argument("--bridge", type=Path, default=Path("data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--sider-dir", type=Path, default=Path("data/raw/sider_4.1"), help="SIDER 4.1 tables: the LLT-to-PT pairing and the raw 2015 label events of the overlap")
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full"))
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full"), help="assembled table whose SIDER-label pairs define the qualified overlap")
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/onsides"))
    parser.add_argument("--non-protein-targets", type=Path, default=NON_PROTEIN_TARGETS_PATH, help="ChEMBL non-protein targets -> Human-GEM metabolites (docs/drug_targets_any_type.md)")
    parser.add_argument("--drugs-acting-as-graph-compounds", type=Path, default=DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH, help="drugs with no mechanism target that are themselves a graph metabolite")
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/onsides_label_slice.md"))
    arguments = parser.parse_args()

    tables = load_onsides_tables(arguments.onsides_parquet_dir)
    preferred_term_to_symptom = read_preferred_term_to_target_symptom(arguments.crosswalk)
    lower_level_to_preferred = sider_lower_level_to_preferred_terms(arguments.sider_dir) if (arguments.sider_dir / "meddra_all_se.tsv.gz").exists() else None
    statements = aggregate_label_statements(tables, preferred_term_to_symptom, lower_level_to_preferred)
    bridges = load_ingredient_identifier_bridge(arguments.bridge)
    mechanisms_by_molecule = mechanisms_by_parent_and_molecule(json.loads((arguments.chembl_dir / "mechanisms.json").read_text()))
    targets = json.loads((arguments.chembl_dir / "targets.json").read_text())
    node_lookup = GraphNodeLookup.from_nodes(pd.read_parquet(arguments.graph_dir / "nodes.parquet"), arguments.non_protein_targets)
    _, compounds_by_chembl_id = load_drugs_acting_as_graph_compounds(arguments.drugs_acting_as_graph_compounds)
    check_graph_compounds_have_no_mechanism(compounds_by_chembl_id, mechanisms_by_molecule)
    reports, counts = onsides_reports(statements, bridges, mechanisms_by_molecule, targets, node_lookup, compounds_by_chembl_id)

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    reports.to_parquet(arguments.output_dir / "onsides_reports.parquet", index=False)
    qualifying = reports[reports.qualifies]
    sider_pairs = qualified_sider_pairs(arguments.evidence_dir)
    raw_sider_pairs = raw_sider_label_pairs(arguments.sider_dir, arguments.crosswalk) if (arguments.sider_dir / "meddra_all_se.tsv.gz").exists() else set()
    overlap = pair_overlap(qualifying, sider_pairs, raw_sider_pairs, bridges)
    statements_by_term_type = statements.meddra_term_type.value_counts().to_dict()
    summary = {"release": ONSIDES_RELEASE, "release_date": ONSIDES_RELEASE_DATE, "counts": counts, "overlap": overlap,
               "statements_by_term_type": statements_by_term_type, "lower_level_terms_paired_through_sider": lower_level_to_preferred is not None,
               "reports_by_symptom": qualifying.symptom.value_counts().to_dict(), "reports_by_country": qualifying.source_country.value_counts().to_dict(),
               "reports_by_section": qualifying.label_section.value_counts().to_dict(), "reports_by_term_type": qualifying.meddra_term_type.value_counts().to_dict(),
               "reports_by_unification_method": {method: int(count) for method, count in qualifying.ingredient_id.map(lambda ingredient_id: bridges[ingredient_id].unification_method).value_counts().items()},
               "qualifying_drugs_unified_with_sider": int(qualifying[qualifying.unified_with_sider].perturbation_id.nunique()),
               "qualifying_drugs_onsides_only": int(qualifying[~qualifying.unified_with_sider].perturbation_id.nunique())}
    (arguments.output_dir / "onsides_summary.json").write_text(json.dumps(summary, indent=1))

    lines = [f"# OnSIDES {ONSIDES_RELEASE} later label slice (generated by experiments/build_onsides_reports.py)", "",
             f"Data release {ONSIDES_RELEASE_DATE} (CC BY 4.0). One report per deduplicated label statement (ingredient, target symptom, MedDRA term, label section, country); "
             "the bridge routes, the report schema, the qualification rule and the label-section semantics are specified in docs/onsides_label_slice_spec.md. "
             f"Statements: {counts['statements']} ({statements_by_term_type.get('PT', 0)} PT-coded, {statements_by_term_type.get('LLT', 0)} LLT-coded and matched through the PT their LLT falls under); "
             f"qualifying reports (the SIDER rule: single mechanism target and an ATC N code; target gene in the graph; at least one single-ingredient label): {counts['qualifying_reports']} "
             f"from {counts['qualifying_ingredients']} ingredients and {counts['qualifying_perturbations']} perturbations ({summary['qualifying_drugs_unified_with_sider']} shared with SIDER, {summary['qualifying_drugs_onsides_only']} OnSIDES-only). "
             f"Disqualified statements by reason: {counts['disqualified']}; statements of ingredients without a bridge entry: {counts['ingredients_without_bridge']}.", "",
             "Label sections exist for US labels only (AR, BW, WP); every UK, EU and JP event carries label_section NA, OnSIDES' value for the undivided non-US label, "
             "coded non_us_label_section and counted as adverse-reactions text by rubric_adverse_reactions_section.", "",
             "MedDRA term names are not reproduced here; identifiers and target symptoms only.", "",
             "## Pair-level overlap with the SIDER drug pairs", "",
             f"Against the qualified SIDER-label pairs of the assembled table (relation induces, the pairs of SIDER drugs that passed SIDER's ChEMBL mapping, single-target and ATC N rules): "
             f"pairs in both {overlap['pairs_in_both']}; SIDER only (on 2015 labels, absent from current labels or from the OnSIDES extraction) {overlap['sider_only']}; "
             f"OnSIDES pairs absent from the qualified table {overlap['onsides_absent_from_qualified_sider_table']}. "
             f"Of those, {overlap['onsides_absent_but_on_sider_labels']} are on 2015 SIDER labels (a PT event of SIDER 4.1 for a SIDER identifier linked to the perturbation, independent of SIDER's qualification) "
             f"and {overlap['onsides_later_slice_candidates']} are candidates for the later slice ({overlap['onsides_later_slice_candidates_on_unified_perturbations']} on perturbations SIDER also holds, "
             f"{overlap['onsides_later_slice_candidates_on_onsides_only_perturbations']} on OnSIDES-only perturbations); a candidate may still have been on a 2015 label that SIDER's extraction missed.", "",
             "| symptom | absent from qualified table | on 2015 SIDER labels | later-slice candidates |", "|---|---|---|---|"]
    lines += [f"| {symptom} | {count} | {overlap['onsides_on_sider_labels_by_symptom'].get(symptom, 0)} | {overlap['onsides_only_by_symptom'].get(symptom, 0)} |" for symptom, count in sorted(overlap["onsides_absent_by_symptom"].items())]
    lines += ["", "## Qualifying reports", "", "| symptom | reports |", "|---|---|"] + [f"| {s} | {c} |" for s, c in sorted(summary["reports_by_symptom"].items())]
    lines += ["", "| country | reports |", "|---|---|"] + [f"| {s} | {c} |" for s, c in sorted(summary["reports_by_country"].items())]
    lines += ["", "| label section | reports |", "|---|---|"] + [f"| {s} | {c} |" for s, c in sorted(summary["reports_by_section"].items())]
    lines += ["", "| source term type | reports |", "|---|---|"] + [f"| {s} | {c} |" for s, c in sorted(summary["reports_by_term_type"].items())]
    lines += ["", "| unification with SIDER | reports |", "|---|---|"] + [f"| {s} | {c} |" for s, c in sorted(summary["reports_by_unification_method"].items())]
    lines += ["", "## Design document changes to apply", "",
              "4.2 (E2, both slices): SIDER 4.1 supplies the 2015 label slice and OnSIDES v3.1.1 the current-label slice; both are grade B label evidence under one qualification rule (single mechanism target, ATC N code), OnSIDES without frequencies but with the label section, the extraction score, the label count and the country as rubric features; ingredients are unified with SIDER drugs through the ChEMBL parents of SIDER's stereo-flattened identifiers or by drug name, so both sources report on one perturbation.",
              f"6.1 (pharmacological time split): SIDER 2015 pairs train; the {overlap['onsides_later_slice_candidates']} later-slice candidates test (OnSIDES pairs absent from the qualified SIDER table and from the raw SIDER 4.1 label events of the linked identifiers), with the caveat that a candidate may have existed on a 2015 label that SIDER missed; the {overlap['onsides_absent_but_on_sider_labels']} pairs absent from the qualified table but on 2015 labels are not test positives; OnSIDES rows are dated by the release (an upper bound).",
              "10 (layout): data/processed/onsides/onsides_reports.parquet; mechanistic_pathway_learning/evidence/load_onsides_label_events.py, onsides_identifier_bridge.py; experiments/fetch_onsides_identifier_bridge.py, build_onsides_reports.py; docs/onsides_label_slice_spec.md.",
              "12 (status): the counts above."]
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(json.dumps({"counts": counts, "overlap": {k: v for k, v in overlap.items() if not k.endswith("_by_symptom")}, "statements_by_term_type": statements_by_term_type}, indent=0))


if __name__ == "__main__":
    main()
