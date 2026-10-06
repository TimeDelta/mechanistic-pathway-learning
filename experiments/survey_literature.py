"""Join the filtered literature relations with their documents and the rule-based appraisal into the literature
report table, and write the survey document with the overlap against the assembled evidence.

Usage:
  python experiments/survey_literature.py
Reads data/processed/literature/{pubtator_relations_filtered,ctd_relations_filtered,documents}.parquet; writes
literature_reports.parquet (evidence_reports columns first, then the literature extras, study_design, species and the
rubric columns), literature_summary.json and docs/literature_survey.md. The loader's evidence_code (the extraction
type) is kept; the appraised study design goes to its own column and is prefixed to model_description.

Literature reports are grade E soft priors (design section 4.2). Under the design they would enter training only for
training-fold perturbations and never as evaluation positives; in the current build the experiment loader keeps grades
A and B as labels (and refuses D and E), so literature rows enter the report table and the reliability fit only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS
from mechanistic_pathway_learning.evidence.rule_based_evidence_appraisal import appraise_from_metadata, rubric_features, rubric_weight_from_metadata

RUBRIC_COLUMNS = ("rubric_human_species", "rubric_animal_only", "rubric_randomized_trial", "rubric_case_report", "rubric_cohort_or_case_control",
                  "rubric_genetic_association", "rubric_review_or_secondary", "rubric_symptom_level_descriptor", "rubric_publication_year_known", "rubric_metadata_available")
SOFT_PRIOR_SENTENCE = ("Literature reports are grade E soft priors: under the design they would enter training only for training-fold perturbations and never as "
                       "evaluation positives; in the current build the experiment loader keeps grades A and B as labels and refuses D and E, so literature rows enter the "
                       "report table and the reliability fit only, and the fold-aware exclusion is the rule to implement when the soft prior is wired into training.")


def appraise_reports(reports: pd.DataFrame, documents: pd.DataFrame) -> pd.DataFrame:
    """Attach study_design, species, rubric features, rubric weight and limitations to every report from its document's
    metadata. evidence_code (the loader's extraction type) is kept; a report whose PMID has no document gets the design
    metadata_not_fetched, distinct from a fetched paper whose metadata say nothing (not_reported)."""
    document_by_pmid = documents.set_index(documents.pmid.astype(str)) if len(documents) else pd.DataFrame()
    appraised = reports.copy()
    designs, species_lists, weights, limitations, features, descriptions = [], [], [], [], [], []
    for row in appraised.itertuples(index=False):
        pmid = str(getattr(row, "pmid", ""))
        document = document_by_pmid.loc[pmid] if len(document_by_pmid) and pmid in document_by_pmid.index else None
        publication_types = [value for value in str(document.publication_types).split(";") if value] if document is not None and isinstance(document.publication_types, str) else []
        species_ids = [value for value in str(document.species_ids).split(";") if value] if document is not None and isinstance(document.species_ids, str) else []
        species_names = [value for value in str(document.species_names).split(";") if value] if document is not None and isinstance(document.species_names, str) else []
        rubric = appraise_from_metadata(publication_types, species_ids, species_names, row.perturbation_type, str(getattr(row, "mesh_descriptor_level", "")), bool(row.evidence_date), metadata_available=document is not None)
        designs.append(rubric.study_design)
        species_lists.append(";".join(rubric.species))
        weights.append(rubric_weight_from_metadata(rubric))
        existing = str(row.limitations or "")
        limitations.append((existing + "; " if existing else "") + rubric.limitations)
        features.append(rubric_features(rubric))
        design_words = rubric.study_design.replace("_", " ")
        subjects = "human" if rubric.study_design.startswith("human") else "non-human" if rubric.study_design.startswith("animal") else "subjects unspecified"
        appraisal_text = "; ".join(part for part in (subjects, design_words, ";".join(rubric.species)) if part)
        descriptions.append(f"{row.model_description}; appraisal: {appraisal_text}" if isinstance(row.model_description, str) and row.model_description else appraisal_text)
    appraised["study_design"] = designs
    appraised["model_description"] = descriptions
    appraised["species"] = species_lists
    appraised["rubric_weight"] = weights
    appraised["limitations"] = limitations
    for column in RUBRIC_COLUMNS:
        appraised[column] = [entry[column] for entry in features]
    return appraised


def overlap_with_evidence(reports: pd.DataFrame, evidence_records_path: Path) -> dict:
    if not evidence_records_path.exists():
        return {}
    records = pd.read_parquet(evidence_records_path)
    labelled = records[records.grade.isin(["A", "B"])] if "grade" in records else records
    if "positive_report_count" in labelled.columns:
        labelled = labelled[labelled.positive_report_count > 0]
    evidence_pairs = {(perturbation, symptom, relation) for perturbation, symptom, relation in zip(labelled.perturbation_id, labelled.symptom, labelled.relation)}
    literature_pairs = {(perturbation, symptom, relation) for perturbation, symptom, relation in zip(reports.perturbation_id, reports.symptom, reports.relation)}
    perturbation_types = dict(zip(labelled.perturbation_id, labelled.perturbation_type)) if "perturbation_type" in labelled.columns else {}
    gene_pairs = {pair for pair in evidence_pairs if perturbation_types.get(pair[0], "gene") == "gene"}
    drug_pairs = evidence_pairs - gene_pairs
    return {
        "evidence_table": str(evidence_records_path),
        "evidence_pairs": len(evidence_pairs), "literature_pairs": len(literature_pairs),
        "evidence_pairs_with_literature": len(evidence_pairs & literature_pairs),
        "gene_pairs_with_literature": len(gene_pairs & literature_pairs), "gene_pairs": len(gene_pairs),
        "drug_pairs_with_literature": len(drug_pairs & literature_pairs), "drug_pairs": len(drug_pairs),
        "literature_pairs_new": len(literature_pairs - evidence_pairs),
        "gene_pairs_with_literature_any_relation": len({pair[:2] for pair in gene_pairs} & {pair[:2] for pair in literature_pairs}),
        "drug_pairs_with_literature_any_relation": len({pair[:2] for pair in drug_pairs} & {pair[:2] for pair in literature_pairs}),
        "note": "pairs split by the perturbation_type column of the assembled table (grades A and B with a positive report); same-relation overlap for genes counts the induces reports from Disease-Gene stimulate and inhibit rows, the any-relation overlap also counts association reports",
    }


def descriptor_level_split_by_symptom(appraised: pd.DataFrame) -> list[dict]:
    """Per symptom: reports, the symptom-level and diagnosis-level shares and the directional statements demoted on a mixed-polarity descriptor."""
    rows = []
    demoted_types = {"cause", "treat", "prevent", "stimulate", "inhibit", "marker/mechanism", "therapeutic"}
    for symptom, group in appraised.groupby("symptom"):
        level = group.mesh_descriptor_level.astype(str)
        demoted = int(((group.relation == "associated_with") & group.relation_type.astype(str).isin(demoted_types)).sum())
        rows.append({"symptom": symptom, "reports": int(len(group)), "symptom_level": int((level == "symptom").sum()), "diagnosis_level": int(level.isin(["diagnosis", "mixed_polarity_diagnosis"]).sum()),
                     "mixed_polarity_descriptor": int((level == "mixed_polarity_diagnosis").sum()), "directional_demoted_to_associated_with": demoted})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--literature-dir", type=Path, default=Path("data/processed/literature"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/literature_survey.md"))
    arguments = parser.parse_args()
    tables = [pd.read_parquet(path) for path in (arguments.literature_dir / "pubtator_relations_filtered.parquet", arguments.literature_dir / "ctd_relations_filtered.parquet") if path.exists()]
    reports = pd.concat(tables, ignore_index=True, sort=False)
    documents_path = arguments.literature_dir / "documents.parquet"
    documents = pd.read_parquet(documents_path) if documents_path.exists() else pd.DataFrame(columns=["pmid", "publication_types", "species_ids", "species_names"])
    for column in REPORT_COLUMNS:  # relation tables written before a shared column existed get its default
        if column not in reports.columns:
            reports[column] = 0.0 if column.startswith("rubric_") else ""
    reports["evidence_date_source"] = ["publication" if isinstance(value, str) and value else "" for value in reports.evidence_date]  # literature dates are publication dates
    appraised = appraise_reports(reports, documents)
    ordered = [column for column in REPORT_COLUMNS] + [column for column in appraised.columns if column not in REPORT_COLUMNS]
    appraised = appraised[ordered]
    appraised.to_parquet(arguments.literature_dir / "literature_reports.parquet", index=False)
    overlap = overlap_with_evidence(appraised, arguments.evidence_dir / "evidence_records.parquet")
    years = pd.to_datetime(appraised.evidence_date, errors="coerce").dt.year
    documents_summary_path = arguments.literature_dir / "documents_summary.json"
    documents_summary = json.loads(documents_summary_path.read_text()) if documents_summary_path.exists() else {}
    level_split = descriptor_level_split_by_symptom(appraised)
    summary = {
        "reports": int(len(appraised)), "distinct_pmids": int(appraised.pmid.nunique()) if "pmid" in appraised else None,
        "distinct_pairs": int(appraised[["perturbation_id", "symptom", "relation"]].drop_duplicates().shape[0]),
        "by_source": appraised.source.value_counts().to_dict(), "by_relation": appraised.relation.value_counts().to_dict(),
        "by_symptom": appraised.symptom.value_counts().to_dict(), "by_study_design": appraised.study_design.value_counts().to_dict(),
        "by_evidence_code": appraised.evidence_code.value_counts().to_dict(),
        "by_perturbation_type": appraised.perturbation_type.value_counts().to_dict(),
        "gene_reports_by_relation": appraised.loc[appraised.perturbation_type == "gene", "relation"].value_counts().to_dict(),
        "documents_with_metadata": int(appraised.pmid.astype(str).isin(documents.pmid.astype(str)).sum()) if len(documents) and "pmid" in appraised else 0,
        "reports_without_document_metadata": int((~appraised.pmid.astype(str).isin(documents.pmid.astype(str))).sum()) if len(documents) and "pmid" in appraised else int(len(appraised)),
        "document_fetch": documents_summary,
        "descriptor_level_split_by_symptom": level_split,
        "year_histogram": {str(int(year)): int(count) for year, count in years.value_counts().sort_index().items() if pd.notna(year)},
        "rubric_weight_quantiles": {str(q): float(appraised.rubric_weight.quantile(q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
        "overlap": overlap,
    }
    (arguments.literature_dir / "literature_summary.json").write_text(json.dumps(summary, indent=1))
    sources_text = ", ".join(f"{source} {count}" for source, count in sorted(summary["by_source"].items()))
    fetched, unfetched = documents_summary.get("pmids_fetched"), documents_summary.get("pmids_unfetched")
    cap_text = (f"Document metadata (species, publication types) was fetched for {fetched} of {(fetched or 0) + (unfetched or 0)} papers (cap --max-pmids {documents_summary.get('max_pmids', documents_summary.get('pmids_selected_by_cap'))}, "
                f"papers with a directional relation first); {unfetched} papers have no metadata, so the {summary['reports_without_document_metadata']} reports on them carry the design metadata_not_fetched by construction, "
                f"distinct from not_reported (metadata read, no design inferable)." if fetched is not None else f"{summary['documents_with_metadata']} reports have document metadata; {summary['reports_without_document_metadata']} have none.")
    lines = ["# Literature survey, evidence class E3 (generated by experiments/survey_literature.py)", "",
             f"Reports: {summary['reports']} ({sources_text}) over {summary['distinct_pmids']} papers and {summary['distinct_pairs']} (perturbation, symptom, relation) triples. {cap_text} {SOFT_PRIOR_SENTENCE}", "",
             "## Reports by relation and by study design (rule-based appraisal from publication types and annotated species)", "",
             "| relation | reports |", "|---|---|"] + [f"| {k} | {v} |" for k, v in summary["by_relation"].items()]
    lines += ["", "| gene reports by relation | reports |", "|---|---|"] + [f"| {k} | {v} |" for k, v in summary["gene_reports_by_relation"].items()]
    lines += ["", "| study design (study_design column) | reports |", "|---|---|"] + [f"| {k} | {v} |" for k, v in summary["by_study_design"].items()]
    lines += ["", "| extraction type (evidence_code) | reports |", "|---|---|"] + [f"| {k} | {v} |" for k, v in sorted(summary["by_evidence_code"].items())]
    lines += ["", "| symptom | reports | symptom-level descriptor | diagnosis-level descriptor | mixed-polarity descriptor | directional statements demoted to associated_with |", "|---|---|---|---|---|---|"]
    lines += [f"| {row['symptom']} | {row['reports']} | {row['symptom_level']} | {row['diagnosis_level']} | {row['mixed_polarity_descriptor']} | {row['directional_demoted_to_associated_with']} |" for row in level_split]
    lines += ["", "A directional statement (cause, treat, prevent, stimulate, inhibit, marker/mechanism, therapeutic) on the mixed-polarity descriptor Bipolar Disorder D001714 is recorded as associated_with, because the diagnosis has the opposite mood pole as a cardinal feature and a treat or cause relation on it carries no sign for mania. Depressed mood has no symptom-level channel in PubTator3 or CTD (Depression D003863 is outside their disease vocabulary), so every depressed_mood report is diagnosis-level."]
    lines += ["", "## Overlap with the assembled evidence (grades A and B, positive report)", "", "| quantity | count |", "|---|---|"] + [f"| {k} | {v} |" for k, v in overlap.items()]
    decades = {}
    for year, count in summary["year_histogram"].items():
        decades[str(int(year) // 10 * 10)] = decades.get(str(int(year) // 10 * 10), 0) + count
    lines += ["", "## Publication years", "", "| decade | reports |", "|---|---|"] + [f"| {k}s | {v} |" for k, v in sorted(decades.items())]
    quantile_text = ", ".join(f"{q}: {value:.3f}" for q, value in summary["rubric_weight_quantiles"].items())
    lines += ["", "## Rubric weights", "", f"Quantiles ({quantile_text}). Weights: randomised trial 1.0, cohort or case-control 0.9, genetic association 0.8, case report 0.6, animal 0.5, review or secondary 0.4, not reported (metadata read, no design inferable) 0.35, metadata not fetched 0.3; times 0.7 for a diagnosis-level descriptor. Species decide before publication types: a paper whose only annotated species is non-human is an animal study whatever its publication types."]
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("by_symptom", "year_histogram")}, indent=0)[:2000])


if __name__ == "__main__":
    main()
