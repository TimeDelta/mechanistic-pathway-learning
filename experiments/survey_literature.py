"""Join the filtered literature relations with their documents and the rule-based appraisal into the literature
report table, and write the survey document with the overlap against the assembled evidence.

Usage:
  python experiments/survey_literature.py
Reads data/processed/literature/{pubtator_relations_filtered,ctd_relations_filtered,documents}.parquet; writes
literature_reports.parquet (evidence_reports columns first, then the literature extras and the rubric columns),
literature_summary.json and docs/literature_survey.md. Literature reports are soft priors (grades C to E), enter
training for training-fold perturbations only and are never evaluation positives (design section 4.2).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS
from mechanistic_pathway_learning.evidence.rule_based_evidence_appraisal import appraise_from_metadata, rubric_features, rubric_weight_from_metadata

RUBRIC_COLUMNS = ("rubric_human_species", "rubric_animal_only", "rubric_randomized_trial", "rubric_case_report", "rubric_cohort_or_case_control",
                  "rubric_genetic_association", "rubric_review_or_secondary", "rubric_symptom_level_descriptor", "rubric_publication_year_known")


def appraise_reports(reports: pd.DataFrame, documents: pd.DataFrame) -> pd.DataFrame:
    """Attach study design, species, rubric features, rubric weight and limitations to every report from its document's metadata."""
    document_by_pmid = documents.set_index(documents.pmid.astype(str)) if len(documents) else pd.DataFrame()
    appraised = reports.copy()
    designs, species_lists, weights, limitations, features = [], [], [], [], []
    for row in appraised.itertuples(index=False):
        pmid = str(getattr(row, "pmid", ""))
        document = document_by_pmid.loc[pmid] if len(document_by_pmid) and pmid in document_by_pmid.index else None
        publication_types = [value for value in str(document.publication_types).split(";") if value] if document is not None and isinstance(document.publication_types, str) else []
        species_ids = [value for value in str(document.species_ids).split(";") if value] if document is not None and isinstance(document.species_ids, str) else []
        species_names = [value for value in str(document.species_names).split(";") if value] if document is not None and isinstance(document.species_names, str) else []
        rubric = appraise_from_metadata(publication_types, species_ids, species_names, row.perturbation_type, str(getattr(row, "mesh_descriptor_level", "")), bool(row.evidence_date))
        designs.append(rubric.study_design); species_lists.append(";".join(rubric.species)); weights.append(rubric_weight_from_metadata(rubric))
        existing = str(row.limitations or "")
        limitations.append((existing + "; " if existing else "") + rubric.limitations)
        features.append(rubric_features(rubric))
    appraised["evidence_code"] = designs
    appraised["model_description"] = ["; ".join(part for part in (("human" if design.startswith("human") else "non-human" if design.startswith("animal") else "unspecified subjects"), design.replace("_", " "), species) if part) for design, species in zip(designs, species_lists)]
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
    evidence_pairs = {(perturbation, symptom, relation) for perturbation, symptom, relation in zip(labelled.perturbation_id, labelled.symptom, labelled.relation)}
    literature_pairs = {(perturbation, symptom, relation) for perturbation, symptom, relation in zip(reports.perturbation_id, reports.symptom, reports.relation)}
    gene_pairs = {pair for pair in evidence_pairs if not str(pair[0]).startswith("CID")}
    drug_pairs = evidence_pairs - gene_pairs
    return {
        "evidence_pairs": len(evidence_pairs), "literature_pairs": len(literature_pairs),
        "evidence_pairs_with_literature": len(evidence_pairs & literature_pairs),
        "gene_pairs_with_literature": len(gene_pairs & literature_pairs), "gene_pairs": len(gene_pairs),
        "drug_pairs_with_literature": len(drug_pairs & literature_pairs), "drug_pairs": len(drug_pairs),
        "literature_pairs_new": len(literature_pairs - evidence_pairs),
    }


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
    appraised = appraise_reports(reports, documents)
    ordered = [column for column in REPORT_COLUMNS] + [column for column in appraised.columns if column not in REPORT_COLUMNS]
    appraised = appraised[ordered]
    appraised.to_parquet(arguments.literature_dir / "literature_reports.parquet", index=False)
    overlap = overlap_with_evidence(appraised, arguments.evidence_dir / "evidence_records.parquet")
    years = pd.to_datetime(appraised.evidence_date, errors="coerce").dt.year
    summary = {
        "reports": int(len(appraised)), "distinct_pmids": int(appraised.pmid.nunique()) if "pmid" in appraised else None,
        "distinct_pairs": int(appraised[["perturbation_id", "symptom", "relation"]].drop_duplicates().shape[0]),
        "by_source": appraised.source.value_counts().to_dict(), "by_relation": appraised.relation.value_counts().to_dict(),
        "by_symptom": appraised.symptom.value_counts().to_dict(), "by_study_design": appraised.evidence_code.value_counts().to_dict(),
        "by_perturbation_type": appraised.perturbation_type.value_counts().to_dict(),
        "documents_with_metadata": int(appraised.pmid.astype(str).isin(documents.pmid.astype(str)).sum()) if len(documents) and "pmid" in appraised else 0,
        "year_histogram": {str(int(year)): int(count) for year, count in years.value_counts().sort_index().items() if pd.notna(year)},
        "rubric_weight_quantiles": {str(q): float(appraised.rubric_weight.quantile(q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
        "overlap": overlap,
    }
    (arguments.literature_dir / "literature_summary.json").write_text(json.dumps(summary, indent=1))
    lines = ["# Literature survey, evidence class E3 (generated by experiments/survey_literature.py)", "",
             f"Reports: {summary['reports']} ({summary['by_source']}) over {summary['distinct_pmids']} papers and {summary['distinct_pairs']} (perturbation, symptom, relation) triples; "
             f"{summary['documents_with_metadata']} reports have document metadata (species, publication types). Literature reports are soft priors (grades C to E): they enter training for training-fold perturbations only, are excluded for any perturbation in the evaluation fold (a report on a held-out gene's symptom would leak its label) and are never evaluation positives.", "",
             "## Reports by relation and by study design (rule-based appraisal from publication types and annotated species)", "",
             "| relation | reports |", "|---|---|"] + [f"| {k} | {v} |" for k, v in summary["by_relation"].items()]
    lines += ["", "| study design (evidence_code) | reports |", "|---|---|"] + [f"| {k} | {v} |" for k, v in summary["by_study_design"].items()]
    lines += ["", "| symptom | reports |", "|---|---|"] + [f"| {k} | {v} |" for k, v in sorted(summary["by_symptom"].items())]
    lines += ["", "## Overlap with the assembled evidence (grades A and B)", "", "| quantity | count |", "|---|---|"] + [f"| {k} | {v} |" for k, v in overlap.items()]
    decades = {}
    for year, count in summary["year_histogram"].items():
        decades[str(int(year) // 10 * 10)] = decades.get(str(int(year) // 10 * 10), 0) + count
    lines += ["", "## Publication years", "", "| decade | reports |", "|---|---|"] + [f"| {k}s | {v} |" for k, v in sorted(decades.items())]
    lines += ["", "## Rubric weights", "", f"Quantiles (0.1, 0.25, 0.5, 0.75, 0.9): {summary['rubric_weight_quantiles']}. Weights: randomised trial 1.0, cohort or case-control 0.9, genetic association 0.8, case report 0.6, animal 0.5, review or secondary 0.3, not reported 0.5; times 0.7 for a diagnosis-level descriptor."]
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("by_symptom", "year_histogram")}, indent=0)[:2000])


if __name__ == "__main__":
    main()
