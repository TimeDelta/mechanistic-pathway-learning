"""Map the human-evidence review of the 58 rodent-readout claims onto the study's symptoms, and write the record.

The Open Targets target safety curation records 58 (symptom, direction, target) claims only through rodent behavioural
readouts (locomotor activity, catalepsy, stereotypy), which are not symptoms people report. A sub-agent searched the
human literature for each claim; docs/rodent_readout_human_evidence_prompt.txt holds the prompt it was given verbatim
and docs/rodent_readout_human_evidence.json holds what it returned. This script does three things:

  1. checks the returned references against PubMed metadata (every cited article exists, the DOI matches the PMID and
     each quote is verbatim in the abstract), writing the outcome to docs/rodent_readout_human_evidence_checks.json;
  2. maps each kept or relabelled claim onto the study's symptoms with the same regular expressions that read the
     curation itself (experiments/scope_open_targets_symptom_labels.py), so the mapping is not done by hand;
  3. writes docs/rodent_readout_human_evidence.md.

Pass --pubmed-metadata to redo step 1 against a file of PubMed records (the tool result of get_article_metadata, which
carries the abstracts); without it the committed check outcome is reused and no abstract is read or stored.
"""
import argparse
import importlib.util
import json
import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_FILE = REPOSITORY_ROOT / "docs/rodent_readout_human_evidence.json"
REFERENCE_FILE = REPOSITORY_ROOT / "docs/rodent_readout_human_evidence_references.json"
CHECKS_FILE = REPOSITORY_ROOT / "docs/rodent_readout_human_evidence_checks.json"
PROMPT_FILE = REPOSITORY_ROOT / "docs/rodent_readout_human_evidence_prompt.txt"
DOCUMENT_FILE = REPOSITORY_ROOT / "docs/rodent_readout_human_evidence.md"
SCOPING_SCRIPT = REPOSITORY_ROOT / "experiments/scope_open_targets_symptom_labels.py"
# The sub-agent's verdict sentences name the symptom the human data support and then, in several cases, the symptom they
# do not support ("somnolence and fatigue, not psychomotor retardation"). Matching the study's patterns against the
# whole sentence would read the negated tail as support, so the tail is cut off first.
NEGATED_TAIL = re.compile(r",?\s+not\s+(the\s+)?\b")


def symptom_patterns() -> dict:
    """The study's symptom regular expressions, read from the script that maps the curation onto symptoms."""
    specification = importlib.util.spec_from_file_location("scope_open_targets_symptom_labels", SCOPING_SCRIPT)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module.SYMPTOM_PATTERNS, module.ANIMAL_READOUT_PATTERN, module.SPECIES_AMBIGUOUS_READOUT_PATTERN


def study_symptoms_named(free_text: str, patterns: dict) -> list:
    """Which of the study's symptoms the human data name, by the study's own patterns, with the negated tail cut off."""
    if not free_text:
        return []
    supported = NEGATED_TAIL.split(free_text, maxsplit=1)[0].lower()
    return sorted(symptom for symptom, pattern in patterns.items() if re.search(pattern, supported))


def check_references(claims: list, reference_records: dict, pubmed_metadata_file: Path | None) -> dict:
    """Every cited article is in the reference records, its DOI agrees, and its quote is verbatim in the abstract."""
    citations = [(claim, reference) for claim in claims for reference in claim["references"]]
    outcome = {
        "num_claims": len(claims),
        "num_citations": len(citations),
        "num_unique_articles": len({reference["pmid"] for _, reference in citations}),
        "articles_without_a_metadata_record": sorted(
            {reference["pmid"] for _, reference in citations if reference["pmid"] not in reference_records}
        ),
        "citations_whose_doi_disagrees_with_the_record": [
            [reference["pmid"], reference["doi"], reference_records[reference["pmid"]]["doi"]]
            for _, reference in citations
            if reference["pmid"] in reference_records
            and (reference["doi"] or "").lower() != (reference_records[reference["pmid"]]["doi"] or "").lower()
        ],
        "references_missing_a_field": [
            [claim["target"], claim["symptom"], field]
            for claim, reference in citations
            for field in ("pmid", "doi", "quote")
            if not reference.get(field)
        ],
    }
    if pubmed_metadata_file is None:
        previous = json.loads(CHECKS_FILE.read_text()) if CHECKS_FILE.exists() else {}
        outcome["quote_check"] = previous.get("quote_check", {"ran": False})
        return outcome
    # PubMed renders Greek letters and the middle dot as characters a quote typed as plain text does not carry, so the
    # comparison folds whitespace, case and those two substitutions; nothing else about a quote is allowed to differ.
    records = json.loads(pubmed_metadata_file.read_text())
    abstracts = {
        str((article.get("identifiers") or {}).get("pmid")): article.get("abstract") or ""
        for article in (records.get("articles") if isinstance(records, dict) else records.values())
    } if not isinstance(records, dict) or "articles" in records else {
        pmid: (article.get("abstract") or "") for pmid, article in records.items()
    }

    def folded(text: str) -> str:
        return re.sub(r"\s+", " ", text).lower().replace("β", "beta").replace("·", ".")

    quotes_not_found = []
    for claim, reference in citations:
        abstract = folded(abstracts.get(reference["pmid"], ""))
        fragments = [fragment.strip() for fragment in reference["quote"].split("...") if fragment.strip()]
        if not all(folded(fragment) in abstract for fragment in fragments):
            quotes_not_found.append([reference["pmid"], reference["quote"]])
    outcome["quote_check"] = {
        "ran": True,
        "source": "PubMed get_article_metadata records",
        "num_quotes": len(citations),
        "num_quotes_verbatim_in_the_abstract": len(citations) - len(quotes_not_found),
        "quotes_not_found": quotes_not_found,
        "folding": "whitespace and case, beta for the Greek letter, a full stop for the middle dot",
    }
    return outcome


def claim_rows(claims: list, patterns: dict) -> list:
    rows = []
    for claim in claims:
        if claim["verdict"] == "drop":
            continue
        human_name = claim.get("symptom_the_human_data_support") or ""
        named = study_symptoms_named(human_name, patterns)
        rows.append(
            {
                "claimed_symptom": claim["symptom"],
                "direction": claim["direction"],
                "target": claim["target"],
                "verdict": claim["verdict"],
                "evidence_kind": claim["human_evidence"],
                "selective_human_drug": claim["selective_human_drug"],
                "study_symptoms_the_human_data_name": named,
                "human_name_matched": NEGATED_TAIL.split(human_name, maxsplit=1)[0].lower(),
                "lands_back_on_the_claimed_symptom": claim["symptom"] in named,
                "references": claim["references"],
                "note": claim["note"],
            }
        )
    return rows


def write_document(evidence: dict, rows: list, checks: dict, reference_records: dict, animal_readout_pattern: str,
                   species_ambiguous_pattern: str) -> None:
    method = evidence["method"]
    claims = evidence["claims"]
    verdict_counts = {}
    for claim in claims:
        verdict_counts[claim["verdict"].split(" to ")[0]] = verdict_counts.get(claim["verdict"].split(" to ")[0], 0) + 1
    kept = [row for row in rows if row["verdict"] == "keep as human"]
    relabelled = [row for row in rows if row["verdict"].startswith("keep but relabel")]
    rodent_only = [row for row in rows if row["verdict"] == "keep as rodent only"]
    lines = [
        "# Human evidence for the rodent-only target safety claims",
        "",
        "Generated by `experiments/map_rodent_readout_human_evidence.py`. Do not edit by hand.",
        "",
        "The Open Targets target safety curation supports 58 (symptom, direction, target) claims only through rodent",
        f"behavioural readouts, which `experiments/scope_open_targets_symptom_labels.py` matches with `{animal_readout_pattern}`",
        f"and `{species_ambiguous_pattern}`, the second being wording the rodent assay and the clinic share.",
        "Rodent locomotor activity is not a symptom a patient reports, so each claim needed human evidence or it leaves the",
        "label set. A sub-agent searched the human literature claim by claim; this document records what it was asked, what",
        "it returned, the checks run on the returned references and where each surviving claim lands among the study's",
        "symptoms.",
        "",
        "## What the search returned",
        "",
        "| Verdict | Claims |",
        "| --- | --- |",
    ]
    for verdict in ("keep as human", "keep but relabel", "keep as rodent only", "drop"):
        lines.append(f"| {verdict} | {verdict_counts.get(verdict, 0)} |")
    lines += [
        "",
        "`drop` means no human pharmacological evidence and either no human drug selective enough to test the claim or a",
        "selective human drug whose controlled data show the symptom does not occur. `keep as rodent only` means the symptom",
        "has no human evidence but the mechanism has indirect human support (inverse-direction pharmacology, an obligatory",
        "subunit or human genetics) worth keeping as a rodent-derived hypothesis rather than as a label.",
        "",
        "## Claims with human evidence for the claimed symptom",
        "",
        "| Symptom | Direction | Target | Evidence | Selective human drug | References |",
        "| --- | --- | --- | --- | --- | --- |",
    ]

    def reference_cell(row: dict) -> str:
        return ", ".join(sorted({f"[{reference['pmid']}](https://doi.org/{reference['doi']})" for reference in row["references"]}))

    for row in kept:
        lines.append(
            f"| {row['claimed_symptom']} | {row['direction']} | {row['target']} | {row['evidence_kind']} | "
            f"{'yes' if row['selective_human_drug'] else 'no'} | {reference_cell(row)} |"
        )
    lines += [
        "",
        "## Claims whose human evidence names a different symptom",
        "",
        "The study's own symptom patterns, read from `experiments/scope_open_targets_symptom_labels.py`, decide which",
        "symptom the human data name; the sub-agent's sentence is matched after its negated tail is cut off. A row whose",
        "`Study symptoms` column is empty names a symptom the study does not score.",
        "",
        "| Claimed symptom | Direction | Target | Human data name | Study symptoms | Evidence | References |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in relabelled:
        named = ", ".join(row["study_symptoms_the_human_data_name"]) or "none of the 23"
        relabel = row["verdict"].removeprefix("keep but relabel to ")
        lines.append(
            f"| {row['claimed_symptom']} | {row['direction']} | {row['target']} | {relabel} | {named} | "
            f"{row['evidence_kind']} | {reference_cell(row)} |"
        )
    lines += [
        "",
        "Rows marked `lands back on the claimed symptom` in the JSON are relabels whose new name the study's pattern for the",
        "claimed symptom already matches: "
        + (", ".join(f"{row['target']} {row['direction']} ({row['claimed_symptom']})" for row in relabelled if row["lands_back_on_the_claimed_symptom"]) or "none")
        + ".",
        "",
        "## Kept as a rodent-derived hypothesis, not as a label",
        "",
        "| Symptom | Direction | Target | Why |",
        "| --- | --- | --- | --- |",
    ]
    for row in rodent_only:
        lines.append(f"| {row['claimed_symptom']} | {row['direction']} | {row['target']} | {row['note']} |")
    lines += [
        "",
        "## What the two conflicting relabels became",
        "",
        "Both were settled by the user on 9 October 2026, and the patterns in",
        "`experiments/scope_open_targets_symptom_labels.py` now carry the decisions.",
        "",
        "- **Bradykinesia is not psychomotor slowing.** The pattern for `psychomotor_retardation` matched `\\bbradykine`",
        "  until then, so the DRD2-inhibition claim, whose human evidence names drug-induced parkinsonism, landed back on",
        "  the claimed symptom. The two differ in treatment, so they are kept apart (the user). The pattern now matches",
        "  only the two MedDRA terms `docs/symptom_crosswalk.csv` gives the symptom, `Psychomotor retardation` and",
        "  `Bradyphrenia`. The claim landed on no study symptom for a few hours, until the user asked for parkinsonism",
        "  to be scored in its own right because the side effect is clinically serious enough that failing to predict",
        "  it would be worse than scoring it imperfectly; the crosswalk's 24th symptom now takes it, so the relabel",
        "  lands on `parkinsonism` rather than back on the claimed symptom or on nothing.",
        "- **Stereotypy was a bug.** The same word marked a row as a rodent readout and matched `compulsive_behavior`, so a",
        "  curation row phrased \"punding / stereotyped repetitive behaviour\" was set aside as an animal readout before it",
        "  could become a label, while the symptom pattern was reading the same word as the symptom. Wording cannot decide",
        "  the species when the rodent assay and the clinic share it, so `stereotyp` is now its own flag",
        "  (`SPECIES_AMBIGUOUS_READOUT_PATTERN`) and such rows come to this reading instead of being decided by a regular",
        "  expression.",
        "",
        "## What the mapping leaves open",
        "",
    ]
    for row in relabelled:
        if not re.search(animal_readout_pattern, row["human_name_matched"]):
            continue
        lines.append(
            f"- The human name for {row['target']} {row['direction']}, \"{row['human_name_matched']}\", also matches the"
            f" rodent-readout filter `{animal_readout_pattern}`, so a curation row phrased that way is set aside as an"
            f" animal readout before it can become a label for {', '.join(row['study_symptoms_the_human_data_name'])}."
            " Reading it as a human symptom needs an exception to that filter, not only this evidence."
        )
    lines += [
        "- Every surviving claim is target-level: it says a drug acting at this target by this direction produces the",
        "  symptom in people, not that a particular drug in the study's evidence table does. They enter as target-level",
        "  rows or not at all.",
        "",
        "## Checks on the references",
        "",
        f"- {checks['num_citations']} citations over {checks['num_unique_articles']} articles for {checks['num_claims']} claims.",
        f"- Articles without a PubMed metadata record: {len(checks['articles_without_a_metadata_record'])}.",
        f"- Citations whose DOI disagrees with the PubMed record: {len(checks['citations_whose_doi_disagrees_with_the_record'])}.",
        f"- References missing a PMID, DOI or quote: {len(checks['references_missing_a_field'])}.",
    ]
    quote_check = checks.get("quote_check", {})
    if quote_check.get("ran"):
        lines.append(
            f"- Quotes verbatim in the abstract: {quote_check['num_quotes_verbatim_in_the_abstract']} of "
            f"{quote_check['num_quotes']}, folding {quote_check['folding']}."
        )
        for pmid, quote in quote_check.get("quotes_not_found", []):
            lines.append(f"  - not found: {pmid} — {quote}")
    else:
        lines.append("- Quotes were not checked against the abstracts in this run.")
    lines += [
        "",
        f"The reference records are in `{REFERENCE_FILE.relative_to(REPOSITORY_ROOT)}` (PMID, DOI, title, journal, year;",
        "no abstracts) and the check outcome in `" + str(CHECKS_FILE.relative_to(REPOSITORY_ROOT)) + "`. The articles are in",
        "`docs/references.bib`. All of them were retrieved from PubMed, and PubMed requires that it be cited and that the",
        "article DOIs be given, which the tables above do.",
        "",
        "## Sub-agent record",
        "",
        f"- Prompt, verbatim: `{PROMPT_FILE.relative_to(REPOSITORY_ROOT)}`.",
        "- Agent type: the general-purpose agent of the Claude Code Agent tool, run in the background from the session of",
        "  2026-10-08. The model identifier is in the session record rather than here, because this repository carries no",
        "  model identifiers.",
        f"- Search tools: {', '.join(method['search_tools_used'])}.",
        f"- Calls: {method['pubmed_searches_run']} searches, {method['pubmed_metadata_fetches']} metadata fetches.",
        f"- Fallbacks: {method['fallback_tools_needed']}.",
        f"- Evidence grading it applied: {method['evidence_grading']}.",
        f"- Returned record: `{EVIDENCE_FILE.relative_to(REPOSITORY_ROOT)}`, one object per claim with the verdict rule it",
        "  applied under `method.verdict_rule`.",
        "",
        "A sub-agent's report is model output, so the references were checked as above and the mapping onto symptoms is done",
        "here by the study's own patterns rather than taken from its sentences. What the checks cannot establish is whether",
        "a quote that is verbatim and correctly attributed supports the claim it is cited for; that reading is the",
        "sub-agent's and is open to review in the JSON, which carries its note for every claim.",
        "",
    ]
    DOCUMENT_FILE.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pubmed-metadata", type=Path, default=None, help="PubMed records with abstracts, to recheck the quotes")
    arguments = parser.parse_args()
    evidence = json.loads(EVIDENCE_FILE.read_text())
    reference_records = json.loads(REFERENCE_FILE.read_text())
    patterns, animal_readout_pattern, species_ambiguous_pattern = symptom_patterns()
    checks = check_references(evidence["claims"], reference_records, arguments.pubmed_metadata)
    CHECKS_FILE.write_text(json.dumps(checks, indent=1) + "\n")
    rows = claim_rows(evidence["claims"], patterns)
    write_document(evidence, rows, checks, reference_records, animal_readout_pattern, species_ambiguous_pattern)
    print(f"{len(rows)} claims kept of {len(evidence['claims'])}; wrote {DOCUMENT_FILE.relative_to(REPOSITORY_ROOT)}")
    for row in rows:
        print(f"  {row['verdict']:55s} {row['claimed_symptom']:25s} {row['direction']:10s} {row['target']:8s} -> {row['study_symptoms_the_human_data_name']}")


if __name__ == "__main__":
    main()
