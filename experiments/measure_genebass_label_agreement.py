"""How the monogenic gene labels look in Genebass's loss-of-function burden results (label_source_review.md section 33).

Genebass (Karczewski et al. 2022, doi:10.1016/j.xgen.2022.100168) tests rare-variant burden in 394,841 UK Biobank
exomes against 4,529 phenotypes. Its carriers are mostly heterozygous for one damaged copy, while the monogenic labels
come from patients with a disorder, often biallelic. This script asks how often a gene labelled "induces" for a symptom
construct shows a nominal burden association in the symptom direction on a UK Biobank item that measures that
construct, against the rate among all genes Genebass tested for the item; and the same split by Gene2Phenotype's
allelic requirement, since a gene whose disorders need one damaged copy should show an effect in heterozygous carriers
if the two sources measure the same thing.

It reads pages saved by experiments/read_public_page.py from the Genebass browser's own interface:

  https://main.genebass.org/api/phenotypes
  https://main.genebass.org/api/analysis/<analysis_id>/gene-manhattan?burdenSet=pLoF   (one per item below)

The gene-level table there gives BETA_Burden and Pvalue_Burden per gene and no carrier count or standard error, and the
analyses were run with SAIGE (the phenotype list names SAIGE_0.44.5), whose burden score need not count carriers one
by one, so a beta is not read here as a per-carrier effect; this script counts nominal associations and makes no claim
about effect sizes or measured absences. The item mapping is this script's reading of each item's wording, written for
this measurement and registered nowhere.

Usage:
  python experiments/measure_genebass_label_agreement.py --pages-directory data/raw/public_pages/genebass_2026-10-17 \\
      --evidence data/processed/evidence_full_v5/evidence_records.parquet --g2p data/raw/gene2phenotype/G2P_all_panels.csv
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

# Genebass phenocode -> (construct, the item's wording). Only items asked of every respondent of their questionnaire:
# 20534 ("Sleeping too much") and 20535 ("Waking too early") are asked about the worst episode of depression, so their
# denominators are people with such an episode, and they are left out (EXCLUDED_ITEMS).
ITEM_TO_CONSTRUCT = {
    "20510": ("depressed_mood", "Recent feelings of depression"),
    "20514": ("anhedonia", "Recent lack of interest or pleasure in doing things"),
    "insomnia_1200_binary": ("insomnia", "Insomnia (case_control)"),
    "20519": ("fatigue", "Recent feelings of tiredness or low energy"),
    "20505": ("irritability", "Recent easy annoyance or irritability"),
    "20502": ("irritability", "Ever had period extreme irritability"),
    "20506": ("anxiety", "Recent feelings or nervousness or anxiety"),
    "20509": ("anxiety", "Recent inability to stop or control worrying"),
    "20463": ("hallucinations", "Ever heard an un-real voice"),
    "20471": ("hallucinations", "Ever seen an un-real vision"),
    "20468": ("paranoia", "Ever believed in an un-real conspiracy against self"),
    "20474": ("delusion", "Ever believed in un-real communications or signs"),
    "20501": ("elevated_mood_or_mania", "Ever had period of mania / excitability"),
    "20016": ("cognitive_impairment", "Fluid intelligence score (a lower score is the symptom direction)"),
    "narcolepsy_1220_binary": ("somnolence_or_hypersomnia", "Daytime dozing, which Genebass names narcolepsy_1220_binary"),
}
EXCLUDED_ITEMS = {
    "20534": "asked about the worst episode of depression only",
    "20535": "asked about the worst episode of depression only",
    "20517": "insomnia and hypersomnia in one item",
    "20518": "psychomotor retardation and agitation in one item",
    "20511": "poor appetite and overeating in one item",
    "20508": "trouble concentrating, which has no construct of its own",
}
LOWER_IS_THE_SYMPTOM = frozenset({"20016"})
NOMINAL_THRESHOLD = 0.05


def read_saved_json(page_path: Path):
    """A page the reader saved, decompressed when the server sent it gzip-encoded."""
    raw = Path(page_path).read_bytes()
    return json.loads(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)


def saved_pages(pages_directory: Path) -> dict[str, Path]:
    """address -> saved file, from the reader's metadata files (<saved file>.json)."""
    pages = {}
    for metadata_path in sorted(Path(pages_directory).glob("*.json")):
        metadata = json.loads(metadata_path.read_text())
        if isinstance(metadata, dict) and metadata.get("read") and "address" in metadata:
            pages[metadata["address"]] = metadata_path.with_suffix("")
    return pages


def analysis_id_of_address(address: str) -> str | None:
    if "/api/analysis/" not in address or "/gene-manhattan" not in address:
        return None
    return address.split("/api/analysis/")[1].split("/gene-manhattan")[0]


def mapped_analyses(phenotypes: list[dict]) -> pd.DataFrame:
    """The both-sexes analyses of the mapped items, with their construct and case counts."""
    rows = []
    for phenotype in phenotypes:
        code = str(phenotype.get("phenocode"))
        if code in ITEM_TO_CONSTRUCT and phenotype.get("pheno_sex") == "both_sexes":
            construct, wording = ITEM_TO_CONSTRUCT[code]
            rows.append({"analysis_id": phenotype["analysis_id"], "phenocode": code, "construct": construct, "wording": wording,
                         "trait_type": phenotype.get("trait_type"), "cases": phenotype.get("n_cases"), "controls": phenotype.get("n_controls")})
    return pd.DataFrame(rows, columns=["analysis_id", "phenocode", "construct", "wording", "trait_type", "cases", "controls"])


def nominal_hits_in_symptom_direction(gene_results: pd.DataFrame, lower_is_the_symptom: bool,
                                      threshold: float = NOMINAL_THRESHOLD) -> dict[str, bool]:
    """gene -> whether its pLoF burden association is nominal (Pvalue_Burden < threshold) and in the symptom direction."""
    usable = gene_results.dropna(subset=["Pvalue_Burden", "BETA_Burden"])
    direction = -1.0 if lower_is_the_symptom else 1.0
    hits = (usable["BETA_Burden"] * direction > 0) & (usable["Pvalue_Burden"] < threshold)
    return dict(zip(usable["gene_symbol"], hits))


def requirement_of_gene(g2p: pd.DataFrame) -> dict[str, str]:
    """gene -> "monoallelic", "biallelic" or "both", from Gene2Phenotype's autosomal records, leaving out the records
    whose confidence is "disputed" or "refuted"; X-linked and mitochondrial records do not count."""
    usable = g2p[~g2p["confidence"].astype(str).str.lower().isin(["disputed", "refuted"])]
    kinds_of_gene: dict[str, set[str]] = {}
    for gene, requirement in zip(usable["gene symbol"], usable["allelic requirement"].astype(str)):
        if requirement == "monoallelic_autosomal":
            kinds_of_gene.setdefault(str(gene), set()).add("monoallelic")
        elif requirement == "biallelic_autosomal":
            kinds_of_gene.setdefault(str(gene), set()).add("biallelic")
    return {gene: (next(iter(kinds)) if len(kinds) == 1 else "both") for gene, kinds in kinds_of_gene.items()}


def labelled_genes(evidence: pd.DataFrame) -> dict[str, set[str]]:
    """construct -> the gene perturbations labelled "induces" for it."""
    genes = evidence[(evidence["perturbation_type"] == "gene") & (evidence["relation"] == "induces")]
    return {construct: set(group["perturbation_id"].astype(str)) for construct, group in genes.groupby("symptom")}


def agreement(analyses: pd.DataFrame, gene_results_of_analysis: dict[str, pd.DataFrame], labels: dict[str, set[str]],
              group_of_gene) -> pd.DataFrame:
    """One row per (item, group of labelled genes): how many were tested, how many are nominal hits in the symptom
    direction, the rate among every gene tested for the item, and a one-sided binomial p of the excess."""
    rows = []
    for analysis in analyses.itertuples():
        if analysis.analysis_id not in gene_results_of_analysis:
            continue
        hits = nominal_hits_in_symptom_direction(gene_results_of_analysis[analysis.analysis_id], analysis.phenocode in LOWER_IS_THE_SYMPTOM)
        background = sum(hits.values()) / len(hits) if hits else float("nan")
        genes_by_group: dict[str, list[str]] = {}
        for gene in sorted(labels.get(analysis.construct, set())):
            genes_by_group.setdefault(group_of_gene(gene), []).append(gene)
        for group, genes in sorted(genes_by_group.items()):
            tested = [gene for gene in genes if gene in hits]
            found = [gene for gene in tested if hits[gene]]
            p_value = stats.binomtest(len(found), len(tested), background, alternative="greater").pvalue if tested else float("nan")
            rows.append({"analysis_id": analysis.analysis_id, "construct": analysis.construct, "cases": analysis.cases, "group": group,
                         "labelled": len(genes), "tested": len(tested), "hits": len(found), "background_rate": background,
                         "binomial_p": p_value, "hit_genes": ";".join(found)})
    return pd.DataFrame(rows, columns=["analysis_id", "construct", "cases", "group", "labelled", "tested", "hits", "background_rate", "binomial_p", "hit_genes"])


def pooled(rows: pd.DataFrame) -> pd.DataFrame:
    """Per group, over every item: gene-item pairs tested, hits, the expected hits at each item's background rate, and a
    one-sided binomial p against the mean background (pairs are not independent: a gene enters once per item of its
    construct, and the items of one construct correlate, so the p is descriptive)."""
    result = []
    for group, part in rows.groupby("group"):
        tested, hits = int(part["tested"].sum()), int(part["hits"].sum())
        expected = float((part["tested"] * part["background_rate"]).sum())
        mean_background = expected / tested if tested else float("nan")
        p_value = stats.binomtest(hits, tested, mean_background, alternative="greater").pvalue if tested else float("nan")
        result.append({"group": group, "pairs_tested": tested, "hits": hits, "expected": round(expected, 1),
                       "ratio": round(hits / expected, 2) if expected else None, "binomial_p": p_value})
    return pd.DataFrame(result)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pages-directory", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--g2p", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None, help="a CSV of the per-item rows")
    arguments = parser.parse_args()
    pages = saved_pages(arguments.pages_directory)
    phenotype_page = pages.get("https://main.genebass.org/api/phenotypes")
    if phenotype_page is None:
        print("the phenotype list (https://main.genebass.org/api/phenotypes) is not among the saved pages", file=sys.stderr)
        return 1
    analyses = mapped_analyses(read_saved_json(phenotype_page))
    gene_results = {analysis_id: pd.DataFrame(read_saved_json(path)) for address, path in pages.items()
                    if (analysis_id := analysis_id_of_address(address)) is not None and "burdenSet=pLoF" in address}
    missing = sorted(set(analyses["analysis_id"]) - set(gene_results))
    if missing:
        print(f"no saved pLoF results for {len(missing)} mapped analyses: {', '.join(missing)}", file=sys.stderr)
    evidence = pd.read_parquet(arguments.evidence)
    labels = labelled_genes(evidence)
    allele_class = evidence[evidence["perturbation_type"] == "gene"].groupby("perturbation_id")["allele_class"].agg(lambda values: ";".join(sorted(set(map(str, values)))))
    by_allele_class = agreement(analyses, gene_results, labels,
                                lambda gene: "loss_of_function label" if "loss_of_function" in allele_class.get(gene, "") else "other label")
    requirement = requirement_of_gene(pd.read_csv(arguments.g2p, low_memory=False))
    by_requirement = agreement(analyses, gene_results, labels, lambda gene: requirement.get(gene, "not in G2P"))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 50)
    print(by_allele_class.drop(columns=["hit_genes"]).to_string(index=False))
    print("\npooled by the label's allele class\n" + pooled(by_allele_class).to_string(index=False))
    print("\npooled by Gene2Phenotype's allelic requirement\n" + pooled(by_requirement).to_string(index=False))
    if arguments.output is not None:
        pd.concat([by_allele_class.assign(split="allele_class"), by_requirement.assign(split="allelic_requirement")]).to_csv(arguments.output, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
