# Evidence reports: one row per report, weighted aggregation and a report-level reliability fit

Implementation specification for the evidence layer (design sections 4.2, 4.3, 10, 11 open question 12 and 12). Status: specification, not yet implemented. Every function below names the file it lives in and its signature; behaviour that today's code already has is marked "as today" so the default output stays byte-comparable with the running sweeps.

## 0. What the raw files hold and what today's code misses

Measured on the pinned inputs (HPO release 2026-09-01, SIDER 4.1) and on the monogenic slice in data/processed/evidence (451 genes, 861 pairs):

- genes_to_phenotype.txt has 333,983 rows and every one of them joins to at least one phenotype.hpoa row on (disease_id, hpo_id). The join is one-to-many: 447 (disease, term) keys have two or more phenotype.hpoa rows (357 with two, 65 with three, up to 13), each with its own reference, evidence code, frequency and biocuration date. genes_to_phenotype.txt collapses them to one row and keeps one of the frequencies (OMIM:612567 HP:0025084: phenotype.hpoa rows 2/2 and 1/1, genes_to_phenotype 1/1). On the slice three genes_to_phenotype rows sit on multi-row keys (APOE OMIM:606889 twice, DNMT1 OMIM:614116). 95 (disease, term, gene) triples appear twice in genes_to_phenotype.txt.
- Neither file carries a single Excluded frequency qualifier (HP:0040285) in this release. The negation that does exist is the NOT qualifier of phenotype.hpoa (733 rows, 479 of them on diseases with a gene), and genes_to_phenotype.txt copies those annotations as ordinary rows with frequency "-" and no marker: 826 genes_to_phenotype rows are NOT-qualified in phenotype.hpoa. On the slice three of them map to a target symptom and all three are positive rows in evidence_records.parquet today: CPT1C ORPHA:444099 NOT Dementia (the only annotation behind the pair, so a pair that asserts absence is a grade A positive with weight 1.0), PDE10A ORPHA:494541 NOT Dementia against a positive OMIM:616922 Mental deterioration row, and GLRX5 ORPHA:401866 NOT Cognitive impairment against a positive Short attention span row of the same disease entry.
- The slice's 1,283 genes_to_phenotype rows expand to 1,288 phenotype.hpoa rows: evidence codes TAS 884, PCS 278, IEA 126; onset on 7 rows; sex on none; references ORPHA 809, PMID 291, OMIM 203, URL 2.
- SIDER meddra_freq.tsv.gz has one row per (label, frequency statement): 59,333 (drug, preferred term) pairs carry a treatment-arm frequency, 11,586 of them also a placebo-arm row (placebo_flag "placebo"), and about 27,000 pairs have two or more frequency rows (both arms counted). frequency_text is a percentage ("21%"), a range ("1-10%") or a word with fixed bounds (postmarketing 0 to 0.001, rare 0 to 0.001, very rare 0 to 0.0001, infrequent and uncommon 0.001 to 0.01, common and frequent 0.01 to 1, very common 0.1 to 1). SIDER has no patient denominators. meddra_all_se.tsv.gz has 163,206 preferred-term rows for 145,321 distinct (drug, preferred term) pairs; the duplicates come from several lower-level terms mapping to one preferred term. meddra_all_indications.tsv.gz has 5,110 NLP_indication preferred-term rows.
- experiment_data.py sets outcome 1 for every row of evidence_records.parquet whatever its weight, and run_main_model.py trains positives at max(weight, 1e-3). A weight 0 row is therefore a positive label in evaluation and a near-zero-weight positive in training. That is intended for the version 0.3 ablation table (data/processed/evidence_two_entries: 408 grade C rows with weight 0 that stay in the test sets) and must not change; it is wrong for a pair whose evidence says the symptom is absent.

The design's learned weighting (section 4.3) is tested in evidence_reliability_model.py but reads a {0, 1} report per (item, source). The report table below is what it needs, and the weighted generalization in section 3 is what makes several reports per cell, with different weights and signs, enter the likelihood.

## 1. The report table: data/processed/<evidence-dir>/evidence_reports.parquet

One row per report. A report is one phenotype.hpoa row behind one (gene, target symptom) or one SIDER frequency statement, label mention or indication behind one (drug, target symptom). Column names are final; the set is fixed.

| column | type | HPO-OMIM and HPO-Orphanet reports | SIDER-label and SIDER-indication reports |
|---|---|---|---|
| report_id | str | `"{source}|{source_record_id}|{source_term_id}|{symptom}|{ordinal}"`, ordinal = position of the phenotype.hpoa row within its (disease, term) group in file order, 0 for an unjoined report (section 1.2) | `"{source}|{stitch_flat_id}|{meddra_cui}|{symptom}|{ordinal}"`, ordinal = position of the meddra_freq row within its (drug, preferred term) treatment-arm group in file order, 0 for a label_text or indication report |
| perturbation_id | str | gene symbol (as today) | STITCH flat id (as today) |
| perturbation_type | str | "gene" | "drug" |
| perturbation_label | str | gene symbol | drug name from drug_names.tsv |
| symptom | str | target symptom whose expansion (roots minus excluded subtrees, as today) contains hpo_id; a term in two expansions yields two reports | target symptom of the preferred term from the crosswalk |
| relation | str | "induces" | "induces" for meddra_all_se and meddra_freq, "relieves" for meddra_all_indications |
| evidence_class | str | "monogenic" | "pharmacological" |
| source | str | "HPO-OMIM" when disease_id starts with "OMIM:", "HPO-Orphanet" for "ORPHA:" | "SIDER-label" or "SIDER-indication" |
| report_value | int | 0 when the phenotype.hpoa row's qualifier is NOT, its frequency is HP:0040285 or its frequency is a patient fraction with numerator 0 (0/35: an observed absence; decision 14), else 1 | 1 (SIDER records no absences) |
| source_record_id | str | disease_id | STITCH flat id |
| source_record_label | str | disease_name from the phenotype.hpoa row (genes_to_phenotype has none) | drug name |
| source_term_id | str | hpo_id | UMLS_meddra_cui of the preferred term (column 5 of meddra_all_se, 9 of meddra_freq, 6 of meddra_all_indications) |
| source_term_label | str | hpo_name from genes_to_phenotype | side_effect_name or meddra_concept_name (the preferred term; licensed, so the table stays gitignored and no name is committed) |
| evidence_code | str | phenotype.hpoa evidence column: "PCS", "TAS", "IEA"; any other value is kept verbatim and treated as "other" by the rubric; "unjoined" when no phenotype.hpoa row exists | "placebo_controlled_frequency" for a meddra_freq treatment row whose (drug, preferred term) also has a placebo row; "frequency" for a treatment row without one; "label_text" for a meddra_all_se pair with no treatment-arm frequency row and for every indication |
| model_description | str | `"human loss-of-function; {disease_id} {disease_name}"` | `"human; drug label; ChEMBL targets {target_chembl_id}:{action_type};..."` sorted by target id, "relieves" reports prefixed `"human; drug indication; ..."` |
| frequency | float or null | parse_frequency_qualifier of the phenotype.hpoa frequency column (midpoint of a qualifier, fraction or percentage; null for empty); 0.0 for Excluded | (frequency_lower + frequency_upper) / 2 of the treatment row; null for label_text |
| frequency_denominator | int or null | parse_frequency_denominator of the same column (5 for "3/5") | null |
| placebo_flag | bool | false | true when evidence_code is placebo_controlled_frequency |
| onset | str or null | onset column (an HP term id) or null | null |
| sex | str or null | sex column ("MALE", "FEMALE") or null | null |
| references | str | the reference column with its ";" separators and prefixes kept (PMID:, OMIM:, ORPHA:, ISBN-13:, http) | "" |
| pubmed_reference_count | int | number of "PMID:" entries in references | 0 |
| evidence_date | str or null | ISO date: for an OMIM row the earliest of its biocuration dates and the publication dates of its PMID references found in docs/hpo_reference_publication_dates.json (the per-row form of today's load_hpo_annotation_dates); null for Orphanet rows (their biocuration date is the release import date) and unjoined reports | null (SIDER 4.1 labels are a 2015 slice, not dated per statement) |
| rubric_log_sample_size | float | log1p(frequency_denominator), 0.0 when null | 0.0 |
| rubric_frequency_known | float | 1.0 when frequency is not null | 1.0 for frequency reports, 0.0 for label_text |
| rubric_evidence_code_pcs, rubric_evidence_code_tas, rubric_evidence_code_iea | float | one-hot of evidence_code; all 0.0 for other and unjoined | 0.0 |
| rubric_placebo_controlled | float | 0.0 | 1.0 when placebo_flag |
| rubric_curated_synopsis | float | 1.0 when references contains the disease entry itself (OMIM:{id} or ORPHA:{id}), meaning the annotation was transcribed from the curated synopsis rather than from a named study; 0.0 otherwise | 0.0 |
| limitations | str | section 1.3 | section 1.3 |
| perturbation_nodes | str | JSON `[[node_id, -1.0, 1.0]]` as today | JSON list of [node_id, sign, 1/len(mapped)] per target gene as today |

Rows whose perturbation has no graph node go to unmapped_records.parquet as today and produce no report.

### 1.1 Monogenic derivation

Join: for every genes_to_phenotype row with gene_symbol in the graph and hpo_id in a symptom expansion, after dropping duplicate (gene_symbol, disease_id, hpo_id) triples, look up phenotype.hpoa rows by (disease_id, hpo_id) and emit one report per phenotype.hpoa row per target symptom. The frequency, qualifier, evidence code, references, onset, sex, disease name and date come from the phenotype.hpoa row, not from genes_to_phenotype; the genes_to_phenotype frequency is one of the phenotype.hpoa frequencies (where a single row exists they agree on every one of 333,983 rows, with "-" standing for empty). DECIPHER prefixes do not occur in genes_to_phenotype.txt; a disease_id with any other prefix is dropped and counted in the summary under reports_dropped_unknown_provenance.

Descendant expansion: a disease annotated to Dementia, Memory impairment and Mental deterioration yields three reports for (gene, cognitive_impairment), as it yields three rows today. The report table keeps them separate and the aggregation (section 2) counts them, which is why report_count is a count of annotation rows and not of independent studies; the limitations text says so when the pair rests on one disease entry.

### 1.2 Unjoined rows

When phenotype.hpoa is not supplied (the loader is callable without it, as the existing tests do) or a (disease_id, hpo_id) key has no row, the genes_to_phenotype row itself becomes one report with evidence_code "unjoined", report_value 0 if its frequency is HP:0040285 else 1, frequency and frequency_denominator from its own frequency column, source_record_label equal to disease_id, empty references, null date, rubric one-hots 0.0 and the "other" rubric factor. The summary counts them (reports_without_hpoa_row); on the pinned release the count is 0.

### 1.3 Limitations text

`limitations_text(report) -> str` in evidence_reports.py joins, with "; ", the clauses that apply: "no frequency reported" (frequency null), "frequency from a curated qualifier, no patient count" (frequency known, denominator null, HPO), "reported in 0 of {denominator} patients" (frequency 0 with a known denominator, the absence claim of a 0/N row), "n = {denominator} patients" (denominator known, below 20), "inferred from electronic annotation (IEA)", "traceable author statement, no primary study cited (TAS)", "no PubMed reference", "annotation imported from the Orphanet clinical-sign table, undated", "no phenotype.hpoa row joined", "label text only, no frequency", "frequency not placebo-controlled", "SIDER reports no patient counts", "label frequencies are not comparable across labels". The mechanism cards (mechanism_cards.py) print the pair's limitations when they list supporting evidence; the column is for display and never enters a weight.

### 1.4 Pharmacological derivation

Inclusion filters as today (ChEMBL dominant target within --max-drug-targets, ATC group N). For each included DrugLabelEvent with relation induces: if meddra_freq has treatment-arm rows (placebo_flag empty, concept_type PT, matching STITCH flat id and preferred term) one report per such row, else one label_text report. For each included indication event (NLP_indication only, as today) one label_text report with source SIDER-indication and relation relieves. The event's label_frequency (mean of the treatment midpoints, placebo rows excluded) is unchanged, so the grade weighting is unchanged.

## 2. Aggregation to evidence_records.parquet

One row per (perturbation_id, symptom, relation) with at least one explicit report. Items without any report are absent, as today: unobserved pairs are unlabelled, not negative.

Existing columns keep name and meaning:

- perturbation_id, perturbation_type, perturbation_label, group_id, disease_cluster_id, symptom, relation, evidence_class, in_metabolic_layer, perturbation_nodes: as today. disease_cluster_ids keeps receiving disease_identifiers from all reports of the pair, positive and negative, so clusters do not move.
- grade: assign_evidence_grade on an EvidenceRecord built from the pair's positive reports (omim_entry_count and orpha_entry_count over positive reports' source_record_id, has_omim_clinical_synopsis = any positive OMIM report, independent_case_series_count = distinct disease entries among positive reports, the version 0.3 proxy; cns_penetrant and has_dominant_target true for included drugs as today). A pair with no positive report gets grade "C" for monogenic and "E" for pharmacological: the grades whose weight is 0 or soft prior, which is what "no positive evidence" means in the grade table.
- weight: section 2.1.
- grade weighting rows (decision 14): the grade, the weight and the descriptive columns below that say "positive reports" are computed over `grade_weighting_rows(pair_reports)`, the positive reports plus the 0/N patient-fraction rows, because the version 0.4 loader read a 0/N row as a frequency-0 presence row and the default weights must reproduce its weights; positive_report_count counts report_value 1 only, so a pair whose rows are all 0/N keeps its grade and floor weight but is not a positive label.
- label_frequency: monogenic, the largest frequency over positive reports (as today: the design's "largest reported frequency"; identical on the slice, checked on the three multi-row keys where the hpoa maximum and the genes_to_phenotype value clip to the same weight); pharmacological, the mean of the frequency reports' midpoints (as today's load_sider_frequencies; a maximum here would change the running weights). Null when no positive report has a frequency.
- source: monogenic `"HPO genes_to_phenotype; diseases=" + ";".join(sorted(all disease ids))`; pharmacological `"SIDER 4.1; targets ..."`, as today.
- omim_entry_count, orpha_entry_count, disease_identifiers: over all explicit reports (today they count NOT-qualified rows too, because today cannot see them).
- annotation_row_count: distinct (source_record_id, source_term_id) among the pair's reports, which equals today's genes_to_phenotype row count except where that file duplicates a triple.
- annotation_patient_count: sum of frequency_denominator over positive reports (today sums over genes_to_phenotype rows, so APOE and DNMT1 change from 4 and null to 16 and 18; descriptive only).
- distinct_reference_count, distinct_pubmed_reference_count: distinct reference strings and their PMID subset over positive reports (today's loaders skip NOT rows, so the meaning is unchanged).
- evidence_date: minimum evidence_date over positive reports (as today; NOT rows were already skipped).

New columns:

- report_count, positive_report_count, negative_report_count: explicit reports of the pair, report_value 1 and 0.
- reliability_posterior: section 3, filled for every pair whatever the weighting.
- weighting: "grade" or "reliability", the value of --weighting, the same on every row.

### 2.1 Weight

Under weighting "grade": `loss_weight_for_record(record, grade_a_policy=...)` on the EvidenceRecord of the positive reports, exactly as today (grade base times clip(frequency × scale, floor, 1), unknown frequency keeps the base), and 0.0 when positive_report_count is 0. Regression requirement: on the 2026-09-01 slice the weights equal today's for 860 of 861 pairs; the one change is CPT1C cognitive_impairment, 1.0 to 0.0, with positive_report_count 0. The pair stays in the table so the absence claim is visible.

Under weighting "reliability": weight = reliability_global_scale × reliability_posterior when positive_report_count > 0, else 0.0 (`reliability_weights`, over `observation_weights_from_posterior`); refused with ValueError when the fit lists a degenerate source (section 3, decision 16).

### 2.2 experiment_data.py

Today load_experiment_data makes every row a positive. Smallest change, in `load_experiment_data` right after the relation filter:

```python
if "positive_report_count" in evidence.columns:  # tables written before reports existed have no such column and every row is a positive claim
    evidence = evidence[evidence.positive_report_count > 0]
```

Filtering on weight > 0 is not acceptable: it would drop the 408 grade C rows of data/processed/evidence_two_entries from the test sets, which the version 0.3 ablation keeps as positives by design. The frequency and date maxima then run over positive rows only, as today.

## 3. Report-level reliability fit

File: mechanistic_pathway_learning/evidence/evidence_reliability_model.py.

Items are the aggregated pairs (perturbation_id, symptom, relation). Sources are the four of section 1 in a fixed order; the fit only uses sources with at least one report in the table. coverage[i, s] = 1 when source s has any report for item i's perturbation_id and relation (for any symptom). Per (item, source) the model sees positive_weights[i, s] = the sum over the source's distinct records (source_record_id: a disease entry, a drug) of the largest rubric weight among that record's explicit reports with report_value 1 on the item, negative_weights[i, s] = the same over report_value 0, plus implicit_negative_weight (default 1.0) when coverage[i, s] = 1 and the source has no explicit report for the item. One record contributes at most one report's weight per side: a synopsis annotated to three descendant terms of one symptom, or several labels of one drug mentioning one term, are not independent observations (section 1.1; decision 15). Likelihood of cell (i, s): sensitivity^pos (1 − sensitivity)^neg under z = 1 and (1 − specificity)^pos specificity^neg under z = 0, so today's model is the special case pos = reports × coverage, neg = (1 − reports) × coverage.

Three guards on the fit (decision 16): prevalence is fitted per evidence class (item_groups), because no source covers both gene and drug items and a pooled prevalence would hand the OMIM-Orphanet agreement to every drug pair; a source with no explicit negative report keeps its specificity at the prior mean a / (a + b) = 0.8, because silence is the only negative it has and the specificity M-step then measures the volume of its positive reports, not a specificity; and a source whose estimated sensitivity + specificity ends at or below 1 (the label-switched optimum, where a positive report has likelihood ratio below 1 and silence above) has its specificity held at the prior mean too and EM runs again, until no estimated source is inverted. A source still at or below 1 after the holds is a degenerate source: the fit records it and --weighting reliability is refused on such a fit.

### 3.1 Rubric weight

File: mechanistic_pathway_learning/evidence/evidence_reports.py. A documented fixed function of the rubric columns, in (0, 1]:

```
rubric_weight = clip(evidence_code_factor × sample_size_factor, minimum_rubric_weight, 1.0)
```

- evidence_code_factor: PCS 1.0, TAS 0.8, IEA 0.6, other (including unjoined) 0.7; SIDER placebo_controlled_frequency 1.0, frequency 0.9, label_text 0.8.
- sample_size_factor: 1.0 when frequency_denominator is null (a curated qualifier or a label frequency is not evidence of a small sample); when known, 0.5 + 0.5 × min(1, log1p(n) / log1p(sample_size_reference)) with sample_size_reference 20, so n = 1 gives 0.61, n = 5 gives 0.79 and n ≥ 20 gives 1.0.
- minimum_rubric_weight 0.1; implicit_negative_weight 1.0 (section 3, overridable here because HPO silence inside a synopsis is weaker than a NOT annotation; see 3.4).

```python
@dataclass(frozen=True)
class RubricWeightDefaults:
    evidence_code_factors: dict[str, float]  # keys: PCS, TAS, IEA, other, placebo_controlled_frequency, frequency, label_text
    sample_size_reference: int = 20
    sample_size_unknown_factor: float = 1.0
    minimum_rubric_weight: float = 0.1
    implicit_negative_weight: float = 1.0

def load_rubric_weight_defaults(overrides_path: Path | None) -> RubricWeightDefaults
    # JSON with any subset of the fields above; unknown keys raise ValueError
def rubric_weight(report: EvidenceReport, defaults: RubricWeightDefaults) -> float
def rubric_weights_for_table(reports: pd.DataFrame, defaults: RubricWeightDefaults) -> np.ndarray
```

The weights are written back to the report table as column rubric_weight and the defaults used are recorded in the summary JSON (reliability.rubric_weight_defaults). Later rubric columns from the document appraisal (llm_evidence_appraisal.EvidenceRubric.feature_vector) can be appended with the same rubric_ prefix; the fixed function is the stand-in for the Raykar feature-dependent form until a third source makes that form estimable.

### 3.2 Functions

```python
@dataclass
class ReportReliabilityFit:
    source_names: list[str]
    prevalence: float
    sensitivity: np.ndarray            # [num_sources]
    specificity: np.ndarray            # [num_sources]
    item_keys: list[tuple[str, str, str]]  # (perturbation_id, symptom, relation) in table order
    posterior: np.ndarray              # [num_items] P(z_i = 1 | reports)
    num_iterations: int
    converged: bool
    weakly_identified: bool            # True when fewer than three sources cover any item (max over items of the number of covering sources < 3)
    log_likelihood_trace: list[float]
    per_source_summary: dict[str, dict] # section 3.3
    prevalence_by_group: dict[str, float]      # one prevalence per item group (evidence class); prevalence is the pooled value
    specificity_held_at_prior: np.ndarray      # [num_sources]
    specificity_hold_reasons: dict[str, str]   # held source -> "no_explicit_negatives" or "degenerate_when_estimated"
    degenerate_sources: list[str]              # sensitivity + specificity <= 1 after the holds

def report_count_matrices(
    reports: pd.DataFrame,             # evidence_reports with rubric_weight filled
    item_keys: list[tuple[str, str, str]],
    source_names: list[str],
    implicit_negative_weight: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]
    """positive_weights, negative_weights, coverage, each [num_items, num_sources]; silent covered cells get implicit_negative_weight in negative_weights."""

def weighted_posterior_truth_probability(positive_weights, negative_weights, prevalence, sensitivity, specificity) -> np.ndarray
def weighted_marginal_log_likelihood(positive_weights, negative_weights, prevalence, sensitivity, specificity) -> float

def fit_weighted_dawid_skene(
    positive_weights: np.ndarray,
    negative_weights: np.ndarray,
    source_names: list[str],
    item_keys: list[tuple[str, str, str]] | None = None,
    num_iterations: int = 100,
    sensitivity_prior_counts: tuple[float, float] = (8.0, 2.0),
    specificity_prior_counts: tuple[float, float] = (8.0, 2.0),
    prevalence_prior_counts: tuple[float, float] = (1.0, 1.0),
    convergence_tolerance: float = 1e-6,
    item_groups: list[str] | None = None,          # one label per item; None pools
    estimate_specificity: np.ndarray | None = None, # one boolean per source; False holds the prior mean
) -> ReportReliabilityFit
```

EM as today, with the per-group prevalence in the E-step and the specificity of a held source replaced by the prior mean after every M-step: majority-vote initialization posterior_0[i] = Σ_s pos[i, s] / Σ_s (pos[i, s] + neg[i, s]) clipped to [1e-6, 1 − 1e-6]; MAP M-step with weighted counts, sensitivity_s = (Σ_i post_i pos[i, s] + a − 1) / (Σ_i post_i (pos + neg)[i, s] + a + b − 2) and the mirror for specificity with 1 − post_i and neg first, prevalence = (Σ post + a) / (N + a + b); E-step weighted_posterior_truth_probability; convergence when the maximum absolute posterior change is below the tolerance; the marginal log likelihood appended per iteration. `fit_dawid_skene` keeps its signature and return type and becomes a wrapper that builds pos = reports × coverage and neg = (1 − reports) × coverage, so tests/test_evidence_reliability_and_appraisal.py passes unchanged.

Assembler side (assemble_evidence_table.py):

```python
def fit_report_reliability(reports: pd.DataFrame, observations: pd.DataFrame, defaults: RubricWeightDefaults) -> ReportReliabilityFit
    # item_keys from observations in row order; source_names = the four sources that have reports;
    # item_groups = evidence_class per observation; estimate_specificity = explicit negative weight > 0 per source
def reliability_weights(fit: ReportReliabilityFit, observations: pd.DataFrame, global_scale: float) -> np.ndarray
    # the --weighting reliability weights; ValueError when fit.degenerate_sources is not empty
```

### 3.3 Per-source summary in evidence_summary.json

Under key "reliability": prevalence, prevalence_by_evidence_class, num_iterations, converged, weakly_identified, degenerate_sources, rubric_weight_defaults, posterior_quantiles, and per source: items_covered, items_with_explicit_reports, explicit_positive_weight, explicit_negative_weight, implicit_negative_cells, sensitivity, specificity, specificity_held_at_prior (null, "no_explicit_negatives" or "degenerate_when_estimated"), degenerate. Under key "reports": rows, by_source_and_value, by_evidence_code, reports_without_hpoa_row, reports_dropped_unknown_provenance, pairs_all_negative, pairs_conflicting. posterior_quantiles are given for four pair categories: covered by two sources and positive in both, positive in one and silent in the other, conflicting (explicit positive and explicit negative), all negative.

### 3.4 What the posterior means on the current sources

A prototype of exactly this fit on the slice (two sources, 861 items, 336 covered by both, 1,285 positive and 3 negative reports, TAS 0.8, PCS 1.0, IEA 0.6) converged in 176 iterations to prevalence 0.655, OMIM sensitivity 0.989 and specificity 0.821, Orphanet sensitivity 0.865 and specificity 0.03.

- One positive HPO-Orphanet report, HPO-OMIM silent but covering (the gene has an OMIM annotation to some other target symptom). The cell log odds are w × log(sens_O / (1 − spec_O)) for the Orphanet report and log((1 − sens_M) / spec_M) for the implicit OMIM negative, so OMIM silence counts as one full negative report from the better sensor and pulls the posterior below the prevalence: the 206 such pairs averaged 0.28 (range 0.01 to 0.99). When the gene has no OMIM annotation to any target symptom the OMIM column is uncovered, contributes nothing and the posterior sits above the prevalence (523 single-source pairs averaged 0.74). The posterior therefore encodes whether OMIM happened to curate the gene, not whether the link is real.
- Conflicting reports (PDE10A: OMIM PCS positive, Orphanet TAS NOT): both cells enter with their rubric weights. The posterior lies strictly between the all-positive and the all-negative value for the same coverage only when every covering source has sensitivity + specificity > 1, since the ordering needs log(sens / (1 − spec)) > log((1 − sens) / spec) for each cell; positive rubric weights and clipping alone do not give it. In the prototype that condition held (measured 0.97, and 0.85 for GLRX5 where the positive and the NOT row come from one Orphanet entry); with OMIM sensitivity 0.681 and specificity 0.034, which the free fit on the full graph produced before decision 16, a PDE10A-type cell gave all-positive 0.979, conflicting 0.053 and all-negative 0.429, the conflicting pair below the all-negative one. The implemented fit holds the specificity of such a source at the prior mean and lists any source still at or below 1 as degenerate (section 3), so the ordering holds on every fit that --weighting reliability accepts; test_inverted_source_is_held_at_prior_then_flagged_and_reliability_weighting_refused covers the inverted case.
- All-negative (CPT1C): posterior 0.86 in the prototype, because the free fit drove Orphanet specificity to 0.03, which makes an Orphanet NOT report uninformative. This is the identifiability failure the design predicts for two sources, and it is why the posterior is written as a column and fitted on every run but is not the default weight (section 6, open question 12). With the holds of decision 16 the fit no longer inverts a source but lands on the other weakly identified optimum (decision 18). The weakly_identified flag is True for every table until a third source exists: HPO sources cover genes only, SIDER-label covers drugs with relation induces only and SIDER-indication drugs with relation relieves only, so no item is ever covered by three sources and drug items are covered by exactly one.

## 4. Command line and assembler

assemble_evidence_table.py gains:

- `--weighting {grade,reliability}` default grade
- `--report-rubric-weights <json path>` optional, read by load_rubric_weight_defaults
- `--reliability-global-scale <float>` default 1.0

It always writes evidence_reports.parquet (atomic replace, next to evidence_records.parquet and unmapped_records.parquet), always fits the reliability model and fills reliability_posterior and the three report counts, whatever the weighting, and records weighting, reliability_global_scale and the rubric defaults in evidence_summary.json. configs/evidence_assembly.yaml is not read by the assembler; its weighting block gets a comment saying the command line is authoritative and that mode learned_reliability is the --weighting reliability arm.

Functions to add or change:

- evidence_reports.py (new): `@dataclass EvidenceReport` with one field per column of section 1 (rubric columns as floats, perturbation_nodes as a JSON string), `REPORT_COLUMNS: tuple[str, ...]`, `reports_to_dataframe(reports: list[EvidenceReport]) -> pd.DataFrame`, `limitations_text`, `RubricWeightDefaults`, `load_rubric_weight_defaults`, `rubric_weight`, `rubric_weights_for_table`.
- load_monogenic_phenotype_annotations.py: `@dataclass HpoaAnnotationRow(disease_id, disease_name, qualifier, hpo_id, references: list[str], evidence, onset, frequency, sex, modifier, aspect, biocuration_dates: list[date], ordinal: int)`; `load_hpo_annotation_rows(phenotype_hpoa_path: Path) -> dict[tuple[str, str], list[HpoaAnnotationRow]]`; `hpoa_row_availability_date(row: HpoaAnnotationRow, dated_provenance_prefixes=("OMIM:",), publication_dates_by_pmid=None) -> date | None`; `monogenic_evidence_reports(annotation_rows, target_symptom_to_hpo_ids, parents_by_term, genes_in_graph, excluded_hpo_ids_by_symptom=None, hpoa_rows_by_key=None, publication_dates_by_pmid=None) -> list[EvidenceReport]` (perturbation_nodes left empty, filled by the assembler); `load_hpo_annotation_dates` and `load_hpo_annotation_references` keep their signatures and are reimplemented over load_hpo_annotation_rows and hpoa_row_availability_date. `monogenic_evidence_records` stays for experiments/run_phase1_counts.py and the existing tests.
- load_drug_label_events.py: `@dataclass SiderFrequencyRow(stitch_flat_id, label_cui, placebo_flag, frequency_text, frequency_lower, frequency_upper, meddra_cui, preferred_term, ordinal)`; `load_sider_frequency_rows(sider_directory) -> dict[tuple[str, str], list[SiderFrequencyRow]]` (PT rows, both arms); `load_sider_frequencies` reimplemented over it with the same output; DrugLabelEvent gains `meddra_cui: str | None` and `frequency_rows: list[SiderFrequencyRow]`; `drug_label_reports(event: DrugLabelEvent, model_description: str, perturbation_nodes_json: str) -> list[EvidenceReport]`.
- assemble_evidence_table.py: `@dataclass AssembledEvidence(observations, unmapped, reports, reliability_fit)`; `assemble(...)` gains `weighting: str = "grade"`, `rubric_weight_defaults: RubricWeightDefaults | None = None`, `reliability_global_scale: float = 1.0` and returns AssembledEvidence (the two existing callers, main and any test, are updated); `aggregate_reports_to_observations(reports: pd.DataFrame, grade_a_policy: str, cluster_by_gene: dict[str, str], metabolic_symbols: set[str]) -> pd.DataFrame`; `evidence_record_from_positive_reports(pair_reports: pd.DataFrame) -> EvidenceRecord`; `fit_report_reliability`; `summarize_reports(reports, fit) -> dict`; `summarize_observations` gains counts pairs_all_negative and pairs_conflicting and keeps every existing key.
- experiment_data.py: the two-line filter of section 2.2.
- mechanism_cards.py: `evidence_limitations_for_symptom(reports: pd.DataFrame, symptom: str) -> list[str]`, printed under a "Supporting evidence" heading when the reports table is passed; optional, cards without it render as today.

## 5. Tests

New file tests/test_evidence_reports_and_aggregation.py, synthetic fixtures only:

- test_not_qualified_and_excluded_annotations_become_negative_reports: a phenotype.hpoa fixture with a NOT row and a row with frequency HP:0040285 for a graph gene's disease; both reports have report_value 0, frequency 0.0 for the Excluded row and null for the NOT row.
- test_evidence_code_references_and_date_join_from_phenotype_hpoa: a (disease, term) key with a PCS row citing two PMIDs and a TAS row citing the OMIM entry yields two reports with evidence_code PCS and TAS, pubmed_reference_count 2 and 0, rubric_curated_synopsis 0.0 and 1.0, evidence_date from the earliest publication date for the PCS row and from the biocuration date for the TAS row; an Orphanet row yields a null date.
- test_unjoined_rows_fall_back_to_genes_to_phenotype: without phenotype.hpoa every report has evidence_code "unjoined" and report_value follows the genes_to_phenotype frequency.
- test_default_weighting_reproduces_grade_weights: on the fixture of tests/test_monogenic_loader_frequency_and_exclusions.py, weights and grades from aggregate_reports_to_observations under weighting grade equal loss_weight_for_record and assign_evidence_grade on monogenic_evidence_records for every pair with a positive report (HMBS 0.68, OTC 1.0), label_frequency equal, and the same for a SIDER fixture (tests/test_drug_label_and_target_loaders.py write_synthetic_sider: insomnia 0.055 mean of midpoints, placebo row excluded).
- test_all_negative_pair_is_kept_with_zero_weight: a gene whose only annotation is NOT-qualified appears in the observations with positive_report_count 0, negative_report_count 1, weight 0.0 and grade "C".
- test_conflicting_pair_posterior_lies_strictly_between: three synthetic sources, one item with one positive and one negative report versus the same item all positive and all negative under identical coverage; posterior_conflicting is strictly between the two.
- test_weighted_fit_reduces_to_dawid_skene_on_unit_weights: fit_weighted_dawid_skene with pos = reports, neg = 1 − reports equals fit_dawid_skene on the simulated three-source data of the existing test.
- test_rubric_weights_order_evidence_codes_and_sample_sizes: PCS > TAS > IEA at equal sample size; n = 1 < n = 20 at equal code; placebo_controlled_frequency > label_text; a JSON override changes one factor and an unknown key raises.
- test_weakly_identified_flag: two sources covering every item sets weakly_identified True; three sources covering one item sets it False.
- test_implicit_negative_only_inside_coverage: an item whose perturbation has no report from a source contributes nothing to that source's column; one with a report for another symptom contributes implicit_negative_weight.
- test_report_ids_are_unique_and_deterministic.

tests/test_splits_and_equifinality.py or a new tests/test_experiment_data_labels.py: test_experiment_data_never_labels_an_all_negative_pair_positive, writing a small evidence_records.parquet with positive_report_count 0 and 1 rows and asserting the outcome matrix is 0 for the former, and test_experiment_data_keeps_zero_weight_positives_without_report_counts for the version 0.3 table layout.

tests/test_monogenic_loader_frequency_and_exclusions.py: the existing tests keep passing; add test_annotation_rows_reader_matches_dates_and_references, asserting load_hpo_annotation_dates and load_hpo_annotation_references give the same dictionaries as before over load_hpo_annotation_rows.

## 6. Documentation changes in docs/experiment_design.md

Section 4.2, paragraph "Evidence class E1": replace the sentence "rows qualified as Excluded assert absence and are dropped" with: "Each annotation row is kept as a report (evidence_reports.parquet, one row per phenotype.hpoa annotation behind a gene-symptom pair, with its evidence code, references, onset, sex, frequency and date); rows qualified NOT or with the Excluded frequency assert absence and are reports with value 0. In the 2026-09-01 release no row carries the Excluded qualifier and 826 genes_to_phenotype rows are NOT-qualified in phenotype.hpoa, which genes_to_phenotype.txt does not mark; three of them fall on the monogenic slice, where one pair (CPT1C, cognitive impairment) rests on a NOT annotation alone." After "An observation does not account for the number of publications behind it." add: "The aggregated observation records report_count, positive_report_count and negative_report_count; a pair with no positive report is kept with weight 0 and is never a positive label."

Section 4.3, first paragraph: after "which is why the model is specified this way." add: "As implemented, each explicit report enters the likelihood with a rubric weight in (0, 1] from a fixed function of its evidence code and sample size (PCS 1.0, TAS 0.8, IEA 0.6, other 0.7; placebo-controlled label frequency 1.0, label frequency 0.9, label text 0.8; sample-size factor 0.5 + 0.5 log1p(n)/log1p(20) when a patient count is reported; overridable by a JSON file), a source's silence on a covered perturbation enters as one implicit negative report of weight 1, and the posterior is written to every evidence table as reliability_posterior. It becomes the weight only under --weighting reliability; the default remains the grade table below (section 12, open question 12)."

Section 10, under evidence/, add two lines: `evidence_reports.py (report table, rubric weights, limitations text)` and, under data/processed/, `evidence_reports.parquet` with the comment "(one row per report; evidence_records.parquet aggregates it)".

Section 12, paragraph "Evidence weighting as implemented", replace with: "Evidence weighting as implemented. The assembler writes one row per report (evidence_reports.parquet) and aggregates to one row per pair; the default weights still come from the grade table of section 4.3 and are unchanged for 860 of the 861 slice pairs, the exception being a pair whose only annotation is NOT-qualified (CPT1C, cognitive impairment), which now carries weight 0 and no positive label. The Dawid-Skene fit runs on every table over rubric-weighted reports and its posterior is a column; with two HPO provenances and no source covering both genes and drugs it is weakly identified (Orphanet specificity 0.03 on the slice) and is not the default weight (open question 12)."

Open question 12, replace with: "12. Identifiability of the learned weighting. With the sources at hand the sensor model sees two HPO provenances on gene items and one label source on drug items, and no item is covered by three sources: it learns OMIM-versus-Orphanet agreement and little else, drives Orphanet specificity to 0.03 on the slice and gives a NOT-only pair a posterior of 0.86. The report table and the weighted fit are in place; a third conditionally independent source covering both genes and drugs, the literature class E3, is the prerequisite for making the posterior the default, and the implicit-negative weight for silence inside a curated synopsis (default 1.0) is the parameter to revisit first when it is. Until then --weighting reliability is a pre-registered arm, judged against the label-permutation and degree-stratified controls of section 6.3 like the multiplicity ablation."

## 7. Subtleties the spec settles

- Descendant expansion: several reports per disease per pair, counted in report_count, named in limitations ("rests on one disease entry") and not treated as independent (section 1.1): inside the reliability likelihood one disease entry (or one drug) contributes at most one report's weight per cell and side (section 3, decision 15).
- Version 0.3 proxies: independent_case_series_count and has_omim_clinical_synopsis are computed from positive reports (section 2), so grade_a_policy two_distinct_disease_entries behaves as today for every pair with a positive report.
- Max-weight collapse in experiment_data.py: unchanged; with one row per pair per relation it takes the row's weight, and the new filter drops only rows without a positive report.
- Monogenic frequency is the maximum over positive reports, SIDER frequency the mean over frequency reports: both as today, both stated.
- NOT rows inside one disease entry with a positive descendant row (GLRX5) are a conflicting pair, not an error; the summary counts them.
- The reliability fit never creates an item: unobserved pairs stay absent, and a source that covers an item only contributes an implicit negative when it has some report for that perturbation and relation.

## Decisions during implementation

Implemented 2026-10-03. The open points the specification left and the choices made, each the simplest one consistent with the repository rules; the checks were run on data/processed/graph_full against data/processed/evidence_full.

1. CPT1C cognitive_impairment (1.0 to 0.0). Implemented as the specification says: the pair stays in evidence_records.parquet with grade C, weight 0.0, positive_report_count 0 and negative_report_count 1, and experiment_data.py no longer labels it positive (3,577 to 3,576 positive cells on the full graph; the gene leaves the perturbation list because it has no other pair). Every other monogenic row of the full table keeps its grade and weight bit for bit (3,033 of 3,034). This is the one intended deviation from the running sweeps and is treated as a bug fix; the orchestrator should confirm before sweeps that read a rebuilt table are compared with the running ones.
2. annotation_patient_count now sums the phenotype.hpoa denominators of the positive reports, so it changes where genes_to_phenotype collapsed several rows (full graph: APOE 6 to 18, DNMT1 8 to 26, PSEN2 24 to 36). Descriptive only; accepted. label_frequency likewise becomes the largest frequency over the hpoa rows and changes on one full-graph row (DNMT1 cognitive_impairment 0.75 to 1.0) without changing its weight, since both values clip to the full weight.
3. implicit_negative_weight stays 1.0 and is a field of RubricWeightDefaults, overridable through --report-rubric-weights, as asked.
4. The placebo-arm frequency is not a column; placebo_flag records that the pair has a placebo row. A treatment frequency at or below the placebo frequency is not counted as evidence against (left for the appraisal work).
5. configs/evidence_assembly.yaml is still not read by the assembler; its weighting block now carries a comment saying the command line is authoritative and naming the --weighting reliability arm.
6. Tests were run in this task: the touched tests plus tests/test_end_to_end_toy_pipeline.py and tests/test_evidence_reliability_and_appraisal.py, then the full suite once (see the task report).

Choices the specification did not foresee:

7. report_id carries the perturbation id. The specification's HPO format `source|disease|term|symptom|ordinal` collides whenever one disease entry is annotated to several graph genes, which is the rule for HPO annotations, so the id is `source|perturbation_id|source_record_id|source_term_id|symptom|ordinal` for every source (for a drug the perturbation id and the source record id are both the STITCH flat id). test_report_ids_are_unique_and_deterministic pins it.
8. Several preferred terms per drug pair. The earlier table wrote one row per (drug, preferred term, relation), so 804 drug rows stood for 572 pairs (431 duplicate rows on the full graph), and experiment_data.py took the maximum weight and the maximum frequency over them. The aggregation reproduces that view exactly: each preferred term keeps its own label frequency (mean of its treatment midpoints, as load_sider_frequencies always gave), the pair's weight is the largest grade weight over its terms and label_frequency the largest term frequency (pharmacological_label_frequency). A mean over all frequency reports of the pair would have changed 20 or so weights (a term without a frequency keeps the full base weight 0.6 while a rare term scales it down). All 572 drug pairs match the collapsed earlier table on grade, weight, label_frequency, group_id and source.
9. aggregate_reports_to_observations takes metabolic_node_ids (graph node ids) rather than gene symbols, because the report table carries node ids in perturbation_nodes and no target gene symbols for drugs, and an optional drug_targets_by_perturbation mapping ("CHEMBLid:ACTION;..." sorted by target id) from which a drug's group_id and source line are built. The disease clusters come from disease_cluster_ids_from_reports, which feeds disease_cluster_ids the disease entries of every report of a gene, positive and negative, so clusters do not move.
10. Rows with a 0/N frequency (MME cognitive_impairment, OMIM:617017, 0/10) were first written as reports with value 1 and frequency 0.0, as the specification's rule (value 0 only for the NOT qualifier or HP:0040285) said, with the grade weight at the floor 0.25 as before. Superseded by decision 14: they are absence reports.
11. AssembledEvidence has a fifth field, reports_dropped_unknown_provenance, and summarize_reports takes the rubric defaults, the observations and that count as further arguments, so the summary can report the dropped rows and the posterior quantiles per pair category. A fifth category, covered_by_one_source_positive, is reported next to the four the specification names because it is the largest.
12. The ordinal of a joined report is 1-based within its (disease, term) group or (drug, preferred term) treatment-arm group, so 0 is reserved for unjoined and label_text reports as the specification requires.
13. Before decisions 14 to 16, the free fit on the full graph converged in 97 iterations to prevalence 0.91, HPO-OMIM sensitivity 0.68 and specificity 0.03, HPO-Orphanet 0.998 and 0.90, SIDER-label 0.999 and 0.07, SIDER-indication 0.97 and 0.82: with Orphanet covering 2,386 gene items and OMIM 1,593, the two-source fit made Orphanet the reference sensor and OMIM silence inside Orphanet coverage the uninformative event, the mirror image of the slice prototype in section 3.4. The posterior of the all-negative pair was 0.07 and the four conflicting pairs lay between 0.05 and 0.83. The review found what those numbers mean: OMIM sensitivity + specificity was 0.71, so an explicit OMIM positive report had a log likelihood ratio of −0.35 per unit weight and OMIM silence +2.24, 180 pairs with an OMIM positive report and a silent Orphanet column sat at 0.012 to 0.055, and the SIDER-label specificity 0.07 was estimated from no negative at all. Decisions 14 to 18 answer these; the flag weakly_identified is True, as predicted, and the default weighting is unchanged.

Review of 2026-10-03 (second round), each applied to the working tree and checked on data/processed/evidence_reports_check (full graph) and data/processed/evidence_reports_check/slice:

14. 0/N rows are absence reports. A phenotype.hpoa row whose frequency is a patient fraction with numerator 0 (TYMP and POLG OMIM:603041 cognitive_impairment 0/35, MME OMIM:617017 0/10) is now report_value 0 with frequency 0.0 and the denominator kept, so the sample-size factor applies to the absence claim (0/35 weighs 1.0, 0/10 weighs 0.89); limitations_text adds "reported in 0 of N patients". The full graph has 31 such reports, all HPO-OMIM, and the summary counts them under HPO-OMIM|0; MME cognitive_impairment is now an all-negative pair (0/10 against an Orphanet NOT) rather than a conflicting one, and the full graph counts 13 all-negative and 21 conflicting pairs (was 1 and 4; the 17 new conflicting pairs each have a 0/N row next to a positive row). The grade weighting keeps reading these rows as the version 0.4 loader did (grade_weighting_rows: positive reports plus 0/N rows), so grade and weight are identical to the running tables on all 3,606 full-graph rows and all 1,137 slice rows (the floor weight 0.25 for a pair whose rows are all 0/N). positive_report_count counts report_value 1 only, so a pair whose rows are all 0/N is no longer a positive label in experiment_data.py: 12 pairs on the full graph (CAMK2A irritability_or_aggression, CHCHD2 cognitive_impairment and psychosis, COQ7, KIF1C, LAMP2, MME, PNPT1 and SAMD9L cognitive_impairment, PLD4 fatigue, PPFIA3 depressed_mood, PPP2R5C irritability_or_aggression, weight 0.25 each) and one on the slice (COQ7 cognitive_impairment). This is the same class of change as CPT1C in decision 1 and is flagged for the orchestrator in the same way; the weight column itself is unchanged, and a rebuilt table differs from the running ones only in those labels. 0% (no denominator) stays a frequency-0 presence row, as it has no patient count behind it.
15. One report's worth per source record. report_count_matrices reduces the rubric weights inside an (item, source) cell per source_record_id to the largest weight among the record's positive reports and, separately, among its negative reports, then sums over distinct records; a table without a source_record_id column treats every row as its own record, which keeps the synthetic tests as they were. On the full graph 763 of 4,577 positive cells held several reports of one record (a synopsis annotated to Dementia, Memory impairment and Mental deterioration; several labels of one drug), the summed positive weight per cell reached 60 for SIDER-label and the 90th percentile was 2.4; now the SIDER-label cell weight is at most 1.0 and HPO-Orphanet's at most 6.4 (distinct Orphanet entries on one pair). The explicit_positive_weight of the summary is the reduced weight (948 for HPO-OMIM instead of 1,093, 485 for SIDER-label instead of 1,592).
16. Three guards on the fit, in fit_weighted_dawid_skene and fit_report_reliability (section 3). Prevalence per evidence class (item_groups), because no source covers both gene and drug items. Specificity held at the prior mean 0.8 for a source with no explicit negative report (SIDER-label and SIDER-indication on both tables; the previous SIDER-label value 0.07 was the specificity M-step reading the volume of positive reports). Specificity held at the prior mean for a source whose estimated sensitivity + specificity ends at or below 1, with EM rerun until no estimated source is inverted, and sources still at or below 1 after the holds listed as degenerate_sources, on which reliability_weights refuses --weighting reliability with ValueError. The summary records prevalence_by_evidence_class, degenerate_sources and, per source, specificity_held_at_prior with its reason and degenerate. Tests: test_silence_only_source_keeps_prior_specificity_and_never_counts_a_positive_against, test_inverted_source_is_held_at_prior_then_flagged_and_reliability_weighting_refused, test_prevalence_is_fitted_per_item_group, test_one_source_record_contributes_one_report_per_cell, test_zero_of_n_pair_keeps_the_grade_weight_but_is_not_a_positive_label and the extended test_not_qualified_and_excluded_annotations_become_negative_reports.
17. The ordering claim of section 3.4 (a conflicting pair's posterior lies strictly between the all-positive and the all-negative value) is restated with its condition, sensitivity + specificity > 1 for every covering source, which the guards of decision 16 enforce on every fit that --weighting reliability accepts.
18. What the guarded fit does on the tables written. Full graph (evidence_reports_check): the free pass inverts HPO-OMIM (sensitivity 0.63, specificity 0.04); with OMIM held, the second pass inverts HPO-Orphanet instead (0.88, 0.01); with both held the fit converges in 28 iterations to HPO-OMIM sensitivity 0.685, HPO-Orphanet 0.923, SIDER-label 0.998, SIDER-indication 0.968, every specificity 0.8, prevalence 0.9992 (monogenic) and 0.9977 (pharmacological), no degenerate source, and every posterior between 0.985 and 1.0: all-negative pairs 0.990 to 0.999 (CPT1C 0.995), conflicting 0.998 to 1.0, drug pairs 0.9993 to 0.9995. Slice (evidence_reports_check/slice): the free pass inverts HPO-OMIM; with OMIM held, HPO-Orphanet stays estimated at 0.85 (explicit negative weight 2.4 against 56 silence cells) and the fit converges in 57 iterations to prevalence 0.998, with every posterior between 0.986 and 1.0 (CPT1C 0.986). So once a sign inversion is ruled out, the two-source fit lands on its other weakly identified optimum: items exist only when some source reported them, the prevalence among them is estimated near 1 and the prior log odds of about 7 outweigh one NOT report (0.8 × log(0.077 / 0.8) = −1.9 nats). The posterior column no longer counts a positive report against a pair, and conflicting and all-negative pairs sit below all-positive ones, but it ranks pairs within a band of 0.015 and is not evidence weighting in any useful sense; the same conclusion as open question 12, now with the degenerate optimum named for each table. Whichever of the two HPO provenances is estimated last is inverted by the other's silence, which is the non-identifiability showing up as a sign.
19. The specificity hold is reported, not silent: specificity_held_at_prior in the per-source summary carries the reason, and the design's open question 12 and section 12 cite the numbers of the tables actually written rather than the prototype's.
20. configs/evidence_assembly.yaml is unchanged; the guards have no command-line switch, since turning them off reproduces the inverted fit of decision 13.
