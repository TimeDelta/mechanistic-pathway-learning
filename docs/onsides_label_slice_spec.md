# OnSIDES later label slice: identifier bridge, report schema and qualification rule

Specification of evidence class E2's later slice (design section 4.2), as implemented in mechanistic_pathway_learning/evidence/onsides_identifier_bridge.py, mechanistic_pathway_learning/evidence/load_onsides_label_events.py, experiments/fetch_onsides_identifier_bridge.py and experiments/build_onsides_reports.py. The generated page docs/onsides_label_slice.md carries the per-symptom, per-country and per-section tables of the current build; this page says what the columns and counts mean. Counts here are from the build of 5 October 2026 on OnSIDES v3.1.1 (data release 2026-04-22), ChEMBL 37 and SIDER 4.1. MedDRA term names are never written to a committed file; identifiers and target symptoms only.

## 1. Inputs

OnSIDES v3.1.1 parquet tables under data/raw/onsides/v3.1.1/parquet: product_adverse_effect (6,928,666 label-event mentions with effect_meddra_id, label_section, match_method and pred1), product_label (41,119 labels with source US, UK, EU or JP), product_to_rxnorm, vocab_rxnorm_ingredient_to_product, vocab_rxnorm_ingredient (1,866 ingredients keyed by RxCUI or an OMOP extension id) and vocab_meddra_adverse_effect (6,423 terms, 3,550 PT and 2,873 LLT). ChEMBL caches under data/raw/chembl: mechanisms.json, targets.json (extended in place by the bridge fetch), molecule_parents.json and pubchem_to_chembl_with_parents.json (SIDER's own ChEMBL mapping). SIDER 4.1 under data/raw/sider_4.1: drug_names.tsv, drug_atc.tsv and meddra_all_se.tsv.gz.

## 2. Identifier bridge (data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json)

One record per ingredient. Routes to a ChEMBL parent, in precedence order; the first that yields a parent is recorded as bridge_method:

| bridge_method | route | ingredients |
|---|---|---|
| unii_unichem | RxCUI -> UNII (RxNav property UNII_CODE) -> ChEMBL ids (UniChem, FDA SRS source 14) -> parent (molecule_hierarchy); an OMOP id never takes this route | 1,172 |
| name_pref_name | ingredient name equals a ChEMBL pref_name ignoring case -> parent | 448 |
| none | no parent | 246 |

When a route returns several parents the lowest-numbered id is taken and the rest are kept in chembl_parent_alternatives (1 ingredient). 1,620 ingredients have a parent (1,609 distinct), 1,608 a UNII, 53 are OMOP extension ids.

ATC classes (atc_codes_rxnav) come from the RxClass byRxcui response. Entries whose minConcept term type is IN (the ingredient) or PIN (a precise ingredient: a salt or ester such as olmesartan medoxomil under olmesartan, where RxClass attaches many drugs' classes) are kept; MIN entries (multi-ingredient concepts with their own codes) are dropped; an untyped entry is kept only when its RxCUI is the queried one. atc_codes_sider adds drug_atc.tsv codes of the unified SIDER drug; atc_codes is the sorted union. 1,523 ingredients have a code, 252 an ATC N code.

Unification with SIDER (unification_method), in precedence order:

| unification_method | rule | ingredients |
|---|---|---|
| chembl_parent | the ingredient's parent is among the ChEMBL parents of a SIDER STITCH flat identifier (pubchem_to_chembl_with_parents.json resolved through molecule_parents.json and the mechanism table) | 842 |
| sider_drug_name | otherwise, the ingredient name equals a SIDER drug name ignoring case (drug_names.tsv); this reaches SIDER drugs whose own ChEMBL mapping is empty or resolves to a stereo or racemate variant of the ingredient's parent | 158 |
| none | no SIDER partner; perturbation_id is RXCUI:<id> or OMOP:<id> | 866 |

A unified ingredient's perturbation_id is the STITCH flat identifier of the lowest matching PubChem CID (CID1 plus the CID zero-padded to eight digits), the string SIDER reports carry, so the evidence table sees one perturbation with two label slices; the other matching CIDs are sider_collision_pubchem_cids (24 ingredients, 23 parents listed under sider_collisions). The join is parent-to-flat-CID, not parent-to-parent: a STITCH flat identifier drops stereochemistry, so one SIDER CID can carry several ChEMBL parents. sider_cid_chembl_parents lists every parent of the chosen CID and sider_cid_merges_distinct_parents is true when there is more than one (196 unified ingredients; in the qualifying rows 232 reports on 20 perturbations whose merged parents are enantiomer or racemate records of one drug, such as citalopram and escitalopram, bupivacaine and levobupivacaine, modafinil and armodafinil, eszopiclone and zopiclone or eslicarbazepine and licarbazepine; 4 of those perturbations carry two ingredients each). A merge of distinct drugs under one CID (betamethasone and dexamethasone, doxorubicin and epirubicin, ephedrine and pseudoephedrine) is therefore visible in the record and in the count unified_cid_merges_distinct_parents, and a later rule change that admits such a pair is caught there.

Mechanism targets: mechanisms.json is keyed by molecule_chembl_id and by parent_molecule_chembl_id (mechanisms_by_parent_and_molecule), the one keying the bridge fetch, the missing-target fetch and the report builder share, so a mechanism ChEMBL records on a salt form counts for the parent. mechanism_target_ids lists the human targets of the parent (organism Homo sapiens or unknown), single_mechanism_target applies has_dominant_target with one target and nervous_system_atc applies is_nervous_system_atc to atc_codes. The targets of every mechanism under a bridged parent are fetched into targets.json when missing (64 records added on 5 October 2026; 631 now); a target id still missing after the fetch is listed under errors and is reported by the builder as target_record_missing rather than as a missing graph gene.

Bridge counts of the current file: ingredients 1,866; with_any_mechanism_target 1,134; single_mechanism_target 892; single_target_and_atc_n 150; single_target_and_atc_n_onsides_only 28; unified_with_sider 1,000 over 973 distinct SIDER drugs. The name-bridge comparison (sider_name_bridge_comparison) records that the 861 ingredients whose name equals a SIDER drug name are all unified, 857 to a CID of that name and 4 to another CID through the parent route; 139 ingredients are unified without a name match.

## 3. Statements (aggregate_label_statements)

Term matching: a PT row of vocab_meddra_adverse_effect matches a target symptom when its name is a crosswalk preferred term (docs/symptom_crosswalk.csv). An LLT row matches through the PT it falls under, using the LLT-to-PT pairing SIDER lists per label term in meddra_all_se.tsv.gz (sider_lower_level_to_preferred_terms; 5,805 LLT names paired without redistributing MedDRA); without the pairing an LLT matches only on an exact name, which never happens in the v3.1.1 vocabulary. 21 OnSIDES LLT ids fall under crosswalk PTs this way; a MedDRA LLT-to-PT table kept outside the repository would be the complete solution. The statement keeps the source term type in meddra_term_type.

One statement per (ingredient_id, target_symptom, meddra_id, label_section, source_country) over the events of every label whose RxNorm products contain the ingredient: label_count distinct labels, single_ingredient_label_count of them with exactly one ingredient, labels_of_ingredient_in_country, max_pred1, mean_pred1 and match_methods (the distinct match_method values joined by ";"). Current build: 5,959 statements, 5,649 PT-coded and 310 LLT-coded.

Label sections exist for US labels only: AR (Adverse Reactions), BW (Boxed Warning) and WP (Warnings and Precautions). Every UK, EU and JP event carries label_section NA, OnSIDES' value for the undivided non-US label whose extraction covers the undesirable-effects section; no US event does (crosstab of the statements: EU 1,550, JP 671 and UK 1,288 NA rows, US 2,216 AR, 20 BW and 214 WP rows).

## 4. Report rows (data/processed/onsides/onsides_reports.parquet)

The 35 shared columns of docs/evidence_reports_spec.md come first (REPORT_COLUMNS), then six OnSIDES rubric columns, then ONSIDES_EXTRA_COLUMNS. Every statement produces a row; qualifies marks the rows that enter the evidence table through assemble_evidence_table.py --extra-reports.

| column | value |
|---|---|
| report_id | `OnSIDES-label|{perturbation_id}|{ingredient_id}|{meddra_id}|{symptom}|{label_section}|{source_country}` |
| perturbation_id, perturbation_type, perturbation_label | the bridge's perturbation_id; "drug"; the SIDER drug name when unified, else the ingredient name |
| symptom, relation, evidence_class, source, report_value | target symptom; "induces"; "pharmacological"; "OnSIDES-label"; 1 (OnSIDES records no absences) |
| source_record_id, source_record_label | ingredient_id and ingredient name |
| source_term_id, source_term_label | meddra_id and the term name (blanked by the data release) |
| evidence_code | adverse_reactions_section (AR), boxed_warning_section (BW), warnings_precautions_section (WP), non_us_label_section (NA), other_section (any other value; none in v3.1.1) |
| model_description | `human; drug label (OnSIDES v3.1.1, {country}); ChEMBL targets {target_chembl_id}:{action_type};...` |
| frequency, frequency_denominator, placebo_flag, onset, sex, references, pubmed_reference_count | null, null, false, null, null, "", 0 |
| evidence_date | 2026-04-22, the release date, an upper bound on when the statement appeared |
| rubric_log_sample_size, rubric_frequency_known, rubric_evidence_code_pcs, rubric_evidence_code_tas, rubric_evidence_code_iea, rubric_placebo_controlled, rubric_curated_synopsis, rubric_causal_association | 0.0 |
| association_type, evidence_date_source | "" |
| limitations | section 6 |
| perturbation_nodes | JSON list of [node_id, sign, 1/len(mapped)] per gene node of the single target |
| rubric_log_label_count | log1p(label_count) |
| rubric_boxed_warning | 1.0 for BW (US labels only) |
| rubric_adverse_reactions_section | 1.0 for AR and for NA: the non-US extraction covers the adverse-reactions text, so the feature is not a second copy of rubric_us_label |
| rubric_nlp_score | max_pred1 scaled to [0, 1] over the model-scored statements (match_method PMB); 0.0 for string matches |
| rubric_us_label | 1.0 for source_country US |
| rubric_atc_known | 1.0 when the ingredient has any ATC code (1.0 on every qualifying row under the rule below) |
| ingredient_id, meddra_term_type, label_section, source_country, label_count, single_ingredient_label_count, labels_of_ingredient_in_country, max_pred1, mean_pred1, match_methods, chembl_parent, unified_with_sider | the statement and bridge fields of sections 2 and 3 |
| qualifies, disqualified_reason | section 5 |

## 5. Qualification rule (qualify_ingredient and onsides_reports)

The rule SIDER rows pass in assemble_evidence_table.py, applied per ingredient, then two row-level tests; the first failing test names disqualified_reason:

| disqualified_reason | test | statements (current build) |
|---|---|---|
| no_chembl_parent | the bridge has no ChEMBL parent | 479 |
| no_single_mechanism_target | the parent's human mechanism targets are not exactly one (has_dominant_target, max_targets 1) | 2,125 |
| no_atc_code | the ingredient has no ATC code at all; there is no waiver, since a SIDER drug without a code does not qualify either | 130 |
| atc_not_nervous_system | no ATC N code (is_nervous_system_atc, the version 1 proxy for central action) | 1,985 |
| combination_label_only | every label behind the statement is a combination product (single_ingredient_label_count 0): the label of carbidopa/levodopa says nothing about carbidopa alone | 96 |
| target_record_missing | the target has no record in targets.json, so no gene symbol can be mapped | 0 |
| no_graph_node | the target's gene symbols are not graph nodes | 0 |

Qualifying: 1,144 reports from 127 ingredients and 123 perturbations (102 shared with SIDER, 21 OnSIDES-only); by unification 100 ingredients through the parent route, 6 by drug name and 21 without a SIDER partner.

## 6. Limitations text

Clauses joined with "; ": "label statement without frequency"; "{label_count} label(s) in {country}"; the section ("adverse reactions section", "boxed warning", "warnings and precautions section", "non-US label, section not recorded by OnSIDES" or "unsectioned mention" for any other value); "combination-product labels only" when single_ingredient_label_count is 0; "string match without model score" or "extraction score {score} of 1"; "no ATC code known, central action not tested" when the ingredient has no ATC code; "dated by the OnSIDES v3.1.1 release (2026-04-22), an upper bound".

## 7. Pair-level overlap with SIDER (build_onsides_reports.py)

Two SIDER sides are measured for the qualifying OnSIDES (perturbation_id, symptom) pairs. The qualified side is the SIDER-label pairs of the assembled table (data/processed/evidence_full/evidence_reports.parquet, source SIDER-label, relation induces): the 543 pairs of SIDER drugs that passed SIDER's own ChEMBL mapping, single-target and ATC N rules, which the time split trains on. The raw side is every PT adverse-event row of SIDER 4.1 whose term is in the crosswalk (5,248 (identifier, symptom) pairs), looked up for every STITCH flat identifier the bridge links to the perturbation (the unified CID and its collisions), independent of SIDER's qualification.

Current build: 587 OnSIDES pairs; 283 in both; 260 SIDER only; 304 absent from the qualified table, of which 204 are on 2015 SIDER labels and 100 are later-slice candidates (41 on perturbations SIDER also holds, 59 on OnSIDES-only perturbations). The design's time split (section 6.1) trains on the SIDER 2015 pairs and tests on the later-slice candidates, not on every pair absent from the qualified table; a candidate may still have been on a 2015 label that SIDER's extraction missed, and OnSIDES rows are dated by the release.

## 8. Known limits

The stereo-flattened SIDER identifier can pool stereoisomers of one drug under one perturbation (section 2); a measured blood-brain barrier flag is still to replace the ATC N proxy; the LLT pairing covers only the terms SIDER listed; the assembled table under data/processed/evidence_full was built before the rule changes of 5 October 2026 (ATC waiver removed, combination-only statements excluded, LLT pairing, name unification) and is to be re-assembled with the new onsides_reports.parquet.
