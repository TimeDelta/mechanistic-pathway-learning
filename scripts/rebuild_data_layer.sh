#!/bin/bash
# Rebuild data/processed in a new container, in dependency order, from the pinned downloads.
#
# data/raw, data/processed and runs/ are not in git, so a new container starts with none of them. This is the chain
# that was run on 10 October 2026 to get the confirmatory inputs back, with the commands the documents left out
# (docs/data_layer_rebuild.md gives what each step reproduced and what it did not). Every step is skipped when its
# output exists, so the script can be rerun after a stop. Run it from the repository root:
#
#   bash scripts/rebuild_data_layer.sh            # graphs, evidence tables, label selections, drug entry nodes, and the two
#                                                   # evidence sources that record absences (trial arms, FAERS)
#   bash scripts/rebuild_data_layer.sh slice        # the monogenic slice tables the slice-era tests read (minutes)
#   bash scripts/rebuild_data_layer.sh descriptors  # then the node descriptors and cell-class weights (about 3 hours
#                                                   # of ESM-2 embeddings on 2 cores, and a 1.24 GB download)
#   bash scripts/rebuild_data_layer.sh fine         # the v5 tables on the fine symptom constructs (decision 9.58; an hour)
#
# Library versions matter for the file hashes: the tables of data/releases/v0.4 were written by pandas 2.3.3 and
# pyarrow 25.0.1, and with those two the rebuilt tables that have a recorded hash and no fitted column match it
# (docs/data_layer_rebuild.md). Install them before running:  pip install "pandas==2.3.3" "pyarrow==25.0.1"
set -euo pipefail
repository_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repository_root"
export PYTHONPATH=. OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
stage="${1:-graphs}"
processed="data/processed"
log_directory="$processed/rebuild_logs"
mkdir -p "$log_directory"

step() {
  # step OUTPUT_PATH DESCRIPTION COMMAND...   runs the command unless OUTPUT_PATH exists
  local output_path="$1" description="$2"
  shift 2
  if [ -e "$output_path" ]; then
    echo "present  $description ($output_path)"
    return 0
  fi
  echo "building $description"
  "$@" > "$log_directory/$(basename "$output_path").log" 2>&1 || { echo "FAILED   $description; see $log_directory/$(basename "$output_path").log" >&2; exit 1; }
}

# --- raw sources ----------------------------------------------------------------------------------------------------
bash scripts/fetch_raw_sources.sh
step data/raw/onsides/v3.1.1/parquet/product_adverse_effect.parquet "OnSIDES tables as parquet" python scripts/convert_onsides_to_parquet.py
# UniChem and ChEMBL answer within a time budget per call; the mapped count stops rising after two or three passes
if [ ! -e data/raw/chembl/targets.json ] || [ ! -e data/raw/chembl/.complete ]; then
  for pass in 1 2 3 4; do python experiments/fetch_chembl_drug_targets.py | tail -1; done
  touch data/raw/chembl/.complete
fi
step data/raw/uniprot/uniprot_human_reviewed.tsv.gz "UniProt proteome, GO, hp-base.owl and ChEBI" python experiments/fetch_node_descriptor_sources.py
step data/raw/reactome/v97/pathways.json "Reactome release 97 pathway exports" python experiments/fetch_reactome_pathways.py
step data/raw/uniprot/uniprot_plasma_binders.tsv.gz "UniProt plasma binder annotations" python experiments/fetch_plasma_binder_annotations.py
step data/raw/plasma_protein_binding/fraction_unbound_database.xlsx "fraction-unbound database" python experiments/fetch_fraction_unbound_database.py
# about 45 minutes uncached (RxNav, UniChem and ChEMBL at a few requests a second); it also extends data/raw/chembl/targets.json,
# which the gene and protein split reads, so it comes before the split
step data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json "OnSIDES identifier bridge" python experiments/fetch_onsides_identifier_bridge.py

# --- graphs ---------------------------------------------------------------------------------------------------------
# The base graph is built, not copied from the release (16 October 2026, decision 9.53). The copy kept the release's
# row order but also kept its physiology: the physiology audit found the copied graph without the 15,853 mirror edges'
# removal (item 27, decision 4.10), the 33 GRK and arrestin sign corrections, the receptor-to-G-alpha coupling filter
# (9.23), the currency-tag fixes and the stoichiometric coefficients (item 20), all built after the release of 6 October.
# The builder now orders a reaction's catalysing genes by id, so the build is reproducible across containers.
step "$processed/graph_full/nodes.parquet" "base graph (build_physiology_graph, deterministic order)" \
  python -m mechanistic_pathway_learning.graph.build_physiology_graph --output-dir "$processed/graph_full"
step "$processed/literature/gene_identifier_map.parquet" "gene identifier map" python -m mechanistic_pathway_learning.evidence.gene_identifier_map
step "$processed/graph_full_manganese/nodes.parquet" "manganese variant" python experiments/build_manganese_graph_variant.py --graph-dirs "$processed/graph_full"
step "$processed/graph_full_neuronal/nodes.parquet" "neuronal and oxidative variant" \
  python experiments/build_neuronal_graph_variant.py --base-graph-dir "$processed/graph_full_manganese" --output-dir "$processed/graph_full_neuronal"
if [ "$stage" = "graphs" ]; then
  # without the descriptor tables; the descriptors stage builds the split again with them, and its nodes and edges are the same
  step "$processed/graph_full_neuronal_split/nodes.parquet" "gene and protein split" \
    python experiments/build_gene_protein_split.py --graph-dir "$processed/graph_full_neuronal" --output-dir "$processed/graph_full_neuronal_split"
  step "$processed/graph_full_neuronal_split_binders/nodes.parquet" "plasma binder variants" \
    python experiments/build_plasma_binder_variant.py --graph-dirs "$processed/graph_full_neuronal" "$processed/graph_full_neuronal_split" --markdown-output none
fi

# --- evidence and label selection -----------------------------------------------------------------------------------
# docs/symptom_crosswalk.csv holds parkinsonism, the 24th symptom, so these commands build the table with it. The table
# docs/drug_targets_any_type.md calls v3 is the same build on the crosswalk of commit 821cef5.
evidence_build() {
  # evidence_build SUFFIX CROSSWALK
  local suffix="$1" crosswalk="$2"
  step "$processed/onsides_v3$suffix/onsides_reports.parquet" "OnSIDES reports v3$suffix" \
    python experiments/build_onsides_reports.py --crosswalk "$crosswalk" --output-dir "$processed/onsides_v3$suffix" \
      --markdown-output "$log_directory/onsides_v3$suffix.md" --evidence-dir "$processed/evidence_full_v3$suffix"
  step "$processed/evidence_full_v3$suffix/evidence_records.parquet" "evidence table v3$suffix" \
    python -m mechanistic_pathway_learning.evidence.assemble_evidence_table --crosswalk "$crosswalk" --graph-dir "$processed/graph_full" \
      --output-dir "$processed/evidence_full_v3$suffix" --extra-reports "$processed/onsides_v3$suffix/onsides_reports.parquet" \
      --onsides-bridge data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json
  step "$processed/label_selection/better_v2_full_v3$suffix.parquet" "label selection better_v2 on v3$suffix" \
    python experiments/build_label_selection.py --evidence-dir "$processed/evidence_full_v3$suffix" --mask-grades C \
      --output "$processed/label_selection/better_v2_full_v3$suffix.parquet"
}
git show 821cef5:docs/symptom_crosswalk.csv > "$log_directory/symptom_crosswalk_before_parkinsonism.csv" 2>/dev/null \
  || echo "commit 821cef5 is not in this clone (git fetch --depth=400 origin main); the table without parkinsonism is skipped"
if [ -s "$log_directory/symptom_crosswalk_before_parkinsonism.csv" ]; then
  evidence_build "" "$log_directory/symptom_crosswalk_before_parkinsonism.csv"
fi
evidence_build "_parkinsonism" docs/symptom_crosswalk.csv

# --- evidence table v4: the label changes the user approved on 10 and 11 October 2026 ---------------------------------
# docs/label_source_review.md, sections 8, 9 and 13: a drug with at most three mechanism targets is admitted, a gene whose
# every causal association is a gain of function is seeded with sign +1, and a kept drug pair never more frequent on the
# drug than on placebo is set aside. No held-out set is drawn on this table yet and nothing registered reads it.
step "$processed/onsides_v4/onsides_reports.parquet" "OnSIDES reports v4" \
  python experiments/build_onsides_reports.py --crosswalk docs/symptom_crosswalk.csv --output-dir "$processed/onsides_v4" \
    --markdown-output "$log_directory/onsides_v4.md" --evidence-dir "$processed/evidence_full_v4" --max-drug-targets 3
step "$processed/evidence_full_v4/evidence_records.parquet" "evidence table v4" \
  python -m mechanistic_pathway_learning.evidence.assemble_evidence_table --crosswalk docs/symptom_crosswalk.csv --graph-dir "$processed/graph_full" \
    --output-dir "$processed/evidence_full_v4" --extra-reports "$processed/onsides_v4/onsides_reports.parquet" \
    --onsides-bridge data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json --max-drug-targets 3 --seed-gain-of-function-genes-positive
step "$processed/label_selection/better_v2_full_v4.parquet" "label selection better_v2 on v4, never-above-placebo pairs set aside" \
  python experiments/build_label_selection.py --evidence-dir "$processed/evidence_full_v4" --mask-grades C --mask-never-above-placebo \
    --output "$processed/label_selection/better_v2_full_v4.parquet"
# the alternative the user is deciding on (11 October 2026): the same pairs counted as negatives in place of being set aside
step "$processed/label_selection/better_v2_full_v4_placebo_negative.parquet" "label selection better_v2 on v4, never-above-placebo pairs as negatives" \
  python experiments/build_label_selection.py --evidence-dir "$processed/evidence_full_v4" --mask-grades C --never-above-placebo-as-negative \
    --output "$processed/label_selection/better_v2_full_v4_placebo_negative.parquet"

# --- drug entry nodes -----------------------------------------------------------------------------------------------
step data/raw/fda_labels/.complete "FDA label protein-binding sentences" bash -c "python experiments/fetch_fda_label_protein_binding.py --evidence-dir $processed/evidence_full_v3 && touch data/raw/fda_labels/.complete"
step data/raw/fda_labels/.complete_v4 "FDA label protein-binding sentences for the drugs v4 adds" bash -c "python experiments/fetch_fda_label_protein_binding.py --evidence-dir $processed/evidence_full_v4 && touch data/raw/fda_labels/.complete_v4"
# configs/drug_plasma_carriers.csv and configs/drug_plasma_carriers_v4.csv are committed. experiments/scope_label_plasma_binder_coverage.py rewrites
# each from the label sentences and the fraction-unbound database, with --evidence-dir naming the table whose drugs are read (the label cache holds
# more drugs than any one table). The variant holds the drug nodes, their mechanism edges and the carriers' sequestration edges;
# both arms of the ablation are read from it (run_main_model.py --drug-entry nodes or targets).
if [ "$stage" = "graphs" ]; then
  step "$processed/graph_full_neuronal_split_binders_drugs/nodes.parquet" "drug entry nodes on the confirmatory graph" \
    python experiments/build_drug_entry_node_variant.py --graph-dir "$processed/graph_full_neuronal_split_binders" --evidence-dir "$processed/evidence_full_v3_parkinsonism"
  step "$processed/graph_full_neuronal_split_binders_drugs_v4/nodes.parquet" "drug entry nodes on the confirmatory graph, for the v4 table" \
    python experiments/build_drug_entry_node_variant.py --graph-dir "$processed/graph_full_neuronal_split_binders" --evidence-dir "$processed/evidence_full_v4" \
      --carrier-table configs/drug_plasma_carriers_v4.csv --output-dir "$processed/graph_full_neuronal_split_binders_drugs_v4"
  # the same graph with a clearance edge from each drug's clearing enzyme (configs/drug_clearing_enzymes.csv); read by the comparative checks only
  step "$processed/graph_full_neuronal_split_binders_drugs_v4_enzymes/nodes.parquet" "drug entry nodes with clearance edges, for the v4 table" \
    python experiments/build_drug_entry_node_variant.py --graph-dir "$processed/graph_full_neuronal_split_binders" --evidence-dir "$processed/evidence_full_v4" \
      --carrier-table configs/drug_plasma_carriers_v4.csv --clearing-enzymes configs/drug_clearing_enzymes.csv \
      --output-dir "$processed/graph_full_neuronal_split_binders_drugs_v4_enzymes"
fi

# --- evidence sources that record absences (15 October 2026) --------------------------------------------------------
# Placebo-controlled trial arms (ClinicalTrials.gov) and FAERS disproportionality (openFDA) for the drugs of the v4
# table; every response is cached under data/raw and the scripts resume. They join the table through --extra-reports.
if [ "$stage" = "graphs" ]; then
  step "$processed/clinical_trials/clinical_trial_reports.parquet" "placebo-controlled trial arms (ClinicalTrials.gov)" \
    python experiments/fetch_clinicaltrials_adverse_events.py --evidence-dir "$processed/evidence_full_v4"
  step "$processed/faers/faers_reports.parquet" "FAERS disproportionality (openFDA)" \
    python experiments/fetch_faers_disproportionality.py --evidence-dir "$processed/evidence_full_v4"
  # IMPC knockouts, mapped to human genes through HGNC and to the graph's gene nodes; positives from the assertions,
  # negatives from the statistical results (tested, no significant call); evidence class animal, grade D
  step "$processed/impc/impc_reports.parquet" "IMPC knockout phenotypes (release 24.0)" \
    python experiments/build_impc_reports.py --graph-dir "$processed/graph_full"
fi

# --- the monogenic slice (pilot) tables the slice-era tests read ---------------------------------------------------
# data/processed/graph and data/processed/evidence are copied from the release, as graph_full is; the manganese,
# neuronal and gene/protein-split variants are built on the slice graph the same way as on the full one.
if [ "$stage" = "slice" ]; then
  step "$processed/graph/nodes.parquet" "slice graph, from data/releases/v0.4" cp -r data/releases/v0.4/graph "$processed/graph"
  step "$processed/evidence/evidence_records.parquet" "slice evidence table, from data/releases/v0.4" cp -r data/releases/v0.4/evidence_slice_sweeps "$processed/evidence"
  step "$processed/graph_manganese/nodes.parquet" "manganese variant of the slice" python experiments/build_manganese_graph_variant.py --graph-dirs "$processed/graph"
  step "$processed/graph_neuronal/nodes.parquet" "neuronal and oxidative variant of the slice" \
    python experiments/build_neuronal_graph_variant.py --base-graph-dir "$processed/graph_manganese" --output-dir "$processed/graph_neuronal"
  step "$processed/graph_split/nodes.parquet" "gene and protein split of the slice" \
    python experiments/build_gene_protein_split.py --graph-dir "$processed/graph" --output-dir "$processed/graph_split"
fi

if [ "$stage" = "descriptors" ]; then
  # --- node descriptors and cell-class weights ------------------------------------------------------------------------
  bash scripts/fetch_raw_sources.sh --with-cellxgene
  if [ ! -e "$processed/protein_embeddings/.complete" ]; then
    num_shards="$(nproc)"
    for shard in $(seq 0 $((num_shards - 1))); do
      OMP_NUM_THREADS=1 python experiments/compute_protein_embeddings.py --shard-index "$shard" --num-shards "$num_shards" > "$log_directory/esm_shard$shard.log" 2>&1 &
    done
    wait
    touch "$processed/protein_embeddings/.complete"
  fi
  descriptors="$processed/node_descriptors"
  step "$descriptors/protein_descriptors.parquet" "protein descriptors" \
    env OMP_NUM_THREADS=2 python experiments/build_protein_descriptors.py --slice-nodes data/releases/v0.4/graph/nodes.parquet --report "$log_directory/protein_descriptor_report.md"
  step "$descriptors/protein_descriptors_per_entry.parquet" "protein descriptors per UniProt entry" \
    env OMP_NUM_THREADS=2 python experiments/build_protein_descriptors.py --per-entry-only --slice-nodes data/releases/v0.4/graph/nodes.parquet --report "$log_directory/protein_descriptor_report.md"
  step "$descriptors/complex_descriptors.parquet" "descriptors of the Reactome protein entities" \
    env OMP_NUM_THREADS=2 python experiments/build_complex_descriptors.py --report "$log_directory/complex_descriptors.md"
  step "$processed/graph_full_neuronal/node_descriptors.parquet" "node descriptors of the merged neuronal graph (84 columns)" \
    python experiments/build_node_descriptors.py --graph-dir "$processed/graph_full_neuronal" --complex-descriptors "$descriptors/complex_descriptors.parquet"
  atlas="data/raw/siletti_cellxgene"
  step "$processed/brain_expression/dopaminergic_siletti_cluster395.parquet" "dopaminergic class from the brain atlas (872 nuclei of cluster 395)" \
    python experiments/build_dopaminergic_expression.py --h5ad "$atlas/dissection_sn_rn.h5ad" "$atlas/dissection_sn.h5ad" "$atlas/dissection_pag_dr.h5ad" "$atlas/dissection_pag.h5ad" \
      --output "$processed/brain_expression/dopaminergic_siletti_cluster395.parquet"
  step "$descriptors/full_neuronal_descriptors_brain_expression.parquet" "brain expression descriptors of the merged neuronal graph" \
    python experiments/build_brain_expression_descriptors.py --graph-dir "$processed/graph_full_neuronal" --output "$descriptors/full_neuronal_descriptors_brain_expression.parquet"
  weights="$processed/cell_class_weights"
  step "$weights/full_neuronal_cell_class_weights.parquet" "cell-class weights of the merged neuronal graph (about 8 minutes)" \
    python experiments/build_cell_class_weights.py --graph-dir "$processed/graph_full_neuronal" --output "$weights/full_neuronal_cell_class_weights.parquet"
  step "$processed/brain_expression/gene_expression_for_descriptors.parquet" "expression tables for the rewired runs" python experiments/write_expression_tables.py
  # the split again, now with its descriptor table and cell-class weights; nodes and edges are the ones the graphs stage wrote
  step "$descriptors/full_neuronal_split_descriptors_brain_expression.parquet" "gene and protein split with its descriptor table and cell-class weights" \
    python experiments/build_gene_protein_split.py --graph-dir "$processed/graph_full_neuronal" --output-dir "$processed/graph_full_neuronal_split" --overwrite \
      --descriptor-table "$descriptors/full_neuronal_descriptors_brain_expression.parquet:$descriptors/full_neuronal_split_descriptors_brain_expression.parquet" \
      --cell-class-weights "$weights/full_neuronal_cell_class_weights.parquet:$weights/full_neuronal_split_cell_class_weights.parquet"
  step "$processed/graph_full_neuronal_split_binders/nodes.parquet" "plasma binder variants" \
    python experiments/build_plasma_binder_variant.py --graph-dirs "$processed/graph_full_neuronal" "$processed/graph_full_neuronal_split" --markdown-output none
  # the drug entry variant with a row per drug node in the descriptor table (zeros) and in the cell-class weights (ones)
  step "$descriptors/full_neuronal_split_drugs_descriptors_brain_expression.parquet" "drug entry nodes with their descriptor table and cell-class weights" \
    python experiments/build_drug_entry_node_variant.py --graph-dir "$processed/graph_full_neuronal_split_binders" --evidence-dir "$processed/evidence_full_v3_parkinsonism" --overwrite \
      --descriptor-table "$descriptors/full_neuronal_split_descriptors_brain_expression.parquet:$descriptors/full_neuronal_split_drugs_descriptors_brain_expression.parquet" \
      --cell-class-weights "$weights/full_neuronal_split_cell_class_weights.parquet:$weights/full_neuronal_split_drugs_cell_class_weights.parquet"
fi
echo "done: $stage"

# --- evidence tables v5: the fine symptom constructs (decision 9.58, 16 October 2026) ---------------------------------
# docs/symptom_crosswalk_fine.csv is derived from the coarse crosswalk by experiments/build_fine_symptom_crosswalk.py
# (committed, held to a fresh build by tests/test_symptom_constructs.py). The same pipeline as v4, pointed at it: the
# label sources, the label selections, the three absence sources (their raw caches are per drug and are reused) and the
# release's literature table re-keyed to the namesake constructs (its raw cache is not on this disk). The drugs and
# their targets are those of v4, so the drug-entry graph of v4 serves the v5 tables.
if [ "$stage" = "fine" ]; then
  step "$processed/onsides_v5/onsides_reports.parquet" "OnSIDES reports v5 (fine constructs)" \
    python experiments/build_onsides_reports.py --crosswalk docs/symptom_crosswalk_fine.csv --output-dir "$processed/onsides_v5" \
      --markdown-output "$log_directory/onsides_v5.md" --evidence-dir "$processed/evidence_full_v5" --max-drug-targets 3
  step "$processed/evidence_full_v5/evidence_records.parquet" "evidence table v5 (fine constructs)" \
    python -m mechanistic_pathway_learning.evidence.assemble_evidence_table --crosswalk docs/symptom_crosswalk_fine.csv --graph-dir "$processed/graph_full" \
      --output-dir "$processed/evidence_full_v5" --extra-reports "$processed/onsides_v5/onsides_reports.parquet" \
      --onsides-bridge data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json --max-drug-targets 3 --seed-gain-of-function-genes-positive
  step "$processed/label_selection/better_v2_full_v5.parquet" "label selection better_v2 on v5, never-above-placebo pairs set aside" \
    python experiments/build_label_selection.py --evidence-dir "$processed/evidence_full_v5" --crosswalk docs/symptom_crosswalk_fine.csv --mask-grades C --mask-never-above-placebo \
      --output "$processed/label_selection/better_v2_full_v5.parquet"
  step "$processed/label_selection/better_v2_full_v5_placebo_negative.parquet" "label selection better_v2 on v5, never-above-placebo pairs as negatives" \
    python experiments/build_label_selection.py --evidence-dir "$processed/evidence_full_v5" --crosswalk docs/symptom_crosswalk_fine.csv --mask-grades C --never-above-placebo-as-negative \
      --output "$processed/label_selection/better_v2_full_v5_placebo_negative.parquet"
  step "$processed/clinical_trials_v5/clinical_trial_reports.parquet" "placebo-controlled trial arms at the fine level" \
    python experiments/fetch_clinicaltrials_adverse_events.py --crosswalk docs/symptom_crosswalk_fine.csv --evidence-dir "$processed/evidence_full_v5" \
      --output-dir "$processed/clinical_trials_v5" --markdown-output "$log_directory/clinical_trial_arm_evidence_v5.md"
  step "$processed/faers_v5/faers_reports.parquet" "FAERS disproportionality at the fine level" \
    python experiments/fetch_faers_disproportionality.py --crosswalk docs/symptom_crosswalk_fine.csv --evidence-dir "$processed/evidence_full_v5" \
      --output-dir "$processed/faers_v5" --markdown-output "$log_directory/faers_disproportionality_evidence_v5.md"
  step "$processed/impc_v5/impc_reports.parquet" "IMPC knockout phenotypes at the fine level" \
    python experiments/build_impc_reports.py --graph-dir "$processed/graph_full" --crosswalk configs/impc_symptom_crosswalk_fine.csv \
      --output-dir "$processed/impc_v5" --markdown-output "$log_directory/impc_knockout_evidence_v5.md"
  step "$processed/literature_v5/literature_reports.parquet" "literature reports re-keyed to the fine namesakes" \
    python experiments/rekey_literature_to_fine.py --input data/releases/v0.4/literature/literature_reports.parquet --output "$processed/literature_v5/literature_reports.parquet"
  # trial outcome measures (decision 9.60): the per-drug cache under data/raw/clinicaltrials_outcomes is reused; fine constructs
  step "$processed/clinical_trial_outcomes/clinical_trial_outcome_reports.parquet" "trial outcome measures on crosswalked scales" \
    python experiments/fetch_clinicaltrials_outcomes.py --evidence-dir "$processed/evidence_full_v4"
fi
