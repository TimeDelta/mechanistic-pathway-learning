#!/bin/bash
# Rebuild data/processed in a new container, in dependency order, from the pinned downloads.
#
# data/raw, data/processed and runs/ are not in git, so a new container starts with none of them. This is the chain
# that was run on 10 October 2026 to get the confirmatory inputs back, with the commands the documents left out
# (docs/data_layer_rebuild.md gives what each step reproduced and what it did not). Every step is skipped when its
# output exists, so the script can be rerun after a stop. Run it from the repository root:
#
#   bash scripts/rebuild_data_layer.sh            # graphs, evidence tables, label selections, drug entry nodes
#   bash scripts/rebuild_data_layer.sh descriptors  # then the node descriptors and cell-class weights (about 3 hours
#                                                   # of ESM-2 embeddings on 2 cores, and a 1.24 GB download)
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
# The base graph is copied from the release, not rebuilt: build_physiology_graph.py orders a reaction's catalysing genes
# by set iteration, so a rebuild has the release's rows in another order (16,056 catalyzed_by rows moved on 10 October).
step "$processed/graph_full/nodes.parquet" "base graph, from data/releases/v0.4" cp -r data/releases/v0.4/graph_full "$processed/graph_full"
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

# --- drug entry nodes -----------------------------------------------------------------------------------------------
step data/raw/fda_labels/.complete "FDA label protein-binding sentences" bash -c "python experiments/fetch_fda_label_protein_binding.py --evidence-dir $processed/evidence_full_v3 && touch data/raw/fda_labels/.complete"
# configs/drug_plasma_carriers.csv is committed; experiments/scope_label_plasma_binder_coverage.py --carrier-table rewrites it from the label
# sentences and the fraction-unbound database. The variant holds the drug nodes, their mechanism edges and the carriers' sequestration edges;
# both arms of the ablation are read from it (run_main_model.py --drug-entry nodes or targets).
if [ "$stage" = "graphs" ]; then
  step "$processed/graph_full_neuronal_split_binders_drugs/nodes.parquet" "drug entry nodes on the confirmatory graph" \
    python experiments/build_drug_entry_node_variant.py --graph-dir "$processed/graph_full_neuronal_split_binders" --evidence-dir "$processed/evidence_full_v3_parkinsonism"
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
