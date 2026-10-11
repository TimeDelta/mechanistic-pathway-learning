#!/bin/bash
# Counting builds of the label changes the user approved on 10 and 11 October 2026 (docs/label_source_review.md,
# sections 9 and 12). Written under data/processed/label_review, which is not in git; none is a study table.
#
#   approved                    gain-of-function sign, placebo mask, at most three mechanism targets per drug
#   approved_no_cap             the same with no limit on the number of targets (the user's statement of 11 October)
#   approved_with_list          approved, with the draft list of admitted ATC groups in place of the nervous-system rule
#   approved_with_list_no_cap   approved_no_cap, with that list
#
# About one minute per build on one core.
#   bash experiments/label_review/build_approved_variants.sh
set -uo pipefail
cd "$(cd "$(dirname "$0")/../.." && pwd)"
export PYTHONPATH=. OMP_NUM_THREADS=1
out=data/processed/label_review
admitted_list=configs/admitted_atc_prefixes_draft.csv
mkdir -p $out/logs
for variant in "approved 3 no_list" "approved_no_cap 0 no_list" "approved_with_list 3 list" "approved_with_list_no_cap 0 list"; do
  set -- $variant; name=$1; cap=$2; list_choice=$3
  list_arguments=""
  if [ "$list_choice" = "list" ]; then list_arguments="--admitted-atc-prefixes $admitted_list"; fi
  echo "== $name"
  python experiments/build_onsides_reports.py --crosswalk docs/symptom_crosswalk.csv --output-dir $out/onsides_$name --markdown-output $out/logs/onsides_$name.md \
      --evidence-dir $out/evidence_$name --max-drug-targets $cap $list_arguments > $out/logs/onsides_$name.log 2>&1 || { echo "onsides failed"; tail -3 $out/logs/onsides_$name.log; continue; }
  python -m mechanistic_pathway_learning.evidence.assemble_evidence_table --crosswalk docs/symptom_crosswalk.csv --graph-dir data/processed/graph_full \
      --output-dir $out/evidence_$name --extra-reports $out/onsides_$name/onsides_reports.parquet --onsides-bridge data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json \
      --max-drug-targets $cap --seed-gain-of-function-genes-positive $list_arguments > $out/logs/evidence_$name.log 2>&1 || { echo "assembly failed"; tail -3 $out/logs/evidence_$name.log; continue; }
  if [ -f $out/selection_$name.parquet ]; then
    # a selection is written once (build_label_selection.py refuses to overwrite one); delete it by hand to rebuild
    echo "selection_$name.parquet exists and is kept"
  else
    python experiments/build_label_selection.py --evidence-dir $out/evidence_$name --mask-grades C --mask-never-above-placebo --output $out/selection_$name.parquet \
        > $out/logs/selection_$name.log 2>&1 || { echo "selection failed"; tail -3 $out/logs/selection_$name.log; continue; }
  fi
  echo finished > $out/logs/$name.DONE
done
echo finished
