#!/bin/bash
# Evidence tables under other drug rules, for counting only (docs/label_source_review.md). Written under
# data/processed/label_review, which is not in git. Six variants: the nervous-system ATC rule kept or lifted, and at most
# one, at most three or any number of mechanism targets per drug. About 25 minutes on two cores.
#   bash experiments/label_review/build_label_variants.sh
set -uo pipefail
cd "$(cd "$(dirname "$0")/../.." && pwd)"
export PYTHONPATH=. OMP_NUM_THREADS=1
runner=experiments/label_review/run_with_any_atc.py
out=data/processed/label_review
mkdir -p $out/logs
for variant in "n_only 1" "any_atc 1" "n_only 3" "any_atc 3" "n_only 0" "any_atc 0"; do
  set -- $variant; atc=$1; cap=$2; name="${atc}_cap${cap}"
  echo "== $name"
  python $runner $atc experiments/build_onsides_reports.py --crosswalk docs/symptom_crosswalk.csv --output-dir $out/onsides_$name --markdown-output $out/logs/onsides_$name.md \
      --evidence-dir $out/evidence_$name --max-drug-targets $cap > $out/logs/onsides_$name.log 2>&1 || { echo "onsides failed"; tail -3 $out/logs/onsides_$name.log; continue; }
  python $runner $atc mechanistic_pathway_learning.evidence.assemble_evidence_table --crosswalk docs/symptom_crosswalk.csv --graph-dir data/processed/graph_full \
      --output-dir $out/evidence_$name --extra-reports $out/onsides_$name/onsides_reports.parquet --onsides-bridge data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json \
      --max-drug-targets $cap > $out/logs/evidence_$name.log 2>&1 || { echo "assembly failed"; tail -3 $out/logs/evidence_$name.log; continue; }
  python experiments/build_label_selection.py --evidence-dir $out/evidence_$name --mask-grades C --output $out/selection_$name.parquet > $out/logs/selection_$name.log 2>&1 || { echo "selection failed"; tail -3 $out/logs/selection_$name.log; continue; }
  tail -1 $out/logs/selection_$name.log | cut -c1-200
done
echo finished
