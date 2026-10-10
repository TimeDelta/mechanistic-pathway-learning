"""Run a repository entry point with the nervous-system ATC rule lifted, without editing the repository: the rule's
function is replaced before the entry point's modules import it. For counting what the rule excludes.

Part of the label source review of 10 October 2026 (docs/label_source_review.md). Counting only: no label, table or
configuration of the study is changed.

Usage:
  python experiments/label_review/run_with_any_atc.py any_atc|n_only <entry point> <its arguments>
"""
import runpy
import sys

sys.path.insert(0, ".")
lift_rule = sys.argv[1] == "any_atc"
target = sys.argv[2]
sys.argv = [target, *sys.argv[3:]]
if lift_rule:
    import mechanistic_pathway_learning.evidence.load_drug_label_events as events
    events.is_nervous_system_atc = lambda atc_codes: True
if target.endswith(".py"):
    runpy.run_path(target, run_name="__main__")
else:
    runpy.run_module(target, run_name="__main__")
