"""The training script records the commit it runs from, so a result names the code that produced it."""
import importlib.util
import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def load_training_script():
    specification = importlib.util.spec_from_file_location("run_main_model", REPOSITORY_ROOT / "experiments" / "run_main_model.py")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_git_provenance_names_the_head_commit_and_lists_tracked_changes() -> None:
    provenance = load_training_script().git_provenance()
    assert re.fullmatch(r"[0-9a-f]{40}", provenance["commit"])
    assert isinstance(provenance["tracked_changes"], list)
