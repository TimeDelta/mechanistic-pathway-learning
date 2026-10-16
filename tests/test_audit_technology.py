"""The technology audit's scans and parsers on fixtures (experiments/audit_technology.py)."""
from __future__ import annotations

from pathlib import Path

from experiments.audit_technology import (
    DEVICE_PATTERNS,
    HAZARD_PATTERNS,
    LONG_FUNCTION_LINES,
    NONDETERMINISM_PATTERN,
    PATH_PATTERN,
    SEED_PATTERN,
    code_shape,
    declared_dependencies,
    dependency_hygiene,
    parse_pytest_output,
    parse_ruff_output,
    python_files,
    render,
    scan,
    test_inventory as inventory_of_tests,
)


def write(root: Path, relative_path: str, text: str) -> Path:
    path = root / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def matches(name: str, line: str, patterns=HAZARD_PATTERNS) -> bool:
    return patterns[name].search(line) is not None


def test_a_log_is_a_hazard_only_without_a_floor():
    assert matches("log_without_floor", "value = torch.log(probability)")
    assert matches("log_without_floor", "value = np.log(counts)")
    assert not matches("log_without_floor", "value = torch.log(probability.clamp(min=1e-6))")
    assert not matches("log_without_floor", "value = torch.log(probability + eps)")
    assert not matches("log_without_floor", "value = torch.log(counts + 1)")
    assert not matches("log_without_floor", "value = np.log1p(counts)"), "log1p is safe at zero and is not the pattern"


def test_a_division_by_a_count_is_a_hazard_only_without_a_floor():
    assert matches("division_by_sum_or_count_without_floor", "share = weight / total_in_degree[target]")
    assert matches("division_by_sum_or_count_without_floor", "mean = value / report_count")
    assert not matches("division_by_sum_or_count_without_floor", "share = weight / total_in_degree.clamp(min=1)")
    assert not matches("division_by_sum_or_count_without_floor", "mean = value / max(report_count, 1)")
    assert not matches("division_by_sum_or_count_without_floor", "mean = value / (report_count + 1e-9)")
    assert not matches("division_by_sum_or_count_without_floor", "half = value / 2")


def test_an_exponential_is_a_hazard_only_without_a_clamp():
    assert matches("exp_without_clamp", "gain = torch.exp(logit)")
    assert not matches("exp_without_clamp", "gain = torch.exp(logit.clamp(max=20.0))")
    assert matches("sqrt_of_difference", "spread = torch.sqrt(second_moment - mean ** 2)")
    assert not matches("sqrt_of_difference", "spread = torch.sqrt(variance)")


def test_device_patterns_find_pinned_devices_and_synchronisation_points():
    assert matches("cuda_call", "tensor = tensor.cuda()", DEVICE_PATTERNS)
    assert matches("device_cpu_literal", "model = Model(device='cpu')", DEVICE_PATTERNS)
    assert matches("to_device", "batch = batch.to(self.device)", DEVICE_PATTERNS)
    assert matches("to_device", "batch = batch.to(target_device)", DEVICE_PATTERNS)
    assert not matches("to_device", "batch = batch.to(torch.float32)", DEVICE_PATTERNS)
    assert matches("numpy_conversion", "scores = logits.cpu().numpy()", DEVICE_PATTERNS)
    assert matches("item_call", "running_loss += loss.item()", DEVICE_PATTERNS)
    assert matches("cuda_is_available", "if torch.cuda.is_available():", DEVICE_PATTERNS)


def test_reproducibility_patterns_find_seeds_set_iteration_and_absolute_paths():
    assert SEED_PATTERN.search("torch.manual_seed(seed)") and SEED_PATTERN.search("generator = np.random.default_rng(20261009)")
    assert not SEED_PATTERN.search("seed_rows = table[table.is_seed]")
    assert NONDETERMINISM_PATTERN.search("for symbol in set(symbols):")
    assert NONDETERMINISM_PATTERN.search("for symbol in frozenset(symbols):")
    assert not NONDETERMINISM_PATTERN.search("for symbol in sorted(set(symbols)):")
    assert not NONDETERMINISM_PATTERN.search("for key in mapping.keys():"), "a dict keeps insertion order"
    assert PATH_PATTERN.search('root = Path("/home/claude/repo")') and PATH_PATTERN.search("cache = '/tmp/cache'")
    assert not PATH_PATTERN.search('root = Path("data/processed")')


def test_scan_counts_matching_lines_skips_comments_and_caps_the_examples(tmp_path):
    write(tmp_path, "package/module.py", "import torch\n# torch.log(commented)\nfirst = torch.log(x)\nsecond = torch.log(y)\nsafe = torch.log(z.clamp(min=1e-6))\n")
    write(tmp_path, "package/audit_me.py", "value = torch.log(x)\n")
    write(tmp_path, "package/__pycache__/module.cpython-312.py", "value = torch.log(x)\n")
    result = scan(HAZARD_PATTERNS["log_without_floor"], ("package",), limit=1, repository_root=tmp_path)
    assert result["count"] == 2, "the comment, the audit script and the cache are left out"
    assert result["examples"] == ["package/module.py:3: first = torch.log(x)"]
    assert [path.name for path in python_files(("package",), tmp_path)] == ["module.py"]


def test_the_inventory_names_the_package_modules_no_test_mentions(tmp_path):
    write(tmp_path, "mechanistic_pathway_learning/__init__.py", "")
    write(tmp_path, "mechanistic_pathway_learning/graph/build.py", "def build():\n    return 1\n")
    write(tmp_path, "mechanistic_pathway_learning/graph/orphan.py", "def orphan():\n    return 2\n")
    write(tmp_path, "mechanistic_pathway_learning/evidence/assemble.py", "def assemble():\n    return 3\n")
    write(tmp_path, "tests/test_build.py", "from mechanistic_pathway_learning.graph.build import build\n\ndef test_one():\n    assert build() == 1\n\ndef test_two():\n    assert True\n")
    write(tmp_path, "tests/test_assemble.py", "import mechanistic_pathway_learning.evidence.assemble as assemble\n\ndef test_three():\n    assert assemble.assemble() == 3\n")
    result = inventory_of_tests(tmp_path)
    assert result == {"test_files": 2, "test_functions": 3, "package_modules": 3, "package_modules_named_in_no_test": ["mechanistic_pathway_learning/graph/orphan.py"]}


def test_pytest_output_parsing_reads_the_summary_and_groups_the_skips():
    stdout = "\n".join([
        "........s.....s",
        "=========================== short test summary info ============================",
        "SKIPPED [4] tests/test_label_review_changes.py:14: data/processed/label_review/drug_names.csv is not built",
        "SKIPPED [1] tests/test_gtopdb.py:9: GtoPdb export not fetched",
        "SKIPPED [2] tests/test_label_review_changes.py:40: data/processed/label_review/drug_names.csv is not built",
        "1012 passed, 7 skipped in 61.52s",
    ])
    result = parse_pytest_output(stdout, 0)
    assert result["counts"] == {"passed": 1012, "skipped": 7} and result["returncode"] == 0 and result["summary"].startswith("1012 passed")
    assert result["skipped"] == [
        {"file": "tests/test_gtopdb.py", "reason": "GtoPdb export not fetched", "count": 1},
        {"file": "tests/test_label_review_changes.py", "reason": "data/processed/label_review/drug_names.csv is not built", "count": 6},
    ]
    assert parse_pytest_output("", 1) == {"summary": "", "counts": {}, "returncode": 1, "skipped": []}


def test_ruff_output_parsing_counts_findings_by_rule():
    stdout = "experiments/a.py:12:5: F841 Local variable `x` is assigned to but never used\nexperiments/b.py:3:1: F401 [*] `os` imported but unused\nexperiments/c.py:9:1: F401 [*] `re` imported but unused\nFound 3 errors.\n"
    result = parse_ruff_output(stdout)
    assert result["findings"] == 3 and result["by_rule"] == {"F401": 2, "F841": 1} and result["examples"][0].startswith("experiments/a.py:12:5: F841")
    assert parse_ruff_output("All checks passed!\n") == {"findings": 0, "by_rule": {}, "examples": []}


PYPROJECT = """[build-system]
requires = ["setuptools>=61"]

[project]
name = "demo"
dependencies = ["numpy>=1.26", "pandas==2.3.3", "scikit-learn", "PyYAML>=6"]

[project.optional-dependencies]
dev = ["pytest", "ruff"]

[tool.ruff]
line-length = 140
"""


def test_declared_dependencies_normalise_names_and_stop_at_the_tool_tables():
    declared = declared_dependencies(PYPROJECT)
    assert declared == {"numpy", "pandas", "scikit_learn", "pyyaml", "pytest", "ruff"}, "the build-system requirement and the tool tables are not dependencies"


def test_dependency_hygiene_names_imports_no_extra_declares(tmp_path):
    write(tmp_path, "pyproject.toml", PYPROJECT)
    write(tmp_path, "mechanistic_pathway_learning/model.py", "import numpy as np\nimport sklearn.linear_model\nimport yaml\nimport torch\nfrom pathlib import Path\nfrom . import sibling\n")
    write(tmp_path, "experiments/run.py", "import pandas as pd\nimport cobra\nimport helper\nfrom mechanistic_pathway_learning.model import fit\n")
    write(tmp_path, "experiments/helper.py", "def help():\n    return 1\n")
    result = dependency_hygiene(tmp_path)
    assert result["imported_but_not_declared"] == {"torch": ["mechanistic_pathway_learning/model.py"], "cobra": ["experiments/run.py"]}
    assert "helper" not in result["third_party_imported"], "a sibling script imported by bare name is local"
    assert "pathlib" not in result["third_party_imported"] and "mechanistic_pathway_learning" not in result["third_party_imported"]
    assert "sklearn" in result["third_party_imported"] and "yaml" in result["third_party_imported"], "declared under their distribution names"


def test_code_shape_finds_long_functions_and_missing_docstrings(tmp_path):
    body = "\n".join("    value += 1" for _ in range(LONG_FUNCTION_LINES))
    write(tmp_path, "mechanistic_pathway_learning/long.py", f"def documented():\n    \"\"\"Says what it does.\"\"\"\n    value = 0\n{body}\n    return value\n\n\ndef _private():\n    return 1\n\n\ndef bare():\n    return 2\n")
    write(tmp_path, "experiments/short.py", "def short():\n    \"\"\"Short.\"\"\"\n    return 3\n")
    result = code_shape(tmp_path)
    assert result["files"] == 2 and result["functions"] == 4
    assert result["functions_over_120_lines"] == [("mechanistic_pathway_learning/long.py", "documented", LONG_FUNCTION_LINES + 4)], "def, docstring, assignment, the body and the return"
    assert result["public_functions_without_docstring"] == 1, "bare() only; the private function is not counted"
    assert result["largest_files"][0][0] == "mechanistic_pathway_learning/long.py"


def test_render_writes_every_section_and_says_when_the_suite_was_not_run():
    empty = {"count": 0, "examples": []}
    sections = {
        "tests": {"skipped_run": True},
        "inventory": {"test_files": 2, "test_functions": 3, "package_modules": 3, "package_modules_named_in_no_test": ["mechanistic_pathway_learning/graph/orphan.py"]},
        "hazards": {name: empty for name in HAZARD_PATTERNS},
        "device": {**{name: empty for name in DEVICE_PATTERNS}, "torch_version": "2.5.1", "cuda_available_here": False, "cuda_build": None},
        "reproducibility": {"seed_calls": {"count": 1, "examples": ["experiments/run.py:3: torch.manual_seed(0)"]}, "iteration_over_a_set_in_a_builder": empty,
                            "absolute_paths_in_code": empty, "pinned_libraries": {"pandas": "2.3.3", "pyarrow": "25.0.1", "numpy": "2.0.0", "torch": "2.5.1", "scikit-learn": "1.5", "cobra": "0.29"},
                            "pins_the_rebuild_notes_name": {"pandas": "2.3.3", "pyarrow": "25.0.1"}},
        "dependencies": {"declared": ["numpy"], "third_party_imported": ["numpy", "torch"], "imported_but_not_declared": {"torch": ["mechanistic_pathway_learning/model.py"]}},
        "shape": {"files": 2, "lines": 140, "functions": 4, "largest_files": [("mechanistic_pathway_learning/long.py", 130)], "functions_over_120_lines": [("mechanistic_pathway_learning/long.py", "documented", 123)],
                  "public_functions_without_docstring": 1},
        "lint": {"ruff_version": "ruff 0.16.8", "classic_rules": "E4,E7,E9,F", "classic": {"findings": 0, "by_rule": {}, "examples": []},
                 "configured": {"findings": 2, "by_rule": {"B023": 2}, "examples": []}, "pyproject_selects_rules": False},
        "environment": {"python": "3.12.3", "platform": "Linux", "git_commit": "abc1234"},
    }
    text = render("Generated by `experiments/audit_technology.py` on 2026-10-16 at commit abc1234 from the working tree; do not edit by hand.", sections)
    assert text.startswith("# Technology audit (generated by experiments/audit_technology.py — do not edit by hand)")
    assert "The suite was not run (--skip-tests)." in text
    assert "**1 modules are named in no test**: `mechanistic_pathway_learning/graph/orphan.py`" in text
    assert "(both match)" in text and "`torch` (mechanistic_pathway_learning/model.py)" in text
    assert "| `mechanistic_pathway_learning/long.py::documented` | 123 |" in text
    assert "in isolation: **0** findings." in text and "under ruff 0.16.8: **2** findings: B023 2. pyproject.toml selects no rules" in text
    for heading in ("## 1. Tests", "## 2. Numerical hazards", "## 3. Device readiness", "## 4. Reproducibility", "## 5. Dependencies", "## 6. Lint", "## 7. Shape of the code"):
        assert heading in text


def test_render_reports_the_suite_and_a_pin_mismatch():
    empty = {"count": 0, "examples": []}
    sections = {
        "tests": {"summary": "10 passed, 1 skipped in 2.0s", "counts": {"passed": 10, "skipped": 1}, "returncode": 0,
                  "skipped": [{"file": "tests/test_x.py", "reason": "not built", "count": 1}]},
        "inventory": {"test_files": 1, "test_functions": 11, "package_modules": 1, "package_modules_named_in_no_test": []},
        "hazards": {name: empty for name in HAZARD_PATTERNS},
        "device": {**{name: empty for name in DEVICE_PATTERNS}, "torch_version": "2.5.1", "cuda_available_here": True, "cuda_build": "12.4"},
        "reproducibility": {"seed_calls": empty, "iteration_over_a_set_in_a_builder": empty, "absolute_paths_in_code": empty,
                            "pinned_libraries": {"pandas": "2.2.0", "pyarrow": "25.0.1", "numpy": "2.0.0", "torch": "2.5.1", "scikit-learn": "1.5", "cobra": "0.29"},
                            "pins_the_rebuild_notes_name": {"pandas": "2.3.3", "pyarrow": "25.0.1"}},
        "dependencies": {"declared": ["numpy"], "third_party_imported": ["numpy"], "imported_but_not_declared": {}},
        "shape": {"files": 1, "lines": 10, "functions": 1, "largest_files": [("mechanistic_pathway_learning/a.py", 10)], "functions_over_120_lines": [], "public_functions_without_docstring": 0},
        "lint": {"ruff_version": "ruff 0.16.8", "classic_rules": "E4,E7,E9,F", "classic": {"findings": 2, "by_rule": {"F401": 2}, "examples": ["experiments/b.py:3:1: F401 `os` imported but unused"]},
                 "configured": {"findings": 2, "by_rule": {"F401": 2}, "examples": []}, "pyproject_selects_rules": True},
        "environment": {"python": "3.12.3", "platform": "Linux", "git_commit": "abc1234"},
    }
    text = render("stamp", sections)
    assert "`pytest tests/`: **10 passed, 1 skipped in 2.0s** (exit 0)." in text and "- 1 in `tests/test_x.py`: not built" in text
    assert "(**mismatch**)" in text and "**0 modules are named in no test**: none." in text
    assert "in isolation: **2** findings: F401 2. First: `experiments/b.py:3:1: F401 `os` imported but unused`" in text and "selects no rules" not in text
    assert "Imported but declared in no extra: none." in text
