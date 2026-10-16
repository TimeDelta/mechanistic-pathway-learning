"""The technology audit: is the code in a state the research can rest on? Correctness as the tests measure it,
numerical hazards, reproducibility, device readiness for the runtime-limited GPU cluster, dependency hygiene and the
shape of the code base. Generated to docs/technology_audit.md (stamped) and data/processed/audit/technology_audit.json.

Every finding is a measured count with file and line examples; the judgement is left to the reader, which is why a
heuristic that lists lines (a bare log, a division without a floor) says so.

Usage:
  python experiments/audit_technology.py            # runs the test suite (about a minute) and the static passes
  python experiments/audit_technology.py --skip-tests
"""
from __future__ import annotations

import argparse
import ast
import importlib.metadata
import json
import platform
import re
import subprocess
import sys
import tomllib
from collections import Counter, defaultdict
from pathlib import Path

from mechanistic_pathway_learning.documentation import REPOSITORY_ROOT, provenance_stamp

PACKAGES = ("mechanistic_pathway_learning", "experiments", "scripts", "tests")
DECLARED_EXTRAS_FILE = REPOSITORY_ROOT / "pyproject.toml"
STDLIB = set(sys.stdlib_module_names)


def python_files(roots=PACKAGES, repository_root: Path = REPOSITORY_ROOT) -> list[Path]:
    """Every .py under the named roots, in path order; the two audit scripts are left out because they hold the
    patterns scanned for and would report their own regular expressions."""
    files = []
    for root in roots:
        files.extend(sorted((repository_root / root).rglob("*.py")))
    return [path for path in files if "__pycache__" not in path.parts and not path.name.startswith("audit_")]


def relative(path: Path, repository_root: Path = REPOSITORY_ROOT) -> str:
    return str(path.relative_to(repository_root))


# ----------------------------------------------------------------------------------------------------------------------
# tests


def parse_pytest_output(stdout: str, returncode: int) -> dict:
    """The last line's counts and every `-rs` skip line, grouped by file and reason."""
    summary = stdout.strip().splitlines()[-1] if stdout.strip() else ""
    counts = {key: int(value) for value, key in re.findall(r"(\d+) (passed|failed|skipped|error|errors)", summary)}
    skipped = Counter()
    for line in stdout.splitlines():
        match = re.match(r"SKIPPED \[(\d+)\] ([^:]+):\d+: (.*)", line)
        if match:
            skipped[(match.group(2), match.group(3))] += int(match.group(1))
    return {"summary": summary, "counts": counts, "returncode": returncode,
            "skipped": [{"file": file, "reason": reason, "count": count} for (file, reason), count in sorted(skipped.items())]}


def run_tests() -> dict:
    result = subprocess.run([sys.executable, "-m", "pytest", "tests/", "-q", "-p", "no:cacheprovider", "-rs"], cwd=REPOSITORY_ROOT, capture_output=True, text=True,
                            env={"PYTHONPATH": ".", "PATH": "/usr/bin:/bin:/usr/local/bin", "OMP_NUM_THREADS": "1"})
    return parse_pytest_output(result.stdout, result.returncode)


def test_inventory(repository_root: Path = REPOSITORY_ROOT, package: str = "mechanistic_pathway_learning") -> dict:
    """How many tests there are and which package modules no test names, by dotted path or by stem."""
    test_files = python_files(("tests",), repository_root)
    test_functions = 0
    for path in test_files:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        test_functions += sum(1 for node in ast.walk(tree) if isinstance(node, ast.FunctionDef) and node.name.startswith("test"))
    modules = [path for path in python_files((package,), repository_root) if path.name != "__init__.py"]
    tested = set()
    for path in test_files:
        text = path.read_text(encoding="utf-8")
        for module in modules:
            dotted = ".".join(module.relative_to(repository_root).with_suffix("").parts)
            if dotted in text or module.stem in text:
                tested.add(relative(module, repository_root))
    untested = sorted(relative(module, repository_root) for module in modules if relative(module, repository_root) not in tested)
    return {"test_files": len(test_files), "test_functions": test_functions, "package_modules": len(modules),
            "package_modules_named_in_no_test": untested}


# ----------------------------------------------------------------------------------------------------------------------
# static passes


HAZARD_PATTERNS = {
    "log_without_floor": re.compile(r"\b(?:torch|np|math)\.log\((?![^)]*(?:clamp|eps|1e-|\+ ?\d))"),
    "division_by_sum_or_count_without_floor": re.compile(r"/ ?\(?[a-z_]*(?:sum|count|total|norm|degree)[a-z_]*(?![^#\n]*(?:clamp|eps|1e-|max\(|\+ ?\d))"),
    "exp_without_clamp": re.compile(r"\btorch\.exp\((?![^)]*clamp)"),
    "sqrt_of_difference": re.compile(r"\bsqrt\([^)]*-[^)]*\)"),
}
DEVICE_PATTERNS = {
    "cuda_call": re.compile(r"\.cuda\(\)"),
    "device_cpu_literal": re.compile(r"device\s*=\s*[\"']cpu[\"']"),
    "to_device": re.compile(r"\.to\((?:device|self\.device|[a-z_]*device)\)"),
    "cuda_is_available": re.compile(r"torch\.cuda\.is_available\(\)"),
    "numpy_conversion": re.compile(r"\.(?:cpu\(\)\.)?numpy\(\)"),
    "item_call": re.compile(r"\.item\(\)"),
}
PATH_PATTERN = re.compile(r"[\"'](/home/|/root/|/tmp/|/mnt/)")
SEED_PATTERN = re.compile(r"(?:torch\.manual_seed|np\.random\.seed|random\.seed|Generator\(|default_rng)")
# iteration over a set orders rows by the hash seed; a dict keeps insertion order and is not a hazard
NONDETERMINISM_PATTERN = re.compile(r"for\s+\w+\s+in\s+(?:set\(|frozenset\()")
LONG_FUNCTION_LINES = 120


def scan(pattern: re.Pattern, roots, limit: int = 12, repository_root: Path = REPOSITORY_ROOT) -> dict:
    """Lines under the roots that match the pattern, comment lines excluded, with the first `limit` as examples."""
    hits, examples = 0, []
    for path in python_files(roots, repository_root):
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            if pattern.search(line) and not line.lstrip().startswith("#"):
                hits += 1
                if len(examples) < limit:
                    examples.append(f"{relative(path, repository_root)}:{number}: {line.strip()[:120]}")
    return {"count": hits, "examples": examples}


def numerical_hazards(repository_root: Path = REPOSITORY_ROOT) -> dict:
    return {name: scan(pattern, ("mechanistic_pathway_learning",), repository_root=repository_root) for name, pattern in HAZARD_PATTERNS.items()}


def device_readiness(repository_root: Path = REPOSITORY_ROOT) -> dict:
    import torch

    results = {name: scan(pattern, ("mechanistic_pathway_learning", "experiments"), repository_root=repository_root) for name, pattern in DEVICE_PATTERNS.items()}
    results["torch_version"] = torch.__version__
    results["cuda_available_here"] = bool(torch.cuda.is_available())
    results["cuda_build"] = torch.version.cuda
    return results


def reproducibility(repository_root: Path = REPOSITORY_ROOT) -> dict:
    return {"seed_calls": scan(SEED_PATTERN, ("mechanistic_pathway_learning", "experiments"), repository_root=repository_root),
            "iteration_over_a_set_in_a_builder": scan(NONDETERMINISM_PATTERN, ("mechanistic_pathway_learning/graph", "mechanistic_pathway_learning/evidence", "experiments"),
                                                      repository_root=repository_root),
            "absolute_paths_in_code": scan(PATH_PATTERN, ("mechanistic_pathway_learning", "experiments", "scripts"), repository_root=repository_root),
            "pinned_libraries": {name: importlib.metadata.version(name) for name in ("pandas", "pyarrow", "numpy", "torch", "scikit-learn", "cobra")},
            "pins_the_rebuild_notes_name": {"pandas": "2.3.3", "pyarrow": "25.0.1"}}


IMPORT_NAME_TO_DISTRIBUTION = {"sklearn": "scikit_learn", "yaml": "pyyaml", "esm": "fair_esm", "PIL": "pillow", "bs4": "beautifulsoup4"}


def declared_dependencies(pyproject_text: str) -> set[str]:
    """The distribution names in [project] dependencies and every optional extra, normalised to lower snake case."""
    project = tomllib.loads(pyproject_text).get("project", {})
    requirements = list(project.get("dependencies", []))
    for extra in project.get("optional-dependencies", {}).values():
        requirements.extend(extra)
    names = set()
    for requirement in requirements:
        match = re.match(r"\s*([A-Za-z0-9_.-]+)", requirement)
        if match:
            names.add(match.group(1).lower().replace("-", "_"))
    return names


def dependency_hygiene(repository_root: Path = REPOSITORY_ROOT, pyproject: Path | None = None) -> dict:
    """Third-party imports across the packages against the names pyproject.toml declares."""
    declared = declared_dependencies((pyproject or repository_root / "pyproject.toml").read_text(encoding="utf-8"))
    repository_files = python_files(("mechanistic_pathway_learning", "experiments", "scripts", "tests"), repository_root)
    # the experiment scripts import their siblings by bare name after putting their directory on sys.path
    local_module_names = {path.stem for path in repository_files} | {"mechanistic_pathway_learning", "experiments", "tests", "scripts"}
    imported: dict[str, set[str]] = defaultdict(set)
    for path in repository_files:
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [alias.name.split(".")[0] for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                names = [node.module.split(".")[0]]
            for name in names:
                if name in STDLIB or name in local_module_names:
                    continue
                imported[name].add(relative(path, repository_root))
    undeclared = {name: sorted(files)[:5] for name, files in imported.items() if IMPORT_NAME_TO_DISTRIBUTION.get(name, name).lower().replace("-", "_") not in declared}
    return {"declared": sorted(declared), "third_party_imported": sorted(imported), "imported_but_not_declared": undeclared}


def code_shape(repository_root: Path = REPOSITORY_ROOT) -> dict:
    """File sizes, function counts, functions over LONG_FUNCTION_LINES and public functions without a docstring."""
    sizes = []
    long_functions = []
    missing_docstrings = 0
    functions = 0
    for path in python_files(("mechanistic_pathway_learning", "experiments"), repository_root):
        text = path.read_text(encoding="utf-8")
        lines = text.count("\n") + 1
        sizes.append((relative(path, repository_root), lines))
        try:
            tree = ast.parse(text)
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                functions += 1
                length = (node.end_lineno or node.lineno) - node.lineno + 1
                if length > LONG_FUNCTION_LINES:
                    long_functions.append((relative(path, repository_root), node.name, length))
                if not ast.get_docstring(node) and not node.name.startswith("_") and not path.parts[-2] == "tests":
                    missing_docstrings += 1
    sizes.sort(key=lambda item: -item[1])
    long_functions.sort(key=lambda item: -item[2])
    return {"files": len(sizes), "lines": sum(count for _, count in sizes), "largest_files": sizes[:10], "functions": functions,
            "functions_over_120_lines": long_functions[:15], "public_functions_without_docstring": missing_docstrings}


def parse_ruff_output(stdout: str) -> dict:
    findings = [line for line in stdout.splitlines() if re.match(r"^[^:]+:\d+:\d+: [A-Z]+\d+", line)]
    by_rule = Counter(re.search(r": ([A-Z]+\d+)", line).group(1) for line in findings)
    return {"findings": len(findings), "by_rule": dict(by_rule.most_common()), "examples": findings[:10]}


CLASSIC_RULES = "E4,E7,E9,F"


def lint(repository_root: Path = REPOSITORY_ROOT) -> dict:
    """Two runs: the classic rule set in isolation (stable across ruff versions) and the repository's configured run,
    whose rule set follows the installed ruff's default when pyproject.toml selects none."""
    classic = subprocess.run(["ruff", "check", "--isolated", "--select", CLASSIC_RULES, "--line-length", "140", "--output-format", "concise", "."], cwd=repository_root,
                             capture_output=True, text=True)
    configured = subprocess.run(["ruff", "check", "--output-format", "concise", "."], cwd=repository_root, capture_output=True, text=True)
    version = subprocess.run(["ruff", "--version"], capture_output=True, text=True).stdout.strip()
    pyproject = (repository_root / "pyproject.toml").read_text(encoding="utf-8") if (repository_root / "pyproject.toml").exists() else ""
    selects_rules = "select" in tomllib.loads(pyproject).get("tool", {}).get("ruff", {}).get("lint", {}) if pyproject else False
    return {"ruff_version": version, "classic_rules": CLASSIC_RULES, "classic": parse_ruff_output(classic.stdout), "configured": parse_ruff_output(configured.stdout),
            "pyproject_selects_rules": selects_rules}


def environment() -> dict:
    return {"python": platform.python_version(), "platform": platform.platform(), "git_commit": subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPOSITORY_ROOT,
                                                                                                              capture_output=True, text=True).stdout.strip()}


# ----------------------------------------------------------------------------------------------------------------------
# rendering


def render(stamp: str, sections: dict) -> str:
    tests, inventory, hazards, device, repro, deps, shape, lint_result, env = (sections[name] for name in ("tests", "inventory", "hazards", "device", "reproducibility", "dependencies", "shape", "lint", "environment"))
    lines = ["# Technology audit (generated by experiments/audit_technology.py — do not edit by hand)", "", stamp, "",
             f"Python {env['python']} at commit {env['git_commit']}. Each finding is a measured count with examples; the heuristics that list lines say so and leave the judgement to the reader.", "",
             "## 1. Tests", ""]
    if tests.get("skipped_run"):
        lines.append("The suite was not run (--skip-tests).")
    else:
        lines.append(f"`pytest tests/`: **{tests['summary']}**" + (" (exit 0)" if tests["returncode"] == 0 else f" (exit {tests['returncode']})") + ".")
        if tests["skipped"]:
            lines += ["", "Skipped, with the reason each gives:", ""]
            for item in tests["skipped"]:
                lines.append(f"- {item['count']} in `{item['file']}`: {item['reason']}")
    lines += ["", f"{inventory['test_files']} test files hold {inventory['test_functions']} test functions over {inventory['package_modules']} package modules; "
              f"**{len(inventory['package_modules_named_in_no_test'])} modules are named in no test**:" + (" " + ", ".join(f"`{m}`" for m in inventory["package_modules_named_in_no_test"][:20]) if inventory["package_modules_named_in_no_test"] else " none") + ".", "",
              "## 2. Numerical hazards (heuristic line scans of the package)", ""]
    for name, result in hazards.items():
        lines.append(f"- {name.replace('_', ' ')}: **{result['count']}**" + ("; " + "; ".join(f"`{e}`" for e in result["examples"][:4]) if result["examples"] else ""))
    lines += ["", "## 3. Device readiness (the runtime-limited GPU cluster)", "",
              f"torch {device['torch_version']} (CUDA build {device['cuda_build']}; CUDA available on the machine that ran this audit: {device['cuda_available_here']}).", ""]
    for name in ("cuda_call", "device_cpu_literal", "to_device", "cuda_is_available", "numpy_conversion", "item_call"):
        result = device[name]
        lines.append(f"- {name.replace('_', ' ')}: **{result['count']}**" + ("; " + "; ".join(f"`{e}`" for e in result["examples"][:3]) if result["examples"] else ""))
    lines += ["", "A `.cuda()` call or a literal cpu device pins the code to one device; `.item()` and `.numpy()` inside a training step force a synchronisation. "
              "The counts above locate them; whether each is in a hot loop is for the reader.", "",
              "## 4. Reproducibility", "",
              f"- seed calls: **{repro['seed_calls']['count']}**" + ("; " + "; ".join(f"`{e}`" for e in repro["seed_calls"]["examples"][:4]) if repro["seed_calls"]["examples"] else ""),
              f"- iteration over a set in a builder (row order then depends on the hash seed; docs/data_layer_rebuild.md): **{repro['iteration_over_a_set_in_a_builder']['count']}**"
              + ("; " + "; ".join(f"`{e}`" for e in repro["iteration_over_a_set_in_a_builder"]["examples"][:6]) if repro["iteration_over_a_set_in_a_builder"]["examples"] else ""),
              f"- absolute paths in code: **{repro['absolute_paths_in_code']['count']}**" + ("; " + "; ".join(f"`{e}`" for e in repro["absolute_paths_in_code"]["examples"][:4]) if repro["absolute_paths_in_code"]["examples"] else ""),
              f"- installed: {', '.join(f'{k} {v}' for k, v in repro['pinned_libraries'].items())}; the rebuild notes pin pandas {repro['pins_the_rebuild_notes_name']['pandas']} and pyarrow {repro['pins_the_rebuild_notes_name']['pyarrow']} for hash-comparable tables"
              + (" (both match)" if repro["pinned_libraries"]["pandas"] == repro["pins_the_rebuild_notes_name"]["pandas"] and repro["pinned_libraries"]["pyarrow"] == repro["pins_the_rebuild_notes_name"]["pyarrow"] else " (**mismatch**)") + ".",
              "", "## 5. Dependencies", "",
              f"Declared in pyproject.toml: {', '.join(deps['declared'])}. Third-party modules imported: {', '.join(deps['third_party_imported'])}.",
              f"**Imported but declared in no extra: {', '.join(f'`{name}` ({files[0]})' for name, files in deps['imported_but_not_declared'].items()) or 'none'}.**", "",
              "## 6. Lint", "",
              f"The classic rule set ({lint_result['classic_rules']}) in isolation: **{lint_result['classic']['findings']}** findings"
              + (": " + ", ".join(f"{rule} {count}" for rule, count in lint_result["classic"]["by_rule"].items()) if lint_result["classic"]["by_rule"] else "") + "."
              + (" First: " + "; ".join(f"`{e}`" for e in lint_result["classic"]["examples"][:5]) if lint_result["classic"]["examples"] else ""), "",
              f"The repository's configured run under {lint_result['ruff_version']}: **{lint_result['configured']['findings']}** findings"
              + (": " + ", ".join(f"{rule} {count}" for rule, count in list(lint_result["configured"]["by_rule"].items())[:12]) if lint_result["configured"]["by_rule"] else "") + "."
              + ("" if lint_result["pyproject_selects_rules"] else " pyproject.toml selects no rules, so this set is the installed ruff's default and changes with its version;"
                 " the classic count above is the stable measure."), "",
              "## 7. Shape of the code", "",
              f"{shape['files']} files, {shape['lines']:,} lines, {shape['functions']:,} functions; **{len(shape['functions_over_120_lines'])} functions over 120 lines**; "
              f"{shape['public_functions_without_docstring']} public functions without a docstring.", "",
              "| file | lines |", "| --- | --- |"]
    for path, count in shape["largest_files"]:
        lines.append(f"| `{path}` | {count:,} |")
    if shape["functions_over_120_lines"]:
        lines += ["", "| function | lines |", "| --- | --- |"]
        for path, name, length in shape["functions_over_120_lines"]:
            lines.append(f"| `{path}::{name}` | {length} |")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/audit"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/technology_audit.md"))
    parser.add_argument("--skip-tests", action="store_true")
    arguments = parser.parse_args()
    sections = {"tests": {"skipped_run": True} if arguments.skip_tests else run_tests(), "inventory": test_inventory(), "hazards": numerical_hazards(),
                "device": device_readiness(), "reproducibility": reproducibility(), "dependencies": dependency_hygiene(), "shape": code_shape(), "lint": lint(),
                "environment": environment()}
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    (arguments.output_dir / "technology_audit.json").write_text(json.dumps(sections, indent=1, default=str) + "\n", encoding="utf-8")
    arguments.markdown_output.write_text(render(provenance_stamp("experiments/audit_technology.py", ["the working tree", "pytest tests/", "ruff"]), sections), encoding="utf-8")
    print(f"wrote {arguments.markdown_output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
