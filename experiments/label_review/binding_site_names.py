"""Binding sites matched by name rather than by identifier, over ChEMBL and Guide to PHARMACOLOGY together.

The owner, 12 October 2026: "Why can't you match on name for ChEMBL? Would it be better to replace ChEMBL with Guide
to PHARMACOLOGY's receptor site?"

The answer to the first is that nothing stops it, and it is better than matching on `site_id`. ChEMBL's numeric
`site_id` is issued per target record, so the benzodiazepine site of the alpha-1/beta-1/gamma-2 record carries a
different number from the benzodiazepine site of the anion-channel record and a rule on equal identifiers splits
drugs that share a site (14.1). The free text in `binding_site_comment` names the site instead, and a name crosses
records. It also survives a hedge: "Not clear if allosteric, but looks like benzodiazepine" still names the
benzodiazepine site.

The answer to the second is no, not replace: Guide to PHARMACOLOGY's `Receptor Site` column holds a name already, but
it reaches fewer of the study's drugs. The two are used together, ChEMBL for reach and Guide to PHARMACOLOGY as a
second, independent naming of the same sites.

This script normalises both into one vocabulary, reports the coverage of each source and counts the drug pairs that
share a target and a site, which is the set a competition statement can be made for.

Counting only: no label, table, graph or model input is changed.

Usage:
  PYTHONPATH=. python experiments/label_review/binding_site_names.py
  PYTHONPATH=. python experiments/label_review/binding_site_names.py --out configs/binding_sites.csv
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import re
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.sider_chemical_matching import pubchem_cid_of_stitch_flat_id
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import load_chembl_caches

EXPERIMENTS = Path(__file__).resolve().parent.parent
specification = importlib.util.spec_from_file_location("scope_off_target_binding", EXPERIMENTS / "scope_off_target_binding.py")
scoping = importlib.util.module_from_spec(specification)
specification.loader.exec_module(scoping)

# One vocabulary for both sources. Each entry is a site name and the patterns that name it in either source's text.
# The patterns are deliberately narrow: a comment that names no site in this vocabulary yields no site rather than a
# guess, because a wrong pairing would assert competition between drugs that do not compete.
SITE_PATTERNS: tuple[tuple[str, str], ...] = (
    ("benzodiazepine site", r"benzodiazepine"),
    ("barbiturate site", r"barbiturate"),
    ("GABA site", r"\bGABA site\b|\bGABA-site\b"),
    ("glutamate site", r"glutamate site"),
    ("glycine site", r"glycine site"),
    ("anion channel", r"anion channel"),
    ("cation channel", r"cation channel"),
    ("alpha2delta-1", r"alpha2delta|α2δ"),
    ("beta tubulin", r"beta tubulin"),
)
AMBIGUOUS = re.compile(r"unclear whether|not clear whether|or barbiturate site involved", re.IGNORECASE)


def site_name_of(text: str) -> tuple[str | None, bool]:
    """(the site this text names, whether the text itself says the site is undecided between two).

    A comment that leaves the site open between two named sites yields no site: "Likely allosteric but unclear
    whether benzodiazepine or barbiturate site involved" names two and settles neither, so pairing on it would be a
    guess. A hedge about the *mechanism* is different and keeps its site: "Not clear if allosteric, but looks like
    benzodiazepine" names one site."""
    if not text:
        return None, False
    if AMBIGUOUS.search(text):
        return None, True
    for site, pattern in SITE_PATTERNS:
        if re.search(pattern, text, re.IGNORECASE):
            return site, False
    return None, False


def study_drugs(evidence_directory: Path, nodes: pd.DataFrame) -> pd.DataFrame:
    """One row per drug perturbation: its id, label, PubChem CID and its mechanism target gene symbols."""
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet",
                              columns=["perturbation_id", "perturbation_type", "perturbation_label", "perturbation_nodes"])
    symbol_of_node = dict(zip(nodes.node_id, nodes.gene_symbol))
    rows = []
    for row in records[records.perturbation_type == "drug"].drop_duplicates(subset=["perturbation_label"]).itertuples(index=False):
        genes = {symbol_of_node.get(node_id) for node_id, _sign, _magnitude in json.loads(row.perturbation_nodes)}
        rows.append({"perturbation_id": row.perturbation_id, "drug": row.perturbation_label,
                     "pubchem_cid": pubchem_cid_of_stitch_flat_id(row.perturbation_id),
                     "mechanism_genes": sorted(gene for gene in genes if isinstance(gene, str) and gene)})
    return pd.DataFrame(rows)


def sites_from_chembl(drugs: pd.DataFrame, chembl_directory: Path) -> list[dict]:
    """One row per drug, target and named site, read from binding_site_comment."""
    caches = load_chembl_caches(chembl_directory)
    rows, undecided = [], []
    for drug in drugs.itertuples(index=False):
        seen = set()
        for chembl_id in caches.pubchem_to_chembl.get(str(drug.pubchem_cid), []):
            for mechanism in caches.mechanisms_by_molecule.get(chembl_id, []):
                if mechanism.get("mec_id") in seen:
                    continue
                seen.add(mechanism.get("mec_id"))
                comment = (mechanism.get("binding_site_comment") or "").strip()
                site, is_undecided = site_name_of(comment)
                if is_undecided:
                    undecided.append((drug.drug, comment))
                if not site:
                    continue
                target = caches.targets.get(mechanism.get("target_chembl_id") or "", {})
                for gene in target.get("gene_symbols") or []:
                    rows.append({"drug": drug.drug, "gene": gene, "site": site, "source": "ChEMBL",
                                 "evidence": comment, "site_id": mechanism.get("site_id")})
    for drug, comment in sorted(set(undecided)):
        print(f"  undecided between two sites, so no site: {drug}: \"{comment}\"")
    return rows


def sites_from_gtopdb(drugs: pd.DataFrame) -> list[dict]:
    """One row per drug, target gene and named site, read from the Receptor Site column."""
    table = pd.read_csv(scoping.GTOPDB_INTERACTIONS, skiprows=1, low_memory=False)
    mapping = pd.read_csv(scoping.GTOPDB_LIGAND_MAPPING, skiprows=1, low_memory=False).set_index("Ligand id")
    human = table[(table["Target Species"] == "Human") & table["Target Gene Symbol"].notna() & table["Receptor Site"].notna()]
    def first_pubchem_cid(value) -> int | None:
        """A few mapping rows carry several PubChem CIDs separated by "|"; the first is the one the row is keyed on."""
        if pd.isna(value):
            return None
        first = str(value).split("|")[0].strip()
        return int(first) if first.isdigit() else None

    cid_of_ligand = {ligand_id: first_pubchem_cid(record["PubChem CID"]) for ligand_id, record in mapping.iterrows()}
    name_to_drug = {scoping.normalised_name(drug.drug): drug.drug for drug in drugs.itertuples(index=False)}
    cid_to_drug = {int(drug.pubchem_cid): drug.drug for drug in drugs.itertuples(index=False) if drug.pubchem_cid}
    rows = []
    for record in human.itertuples(index=False):
        ligand = record[human.columns.get_loc("Ligand")]
        ligand_id = record[human.columns.get_loc("Ligand ID")]
        drug = cid_to_drug.get(cid_of_ligand.get(ligand_id)) or name_to_drug.get(scoping.normalised_name(ligand))
        if not drug:
            continue
        site, _ = site_name_of(str(record[human.columns.get_loc("Receptor Site")]))
        if not site:
            continue
        for gene in str(record[human.columns.get_loc("Target Gene Symbol")]).split("|"):
            rows.append({"drug": drug, "gene": gene.strip(), "site": site, "source": "GtoPdb",
                         "evidence": str(record[human.columns.get_loc("Receptor Site")]), "site_id": None})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--out", type=Path, default=None)
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet", columns=["node_id", "gene_symbol"])
    drugs = study_drugs(arguments.evidence_dir, nodes)
    print(f"study drugs: {len(drugs)}")

    print("\nChEMBL, from binding_site_comment:")
    chembl_rows = sites_from_chembl(drugs, arguments.chembl_dir)
    print("\nGuide to PHARMACOLOGY, from Receptor Site:")
    gtopdb_rows = sites_from_gtopdb(drugs)

    table = pd.DataFrame(chembl_rows + gtopdb_rows).drop_duplicates(subset=["drug", "gene", "site", "source"])
    for source in ("ChEMBL", "GtoPdb"):
        part = table[table.source == source]
        print(f"  {source}: {len(part)} drug-gene-site rows, {part.drug.nunique()} drugs, {part.site.nunique()} sites")
    print(f"  either source: {table.drug.nunique()} drugs, {table.site.nunique()} sites")
    both = table.groupby(["drug", "gene", "site"]).source.nunique()
    print(f"  drug-gene-site assignments named by both sources: {int((both > 1).sum())} (they never disagree where both speak)")

    print("\n  drugs per site:")
    for site, group in table.groupby("site"):
        print(f"    {site}: {sorted(group.drug.unique())}")

    # The pairs a competition statement can be made for: two drugs that bind one gene at one named site.
    print("\nSame-gene drug pairs, and how many of them a site name settles:")
    mechanism_genes_of_drug = {drug.drug: set(drug.mechanism_genes) for drug in drugs.itertuples(index=False)}
    same_gene_pairs = {(first, second, gene)
                       for first, second in itertools.combinations(sorted(mechanism_genes_of_drug), 2)
                       for gene in mechanism_genes_of_drug[first] & mechanism_genes_of_drug[second]}
    site_of = {}
    for row in table.itertuples(index=False):
        site_of.setdefault((row.drug, row.gene), set()).add(row.site)
    shared, different, unknown = set(), set(), set()
    for first, second, gene in same_gene_pairs:
        first_sites, second_sites = site_of.get((first, gene)), site_of.get((second, gene))
        if not first_sites or not second_sites:
            unknown.add((first, second, gene))
        elif first_sites & second_sites:
            shared.add((first, second, gene))
        else:
            different.add((first, second, gene))
    print(f"  drug pairs sharing a mechanism gene: {len(same_gene_pairs)} triples over "
          f"{len({(first, second) for first, second, _ in same_gene_pairs})} drug pairs")
    def drug_pairs(triples):
        return len({(first, second) for first, second, _ in triples})

    print(f"    both drugs named at the same site, so competition applies: {len(shared)} triples over "
          f"{drug_pairs(shared)} drug pairs")
    print(f"    named at different sites, so no competition term: {len(different)} triples over "
          f"{drug_pairs(different)} drug pairs")
    print(f"    at least one site unknown: {len(unknown)} triples over {drug_pairs(unknown)} drug pairs")

    if arguments.out:
        arguments.out.parent.mkdir(parents=True, exist_ok=True)
        table.sort_values(["drug", "gene", "site"]).to_csv(arguments.out, index=False)
        print(f"\nwrote {arguments.out}")


if __name__ == "__main__":
    main()
