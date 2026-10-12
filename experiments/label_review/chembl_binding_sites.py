"""What ChEMBL says about where a study drug binds: the site_id rows and the free-text binding-site comments.

Part of the label source review (docs/label_source_review.md, section 13.4, the question of ligand binding
competition). The question this answers is whether two drugs that act on one gene act at the same site, because only
then does one displace the other. ChEMBL's drug_mechanism table carries two fields for that: site_id, a curated
binding-site identifier, and binding_site_comment, free text. The text is where the curator's uncertainty is written
down, which is why this script prints every comment verbatim and marks the hedged ones by the word that hedges them.

Counting and quoting only: no label, table, graph or configuration of the study is changed.

Usage:
  python experiments/label_review/chembl_binding_sites.py
  python experiments/label_review/chembl_binding_sites.py --evidence-dir data/processed/evidence_full_v4 \
      --chembl-dir data/raw/chembl --out data/processed/label_review/chembl_binding_sites.csv
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.sider_chemical_matching import pubchem_cid_of_stitch_flat_id
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import load_chembl_caches

# A word in a binding-site comment that marks the curator's own uncertainty about the site, as distinct from a comment
# that states a site plainly ("Beta subunit is binding site"). Matched on word boundaries, lower-cased.
HEDGING_WORDS = (
    "likely", "unclear", "not clear", "unknown", "probably", "probable", "possibly", "possible", "putative",
    "presumed", "presumably", "may be", "might", "thought to", "appears", "looks like", "suggested", "uncertain",
)


def study_drugs(evidence_directory: Path) -> pd.DataFrame:
    """One row per drug perturbation of the evidence table: its perturbation id, its label and its PubChem CID."""
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet",
                              columns=["perturbation_id", "perturbation_type", "perturbation_label"])
    drugs = (records[records.perturbation_type == "drug"][["perturbation_id", "perturbation_label"]]
             .drop_duplicates().reset_index(drop=True))
    drugs["pubchem_cid"] = [pubchem_cid_of_stitch_flat_id(identifier) for identifier in drugs.perturbation_id]
    return drugs


def hedging_words_in(comment: str) -> list[str]:
    lowered = comment.lower()
    return [word for word in HEDGING_WORDS if re.search(rf"(?<!\w){re.escape(word)}(?!\w)", lowered)]


def binding_site_rows(evidence_directory: Path, chembl_directory: Path) -> pd.DataFrame:
    """Every ChEMBL mechanism row of a study drug, with its site identifier, its comment and the gene of its target."""
    drugs = study_drugs(evidence_directory)
    caches = load_chembl_caches(chembl_directory)

    rows = []
    for drug in drugs.itertuples(index=False):
        chembl_ids = caches.pubchem_to_chembl.get(str(drug.pubchem_cid), [])
        seen_mechanism_ids = set()
        for chembl_id in chembl_ids:
            for mechanism in caches.mechanisms_by_molecule.get(chembl_id, []):
                if mechanism.get("mec_id") in seen_mechanism_ids:
                    continue
                seen_mechanism_ids.add(mechanism.get("mec_id"))
                target = caches.targets.get(mechanism.get("target_chembl_id") or "", {})
                comment = (mechanism.get("binding_site_comment") or "").strip()
                rows.append({
                    "drug": drug.perturbation_label,
                    "perturbation_id": drug.perturbation_id,
                    "molecule_chembl_id": mechanism.get("molecule_chembl_id"),
                    "target_chembl_id": mechanism.get("target_chembl_id"),
                    "target_pref_name": target.get("pref_name"),
                    "gene_symbols": "|".join(target.get("gene_symbols") or []),
                    "action_type": mechanism.get("action_type"),
                    "site_id": mechanism.get("site_id"),
                    "binding_site_comment": comment,
                    "hedging_words": "|".join(hedging_words_in(comment)),
                    "mechanism_of_action": mechanism.get("mechanism_of_action"),
                })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--out", type=Path, default=None)
    arguments = parser.parse_args()

    table = binding_site_rows(arguments.evidence_dir, arguments.chembl_dir)
    drugs_in_table = table.drug.nunique()
    with_site_id = table[table.site_id.notna()]
    with_comment = table[table.binding_site_comment.astype(bool)]
    hedged = with_comment[with_comment.hedging_words.astype(bool)]

    print(f"mechanism rows for study drugs: {len(table)} over {drugs_in_table} drugs")
    print(f"  with a curated site_id:        {len(with_site_id)} over {with_site_id.drug.nunique()} drugs")
    print(f"  with a binding-site comment:   {len(with_comment)} over {with_comment.drug.nunique()} drugs")
    print(f"    of those, hedged:            {len(hedged)} over {hedged.drug.nunique()} drugs")

    print("\nEvery binding-site comment, verbatim:")
    for row in with_comment.sort_values(["drug", "target_pref_name"]).itertuples(index=False):
        marker = f"  [hedged: {row.hedging_words}]" if row.hedging_words else ""
        print(f"  {row.drug} -> {row.gene_symbols or row.target_pref_name} ({row.action_type}), site_id {row.site_id}")
        print(f"      \"{row.binding_site_comment}\"{marker}")

    print("\nThe targets that a site_id splits, which is what a competition statement needs:")
    for target_name, group in with_site_id.groupby("target_pref_name"):
        sites = sorted({int(site) for site in group.site_id})
        print(f"  {target_name}: {len(group)} rows over {group.drug.nunique()} drugs, site ids {sites}")

    if arguments.out:
        arguments.out.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(arguments.out, index=False)
        print(f"\nwrote {arguments.out}")


if __name__ == "__main__":
    main()
