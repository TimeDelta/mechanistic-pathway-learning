"""How many study drugs have published data on their concentration at the site of action in the brain?

The question of 10 October 2026: can a prediction at the peak of a drug's effect-site concentration be supplied with
data? Three kinds of data would be needed, and this script counts the study drugs that have each:

  extent     a measured unbound brain-to-plasma partition coefficient (Kp,uu), which scales the unbound plasma
             concentration to the unbound brain concentration at steady state;
  rate       what sets how fast the brain follows plasma: a blood-brain barrier permeability together with a brain tissue
             binding (the inputs of the intrinsic brain equilibration half-life of Liu et al. 2005,
             doi:10.1124/jpet.104.079319), a brain equilibration half-life itself, or a human effect-compartment rate
             constant (ke0);
  occupancy  a human target occupancy by PET or SPECT measured at more than one time after a dose, which is the direct
             measurement of what the first two would be used to calculate.

Inputs. The tables under --tables-dir were copied on 10 and 11 October 2026 by three research passes from the sources
listed in the sources.md of each subdirectory; a row was written only for a compound and value read in a retrieved
file. They are transcriptions, not downloads, and sit under data/raw, which git ignores. Many rows of different files
repeat one original measurement, so the count here is of drugs with at least one row, never of measurements.
--keyword-counts is the output of effect_site_keyword_counts.py: PubMed keyword counts, which are upper bounds.

A table name is matched to a study drug when, after lower-casing and removing salt names, stereo prefixes, isotope
labels and footnote marks, it equals the drug's study label, one of its RxNorm names or a listed alias. No partial
match is made, so a metabolite or an ester is never counted for its parent.

Writes one row per drug to --output and the summary to --markdown-output. Counting only: no label, table or
configuration of the study is changed.

Usage:
  python experiments/label_review/effect_site_keyword_counts.py
  python experiments/label_review/effect_site_coverage.py
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

SETS = (("in_current", "current rule"), ("in_approved", "with the three approved label changes"), ("in_approved_with_list", "and the draft list of admitted ATC groups"))
SELECTIONS = {"in_current": Path("data/processed/label_selection/better_v2_full_v3_parkinsonism.parquet"), "in_approved": Path("data/processed/label_review/selection_approved.parquet"),
              "in_approved_with_list": Path("data/processed/label_review/selection_approved_with_list.parquet")}
SALT_WORDS = ("hydrochloride", "dihydrochloride", "hcl", "hydrobromide", "hbr", "maleate", "fumarate", "hemifumarate", "tartrate", "bitartrate", "citrate", "succinate", "mesylate",
              "mesilate", "besylate", "sulfate", "sulphate", "phosphate", "acetate", "sodium", "potassium", "calcium", "oxalate", "lactate",
              "nitrate", "bromide", "chloride", "monohydrate", "dihydrate", "anhydrous", "free base", "base")
# spellings of the sources and names under which a study drug is printed elsewhere; left side as normalised, right side a study name
ALIASES = {
    "midasolam": "midazolam", "alprasolam": "alprazolam", "triasolam": "triazolam", "carbamazapine": "carbamazepine", "galanthamine": "galantamine", "paroxetin": "paroxetine",
    "apomorphine-r": "apomorphine", "o-desmethyl-venlafaxine": "desvenlafaxine", "9-oh risperidone": "paliperidone", "buproprion": "bupropion", "selegeline": "selegiline", "l-deprenyl": "selegiline", "deprenyl": "selegiline", "paracetamol": "acetaminophen",
    "valproic acid": "valproate", "pethidine": "meperidine", "dosulepin": "dothiepin", "lignocaine": "lidocaine", "l-dopa": "levodopa", "hydrocortisone": "cortisol",
    "9-hydroxyrisperidone": "paliperidone", "9-oh-risperidone": "paliperidone", "9-hydroxy-risperidone": "paliperidone", "o-desmethylvenlafaxine": "desvenlafaxine",
    "gamma-hydroxybutyrate": "hydroxybutyrate", "ghb": "hydroxybutyrate", "oxybate": "hydroxybutyrate", "tetrahydrocannabinol": "dronabinol", "delta-9-tetrahydrocannabinol": "dronabinol",
    "thc": "dronabinol", "diphenylhydantoin": "phenytoin", "amfetamine": "amphetamine", "d-amphetamine": "dextroamphetamine", "dexamfetamine": "dextroamphetamine",
    "metamfetamine": "methamphetamine", "oxcarbamazepine": "oxcarbazepine", "levetiracetame": "levetiracetam", "escitalopram oxalate": "escitalopram", "noradrenaline": "norepinephrine",
    "adrenaline": "epinephrine", "glibenclamide": "glyburide", "salbutamol": "albuterol", "esketamine": "esketamine", "s-ketamine": "esketamine", "(s)-ketamine": "esketamine",
    "levothyroxine": "thyroxine", "norethindrone": "norethisterone", "liothyronine": "triiodothyronine", "ezogabine": "retigabine", "vortioxetine": "lu", "etonogestrel": "implanon", "l-t3": "triiodothyronine",
}
PERMEABILITY = re.compile(r"\bPS\b|logPS|PcS|PoS|\bPc\b|\bPo\b|\bKin\b|CLin|\bK1\b|Q_PL_ECF|permeab", re.I)
BRAIN_BINDING = re.compile(r"fu,b\b|fu,brain|fu\(brain\)|Vu,br|Vu\(brain\)|Vu,brain|free fraction in brain", re.I)
EQUILIBRATION = re.compile(r"equilibration half-life|equilibrium half-life|brain entry time|Tmax brain|displacement half-time", re.I)
EFFECT_COMPARTMENT = re.compile(r"ke0|keo\b|tPeak|effect-site|effect site", re.I)
# Depot esters (decanoate, enanthate, pamoate) are not in SALT_WORDS: an occupancy followed for weeks after a depot injection says
# nothing about the peak and trough of an oral dose, so "haloperidol decanoate" is not counted for haloperidol.
# Occupancy time courses that are not a course within hours to days of a dose, by my reading of the times each study states:
OCCUPANCY_COURSE_NOT_WITHIN_DAYS_OF_A_DOSE = {
    "cariprazine": "occupancy 4 h after the dose on days 1, 4 and 15 of daily dosing: accumulation across days, one time per dose",
    "rasagiline": "enzyme recovery 2 to 6 weeks after the last dose of an irreversible inhibitor",
    "selegiline": "1 h after a dose and 9 to 28 days later, an irreversible inhibitor",
}
SOFTWARE_DEFAULT_FILES = {"shafer2026.csv"}  # a dosing-simulation program's defaults, several marked provisional in the source; counted apart
WITHIN_DAYS = "occupancy: at more than one time within hours to days of a dose"
UNCHECKED_FILES = {"chembl_kpuu_typed_other_documents.csv"}  # ChEMBL activity records typed Kp,uu whose assay descriptions were not read against the papers


def normalise(name: str) -> str:
    text = str(name).lower().replace("‐", "-").replace("‑", "-").replace("–", "-")
    text = re.sub(r"\[[^\]]*\]", " ", text)  # isotope labels such as [11C]
    text = re.sub(r"[*#†‡§¶]+", " ", text)
    text = re.sub(r"\((?:r|s|rs|r,s|\+|-|\+/-|±|\+-)\)-?", " ", text)
    text = re.sub(r"^(?:dl|d,l|rac|racemic|\(\+\)|\(-\))[- ]", " ", text.strip())
    text = re.sub(r"\([^)]*\)", " ", text)
    for salt in SALT_WORDS:
        text = re.sub(rf"\b{re.escape(salt)}\b", " ", text)
    text = re.sub(r"\b(?:extended|sustained|immediate|controlled|prolonged)[- ]release\b|\b(?:er|xr|sr|ir|xl)\b", " ", text)  # oral formulations of one drug
    text = re.sub(r"\binhaler\b|\band\s*$", " ", text)
    text = re.sub(r"\s+", " ", text).strip(" -,;")
    return ALIASES.get(text, text)


def study_keys(drugs: pd.DataFrame) -> dict[str, str]:
    """normalised name -> perturbation_id, over the study label and the RxNorm names of every drug."""
    key_to_drug: dict[str, str] = {}
    for drug in drugs.itertuples(index=False):
        names = [drug.perturbation_label, *[name for name in str(drug.rxnorm_names).split(";") if name]]
        for name in names:
            key = normalise(name)
            if key and key_to_drug.setdefault(key, drug.perturbation_id) != drug.perturbation_id:
                raise ValueError(f"the name {key!r} belongs to two study drugs")
    for alias, study_name in ALIASES.items():
        if normalise(study_name) in key_to_drug or study_name in key_to_drug:
            key_to_drug.setdefault(alias, key_to_drug.get(study_name, key_to_drug.get(normalise(study_name))))
    return key_to_drug


def read_rows(directory: Path) -> list[dict]:
    rows = []
    for path in sorted(directory.glob("*.csv")):
        with open(path, newline="", encoding="utf-8") as table_file:
            for row in csv.DictReader(table_file):
                rows.append({**row, "file": path.name})
    return rows


def species_class(species: str) -> str:
    text = species.lower()
    if "human" in text and "corrected for human" not in text and "scaling to human" not in text:
        return "human"
    if any(word in text for word in ("monkey", "primate")):
        return "non-human primate"
    return "other animal" if any(word in text for word in ("rat", "mouse", "mice", "rodent", "dog", "pig", "guinea")) else "not stated"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--drug-names", type=Path, default=Path("data/processed/label_review/drug_names.csv"))
    parser.add_argument("--tables-dir", type=Path, default=Path("data/raw/label_review/effect_site_tables"))
    parser.add_argument("--keyword-counts", type=Path, default=Path("data/processed/label_review/effect_site_keyword_counts.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_review/effect_site_coverage.csv"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/effect_site_data_coverage.md"))
    arguments = parser.parse_args()

    drugs = pd.read_csv(arguments.drug_names).fillna("")
    key_to_drug = study_keys(drugs)
    flags: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))  # perturbation_id -> flag -> source files

    def match(name: str) -> str | None:
        """The study drug a printed name stands for; "brand/drug" and "drug/synonym" are matched by either part."""
        for candidate in [name, *(name.split("/") if "/" in name else [])]:
            drug = key_to_drug.get(normalise(candidate))
            if drug is not None:
                return drug
        return None

    unmatched: dict[str, set[str]] = defaultdict(set)
    for row in read_rows(arguments.tables_dir / "kpuu"):
        drug = match(row["compound_as_printed"])
        if drug is None:
            unmatched["kpuu"].add(normalise(row["compound_as_printed"]))
            continue
        if row["file"] in UNCHECKED_FILES:
            flags[drug]["extent, unchecked ChEMBL record only"].add(row["file"])
            continue
        in_cerebrospinal_fluid = "csf" in row["quantity"].lower() or row["method"].strip().lower().startswith("csf")
        kind = species_class(row["species"])
        flags[drug]["extent: any measured Kp,uu"].add(row["file"])
        if not in_cerebrospinal_fluid:
            flags[drug]["extent: brain tissue or brain fluid, not cerebrospinal fluid"].add(row["file"])
        if "microdialysis" in row["method"].lower():
            flags[drug]["extent: by microdialysis"].add(row["file"])
        if kind == "human":
            flags[drug]["extent: a human value"].add(row["file"])
        if kind == "non-human primate":
            flags[drug]["extent: a non-human primate value"].add(row["file"])
    for row in read_rows(arguments.tables_dir / "rate"):
        drug = match(row["compound_as_printed"])
        if drug is None:
            unmatched["rate"].add(normalise(row["compound_as_printed"]))
            continue
        quantity = row["quantity"]
        if EFFECT_COMPARTMENT.search(quantity):
            name = "rate: effect-compartment rate, simulation program default" if row["file"] in SOFTWARE_DEFAULT_FILES else "rate: human effect-compartment rate (ke0) from a journal"
            flags[drug][name].add(row["file"])
        elif EQUILIBRATION.search(quantity):
            flags[drug]["rate: a brain equilibration half-life or time, measured or calculated in the source"].add(row["file"])
        elif BRAIN_BINDING.search(quantity):
            flags[drug]["rate: brain tissue binding"].add(row["file"])
        elif PERMEABILITY.search(quantity) and row["file"] != "liu2005.csv":  # liu2005.csv holds the seven names without values
            flags[drug]["rate: blood-brain barrier permeability"].add(row["file"])
    for row in read_rows(arguments.tables_dir / "occupancy"):
        drug = match(row["drug_as_printed"])
        if drug is None:
            unmatched["occupancy"].add(normalise(row["drug_as_printed"]))
            continue
        flags[drug]["occupancy: any human PET or SPECT study read"].add(row["pmid"])
        if row["time_course"].strip().lower() == "yes":
            flags[drug]["occupancy: measured at more than one time"].add(row["pmid"])
            if normalise(row["drug_as_printed"]) not in OCCUPANCY_COURSE_NOT_WITHIN_DAYS_OF_A_DOSE:
                flags[drug][WITHIN_DAYS].add(row["pmid"])

    both_inputs = "rate: permeability and brain binding both (the inputs of the equilibration half-life)"
    for drug in drugs.perturbation_id:
        flags[drug].setdefault(WITHIN_DAYS, set())
    any_rate = "rate: any of the above except a program default"
    all_three = "extent, rate and an occupancy time course"
    for drug in drugs.perturbation_id:
        if flags[drug].get("rate: blood-brain barrier permeability") and flags[drug].get("rate: brain tissue binding"):
            flags[drug][both_inputs] = flags[drug]["rate: blood-brain barrier permeability"] | flags[drug]["rate: brain tissue binding"]
        rate_sources = set().union(*[flags[drug].get(name, set()) for name in (both_inputs, "rate: a brain equilibration half-life or time, measured or calculated in the source",
                                                                              "rate: human effect-compartment rate (ke0) from a journal")])
        if rate_sources:
            flags[drug][any_rate] = rate_sources
        if flags[drug].get("extent: any measured Kp,uu") and rate_sources and flags[drug].get(WITHIN_DAYS):
            flags[drug][all_three] = {"all"}
        if flags[drug].get("extent: any measured Kp,uu") and rate_sources:
            flags[drug]["extent and rate"] = {"both"}
    order = ["extent: any measured Kp,uu", "extent: brain tissue or brain fluid, not cerebrospinal fluid", "extent: by microdialysis", "extent: a non-human primate value", "extent: a human value",
             "extent, unchecked ChEMBL record only",
             "rate: blood-brain barrier permeability", "rate: brain tissue binding", both_inputs, "rate: a brain equilibration half-life or time, measured or calculated in the source",
             "rate: human effect-compartment rate (ke0) from a journal", "rate: effect-compartment rate, simulation program default", any_rate,
             "occupancy: any human PET or SPECT study read", "occupancy: measured at more than one time", WITHIN_DAYS, "extent and rate", all_three]
    for drug in drugs.perturbation_id:  # "unchecked only" means no checked source has the drug
        if flags[drug].get("extent: any measured Kp,uu"):
            flags[drug].pop("extent, unchecked ChEMBL record only", None)

    table = drugs[["perturbation_id", "perturbation_label", "rxnorm_names", "in_current", "in_approved", "in_approved_with_list"]].copy()
    for name in order:
        table[name] = [bool(flags[drug].get(name)) for drug in table.perturbation_id]
        table[name + " | sources"] = [";".join(sorted(flags[drug].get(name, set()))) for drug in table.perturbation_id]
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(arguments.output, index=False)

    keyword = pd.read_csv(arguments.keyword_counts) if arguments.keyword_counts.exists() else None
    lines = [
        "# Data on a drug's concentration at its site of action: how many study drugs have any",
        "",
        "Generated by `experiments/label_review/effect_site_coverage.py`. Do not edit by hand.",
        "",
        "Why this table exists. A prediction at the peak of a drug's concentration at its site of action in the brain needs, per drug, the extent of unbound brain "
        "exposure (the unbound brain-to-plasma partition coefficient, Kp,uu), the rate at which the brain follows plasma, and for a check a measured target occupancy "
        "over time. This counts the study drugs for which each was found in a retrieved source.",
        "",
        "How to read it. A drug is counted when a retrieved table names it; the count is of drugs, not of measurements, and many tables repeat one original value. "
        "The tables were copied by three research passes on 10 and 11 October 2026 from the sources listed beside them (`data/raw/label_review/effect_site_tables/*/sources.md`, "
        "not in git); a row was written only for a compound and value read in a retrieved file. Sources that could not be retrieved are listed there too, "
        "so every count is a lower bound on what is published and the keyword counts at the end are upper bounds. Most extent and rate values are from rats or mice.",
        "",
        "## Counts",
        "",
        "| | " + " | ".join(f"{label} ({int(table[column].astype(bool).sum())} drugs)" for column, label in SETS) + " |",
        "| --- | " + " | ".join("---" for _ in SETS) + " |",
    ]
    for name in order:
        lines.append(f"| {name} | " + " | ".join(str(int((table[column].astype(bool) & table[name]).sum())) for column, _ in SETS) + " |")
    kept_pairs = {}
    for column, _ in SETS:
        if SELECTIONS[column].exists():
            selection = pd.read_parquet(SELECTIONS[column])
            kept_pairs[column] = selection[selection.keep.astype(bool) & (selection.perturbation_type == "drug")].groupby("perturbation_id").size()
    if len(kept_pairs) == len(SETS):
        lines += ["", "## The same counts in kept drug pairs", "", "Kept positive drug pairs of each label selection that sit on a drug with the data, and all kept drug pairs in the first row.", "",
                  "| | " + " | ".join(label for _, label in SETS) + " |", "| --- | " + " | ".join("---" for _ in SETS) + " |",
                  "| all kept drug pairs | " + " | ".join(str(int(kept_pairs[column].sum())) for column, _ in SETS) + " |"]
        for name in ("extent: any measured Kp,uu", "extent: a human value", any_rate, "extent and rate", WITHIN_DAYS, all_three):
            lines.append(f"| {name} | " + " | ".join(str(int(kept_pairs[column].reindex(table[table[column].astype(bool) & table[name]].perturbation_id).fillna(0).sum())) for column, _ in SETS) + " |")
    if keyword is not None:
        merged = table.merge(keyword[["perturbation_id", "occupancy_by_pet_or_spect", "occupancy_with_time_words", "effect_compartment_rate", "unbound_brain_concentration"]], on="perturbation_id")
        lines += ["", "## PubMed keyword counts (upper bounds)", "",
                  "Drugs with at least one abstract that names the drug with the words of the query (`experiments/label_review/effect_site_keyword_counts.py`). "
                  "An abstract that names a drug beside the words need not have measured the quantity for it.", "",
                  "| query | " + " | ".join(label for _, label in SETS) + " |", "| --- | " + " | ".join("---" for _ in SETS) + " |"]
        for column, label in (("unbound_brain_concentration", "an unbound brain concentration or Kp,uu"), ("effect_compartment_rate", "an effect compartment or ke0"),
                              ("occupancy_by_pet_or_spect", "occupancy by PET or SPECT"), ("occupancy_with_time_words", "occupancy by PET or SPECT with words for time")):
            lines.append(f"| {label} | " + " | ".join(str(int((merged[set_column].astype(bool) & (merged[column] > 0)).sum())) for set_column, _ in SETS) + " |")
    approved = table[table.in_approved.astype(bool)]
    lines += ["", "## The drugs, under the three approved label changes", ""]
    for name in (WITHIN_DAYS, "extent: a human value", "rate: human effect-compartment rate (ke0) from a journal", all_three):
        lines.append(f"- **{name}** ({int(approved[name].sum())}): " + ", ".join(sorted(approved[approved[name]].perturbation_label, key=str.lower)) + ".")
    without = approved[~approved["extent: any measured Kp,uu"] & ~approved[any_rate] & ~approved["occupancy: any human PET or SPECT study read"]]
    lines.append(f"- **none of the three in any retrieved source** ({len(without)}): " + ", ".join(sorted(without.perturbation_label, key=str.lower)) + ".")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    for column, label in SETS:
        print(label, int(table[column].astype(bool).sum()), {name: int((table[column].astype(bool) & table[name]).sum()) for name in order})
    for kind, names in unmatched.items():
        print(f"{kind}: {len(names)} table names matched no study drug")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
