"""One row per study drug: clearance, extraction ratio and the concentration the drug reaches at its site of action.

Why. Three open points of the study need a number per drug that no table of the repository held:

  - the sequestration edge of a plasma carrier means what it says only for a drug of high hepatic extraction given past
    the first pass (docs/drug_entry_nodes.md); the route half of that condition is tabulated
    (configs/drug_administration_routes.csv) and the extraction half needed a clearance per drug;
  - an occupancy at a target, theta = C / (C + K), needs the free concentration C at the target
    (docs/off_target_scoping.md, decision 3: the 90-percent placeholder or a measured concentration);
  - a reading at the peak of a drug's concentration at its site of action (docs/label_source_review.md, section 10).

What is joined, each from a published table kept under data/raw/exposure (not in git; docs/exposure_sources/ records
where each was read and experiments/exposure_sources/ holds the scripts that transcribed it):

  Lombardo et al. 2018     intravenous plasma clearance, volume of distribution, half-life, fraction unbound in plasma
  FDA drug labels          percent bound to plasma proteins, where Lombardo's table has no fraction unbound
                           (configs/drug_label_protein_binding.csv, read by experiments/read_label_protein_binding.py)
  Benet et al. 2011        percent of the dose excreted unchanged in urine (four values corrected by Hosey et al. 2016)
  Paixao 2012, Uchimura 2010 and the EPA httk table    blood-to-plasma concentration ratio
  Schulz et al. 2020       therapeutic total plasma concentration (mg/L) and half-life
  the brain-exposure tables of docs/effect_site_data_coverage.md    unbound brain-to-plasma ratio, Kp,uu,brain
  PubChem                  molecular weight of the study compound, for the molar conversion (by compound identifier, or
                           by name for a drug whose identifier is an RxNorm concept)

What is computed, with every assumption named in the output:

  hepatic extraction ratio   E = CL_plasma * (1 - fraction excreted unchanged in urine) / (blood-to-plasma ratio) / Q_H.
                             Assumptions: clearance that is not renal excretion of unchanged drug is hepatic; a missing
                             blood-to-plasma ratio is taken as 1 and flagged; Q_H = 20.7 mL/min/kg (1,450 mL/min in a 70 kg
                             adult). A value above 1 means one of the assumptions fails for that drug and is flagged.
  unbound plasma concentration   therapeutic total concentration * fraction unbound, in micromolar.
  unbound brain concentration    unbound plasma concentration * Kp,uu,brain. The Kp,uu is the median of the retrieved
                                 values of the best species class available (human, then non-human primate, then other
                                 animals), brain tissue or brain fluid only, wild-type animals only.

Nothing here changes a label, a graph or a model input. The joined table is written under data/processed/exposure,
which is not in git, because the licences of two of the sources do not cover redistribution of their tables
(docs/drug_exposure_table.md states each). The generated document carries counts and the method.

Usage:
  PYTHONPATH=. python experiments/build_drug_exposure_table.py
  PYTHONPATH=. python experiments/build_drug_exposure_table.py --evidence-dir data/processed/label_review/evidence_approved_no_cap --name approved_no_cap
"""
from __future__ import annotations

import argparse
import json
import math
import re
import statistics
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "label_review"))
from effect_site_coverage import UNCHECKED_FILES, normalise, species_class, study_keys

EXPOSURE_DIRECTORY = Path("data/raw/exposure")
BRAIN_EXPOSURE_DIRECTORY = Path("data/raw/label_review/effect_site_tables/kpuu")
DRUG_NAMES = Path("data/processed/label_review/drug_names.csv")
ROUTES = Path("configs/drug_administration_routes_v4.csv")  # experiments/tabulate_drug_routes.py on the v4 table
CARRIERS = Path("configs/drug_plasma_carriers_v4.csv")  # experiments/scope_label_plasma_binder_coverage.py on the v4 table
LABEL_PROTEIN_BINDING = Path("configs/drug_label_protein_binding.csv")  # experiments/read_label_protein_binding.py
OUTPUT_DIRECTORY = Path("data/processed/exposure")
OUTPUT_DOCUMENT = Path("docs/drug_exposure_table.md")
HEPATIC_BLOOD_FLOW_ML_PER_MIN_PER_KG = 20.7  # 1,450 mL/min in a 70 kg adult
LOW_EXTRACTION_BELOW, HIGH_EXTRACTION_ABOVE = 0.3, 0.7  # the usual cut-offs for low and high hepatic extraction
SPECIES_PREFERENCE = ("human", "non-human primate", "other animal", "not stated")
NOT_WILD_TYPE = re.compile(r"knock|deficient|mdr1a|bcrp\(|\(-/-\)|-/-", re.IGNORECASE)
BRAIN_NOT_CEREBROSPINAL_FLUID = re.compile(r"kp,uu,(brain|bbb)|unbound brain", re.IGNORECASE)
REGIONAL_VALUE = re.compile(r"kp,uu,roi|region", re.IGNORECASE)
USER_AGENT = "mechanistic-pathway-learning/0.1 (drug exposure table)"


def number_or_none(value) -> float | None:
    try:
        number = float(str(value).replace(",", ".").replace("<", "").replace(">", "").replace("~", "").strip())
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def unbracketed_range(printed_cell: str) -> tuple[float, float] | None:
    """The range a Schulz cell prints outside its brackets: "(5-15) 10-20" and "10-20 (5-15)" both give (10, 20).

    Bracketed parts are the wider or older ranges the table adds beside its main one. A cell that prints a list, another
    unit or no number gives None."""
    text = re.sub(r"\^\{[^}]*\}", "", str(printed_cell))
    text = re.sub(r"\([^)]*\)", " ", text)
    text = re.sub(r"\b(appr|approx|ca)\.?\s*", "", text, flags=re.IGNORECASE).strip()
    if not text or ";" in text or re.search(r"[a-zA-Z%/]", text):
        return None
    match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*[-–]\s*(\d+(?:\.\d+)?)", text)
    if match:
        return float(match.group(1)), float(match.group(2))
    match = re.fullmatch(r"(\d+(?:\.\d+)?)", text)
    return (float(match.group(1)), float(match.group(1))) if match else None


def study_drugs(evidence_directory: Path) -> pd.DataFrame:
    """perturbation_id, perturbation_label and the other names a drug is printed under (RxNorm names, when the name table holds the drug)."""
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet", columns=["perturbation_id", "perturbation_type", "perturbation_label"])
    drugs = records[records.perturbation_type == "drug"].drop_duplicates("perturbation_id")[["perturbation_id", "perturbation_label"]].reset_index(drop=True)
    drugs["rxnorm_names"] = ""
    if DRUG_NAMES.exists():
        names = pd.read_csv(DRUG_NAMES)
        rxnorm_names_of = dict(zip(names.perturbation_id, names.rxnorm_names.fillna("")))
        drugs["rxnorm_names"] = drugs.perturbation_id.map(rxnorm_names_of).fillna("")
    return drugs


def pubchem_compound_id(perturbation_id: str) -> int | None:
    match = re.fullmatch(r"CID[01](\d+)", perturbation_id)  # CID1 marks a flat (stereo-merged) compound and CID0 a stereo-specific one in the source of the labels
    return int(match.group(1)) if match else None


def molecular_weights(perturbation_ids: list[str], cache_path: Path) -> dict[str, float]:
    """Molecular weight of each study compound from PubChem's property service, cached beside the raw tables."""
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    wanted = {perturbation_id: pubchem_compound_id(perturbation_id) for perturbation_id in perturbation_ids}
    missing = sorted({compound_id for perturbation_id, compound_id in wanted.items() if compound_id is not None and str(compound_id) not in cache})
    for start in range(0, len(missing), 100):
        batch = missing[start:start + 100]
        address = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{','.join(str(compound_id) for compound_id in batch)}/property/MolecularWeight/JSON"
        with urllib.request.urlopen(urllib.request.Request(address, headers={"User-Agent": USER_AGENT}), timeout=120) as response:
            for entry in json.load(response)["PropertyTable"]["Properties"]:
                cache[str(entry["CID"])] = float(entry["MolecularWeight"])
        time.sleep(0.4)
    if missing:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, indent=1, sort_keys=True))
    return {perturbation_id: cache[str(compound_id)] for perturbation_id, compound_id in wanted.items() if compound_id is not None and str(compound_id) in cache}


def molecular_weights_by_name(drug_names: list[str], cache_path: Path) -> dict[str, float]:
    """Molecular weight from PubChem's name lookup, for a study drug whose identifier is an RxNorm concept and not a compound.

    Only the drug's name goes into the request. A name PubChem does not resolve (an antibody, for one) is cached as absent."""
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    looked_up = False
    for drug_name in drug_names:
        if drug_name in cache:
            continue
        address = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(drug_name)}/property/MolecularWeight/JSON"
        try:
            with urllib.request.urlopen(urllib.request.Request(address, headers={"User-Agent": USER_AGENT}), timeout=60) as response:
                cache[drug_name] = float(json.load(response)["PropertyTable"]["Properties"][0]["MolecularWeight"])
        except urllib.error.HTTPError as error:
            if error.code != 404:
                raise
            cache[drug_name] = None
        looked_up = True
        time.sleep(0.4)
    if looked_up:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(cache, indent=1, sort_keys=True))
    return {drug_name: cache[drug_name] for drug_name in drug_names if cache.get(drug_name) is not None}


def rows_by_drug(table: pd.DataFrame, key_to_drug: dict[str, str], name_column: str = "name") -> dict[str, list[dict]]:
    """Rows of a source table grouped by the study drug their printed name matches; rows naming two analytes are left out."""
    matched: dict[str, list[dict]] = {}
    for row in table.to_dict("records"):
        printed_name = str(row[name_column])
        if " + " in printed_name:
            continue
        drug = key_to_drug.get(normalise(printed_name))
        if drug is not None:
            matched.setdefault(drug, []).append(row)
    return matched


def brain_exposure_summary(key_to_drug: dict[str, str]) -> dict[str, dict]:
    """Per drug: the median Kp,uu,brain of the best species class available, with how many values and from which files."""
    values: dict[str, dict[str, list[tuple[float, str]]]] = {}
    for path in sorted(BRAIN_EXPOSURE_DIRECTORY.glob("*.csv")):
        if path.name in UNCHECKED_FILES:
            continue
        for row in pd.read_csv(path, dtype=str).fillna("").to_dict("records"):
            quantity, species = row.get("quantity", ""), row.get("species", "")
            if not BRAIN_NOT_CEREBROSPINAL_FLUID.search(quantity) or REGIONAL_VALUE.search(quantity) or NOT_WILD_TYPE.search(species + " " + row.get("note", "")[:120]):
                continue
            drug = key_to_drug.get(normalise(row.get("compound_as_printed", "")))
            value = number_or_none(row.get("value_as_printed"))
            if drug is None or value is None or value <= 0:
                continue
            values.setdefault(drug, {}).setdefault(species_class(species), []).append((value, path.name))
    summary = {}
    for drug, by_species in values.items():
        species = next(name for name in SPECIES_PREFERENCE if name in by_species)
        numbers = [value for value, _ in by_species[species]]
        summary[drug] = {"kpuu_brain": statistics.median(numbers), "kpuu_brain_species": species, "kpuu_brain_values": len(numbers),
                         "kpuu_brain_min": min(numbers), "kpuu_brain_max": max(numbers), "kpuu_brain_files": ";".join(sorted({name for _, name in by_species[species]}))}
    return summary


def build_table(evidence_directory: Path, routes_path: Path = ROUTES, carriers_path: Path = CARRIERS) -> pd.DataFrame:
    drugs = study_drugs(evidence_directory)
    key_to_drug = study_keys(drugs)
    lombardo = rows_by_drug(pd.read_csv(EXPOSURE_DIRECTORY / "lombardo2018/lombardo2018.csv"), key_to_drug)
    benet = rows_by_drug(pd.read_csv(EXPOSURE_DIRECTORY / "benet2011_bddcs/benet2011_bddcs.csv"), key_to_drug)
    corrections = pd.read_csv(EXPOSURE_DIRECTORY / "hosey2016_bddcs/hosey2016_table2_parameter_changes.csv")
    corrected_percent = {key_to_drug.get(normalise(row.name)): number_or_none(row.updated_parameter_value)
                         for row in corrections.itertuples() if "excreted" in str(row.parameter).lower()}
    schulz = rows_by_drug(pd.read_csv(EXPOSURE_DIRECTORY / "schulz2020/schulz2020.csv"), key_to_drug)
    ratio_sources = (("Paixao 2012", "blood_to_plasma/paixao2012_table2.csv", "blood_to_plasma_ratio_rb"),
                     ("Uchimura 2010", "blood_to_plasma/uchimura2010.csv", "blood_to_plasma_ratio_human"),
                     ("EPA httk", "blood_to_plasma/httk_human_rblood2plasma.csv", "blood_to_plasma_ratio_human"))
    ratios = [(label, rows_by_drug(pd.read_csv(EXPOSURE_DIRECTORY / file_name), key_to_drug), column) for label, file_name, column in ratio_sources]
    brain = brain_exposure_summary(key_to_drug)
    weights = molecular_weights(list(drugs.perturbation_id), EXPOSURE_DIRECTORY / "pubchem_molecular_weight.json")
    drugs_without_compound_id = [drug.perturbation_label for drug in drugs.itertuples(index=False) if drug.perturbation_id not in weights]
    weights_by_name = molecular_weights_by_name(drugs_without_compound_id, EXPOSURE_DIRECTORY / "pubchem_molecular_weight_by_name.json")
    routes = pd.read_csv(routes_path).set_index("perturbation_id").route_class.to_dict() if routes_path.exists() else {}
    carriers = pd.read_csv(carriers_path).set_index("perturbation_id").carrier.to_dict() if carriers_path.exists() else {}
    label_binding = pd.read_csv(LABEL_PROTEIN_BINDING).set_index("perturbation_id").to_dict("index") if LABEL_PROTEIN_BINDING.exists() else {}

    rows = []
    for drug in drugs.itertuples(index=False):
        row: dict = {"perturbation_id": drug.perturbation_id, "drug": drug.perturbation_label, "molecular_weight": weights.get(drug.perturbation_id),
                     "molecular_weight_source": "PubChem, by compound identifier" if drug.perturbation_id in weights else None,
                     "route_class": routes.get(drug.perturbation_id), "plasma_carrier": carriers.get(drug.perturbation_id)}
        if row["molecular_weight"] is None and drug.perturbation_label in weights_by_name:
            row |= {"molecular_weight": weights_by_name[drug.perturbation_label], "molecular_weight_source": "PubChem, by name"}
        pharmacokinetics = lombardo.get(drug.perturbation_id, [])
        if pharmacokinetics:
            first = pharmacokinetics[0]
            row |= {"clearance_ml_per_min_per_kg": number_or_none(first["clearance_ml_per_min_per_kg"]), "vdss_l_per_kg": number_or_none(first["vdss_l_per_kg"]),
                    "half_life_h_intravenous": number_or_none(first["half_life_h"]), "fraction_unbound_lombardo": number_or_none(first["fraction_unbound_plasma"]),
                    "lombardo_name": first["name"], "lombardo_rows": len(pharmacokinetics)}
        from_label = label_binding.get(drug.perturbation_id)
        if from_label:
            row |= {"fraction_unbound_label": from_label["fraction_unbound"], "fraction_unbound_label_is": from_label["fraction_unbound_is"],
                    "label_bound_99_percent_or_more": bool(from_label["bound_99_percent_or_more"])}
        # the measured table first; the label where the table has no value
        if row.get("fraction_unbound_lombardo") is not None:
            row |= {"fraction_unbound_plasma": row["fraction_unbound_lombardo"], "fraction_unbound_source": "Lombardo et al. 2018"}
        elif from_label:
            row |= {"fraction_unbound_plasma": from_label["fraction_unbound"], "fraction_unbound_source": "FDA label, " + from_label["fraction_unbound_is"]}
        if drug.perturbation_id in benet:
            percent = number_or_none(benet[drug.perturbation_id][0]["percent_excreted_unchanged_urine"])
            if corrected_percent.get(drug.perturbation_id) is not None:
                percent, row["percent_excreted_corrected_by_hosey_2016"] = corrected_percent[drug.perturbation_id], True
            row |= {"percent_excreted_unchanged_urine": percent, "bddcs_class": benet[drug.perturbation_id][0]["bddcs_class"]}
        for label, matched, column in ratios:
            value = number_or_none(matched[drug.perturbation_id][0][column]) if drug.perturbation_id in matched else None
            if value is not None and "blood_to_plasma_ratio" not in row:
                row |= {"blood_to_plasma_ratio": value, "blood_to_plasma_source": label}
        for entry in schulz.get(drug.perturbation_id, []):
            low, high = number_or_none(entry["therapeutic_low_mg_per_l"]), number_or_none(entry["therapeutic_high_mg_per_l"])
            if low is None or high is None:
                parsed = unbracketed_range(entry["therapeutic_printed"])
                low, high = parsed if parsed else (None, None)
            if low is not None and "therapeutic_low_mg_per_l" not in row:
                row |= {"therapeutic_low_mg_per_l": low, "therapeutic_high_mg_per_l": high, "therapeutic_printed": entry["therapeutic_printed"], "schulz_name": entry["name"],
                        "half_life_printed_schulz": entry["half_life_printed"]}
        row |= brain.get(drug.perturbation_id, {})

        clearance, percent_unchanged = row.get("clearance_ml_per_min_per_kg"), row.get("percent_excreted_unchanged_urine")
        if clearance is not None and percent_unchanged is not None:
            ratio = row.get("blood_to_plasma_ratio")
            row["blood_to_plasma_ratio_assumed_one"] = ratio is None
            extraction = clearance * (1.0 - min(percent_unchanged, 100.0) / 100.0) / (ratio or 1.0) / HEPATIC_BLOOD_FLOW_ML_PER_MIN_PER_KG
            row |= {"hepatic_extraction_ratio": extraction, "extraction_ratio_above_one": extraction > 1.0,
                    "extraction_class": "high" if extraction > HIGH_EXTRACTION_ABOVE else ("low" if extraction < LOW_EXTRACTION_BELOW else "intermediate")}
        unbound_fraction, weight = row.get("fraction_unbound_plasma"), row.get("molecular_weight")
        if row.get("therapeutic_low_mg_per_l") is not None and unbound_fraction is not None and weight:
            for bound in ("low", "high"):
                row[f"unbound_plasma_{bound}_micromolar"] = row[f"therapeutic_{bound}_mg_per_l"] / weight * 1000.0 * unbound_fraction
                if row.get("kpuu_brain") is not None:
                    row[f"unbound_brain_{bound}_micromolar"] = row[f"unbound_plasma_{bound}_micromolar"] * row["kpuu_brain"]
        rows.append(row)
    return pd.DataFrame(rows)


def count(table: pd.DataFrame, column: str) -> int:
    return int(table[column].notna().sum()) if column in table.columns else 0


def write_document(table: pd.DataFrame, name: str, evidence_directory: Path, output_path: Path) -> None:
    drugs = len(table)
    with_extraction = table[table.hepatic_extraction_ratio.notna()] if "hepatic_extraction_ratio" in table.columns else table.iloc[:0]
    lines = [
        "# Clearance, extraction and concentration at the site of action, per study drug",
        "",
        "Generated by `experiments/build_drug_exposure_table.py`. Do not edit by hand.",
        "",
        f"Drugs: the {drugs} drugs of `{evidence_directory}` (build name `{name}`). The joined table is `{OUTPUT_DIRECTORY}/drug_exposure_{name}.csv`, which is not in git.",
        "Nothing here changes a label, a graph or a model input.",
        "",
        "## Coverage",
        "",
        "| quantity | source | drugs with a value |",
        "| --- | --- | --- |",
        f"| intravenous plasma clearance | Lombardo et al. 2018 | {count(table, 'clearance_ml_per_min_per_kg')} |",
        f"| fraction unbound in plasma | Lombardo et al. 2018 ({count(table, 'fraction_unbound_lombardo')}); the drug's FDA label where that table has none ({int(table.fraction_unbound_source.astype(str).str.startswith('FDA label').sum())}) | {count(table, 'fraction_unbound_plasma')} |",
        f"| half-life after an intravenous dose | Lombardo et al. 2018 | {count(table, 'half_life_h_intravenous')} |",
        f"| percent excreted unchanged in urine | Benet et al. 2011 | {count(table, 'percent_excreted_unchanged_urine')} |",
        f"| blood-to-plasma ratio | Paixão 2012, Uchimura 2010, EPA httk | {count(table, 'blood_to_plasma_ratio')} |",
        f"| therapeutic total plasma concentration | Schulz et al. 2020 | {count(table, 'therapeutic_low_mg_per_l')} |",
        f"| unbound brain-to-plasma ratio, Kp,uu,brain | the tables of `docs/effect_site_data_coverage.md` | {count(table, 'kpuu_brain')} |",
        f"| molecular weight | PubChem, by compound identifier ({int((table.molecular_weight_source == 'PubChem, by compound identifier').sum())}) or by name ({int((table.molecular_weight_source == 'PubChem, by name').sum())}) | {count(table, 'molecular_weight')} |",
        f"| **hepatic extraction ratio** (computed) | | **{len(with_extraction)}** |",
        f"| **unbound plasma concentration at therapeutic levels** (computed) | | **{count(table, 'unbound_plasma_low_micromolar')}** |",
        f"| **unbound brain concentration at therapeutic levels** (computed) | | **{count(table, 'unbound_brain_low_micromolar')}** |",
        "",
    ]
    if "kpuu_brain_species" in table.columns:
        species_counts = Counter(table.kpuu_brain_species.dropna())
        with_brain = table[table.unbound_brain_low_micromolar.notna()] if "unbound_brain_low_micromolar" in table.columns else table.iloc[:0]
        lines += ["The Kp,uu,brain used is from: " + ", ".join(f"{species} {species_counts[species]}" for species in SPECIES_PREFERENCE if species_counts.get(species)) + ".",
                  "Of the drugs with an unbound brain concentration, the Kp,uu is from: "
                  + ", ".join(f"{species} {int((with_brain.kpuu_brain_species == species).sum())}" for species in SPECIES_PREFERENCE if (with_brain.kpuu_brain_species == species).any()) + ".", ""]
    both = table[table.fraction_unbound_lombardo.notna() & table.fraction_unbound_label.notna() & (table.fraction_unbound_label_is == "a value")] if "fraction_unbound_label" in table.columns else table.iloc[:0]
    if len(both):
        fold = (both.fraction_unbound_label / both.fraction_unbound_lombardo).map(lambda ratio: max(ratio, 1.0 / ratio))
        from_label = table[table.fraction_unbound_source.astype(str).str.startswith("FDA label")]
        lines += ["## Fraction unbound: the two sources compared", "",
                  f"{len(both)} drugs have a fraction unbound in Lombardo's table and a percentage bound in their label (`{LABEL_PROTEIN_BINDING}`, each with its sentence and label version). "
                  f"The two differ by a median factor of {fold.median():.2f}; by less than 1.5 for {int((fold < 1.5).sum())} drugs, by 2 or more for {int((fold >= 2).sum())} and by 3 or more for {int((fold >= 3).sum())}. "
                  f"The largest differences: {', '.join(f'{drug} ({label_value:.3g} from the label, {table_value:.3g} in the table)' for drug, label_value, table_value in both.assign(fold=fold).sort_values('fold', ascending=False).head(5)[['drug', 'fraction_unbound_label', 'fraction_unbound_lombardo']].itertuples(index=False))}.", "",
                  f"Of the {len(from_label)} drugs whose fraction unbound comes from the label, it is a value for {int((from_label.fraction_unbound_label_is == 'a value').sum())}, an upper limit for {int((from_label.fraction_unbound_label_is == 'an upper limit').sum())} "
                  f"(the label says \"greater than\") and a lower limit for {int((from_label.fraction_unbound_label_is == 'a lower limit').sum())}; {int((from_label.label_bound_99_percent_or_more == True).sum())} are bound 99 percent or more, "
                  "where the rounding of a label moves the fraction unbound by a factor of two or more.", ""]
    lines += ["## Hepatic extraction ratio", "",
              f"E = CL × (1 − fraction excreted unchanged in urine) / (blood-to-plasma ratio) / Q_H, with Q_H = {HEPATIC_BLOOD_FLOW_ML_PER_MIN_PER_KG} mL/min/kg.",
              f"Low is below {LOW_EXTRACTION_BELOW} and high above {HIGH_EXTRACTION_ABOVE}, the usual cut-offs.", ""]
    if len(with_extraction):
        class_counts = Counter(with_extraction.extraction_class)
        lines += ["| | drugs |", "| --- | --- |"]
        lines += [f"| {label} extraction | {class_counts.get(label, 0)} |" for label in ("low", "intermediate", "high")]
        lines += [f"| blood-to-plasma ratio missing and taken as 1 | {int(with_extraction.blood_to_plasma_ratio_assumed_one.sum())} |",
                  f"| computed value above 1 (an assumption fails for the drug) | {int(with_extraction.extraction_ratio_above_one.sum())} |", ""]
        with_carrier = with_extraction[with_extraction.plasma_carrier.notna()]
        if len(with_carrier):
            lines += ["### The drugs with a sequestration edge", "",
                      f"{int(table.plasma_carrier.notna().sum())} drugs have a known plasma carrier (`{CARRIERS}`); {len(with_carrier)} of them have an extraction ratio here.",
                      "The edge states that more carrier lowers the drug's effect. For the unbound exposure that holds for a drug of high extraction given past the first pass (`docs/drug_entry_nodes.md`).", "",
                      "| route class | low | intermediate | high |", "| --- | --- | --- | --- |"]
            for route_class in sorted(with_carrier.route_class.dropna().unique()):
                subset = with_carrier[with_carrier.route_class == route_class]
                lines.append(f"| {route_class} | " + " | ".join(str(int((subset.extraction_class == label).sum())) for label in ("low", "intermediate", "high")) + " |")
            inside = with_carrier[(with_carrier.extraction_class == "high") & with_carrier.route_class.isin(["intravenous", "other route past the first pass"])]
            lines += ["", f"High extraction and a route past the first pass: {len(inside)} drugs ({', '.join(sorted(inside.drug)) or 'none'}).", ""]
    lines += ["## What the numbers assume", "",
              "- Clearance that is not renal excretion of unchanged drug is taken as hepatic. A drug cleared by plasma esterases or in other tissues gets too high a value.",
              "- Clearance is of plasma after an intravenous dose; it is divided by the blood-to-plasma ratio to give a clearance of blood, and by 1 where no ratio was found.",
              "- The therapeutic concentration is a total plasma or serum concentration at steady state as compiled by Schulz et al.; the range is the one the table prints outside brackets.",
              "- The Kp,uu,brain is mostly from rats and mice. A human value is used where one exists. Values from knockout animals, from brain regions and from cerebrospinal fluid are left out.",
              "- The unbound brain concentration assumes the steady-state ratio holds at the therapeutic concentration. It says nothing about the time course after a dose.", "",
              "## Sources and their licences", "",
              "| source | licence, as read | where |",
              "| --- | --- | --- |",
              "| Schulz et al. 2020, [doi:10.1186/s13054-020-02915-5](https://doi.org/10.1186/s13054-020-02915-5) | \"This article is licensed under a Creative Commons Attribution 4.0 International License\"; \"The Creative Commons Public Domain Dedication waiver (http://creativecommons.org/publicdomain/zero/1.0/) applies to the data made available in this article, unless otherwise stated in a credit line to the data.\" | the article's permissions in PubMed Central |",
              "| Lombardo et al. 2018, [doi:10.1124/dmd.118.082966](https://doi.org/10.1124/dmd.118.082966) | **permission required.** The copyright holder is the journal's society; the article is free to read under no licence; the only licence its Crossref record carries is for text and data mining and states that a user may not \"Substantially or systematically reproduce, retain or redistribute the Dataset\". The US EPA repository the table was read from carries no permission statement for it. The table is used and not redistributed, which is what those terms permit | `docs/exposure_sources/lombardo2018.md`, `docs/label_source_review.md` 13.9 |",
              "| Benet et al. 2011, [doi:10.1208/s12248-011-9290-9](https://doi.org/10.1208/s12248-011-9290-9); Hosey et al. 2016, [doi:10.1208/s12248-015-9845-2](https://doi.org/10.1208/s12248-015-9845-2) | not open access; copyright of the journal's society. Read from the free copies in PubMed Central | `docs/exposure_sources/benet2011_bddcs.md` and `hosey2016_bddcs.md` |",
              "| blood-to-plasma ratios | the EPA `httk` package (\"License: MIT\"); the EPA's typed copies of Paixão et al. 2012 and Uchimura et al. 2010, not checked against the articles | `docs/exposure_sources/blood_to_plasma.md` |",
              "| FDA drug labels through openFDA | \"Unless otherwise noted, the content, data, documentation, code, and related materials on openFDA is public domain\" (open.fda.gov/license, read 11 October 2026); \"we ask that proper credit be given\", with the example \"Data provided by the U.S. Food and Drug Administration (https://open.fda.gov)\" (open.fda.gov/terms, read 11 October 2026) | `configs/drug_label_protein_binding.csv`, with the label's set_id and effective_time per sentence |",
              "| Kp,uu,brain | per table, in `docs/effect_site_data_coverage.md` | |", "",
              "Because two of these tables are not licensed for redistribution, the joined table stays out of git. The transcription scripts (`experiments/exposure_sources/`), "
              "the notes on each source with its addresses and checksums (`docs/exposure_sources/`) and this script reproduce it.", "",
              "## Not retrieved", "",
              "The consensus guideline for therapeutic drug monitoring in psychiatry (Hiemke et al. 2018, [doi:10.1055/s-0043-116492](https://doi.org/10.1055/s-0043-116492)) gives therapeutic reference ranges for the psychotropic drugs, "
              "and would be the better source for them. It was not read: the publisher's robots file disallows automated readers, the article is not in PubMed Central, and no open licence was found "
              "(`docs/exposure_sources/hiemke2018_agnp.md` lists the 22 routes tried). The table of Schulz et al. already cites the guideline for many of its ranges. "
              "For the owner to open: <https://www.thieme-connect.com/products/ejournals/pdf/10.1055/s-0043-116492.pdf>, and the repository copy that one index labels CC BY, <https://resolver.sub.uni-goettingen.de/purl?gro-2/74377>. A PDF saved under `data/raw/exposure/hiemke2018_agnp/raw/` can be transcribed as the other tables were.", ""]
    output_path.write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--name", default="v4", help="names the output table")
    parser.add_argument("--routes", type=Path, default=ROUTES, help="route class per drug")
    parser.add_argument("--carriers", type=Path, default=CARRIERS, help="plasma carrier per drug, which is the drugs with a sequestration edge")
    parser.add_argument("--document", type=Path, default=OUTPUT_DOCUMENT)
    arguments = parser.parse_args()
    table = build_table(arguments.evidence_dir, arguments.routes, arguments.carriers)
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIRECTORY / f"drug_exposure_{arguments.name}.csv"
    table.to_csv(output_path, index=False)
    write_document(table, arguments.name, arguments.evidence_dir, arguments.document)
    print(f"wrote {output_path} ({len(table)} drugs) and {arguments.document}")
    for column in ("clearance_ml_per_min_per_kg", "fraction_unbound_plasma", "percent_excreted_unchanged_urine", "blood_to_plasma_ratio", "therapeutic_low_mg_per_l", "kpuu_brain",
                   "hepatic_extraction_ratio", "unbound_plasma_low_micromolar", "unbound_brain_low_micromolar"):
        print(f"  {column}: {count(table, column)}")


if __name__ == "__main__":
    main()
