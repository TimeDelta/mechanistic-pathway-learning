"""A drug's targets per subunit, from the IUPHAR/BPS Guide to PHARMACOLOGY, in place of ChEMBL's complex targets.

Why (measured 13 October 2026, and the owner's decision of that date: "agreed on both 1 and 2"). ChEMBL's
`drug_mechanism` rows often name a protein *complex*, whose `gene_symbols` field lists every subunit, and
map_drug_targets_to_graph_nodes expands a target to all of them and splits its magnitude 1 / n across them. On the
study's drugs that is 114 of 326 mechanism rows over 91 drugs, and the result is not what the drugs do:

- gabapentin on 26 calcium-channel subunits at 1/22 each, including the skeletal-muscle CACNA1S and the T-type
  CACNA1G, CACNA1H and CACNA1I, where it binds the alpha-2-delta subunits CACNA2D1 and CACNA2D2;
- every benzodiazepine on 16 GABA-A subunits, among them GABRA4, GABRA6 and GABRD, which form
  benzodiazepine-insensitive receptors, since the site needs a gamma subunit and an alpha-1, 2, 3 or 5;
- propofol and the barbiturates on the *same* 16 subunits as the benzodiazepines, although they act at a different
  pocket, which erases the one distinction `configs/binding_sites.csv` exists to record;
- amisulpride on 17 serotonin receptors.

GtoPdb names one gene per interaction: of 357 human rows for the drugs matched here, 356 carry a single
`Target Gene Symbol`. Where both sources describe a drug it gives 3 or 4 subunits against ChEMBL's 16 (alprazolam 4,
clonazepam 3, diazepam 3, flunitrazepam 4, zolpidem 1). What replaces what, and the unit is a ChEMBL **row** rather than a drug. The first version of this module replaced a
drug's whole ChEMBL target set wherever GtoPdb covered the drug, and the owner caught what that does (13 October 2026:
"cocaine should NOT narrow to a single target ... it blocks both DAT and NET at minimum ... double check your
assumptions and data availability for these drugs"). GtoPdb holds two human rows for cocaine, both at HTR3A, so
replacing wholesale would have left cocaine with HTR3A alone. ChEMBL's rows come in two kinds and only one is broken:

- a **complex row**, whose `gene_symbols` lists several genes, is the expansion this module exists to undo, and
  GtoPdb's rows naming genes of that complex replace it;
- a **single-gene row** names one target specifically and is kept, always;
- a GtoPdb gene belonging to no ChEMBL row of that drug is added.

So a drug can only gain specificity within a complex GtoPdb describes, never lose a target to a gap in GtoPdb. Under
this rule cocaine goes from 10 targets to 11 - its sodium-channel complex untouched, because GtoPdb names none of its
subunits, and HTR3A added - rather than to 1.

**Where it runs.** Until 17 October 2026 only measurement scripts called this module, and the evidence build kept
seeding every subunit of a ChEMBL complex (decision 8.11). ComplexTargetNarrowing, at the end of the module, is the
rule wired into both label routes: the SIDER route of evidence/assemble_evidence_table.py and the OnSIDES route of
evidence/load_onsides_label_events.py, behind --narrow-complex-targets (on by default) in each of their scripts.

**What this does not fix, and is not for.** Cocaine's sodium-channel rows are themselves the wrong pharmacology for
this study: ChEMBL curates mechanism by approved indication, cocaine is approved as a topical anaesthetic, and its
monoamine transport sits in the off-target table at 300 times the occupancy of its mechanism targets. That is a
question about which table a target is read from, not about subunit expansion, and it is answered elsewhere.

Licence, already verified and recorded in docs/data_sources.md: "The Guide to PHARMACOLOGY database is licensed under
the Open Data Commons Open Database License (ODbL)."; "Its contents are licensed under a Creative Commons
Attribution-ShareAlike 4.0 International License" (guidetopharmacology.org/about.jsp, read 11 October 2026). Release
2026.3, published 16 September 2026. Not redistributed.
"""
from __future__ import annotations

import dataclasses
import hashlib
import re
from pathlib import Path

import pandas as pd

GTOPDB_DIRECTORY = Path("data/raw/gtopdb_2026")
INTERACTIONS_FILE = "interactions.csv"
LIGAND_MAPPING_FILE = "ligand_id_mapping.csv"
HUMAN_SPECIES = "Human"
VERSION_HEADER_ROWS = 1  # the first line is "# GtoPdb Version: ..."

# GtoPdb's Action vocabulary to the same signs ChEMBL's action types take
# (map_drug_targets_to_graph_nodes.ACTION_TYPE_SIGN), so the two routes put a drug on one scale.
GTOPDB_ACTION_SIGN: dict[str, float] = {
    "FULL AGONIST": 1.0,
    "AGONIST": 1.0,
    "ACTIVATION": 1.0,
    "POTENTIATION": 0.5,
    "PARTIAL AGONIST": 0.5,
    "POSITIVE": 0.5,  # positive allosteric modulation
    "BIASED AGONIST": 0.5,  # agonist at some of a receptor's pathways only, so not a full agonist on any one of them
    "NEGATIVE": -0.5,  # negative allosteric modulation
    "ANTAGONIST": -1.0,
    "INHIBITION": -1.0,
    "INVERSE AGONIST": -1.0,
    "PORE BLOCKER": -1.0,
    "IRREVERSIBLE INHIBITION": -1.0,
    "BINDING": 0.0,
    "NEUTRAL": 0.0,
    "MIXED": 0.0,  # opposite directions at one target; a direction-free row rather than a guessed one
    # the rest of the release's vocabulary, which no row of this study's drugs carries (measured 13 October 2026) but
    # which a release or a drug set could bring in, and which a guard test holds this table to
    "COMPETITIVE": -1.0,  # Type Inhibitor throughout: competitive inhibition lowers the target's activity
    "NON-COMPETITIVE": -1.0,  # Type Antagonist or Inhibitor; the site differs, the direction does not
    "FEEDBACK INHIBITION": -1.0,
    "VOLTAGE-DEPENDENT INHIBITION": -1.0,  # Type Gating inhibitor
    "IRREVERSIBLE AGONIST": 1.0,
    "SLOWS INACTIVATION": 1.0,  # a channel that inactivates more slowly carries more current, so the flux rises
    "UNKNOWN": 0.0,
}
# GtoPdb's two competition words, which section 13.4 of docs/label_source_review.md counted: 53 human rows read
# "Competitive" and 4 "Non-competitive" across the release, and none of them belongs to a drug of this study. They are
# in the sign table so a different drug set does not meet an unmapped action, not because they reach anything here.
COMPETITION_ACTIONS = ("COMPETITIVE", "NON-COMPETITIVE")
# Action values that modulate an endogenous ligand's response instead of supplying one of their own, which the
# modulator relation reads (preregistration item 24, partial agonism). Kept beside the sign so one table answers both.
MODULATES_THE_ENDOGENOUS_LIGAND = ("PARTIAL AGONIST", "POSITIVE", "NEGATIVE", "POTENTIATION")
ALLOSTERIC_TYPE = "Allosteric modulator"


def read_interactions(gtopdb_directory: Path = GTOPDB_DIRECTORY, human_only: bool = True) -> pd.DataFrame:
    """The interaction table, human rows only by default. The version line is skipped, not parsed, so a release change
    is visible in the file's own hash rather than silently absorbed."""
    table = pd.read_csv(gtopdb_directory / INTERACTIONS_FILE, low_memory=False, skiprows=VERSION_HEADER_ROWS)
    if human_only:
        table = table[table["Target Species"].astype(str).str.strip() == HUMAN_SPECIES]
    return table


def read_ligand_mapping(gtopdb_directory: Path = GTOPDB_DIRECTORY) -> pd.DataFrame:
    return pd.read_csv(gtopdb_directory / LIGAND_MAPPING_FILE, low_memory=False, skiprows=VERSION_HEADER_ROWS)


def ligand_id_of_drug(pubchem_cid_of_drug: dict[str, int], name_of_drug: dict[str, str],
                      gtopdb_directory: Path = GTOPDB_DIRECTORY) -> tuple[dict[str, int], dict[str, str]]:
    """perturbation id -> GtoPdb ligand id, and perturbation id -> how it was matched.

    PubChem CID first, because an identifier is exact; the drug's name second, lowercased and stripped, because GtoPdb
    carries approved names for drugs whose CID it does not record. Matching by name alone reaches more of the study's
    drugs than by CID alone, so both are used and which one matched is reported, never guessed at."""
    mapping = read_ligand_mapping(gtopdb_directory)
    cid_column = pd.to_numeric(mapping["PubChem CID"], errors="coerce")
    ligand_of_cid = {int(cid): int(ligand) for cid, ligand in zip(cid_column, mapping["Ligand id"]) if pd.notna(cid)}
    ligand_of_name: dict[str, int] = {}
    for name, ligand in zip(mapping["Name"], mapping["Ligand id"]):
        key = str(name).strip().lower()
        if key:
            ligand_of_name.setdefault(key, int(ligand))
    ligand_of_drug: dict[str, int] = {}
    route: dict[str, str] = {}
    for perturbation_id in {*pubchem_cid_of_drug, *name_of_drug}:
        cid = pubchem_cid_of_drug.get(perturbation_id)
        name = str(name_of_drug.get(perturbation_id, "")).strip().lower()
        if cid is not None and int(cid) in ligand_of_cid:
            ligand_of_drug[perturbation_id], route[perturbation_id] = ligand_of_cid[int(cid)], "pubchem_cid"
        elif name in ligand_of_name:
            ligand_of_drug[perturbation_id], route[perturbation_id] = ligand_of_name[name], "name"
    return ligand_of_drug, route


def sign_of_action(action: str | float | None) -> float | None:
    """The sign of a GtoPdb Action, or None where the vocabulary does not name one, so an unknown action is reported
    rather than silently taken as zero."""
    if action is None or (isinstance(action, float) and pd.isna(action)):
        return None
    return GTOPDB_ACTION_SIGN.get(str(action).strip().upper())


def targets_per_drug(ligand_of_drug: dict[str, int], interactions: pd.DataFrame | None = None,
                     gtopdb_directory: Path = GTOPDB_DIRECTORY, primary_only: bool = False) -> dict[str, pd.DataFrame]:
    """perturbation id -> its GtoPdb rows, one gene per row, with the sign its Action takes.

    primary_only keeps the rows GtoPdb marks `Primary Target`, which is the narrower reading of a drug's mechanism.
    Rows with no gene symbol, and rows whose Action the vocabulary does not name, are kept in the frame with a null
    sign so a caller counts them rather than losing them."""
    table = read_interactions(gtopdb_directory) if interactions is None else interactions
    if primary_only:
        table = table[table["Primary Target"].astype(str).str.strip().str.lower().isin(("t", "true", "yes"))]
    ligand_column = pd.to_numeric(table["Ligand ID"], errors="coerce")
    per_drug: dict[str, pd.DataFrame] = {}
    for perturbation_id, ligand_id in ligand_of_drug.items():
        rows = table[ligand_column == ligand_id].copy()
        if rows.empty:
            continue
        rows["gene_symbol"] = rows["Target Gene Symbol"].astype(str).str.strip()
        rows["sign"] = [sign_of_action(action) for action in rows["Action"]]
        rows["modulates_the_endogenous_ligand"] = [
            str(action).strip().upper() in MODULATES_THE_ENDOGENOUS_LIGAND for action in rows["Action"]]
        per_drug[perturbation_id] = rows
    return per_drug


def one_gene_per_row(rows: pd.DataFrame) -> pd.DataFrame:
    """The rows that name exactly one gene and take a sign: the ones that can replace a ChEMBL complex target.

    A GtoPdb gene symbol field occasionally carries a subunit list of its own ("GABRA1|GABRB2|GABRG2" for a named
    receptor assembly). Those are left out here rather than expanded, because expanding them is the error this module
    exists to undo; a caller that wants them has them in the full frame."""
    named = rows[rows.gene_symbol.notna() & (rows.gene_symbol != "") & (rows.gene_symbol.str.lower() != "nan")]
    single = named[~named.gene_symbol.str.contains(r"[|,;]", regex=True, na=False)]
    return single[single.sign.notna()]


# Subunits a binding site needs that an affinity table will not name (the owner, 13 October 2026: "ya absolutely add
# back obligate partners. that sounds obvious"). The narrowing drops them systematically, and the reason generalises:
# GtoPdb gives a row per subunit it holds affinity data for, and affinity is measured against the subunit that
# differentiates the pharmacology, so an obligate partner that does not change affinity never earns a row of its own.
# Measured on this study's drugs, the narrowing dropped GABRG2 from all five benzodiazepines while keeping GABRA1 to
# GABRA3 - and the benzodiazepine site IS the alpha/gamma2 interface, so that is half a binding site.
# The sites are those `configs/binding_sites.csv` names, and each entry is a structural fact about the site rather
# than a measurement: a subunit the site cannot form without.
OBLIGATE_SUBUNITS_OF_A_SITE: dict[str, tuple[str, ...]] = {
    # the classical benzodiazepine site sits at the interface between an alpha subunit and gamma2; receptors built on
    # alpha4 or alpha6, and delta-containing receptors, which carry no gamma at all, are insensitive to it
    "benzodiazepine site": ("GABRG2",),
    # the barbiturate and general-anaesthetic site is on the beta subunit's transmembrane domain and needs no gamma,
    # which is the distinction from the benzodiazepine site that the complex expansion erased
    "barbiturate site": ("GABRB1", "GABRB2", "GABRB3"),
    # the orthosteric GABA site is at the beta/alpha interface
    "GABA site": ("GABRB1", "GABRB2", "GABRB3"),
    # gabapentin and pregabalin bind the alpha-2-delta auxiliary subunits, both of them
    "alpha2delta-1": ("CACNA2D1", "CACNA2D2"),
}


def obligate_partners(sites, gene_symbols_of_the_complex: set[str] | None = None) -> set[str]:
    """The subunits the named sites cannot form without, restricted to the complex's own genes when one is given.

    A partner the drug's complex does not contain is left out rather than invented: the restriction is what keeps this
    a correction to the narrowing and not a second source of targets."""
    partners: set[str] = set()
    for site in sites:
        partners.update(OBLIGATE_SUBUNITS_OF_A_SITE.get(str(site).strip(), ()))
    return partners if gene_symbols_of_the_complex is None else partners & gene_symbols_of_the_complex


# How many subunits of a complex GtoPdb must name before its rows narrow that complex. One row is an assay against
# whichever subtype was convenient to test, not a statement of selectivity: GtoPdb holds a single human row for
# lidocaine (SCN5A, the cardiac subtype its antiarrhythmic data comes from), for phenytoin (SCN2A) and for lacosamide
# (SCN9A), and all three block sodium channels broadly. Two or more named subunits is evidence that selectivity was
# looked at, and it is what narrows the benzodiazepines from 16 GABA-A subunits to the 3 or 4 they act through.
SUBUNITS_GTOPDB_MUST_NAME_TO_NARROW_A_COMPLEX = 2


def replace_complex_rows(chembl_rows: pd.DataFrame, gtopdb_rows: pd.DataFrame,
                         gene_symbols_column: str = "gene_symbols",
                         subunits_needed: int = SUBUNITS_GTOPDB_MUST_NAME_TO_NARROW_A_COMPLEX,
                         sites: "tuple[str, ...] | None" = None) -> tuple[dict[str, float], dict[str, str]]:
    """gene -> sign for one drug, and gene -> why it is there, under the row-wise rule in this module's docstring.

    chembl_rows are that drug's rows of the mechanism table, each with a gene_symbols field that may list several genes
    and an action; gtopdb_rows are one_gene_per_row's output for the same drug. A ChEMBL row naming one gene is kept
    whatever GtoPdb says. A ChEMBL row naming several is replaced by the GtoPdb rows whose gene is among them, but only
    where GtoPdb names at least subunits_needed of them: GtoPdb's silence about a complex is not evidence against it,
    and one named subunit is an assay rather than a selectivity claim
    (SUBUNITS_GTOPDB_MUST_NAME_TO_NARROW_A_COMPLEX)."""
    sign_of_gene: dict[str, float] = {}
    reason_of_gene: dict[str, str] = {}
    gtopdb_sign = {str(row.gene_symbol): float(row.sign) for row in gtopdb_rows.itertuples()}
    for chembl_row in chembl_rows.itertuples():
        genes = [gene for gene in str(getattr(chembl_row, gene_symbols_column)).split("|") if gene]
        chembl_sign = float(getattr(chembl_row, "sign", 0.0))
        if len(genes) <= 1:
            for gene in genes:
                sign_of_gene[gene], reason_of_gene[gene] = chembl_sign, "chembl single-gene row"
            continue
        named_by_gtopdb = [gene for gene in genes if gene in gtopdb_sign]
        if len(named_by_gtopdb) >= subunits_needed:
            for gene in named_by_gtopdb:
                sign_of_gene[gene], reason_of_gene[gene] = gtopdb_sign[gene], "gtopdb, narrowing a chembl complex"
            # the subunits the site cannot form without, which an affinity table does not name (see above)
            for gene in sorted(obligate_partners(sites or (), set(genes)) - set(named_by_gtopdb)):
                sign_of_gene[gene] = chembl_sign
                reason_of_gene[gene] = "obligate partner of the site, which the narrowing would drop"
        else:
            reason = ("chembl complex, no gtopdb row names a subunit of it" if not named_by_gtopdb else
                      f"chembl complex, gtopdb names only {len(named_by_gtopdb)} of its subunits")
            for gene in genes:
                sign_of_gene.setdefault(gene, chembl_sign)
                reason_of_gene.setdefault(gene, reason)
    for gene, sign in gtopdb_sign.items():  # a GtoPdb target belonging to no ChEMBL row of this drug
        if gene not in sign_of_gene:
            sign_of_gene[gene], reason_of_gene[gene] = sign, "gtopdb, outside every chembl row"
    return sign_of_gene, reason_of_gene


def pubchem_cid_from_perturbation_id(perturbation_id: str) -> int | None:
    """The PubChem CID of a STITCH flat identifier CID1xxxxxxxx, or None for an identifier of another kind."""
    match = re.fullmatch(r"CID1(\d{8})", str(perturbation_id))
    return int(match.group(1)) if match else None


# The two reasons replace_complex_rows gives a subunit it keeps from a complex row it narrows.
NARROWING_REASON = "gtopdb, narrowing a chembl complex"
OBLIGATE_PARTNER_REASON = "obligate partner of the site, which the narrowing would drop"
DEFAULT_BINDING_SITES_PATH = Path("configs/binding_sites.csv")


def binding_sites_by_drug_name(site_table_path: Path | None = DEFAULT_BINDING_SITES_PATH) -> dict[str, tuple[str, ...]]:
    """drug name, lowercased -> the binding sites configs/binding_sites.csv names for it, which decide the obligate
    subunits a narrowing restores (OBLIGATE_SUBUNITS_OF_A_SITE)."""
    if site_table_path is None or not Path(site_table_path).exists():
        return {}
    table = pd.read_csv(site_table_path)
    return {str(drug).strip().lower(): tuple(sorted(set(group.site.astype(str)))) for drug, group in table.groupby("drug")}


@dataclasses.dataclass
class ComplexTargetNarrowing:
    """A drug's ChEMBL complex targets narrowed to the subunits the Guide to PHARMACOLOGY names for that drug, with the
    binding site's obligate subunits restored, and nothing else changed (decision 8.11; the owner, 13 October 2026:
    "agreed on both 1 and 2", the unit a ChEMBL row after "cocaine should NOT narrow to a single target", and "ya
    absolutely add back obligate partners").

    The unit is one ChEMBL target, as replace_complex_rows has it: a single-gene target is kept whatever GtoPdb says, a
    complex target is narrowed only where GtoPdb names at least SUBUNITS_GTOPDB_MUST_NAME_TO_NARROW_A_COMPLEX of its
    subunits, and the target keeps its ChEMBL identifier and its one unit of magnitude, now spread over the subunits
    kept, each with the sign its own GtoPdb row states (an obligate subunit takes the ChEMBL target's sign). A GtoPdb
    target outside every ChEMBL row of the drug is not added: measured on 13 October the narrowing cost nothing in
    leakage and the additions took the largest leakage block from 114 drugs to 200 (docs/gtopdb_leakage_cost.md), and
    whether they enter is the owner's to decide."""

    ligand_of_pubchem_cid: dict[int, int]
    ligand_of_name: dict[str, int]
    rows_of_ligand: dict[int, pd.DataFrame]
    sites_of_drug_name: dict[str, tuple[str, ...]]
    provenance: dict

    @classmethod
    def from_gtopdb(cls, gtopdb_directory: Path = GTOPDB_DIRECTORY,
                    binding_sites_path: Path | None = DEFAULT_BINDING_SITES_PATH) -> "ComplexTargetNarrowing":
        gtopdb_directory = Path(gtopdb_directory)
        for name in (INTERACTIONS_FILE, LIGAND_MAPPING_FILE):
            if not (gtopdb_directory / name).is_file():
                raise FileNotFoundError(f"{gtopdb_directory / name} is missing: fetch the Guide to PHARMACOLOGY (scripts/fetch_raw_sources.sh), or build "
                                        "without the narrowing of complex targets (--no-narrow-complex-targets)")
        mapping = read_ligand_mapping(gtopdb_directory)
        cid_column = pd.to_numeric(mapping["PubChem CID"], errors="coerce")
        ligand_of_pubchem_cid = {int(cid): int(ligand) for cid, ligand in zip(cid_column, mapping["Ligand id"]) if pd.notna(cid)}
        ligand_of_name: dict[str, int] = {}
        for name, ligand in zip(mapping["Name"], mapping["Ligand id"]):
            key = str(name).strip().lower()
            if key:
                ligand_of_name.setdefault(key, int(ligand))
        interactions = read_interactions(gtopdb_directory)
        ligand_column = pd.to_numeric(interactions["Ligand ID"], errors="coerce")
        rows_of_ligand: dict[int, pd.DataFrame] = {}
        for ligand, rows in interactions.groupby(ligand_column):
            rows = rows.copy()
            rows["gene_symbol"] = rows["Target Gene Symbol"].astype(str).str.strip()
            rows["sign"] = [sign_of_action(action) for action in rows["Action"]]
            usable = one_gene_per_row(rows)
            if not usable.empty:
                rows_of_ligand[int(ligand)] = usable
        provenance = {"gtopdb_directory": str(gtopdb_directory),
                      **{f"{name}_sha256": hashlib.sha256((gtopdb_directory / name).read_bytes()).hexdigest() for name in (INTERACTIONS_FILE, LIGAND_MAPPING_FILE)},
                      "binding_sites": str(binding_sites_path) if binding_sites_path is not None else None,
                      "subunits_gtopdb_must_name": SUBUNITS_GTOPDB_MUST_NAME_TO_NARROW_A_COMPLEX, "gtopdb_additions": "not added"}
        return cls(ligand_of_pubchem_cid, ligand_of_name, rows_of_ligand, binding_sites_by_drug_name(binding_sites_path), provenance)

    def ligand_of(self, pubchem_cid: int | None, drug_name: str | None) -> tuple[int | None, str]:
        """The drug's GtoPdb ligand, by PubChem CID first and by its name second, and which matched."""
        if pubchem_cid is not None and int(pubchem_cid) in self.ligand_of_pubchem_cid:
            return self.ligand_of_pubchem_cid[int(pubchem_cid)], "pubchem_cid"
        name = str(drug_name or "").strip().lower()
        if name in self.ligand_of_name:
            return self.ligand_of_name[name], "name"
        return None, ""

    def narrow(self, drug_targets: list, pubchem_cid: int | None = None, drug_name: str | None = None) -> tuple[list, dict]:
        """(the drug's targets with each complex target narrowed where the rule allows, a record of what changed)."""
        ligand, matched_by = self.ligand_of(pubchem_cid, drug_name)
        record: dict = {"ligand": ligand, "matched_by": matched_by, "narrowed": {}, "kept_whole": {}}
        complex_targets = [target for target in drug_targets if len(target.gene_symbols) > 1]
        usable = self.rows_of_ligand.get(ligand) if ligand is not None else None
        if not complex_targets or usable is None:
            for target in complex_targets:
                record["kept_whole"][target.target_chembl_id] = "no GtoPdb ligand" if ligand is None else "no GtoPdb row names one gene with a sign"
            return list(drug_targets), record
        sites = self.sites_of_drug_name.get(str(drug_name or "").strip().lower(), ())
        narrowed_targets = []
        for target in drug_targets:
            if len(target.gene_symbols) <= 1:
                narrowed_targets.append(target)
                continue
            one_row = pd.DataFrame([{"gene_symbols": "|".join(target.gene_symbols), "sign": target.sign}])
            sign_of_gene, reason_of_gene = replace_complex_rows(one_row, usable, sites=sites)
            kept = [gene for gene in target.gene_symbols if reason_of_gene.get(gene) in (NARROWING_REASON, OBLIGATE_PARTNER_REASON)]
            if not kept:
                narrowed_targets.append(target)
                record["kept_whole"][target.target_chembl_id] = reason_of_gene.get(target.gene_symbols[0], "")
                continue
            narrowed_targets.append(dataclasses.replace(target, gene_symbols=kept, gene_signs={gene: float(sign_of_gene[gene]) for gene in kept}))
            record["narrowed"][target.target_chembl_id] = {"before": list(target.gene_symbols), "after": kept,
                                                           "obligate_partners": [gene for gene in kept if reason_of_gene[gene] == OBLIGATE_PARTNER_REASON],
                                                           "signs": {gene: float(sign_of_gene[gene]) for gene in kept}}
        return narrowed_targets, record
