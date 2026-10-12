"""Percent of each study drug bound to plasma proteins, as its FDA label states it, with the sentence and the label version.

Why. The unbound concentration of a drug needs its fraction unbound in plasma. The table of Lombardo et al. 2018
holds one for drugs with intravenous data only, which leaves out most oral psychotropic drugs (fluoxetine, sertraline,
olanzapine, quetiapine, lamotrigine). A drug label states the binding in its pharmacokinetics section, and
experiments/fetch_fda_label_protein_binding.py already caches those sentences with the label's set_id and
effective_time.

How. The readings below were made by hand, one drug at a time, from the cached sentences: a parser that takes the first
percentage of a sentence reads "bioavailability approximately 60%" as binding and reads the binding of penicillin G as
that of procaine. This script does not trust the hand either. For each reading it finds the cached sentence that holds
the quoted fragment and checks that the percentages of the reading are printed in that sentence; a reading that fails
is reported and left out. The output carries the full sentence, the label's set_id and its effective_time, so every
number can be checked against the label.

What was left out on purpose: a sentence about another ingredient of a combination product (naloxone under
buprenorphine, acetaminophen under hydrocodone, ethinyl estradiol under desogestrel), a sentence about a metabolite
only (the monohydroxy derivative of oxcarbazepine, the active metabolite of deflazacort), and a sentence that names
the binding protein without a percentage.

A reading of "at least" 99 percent bound gives only an upper limit for the fraction unbound. For a drug bound 99 percent
or more, a label's rounding moves the fraction unbound by a factor of two or more; the output flags those rows.

Usage:
  python experiments/read_label_protein_binding.py
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

LABEL_CACHE = Path("data/raw/fda_labels")
OUTPUT = Path("configs/drug_label_protein_binding.csv")
EVIDENCE_DIRECTORY = Path("data/processed/label_review/evidence_approved_with_list_no_cap")
HIGHLY_BOUND_FROM_PERCENT = 99.0

# drug as the study names it | percent bound, low | percent bound, high | "", "at least" or "at most" | fragment of the label sentence | note
READINGS = """
bupropion|84|84||bupropion is 84% bound to human plasma proteins|
alprazolam|80|80||Alprazolam is 80% bound to human serum protein|
amantadine|67|67||Amantadine is approximately 67% bound to plasma proteins|
amisulpride|25|30||Plasma protein binding is 25% to 30%|the sentence does not name the drug; the label has one active ingredient
bromocriptine|90|96||bromocriptine was 90%|
bupivacaine|95|95||Bupivacaine with a high protein binding capacity (95%)|
buspirone|86|86||approximately 86% of buspirone is bound to plasma proteins|
butorphanol|80|80||with a bound fraction of approximately 80%|
cabergoline|40|42||Protein binding of cabergoline was 40% to 42%|
carbamazepine|76|76||Carbamazepine is 76% bound to plasma proteins|
cevimeline|20|20|at most|is <20% bound to human plasma proteins|
clobazam|80|90||plasma protein binding of clobazam and N-desmethylclobazam is approximately 80% to 90% and 70%, respectively|the first range is clobazam's
clomipramine|97|97||The protein binding of CMI is approximately 97%|CMI is the label's abbreviation of clomipramine
clonazepam|85|85||Clonazepam is approximately 85% bound to plasma proteins|
clozapine|97|97||Clozapine is approximately 97% bound to serum proteins|
cocaine|84|92||Cocaine has been described in literature as approximately 84-92% bound to human plasma proteins|
codeine|7|25||Codeine is about 7 to 25% bound to plasma proteins|
diazepam|98|98||highly bound to plasma proteins (diazepam 98%)|
diflunisal|99|99|at least|More than 99% of diflunisal in plasma is bound to proteins|
valproate|81.5|90||the free fraction increases from approximately 10% at 40 mcg/mL to 18.5% at 130 mcg/mL|the label gives the free fraction, 10 to 18.5 percent, which rises with concentration; the bound percentages are 100 minus those
donepezil|96|96||Donepezil is approximately 96% bound to human plasma proteins|
doxepin|80|80||Doxepin is approximately 80% bound to plasma proteins|one other label gives 76%
felbamate|22|25||Binding ranged from 22% to 25%, mostly to albumin|the sentence does not name the drug; the label has one active ingredient
fluoxetine|94.5|94.5||approximately 94.5% of fluoxetine is bound in vitro to human serum proteins|
fluvoxamine|80|80||Approximately 80% of fluvoxamine is bound to plasma protein|
gabapentin|3|3|at most|Less than 3% of gabapentin circulates bound to plasma protein|
galantamine|18|18||The plasma protein binding of galantamine is 18%|
guanfacine|70|70||approximately 70% bound to plasma proteins independent of drug concentration|
hydromorphone|8|19||hydromorphone is approximately 8 to 19% bound to plasma proteins|
lamotrigine|55|55||lamotrigine is approximately 55% bound to human plasma proteins|
lorazepam|85|85||lorazepam is approximately 85% bound to plasma proteins|
memantine|45|45||the plasma protein binding is low (45%)|
mepivacaine|75|75||Mepivacaine is approximately 75% bound to plasma proteins|
methadone|85|90||methadone is predominantly bound to α 1 -acid glycoprotein (85% to 90%)|
methocarbamol|46|50||the plasma protein binding ranges between 46% and 50%|
methylphenidate|10|33||Plasma protein binding is 10% to 33%|the sentence does not name the drug; the label has one active ingredient
metoclopramide|30|30||Metoclopramide is not extensively bound to plasma proteins (about 30%)|
midazolam|97|97||midazolam is approximately 97% bound to plasma protein, principally albumin|
mifepristone|98|98||Mifepristone is 98% bound to plasma proteins|
mirtazapine|85|85||Mirtazapine is approximately 85% bound to plasma proteins|
modafinil|60|60||modafinil is moderately bound to plasma protein (approximately 60%)|
morphine|30|35||morphine is 30 to 35% reversibly bound to plasma proteins|one other label gives 20% to 35%
naltrexone|21|21||show naltrexone to be 21% bound to plasma proteins|
naratriptan|28|31||Plasma protein binding is 28% to 31%|the sentence does not name the drug; the label has one active ingredient
nefazodone|99|99|at least|nefazodone is extensively (> 99%) bound to human plasma proteins|
neostigmine|15|25||to human serum albumin ranges from 15 to 25%|
levonorgestrel|97.5|99||Levonorgestrel is about 97.5 to 99% protein-bound|
olanzapine|93|93||It is 93% bound to plasma proteins over the concentration range of 7 to 1100 ng/mL|the sentence does not name the drug; it stands in the distribution paragraph of the olanzapine label
oxycodone|45|45||Plasma protein binding of oxycodone at 37°C and a pH of 7.4 was about 45%|
paroxetine|93|95||Approximately 95% and 93% of paroxetine is bound to plasma protein at 100 ng/mL and 400 ng/mL, respectively|
pramipexole|15|15||It is about 15% bound to plasma proteins|the sentence does not name the drug; the label has one active ingredient
prilocaine|55|55||Prilocaine is 55% bound to plasma proteins|
quazepam|95|95|at least|plasma protein binding for quazepam and its two major metabolites is greater than 95%|
quetiapine|83|83||It is 83% bound to plasma proteins at therapeutic concentrations|the sentence does not name the drug; the label has one active ingredient
raloxifene|95|95||Raloxifene and the monoglucuronide conjugates are highly (95%) bound to plasma proteins|
riluzole|96|96||Plasma Protein Binding 96%|a row of the label's table of pharmacokinetic parameters
risperidone|90|90||The plasma protein binding of risperidone is approximately 90%|
rivastigmine|40|40||Rivastigmine is weakly bound to plasma proteins (approximately 40%)|
ropinirole|40|40|at most|It is up to 40% bound to plasma proteins|the sentence does not name the drug; the label has one active ingredient
sertraline|98|98||sertraline is highly bound to serum proteins (98%)|
sumatriptan|14|21||is low, approximately 14% to 21%|the sentence does not name the drug
temazepam|96|96||The unchanged drug was 96% bound to plasma proteins|
tiagabine|96|96||Tiagabine is 96% bound to human plasma proteins|
tizanidine|30|30||Tizanidine is approximately 30% bound to plasma proteins|
topiramate|15|41||Topiramate is 15% to 41% bound to human plasma proteins|
trazodone|89|95||Trazodone is 89 to 95% protein bound in vitro|
venlafaxine|27|27||Venlafaxine is 27% and ODV is 30% bound to plasma proteins|
zaleplon|60|60||The in vitro plasma protein binding is approximately 60%|
zolmitriptan|25|25||Plasma protein binding of zolmitriptan is 25%|
zolpidem|92.5|92.5||Total protein binding was found to be 92.5|
zonisamide|40|40||is approximately 40% bound to human plasma proteins|
sufentanil|93|93||was approximately 93% in healthy males|
atomoxetine|98|98||98% of atomoxetine in plasma is bound to protein|
milnacipran|13|13||Plasma protein binding is 13%|the sentence does not name the drug; the label has one active ingredient
fosphenytoin|95|99||Fosphenytoin is extensively bound (95% to 99%) to human plasma proteins|
levetiracetam|10|10|at most|Levetiracetam and its major metabolite are less than 10% bound to plasma proteins|
dexmedetomidine|94|94||The average protein binding was 94%|the sentence does not name the drug; the label has one active ingredient
aripiprazole|99|99|at least|aripiprazole and its major metabolite are greater than 99% bound to serum proteins|
remifentanil|70|70||Remifentanil is approximately 70% bound to plasma proteins|
duloxetine|90|90|at least|Duloxetine is highly bound (>90%) to proteins in human plasma|
ziprasidone|99|99|at least|It is greater than 99% bound to plasma proteins|the sentence does not name the drug; the label has one active ingredient
ropivacaine|94|94||Ropivacaine is 94% protein bound|
eletriptan|85|85||Plasma protein binding is moderate and approximately 85%|the sentence does not name the drug; the label has one active ingredient
paricalcitol|99.8|99.8|at least|Paricalcitol is extensively bound to plasma proteins (≥ 99.8%)|
paliperidone|74|74||The plasma protein binding of racemic paliperidone is 74%|
rasagiline|88|94||Plasma protein binding ranges from 88 to 94%|the sentence does not name the drug; the label has one active ingredient
desvenlafaxine|30|30||Plasma protein binding of desvenlafaxine is 30%|
rufinamide|34|34||Only a small fraction of rufinamide (34%) is bound to human serum proteins|
aprepitant|95|95|at least|Aprepitant is greater than 95% bound to plasma proteins|
asenapine|95|95||Asenapine is highly bound (95%) to plasma proteins|
cinacalcet|93|97||Cinacalcet is approximately 93% to 97% bound to plasma protein|
varenicline|20|20|at most|Plasma protein binding of varenicline is low (≤20%)|
ramelteon|82|82||In vitro protein binding of ramelteon is approximately 82% in human serum|
lacosamide|15|15|at most|Lacosamide is less than 15% bound to plasma proteins|
entacapone|98|98||the plasma protein binding of entacapone is 98%|
tolcapone|99.9|99.9|at least|The plasma protein binding of tolcapone is >99.9%|
tapentadol|20|20||Plasma protein binding of tapentadol is low (approximately 20%)|
ziconotide|50|50||Ziconotide is about 50% bound to human plasma proteins|
Nafarelin|80|80||About 80% of nafarelin acetate is bound to plasma proteins at 4°C|measured at 4°C
ulipristal|94|94|at least|Ulipristal acetate is highly bound (> 94%) to plasma proteins|
perampanel|95|96||perampanel is approximately 95 to 96% bound to plasma proteins|
tasimelteon|90|90||tasimelteon is about 90% bound to proteins|
brexpiprazole|99|99|at least|Brexpiprazole is highly protein bound in plasma (greater than 99%)|
brivaracetam|20|20|at most|Brivaracetam is weakly bound to plasma proteins (≤20%)|
cannabidiol|94|94|at least|Protein binding of the cannabidiol and its metabolites was >94% in vitro|
elagolix|80|80||% Bound to human plasma proteins 80|a row of the label's table of pharmacokinetic properties
pitolisant|91|96||Serum protein binding is approximately 91% to 96%|the sentence does not name the drug; the label has one active ingredient
istradefylline|98|98||The plasma protein binding of istradefylline was approximately 98%|
lasmiditan|55|60||The human plasma protein binding of lasmiditan is approximately 55% to 60%|
cenobamate|60|60||Plasma protein binding of cenobamate is 60%|
lemborexant|88|94||Plasma protein binding of lemborexant is approximately 88% in vitro and 94% in clinical samples|
Lumateperone|97.4|97.4||Protein binding of lumateperone is 97.4% at 5 µM|measured at about 70 times the therapeutic concentration, as the sentence says
osilodrostat|36.4|36.4||Protein binding is low (36.4%)|the sentence does not name the drug; the label has one active ingredient
opicapone|99|99|at least|Opicapone is highly bound to plasma proteins (>99%)|
remimazolam|91|91|at least|Plasma protein binding of BYFAVO was >91%|BYFAVO is the product name of remimazolam
setmelanotide|79.1|79.1||Protein binding of setmelanotide is 79.1%|
daridorexant|99.7|99.7||Daridorexant is 99.7% bound to plasma proteins|
lofexidine|55|55||Lofexidine protein binding is approximately 55%|
nalmefene|45|45||45% (CV 4.1%) is bound to plasma proteins|
ethinyl estradiol|97|97||Ethinyl estradiol is about 97% bound to serum albumin|
fenfluramine|50|50||Fenfluramine is 50% bound to human plasma proteins in vitro|
ketoconazole|99|99||the plasma protein binding is about 99% mainly to the albumin fraction|the sentence does not name the drug
testosterone|98|98||2% remains unbound (free)|the label gives the unbound percentage, 2; the bound percentage is 100 minus it
"""


def normalised(text: str) -> str:
    """Lower case with every run of white space removed, so that "90% -96%" and "90%-96%" compare equal."""
    return re.sub(r"\s+", "", str(text).lower()).replace("–", "-").replace("—", "-")


def parsed_readings(readings_text: str = READINGS) -> list[dict]:
    readings = []
    for line in readings_text.strip().splitlines():
        drug, low, high, qualifier, fragment, note = line.split("|")
        readings.append({"drug": drug, "percent_bound_low": float(low), "percent_bound_high": float(high), "qualifier": qualifier, "fragment": fragment, "note": note})
    return readings


def printed_form(number: float) -> str:
    return str(int(number)) if float(number).is_integer() else str(number)


def numbers_are_printed(reading: dict, sentence: str) -> bool:
    """True when the sentence prints the percentages of the reading, or the unbound percentages they were derived from."""
    numbers_in_sentence = set(re.findall(r"\d+(?:\.\d+)?", sentence))
    bound = {printed_form(reading["percent_bound_low"]), printed_form(reading["percent_bound_high"])}
    unbound = {printed_form(round(100.0 - reading["percent_bound_low"], 3)), printed_form(round(100.0 - reading["percent_bound_high"], 3))}
    return bound <= numbers_in_sentence or unbound <= numbers_in_sentence


def verified_reading(reading: dict, cached_label_file: dict) -> dict | None:
    """The reading with the newest label that holds its sentence, or None when no cached label holds it or its numbers."""
    wanted = normalised(reading["fragment"])
    holding = [(label.get("effective_time", ""), label.get("set_id", ""), sentence)
               for label in cached_label_file.get("labels", []) for sentence in label.get("sentences", []) if wanted in normalised(sentence)]
    holding = [entry for entry in holding if numbers_are_printed(reading, entry[2])]
    if not holding:
        return None
    effective_time, set_id, sentence = max(holding)
    low, high = reading["percent_bound_low"], reading["percent_bound_high"]
    # "at least X percent bound" leaves the fraction unbound at most 1 - X/100; "at most" leaves it at least that
    return {"percent_bound_low": low, "percent_bound_high": high, "qualifier": reading["qualifier"],
            "fraction_unbound": round(1.0 - (low + high) / 200.0, 5),
            "fraction_unbound_is": {"at least": "an upper limit", "at most": "a lower limit"}.get(reading["qualifier"], "a value"),
            "bound_99_percent_or_more": high >= HIGHLY_BOUND_FROM_PERCENT,
            "labels_with_the_sentence": len({entry[1] for entry in holding}), "label_set_id": set_id, "label_effective_time": effective_time,
            "label_sentence": sentence, "note": reading["note"]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=EVIDENCE_DIRECTORY, help="names the study drugs; the widest counting build by default")
    parser.add_argument("--cache-dir", type=Path, default=LABEL_CACHE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet", columns=["perturbation_id", "perturbation_type", "perturbation_label"])
    drugs = records[records.perturbation_type == "drug"].drop_duplicates("perturbation_id")
    perturbation_id_of = dict(zip(drugs.perturbation_label, drugs.perturbation_id))
    rows, not_verified = [], []
    for reading in parsed_readings():
        perturbation_id = perturbation_id_of.get(reading["drug"])
        cache_path = arguments.cache_dir / f"{perturbation_id}.json"
        if perturbation_id is None or not cache_path.exists():
            not_verified.append((reading["drug"], "not a drug of the evidence table" if perturbation_id is None else "no cached label"))
            continue
        verified = verified_reading(reading, json.loads(cache_path.read_text()))
        if verified is None:
            not_verified.append((reading["drug"], "no cached label holds the fragment with these percentages"))
            continue
        rows.append({"perturbation_id": perturbation_id, "drug": reading["drug"], **verified})
    table = pd.DataFrame(rows).sort_values("drug", key=lambda names: names.str.lower())
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(arguments.output, index=False)
    print(f"readings: {len(parsed_readings())}; verified against a cached label: {len(table)}; wrote {arguments.output}")
    print(f"  fraction unbound is a value for {int((table.fraction_unbound_is == 'a value').sum())}, an upper limit for {int((table.fraction_unbound_is == 'an upper limit').sum())}, "
          f"a lower limit for {int((table.fraction_unbound_is == 'a lower limit').sum())}; bound 99 percent or more: {int(table.bound_99_percent_or_more.sum())}")
    for drug, reason in not_verified:
        print(f"  not verified: {drug}: {reason}")
    return 1 if not_verified else 0


if __name__ == "__main__":
    raise SystemExit(main())
