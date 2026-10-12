"""Print spot checks: values read straight from the raw files next to the values in the CSV files.

Usage: python3 -I spot_checks.py EXPOSURE_DIRECTORY
"""
import csv
import pathlib
import re
import sys
import zipfile

from bs4 import BeautifulSoup

base = pathlib.Path(sys.argv[1])


def read_csv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


# Lombardo 2018: read the sheet XML and shared strings directly, without openpyxl
workbook = base / "lombardo2018/raw/usepa_httk/datatables_Lombardo2018-Supplemental_82966_revised_corrected.xlsx"
with zipfile.ZipFile(workbook) as archive:
    shared = re.findall(r"<si>(.*?)</si>", archive.read("xl/sharedStrings.xml").decode("utf-8"), flags=re.S)
    shared = ["".join(re.findall(r"<t[^>]*>(.*?)</t>", item, flags=re.S)) for item in shared]
    sheet_xml = archive.read("xl/worksheets/sheet1.xml").decode("utf-8")
raw_rows = {}
for row_number, row_xml in re.findall(r'<row r="(\d+)"[^>]*>(.*?)</row>', sheet_xml, flags=re.S):
    cells = {}
    for reference, attributes, inner in re.findall(r'<c r="([A-Z]+)\d+"([^>]*?)(?:/>|>(.*?)</c>)', row_xml, flags=re.S):
        value = re.search(r"<v>(.*?)</v>", inner or "")
        if value is None:
            cells[reference] = ""
        elif 't="s"' in attributes:
            cells[reference] = shared[int(value.group(1))]
        else:
            cells[reference] = value.group(1)
    raw_rows[int(row_number)] = cells
ours = {row["name"]: row for row in read_csv(base / "lombardo2018/lombardo2018.csv")}
print("LOMBARDO 2018 (raw sheet XML text | CSV)")
for name in ["Amiodarone", "Diazepam", "Warfarin", "Haloperidol", "Morphine", "Caffeine", "Midazolam", "Propranolol", "Olanzapine", "Hydroxychloroquine", "Alfentanil", "Artesunate", "Lithium", "Citalopram", "Venlafaxine"]:
    if name not in ours:
        print("   not in the table under this name:", name, [n for n in ours if n.lower().startswith(name[:6].lower())][:5])
        continue
    row = ours[name]
    raw = raw_rows[int(row["source_row"])]
    print(f"| {raw['A']!r} | {raw['B']} | VDss {raw['D']} / {row['vdss_l_per_kg']} | CL {raw['E']} / {row['clearance_ml_per_min_per_kg']} | fu {raw['F']} / {row['fraction_unbound_plasma']} | MRT {raw['G']} / {row['mrt_h']} | t1/2 {raw['H']} / {row['half_life_h']} | row {row['source_row']} |")


def html_rows(path):
    soup = BeautifulSoup(open(path, encoding="utf-8").read(), "lxml")
    found = {}
    for table in soup.find_all("table"):
        label = re.sub(r"\s+", " ", table.find_parent("section").find("h4").get_text("", strip=True))
        for row in table.find_all("tr")[1:]:
            cells = [re.sub(r"\s+", " ", cell.get_text("", strip=True)) for cell in row.find_all(["td", "th"])]
            found[(label, cells[0])] = cells
    return found


raw = html_rows(base / "benet2011_bddcs/raw/2f21927cccf25526.html")
ours = {row["name"]: row for row in read_csv(base / "benet2011_bddcs/benet2011_bddcs.csv")}
print("BENET 2011 (raw table cells | CSV)")
for name in ["Amiodarone hydrochloride", "Diazepam", "Warfarin", "Haloperidol", "Fluoxetine", "Lithium carbonate", "Metformin", "Gabapentin", "Acetazolamide", "Amphetamine sulfate", "Atenolol", "Penicillin G; benzylpenicillin"]:
    if name not in ours:
        print("   not found:", name, [n for n in ours if n.lower().startswith(name.split()[0].lower())][:5])
        continue
    row = ours[name]
    cells = raw[(row["source_table"] + ".", name)]
    print(f"| {cells[0]} | {row['source_table']} row {row['source_row']} | class {row['bddcs_class']} | % excreted unchanged in urine: raw {cells[7]!r} / CSV {row['percent_excreted_unchanged_urine']!r} | dose {cells[1]!r} {cells[2]!r} / {row['maximum_strength_dose_value']!r} {row['maximum_strength_dose_unit']!r} | CLogP raw {cells[13]!r} / CSV {row['clogp']!r} |")

raw = html_rows(base / "hosey2016_bddcs/raw/6e9eae3c00ba1e65.html")
ours = {row["name"]: row for row in read_csv(base / "hosey2016_bddcs/hosey2016_bddcs.csv")}
print("HOSEY 2016 Table III (raw | CSV)")
for name in list(ours)[::18]:
    cells = raw[("Table III.", name)]
    print(f"| {cells[0]} | raw class {cells[1]!r} | CSV class {ours[name]['bddcs_class']!r} | row {ours[name]['source_row']} |")
