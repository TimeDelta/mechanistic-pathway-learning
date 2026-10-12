"""Check the transcribed Lombardo 2018 CSV against the raw workbook XML and against independent copies.

Usage: python3 -I verify_lombardo2018.py OUTPUT.csv EPA.xlsx EPA.csv DREADLESSS.xlsx DDUPRED.csv.gz
"""
import csv
import gzip
import re
import sys
import zipfile

import openpyxl

output_csv, epa_xlsx, epa_csv, dreadlesss_xlsx, ddupred_csv = sys.argv[1:6]
with open(output_csv, encoding="utf-8", newline="") as handle:
    ours = list(csv.DictReader(handle))
print("our rows:", len(ours))
PK = ["vdss_l_per_kg", "clearance_ml_per_min_per_kg", "fraction_unbound_plasma", "mrt_h", "half_life_h"]

# 1. raw XML numbers of the EPA workbook: columns D..H, rows 10..1361
with zipfile.ZipFile(epa_xlsx) as archive:
    sheet_xml = archive.read("xl/worksheets/sheet1.xml").decode("utf-8")
cells = dict()
for reference, attributes, inner in re.findall(r'<c r="([A-Z]+\d+)"([^>]*?)(?:/>|>(.*?)</c>)', sheet_xml, flags=re.S):
    value = re.search(r"<v>(.*?)</v>", inner or "")
    cells[reference] = (attributes, value.group(1) if value else None)
different_text = 0
mismatch = 0
examples = []
for row in ours:
    source_row = row["source_row"]
    for letter, column in zip("DEFGH", PK):
        attributes, raw = cells.get(f"{letter}{source_row}", ("", None))
        if 't="s"' in attributes or raw is None:
            if row[column] != "":
                mismatch += 1
            continue
        if float(raw) != float(row[column]):
            mismatch += 1
        if raw != row[column]:
            different_text += 1
            if len(examples) < 8:
                examples.append((row["name"], column, raw, row[column]))
print("numeric PK cells whose value differs from the raw XML number:", mismatch)
print("numeric PK cells whose text differs from the raw XML text (same number):", different_text, examples)

# 2. dreadlesss copy, row by row in order
sheet = openpyxl.load_workbook(dreadlesss_xlsx).worksheets[0]
theirs = [[cell for cell in row] for row in sheet.iter_rows(min_row=2, values_only=True)]
print("dreadlesss rows:", len(theirs))
by_name = {}
differences = []
for row in theirs:
    by_name.setdefault((str(row[0]).strip(), str(row[1]).strip()), []).append(row)
matched = 0
for row in ours:
    key = (row["name"], row["cas_number"])
    candidates = by_name.get(key)
    if not candidates:
        differences.append(("missing in dreadlesss", key))
        continue
    matched += 1
    other = candidates[0]
    for index, column in zip((3, 4, 5, 6, 7), PK):
        other_value = other[index]
        other_text = "" if other_value in (None, "") else other_value
        our_text = row[column]
        same = (other_text == "" and our_text == "") or (other_text != "" and our_text != "" and float(other_text) == float(our_text))
        if not same:
            differences.append((key, column, our_text, other_value))
print("matched on (name, CAS) in dreadlesss copy:", matched, "differences:", len(differences), differences[:20])

# 3. EPA csv in the same repository
with open(epa_csv, encoding="utf-8", errors="replace", newline="") as handle:
    epa_rows = list(csv.reader(handle))
header_index = next(index for index, row in enumerate(epa_rows) if row and row[0] == "Name")
epa_header = epa_rows[header_index]
epa_data = [row for row in epa_rows[header_index + 1:] if any(cell.strip() for cell in row)]
print("EPA csv header:", epa_header[:9], "rows:", len(epa_data))
epa_by_cas = {row[1].strip(): row for row in epa_data}
differences = []
for row in ours:
    other = epa_by_cas.get(row["cas_number"])
    if other is None:
        differences.append(("missing in EPA csv", row["name"], row["cas_number"]))
        continue
    for index, column in zip((3, 4, 5, 6, 7), PK):
        other_text, our_text = other[index].strip(), row[column]
        same = (other_text == "" and our_text == "") or (other_text != "" and our_text != "" and float(other_text) == float(our_text))
        if not same:
            differences.append((row["name"], column, our_text, other_text))
print("EPA csv differences:", len(differences), differences[:20])

# 4. DDUpred copy (subset)
with gzip.open(ddupred_csv, "rt", encoding="utf-8", errors="replace", newline="") as handle:
    ddu = list(csv.DictReader(handle))
ours_by_cas = {row["cas_number"]: row for row in ours}
checked = 0
differences = []
for other in ddu:
    row = ours_by_cas.get(other["CAS #"].strip())
    if row is None:
        differences.append(("not in ours", other["Name"], other["CAS #"]))
        continue
    checked += 1
    for source_column, column in (("human VDss (L/kg)", "vdss_l_per_kg"), ("human CL (mL/min/kg)", "clearance_ml_per_min_per_kg"),
                                  ("UnknownPropertyName1[2]", "fraction_unbound_plasma"), ("MRT (h)", "mrt_h"), ("terminal  t1/2 (h)", "half_life_h")):
        other_text, our_text = other[source_column].strip(), row[column]
        same = (other_text == "" and our_text == "") or (other_text != "" and our_text != "" and float(other_text) == float(our_text))
        if not same:
            differences.append((row["name"], column, our_text, other_text))
print("DDUpred rows:", len(ddu), "matched on CAS:", checked, "differences:", len(differences), differences[:40])
