"""Rows and non-missing values per column of a CSV file given as argument."""
import csv
import sys

with open(sys.argv[1], encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle))
print("rows:", len(rows))
for column in rows[0]:
    non_missing = sum(1 for row in rows if row[column] != "")
    distinct = len({row[column] for row in rows if row[column] != ""})
    print(f"| `{column}` | {non_missing} | {distinct} |")
