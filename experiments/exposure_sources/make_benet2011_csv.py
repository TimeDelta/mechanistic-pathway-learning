"""Transcribe Tables I-V of Benet, Broccatelli and Oprea (2011) from the PubMed Central article page to one CSV.

Usage: python3 -I make_benet2011_csv.py ARTICLE.html OUTPUT.csv

Each table lists the drugs of one BDDCS class; the class is taken from the table caption. Cell text is kept as
printed, with runs of whitespace (line wraps of the page source) collapsed to one space. In numeric columns the
typographic minus sign U+2212 is written as the ASCII hyphen-minus; nothing else is changed and nothing is computed.
"""
import csv
import re
import sys

from bs4 import BeautifulSoup

EXPECTED_HEADER = ["generic name", "maximum strength dose value", "maximum strength dose unit", "formulation", "route",
                   "measured solubility (mg/ml)", "dose number", "% excreted unchanged in urine", "mw drug",
                   "measured logs molar", "measured logp", "measured logd74", "minvslgs 3–7.5", "clogp", "hba", "hbd", "psa"]
OUTPUT_HEADER = ["name", "maximum_strength_dose_value", "maximum_strength_dose_unit", "formulation", "route",
                 "measured_solubility_mg_per_ml", "dose_number", "percent_excreted_unchanged_urine", "mw_drug",
                 "measured_log_s_molar", "measured_log_p", "measured_log_d_7_4", "min_vs_lgs_3_7_5", "clogp", "hba",
                 "hbd", "psa"]
TEXT_COLUMNS = {"name", "maximum_strength_dose_unit", "formulation", "route"}
CAPTION_PATTERN = re.compile(r"Measured and In Silico Data for (\d+) BDDCS Class (\d) Drugs")
NUMBER_PATTERN = re.compile(r"-?\d+(\.\d+)?")


def cell_text(cell) -> str:
    return re.sub(r"\s+", " ", cell.get_text("", strip=False)).strip()


def main() -> int:
    article_path, output_path = sys.argv[1], sys.argv[2]
    soup = BeautifulSoup(open(article_path, encoding="utf-8").read(), "lxml")
    total = 0
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "bddcs_class"] + OUTPUT_HEADER[1:] + ["source_table", "source_row", "note"])
        for table in soup.find_all("table"):
            section = table.find_parent("section")
            label = cell_text(section.find("h4"))
            caption = cell_text(section.find(class_="caption"))
            match = CAPTION_PATTERN.fullmatch(caption)
            if not match:
                raise SystemExit(f"unexpected caption: {caption!r}")
            declared_count, bddcs_class = int(match.group(1)), match.group(2)
            rows = table.find_all("tr")
            header = [re.sub(r"\s+", "", cell_text(cell)).lower() for cell in rows[0].find_all(["th", "td"])]
            if header != [re.sub(r"\s+", "", name) for name in EXPECTED_HEADER]:
                raise SystemExit(f"unexpected header in {label}: {header}")
            data_rows = rows[1:]
            if len(data_rows) != declared_count:
                raise SystemExit(f"{label}: caption says {declared_count} drugs, table has {len(data_rows)} rows")
            for position, row in enumerate(data_rows, start=1):
                cells = row.find_all(["td", "th"])
                if len(cells) != len(OUTPUT_HEADER) or any(cell.get("colspan", "1") != "1" or cell.get("rowspan", "1") != "1" for cell in cells):
                    raise SystemExit(f"{label} row {position}: unexpected cell layout")
                values, notes = {}, []
                for column, cell in zip(OUTPUT_HEADER, cells):
                    text = cell_text(cell)
                    if re.search(r"-\s*\n", cell.get_text("", strip=False).strip()):
                        notes.append(f"{column}: the page source breaks the line after a hyphen; the space after that hyphen comes from the line break")
                    if cell.find(["sup", "sub"]):
                        notes.append(f"{column}: cell carries a superscript or subscript: {text!r}")
                    if column not in TEXT_COLUMNS and text != "":
                        text = text.replace("−", "-")
                        if not NUMBER_PATTERN.fullmatch(text):
                            notes.append(f"{column}: cell is not a plain number, kept as printed: {text!r}")
                    values[column] = text
                writer.writerow([values["name"], bddcs_class] + [values[column] for column in OUTPUT_HEADER[1:]] + [label.rstrip("."), position, "; ".join(notes)])
                total += 1
            print(label, caption, "->", len(data_rows), "rows")
    print("rows written:", total)
    return 0


if __name__ == "__main__":
    sys.exit(main())
