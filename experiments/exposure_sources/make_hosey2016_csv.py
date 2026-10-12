"""Transcribe Tables I-III of Hosey, Chan and Benet (2016) from the PubMed Central article page to three CSV files.

Usage: python3 -I make_hosey2016_csv.py ARTICLE.html OUTPUT_DIRECTORY

Cell text is kept as printed, with runs of whitespace collapsed to one space. Nothing is computed.
"""
import csv
import pathlib
import re
import sys

from bs4 import BeautifulSoup


def cell_text(cell) -> str:
    return re.sub(r"\s+", " ", cell.get_text("", strip=False)).strip()


def main() -> int:
    article_path, output_directory = sys.argv[1], pathlib.Path(sys.argv[2])
    soup = BeautifulSoup(open(article_path, encoding="utf-8").read(), "lxml")
    tables = {cell_text(table.find_parent("section").find("h4")): table for table in soup.find_all("table")}
    captions = {label: cell_text(table.find_parent("section").find(class_="caption")) for label, table in tables.items()}
    if captions != {"Table I.": "BDDCS Class Changes from Initial Publication", "Table II.": "Parameter Changes from Initial Publication",
                    "Table III.": "Newly BDDCS Classified Drugs"}:
        raise SystemExit(f"unexpected tables: {captions}")

    # Table III: newly classified drugs
    rows = tables["Table III."].find_all("tr")
    if [cell_text(cell) for cell in rows[0].find_all(["th", "td"])] != ["Drug", "Class"]:
        raise SystemExit("unexpected header in Table III")
    with open(output_directory / "hosey2016_bddcs.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "bddcs_class", "source_table", "source_row", "note"])
        for position, row in enumerate(rows[1:], start=1):
            cells = row.find_all(["td", "th"])
            if len(cells) != 2:
                raise SystemExit(f"Table III row {position}: unexpected layout")
            notes = [f"cell carries a superscript or subscript: {cell_text(cell)!r}" for cell in cells if cell.find(["sup", "sub"])]
            writer.writerow([cell_text(cells[0]), cell_text(cells[1]), "Table III", position, "; ".join(notes)])
    print("Table III rows:", len(rows) - 1)

    # Table I: class changes
    table = tables["Table I."]
    rows = table.find_all("tr")
    header = [cell_text(cell) for cell in rows[0].find_all(["th", "td"])]
    if header != ["Compound", "Listed classa", "Updated class", "Major route of elimination"]:
        raise SystemExit(f"unexpected header in Table I: {header}")
    footnote = cell_text(table.find_parent("section").find(class_="tw-foot"))
    with open(output_directory / "hosey2016_table1_class_changes.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "listed_bddcs_class", "updated_bddcs_class", "major_route_of_elimination", "source_table", "source_row", "note"])
        for position, row in enumerate(rows[1:], start=1):
            cells = [cell_text(cell) for cell in row.find_all(["td", "th"])]
            if len(cells) != 4:
                raise SystemExit(f"Table I row {position}: unexpected layout")
            writer.writerow(cells + ["Table I", position, f"column header 'Listed class' carries footnote marker a; footnote as printed: {footnote!r}"])
    print("Table I rows:", len(rows) - 1, "| footnote:", footnote)

    # Table II: parameter changes, with sub-headings that span the two value columns
    rows = tables["Table II."].find_all("tr")
    header = [cell_text(cell) for cell in rows[0].find_all(["th", "td"])]
    if header != ["Compound", "Parameter value listed", "Updated parameter value"]:
        raise SystemExit(f"unexpected header in Table II: {header}")
    written = 0
    with open(output_directory / "hosey2016_table2_parameter_changes.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "parameter", "parameter_value_listed", "updated_parameter_value", "source_table", "source_row", "note"])
        parameter = None
        for position, row in enumerate(rows[1:], start=1):
            cells = row.find_all(["td", "th"])
            if len(cells) == 2 and cells[1].get("colspan") == "2" and cell_text(cells[0]) == "":
                parameter = cell_text(cells[1])
                continue
            if len(cells) != 3 or parameter is None:
                raise SystemExit(f"Table II row {position}: unexpected layout")
            writer.writerow([cell_text(cells[0]), parameter, cell_text(cells[1]), cell_text(cells[2]), "Table II", position, ""])
            written += 1
    print("Table II rows:", written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
