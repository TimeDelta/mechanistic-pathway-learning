"""Transcribe the blood-to-plasma tables held in the US EPA httk repository to CSV.

Usage: python3 -I make_blood_to_plasma_csv.py UCHIMURA.xlsx PAIXAO.xlsx HTTK_CHEM_PROPS.txt OUTPUT_DIRECTORY

Nothing is computed. A numeric cell is written with the shortest decimal text that reads back to the stored number;
a text cell is written as it stands; a spreadsheet error marker such as '#N/A' is left empty and quoted in `note`.
"""
import csv
import pathlib
import sys

import openpyxl


def cell_text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        raise ValueError("boolean cell")
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return str(int(value)) if value == int(value) and abs(value) < 1e15 else repr(value)
    return str(value)


def transcribe_sheet(sheet, expected_header, output_header, output_path, name_column) -> int:
    rows = list(sheet.iter_rows(values_only=True))
    if list(rows[0][: len(expected_header)]) != expected_header:
        raise SystemExit(f"unexpected header: {rows[0]}")
    written = 0
    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(output_header + ["source_row", "note"])
        for row_number, row in enumerate(rows[1:], start=2):
            values = list(row[: len(expected_header)])
            if all(value is None or value == "" for value in values):
                continue
            output, notes = [], []
            for source_column, column, value in zip(expected_header, output_header, values):
                text = cell_text(value)
                if isinstance(value, str) and value.startswith("#"):
                    notes.append(f"{source_column!r} cell holds the spreadsheet marker {value!r}, left empty here")
                    text = ""
                elif isinstance(value, str) and value != value.strip():
                    notes.append(f"{source_column!r} cell has outer whitespace in the source ({value!r}), stripped here")
                    text = value.strip()
                output.append(text)
            name_position = output_header.index(name_column)
            writer.writerow([output[name_position]] + output[:name_position] + output[name_position + 1:] + [row_number, "; ".join(notes)])
            written += 1
    return written


def main() -> int:
    uchimura_path, paixao_path, httk_path, output_directory = sys.argv[1], sys.argv[2], sys.argv[3], pathlib.Path(sys.argv[4])

    sheet = openpyxl.load_workbook(uchimura_path).worksheets[0]
    expected = ["cas", "Name", "Polarity", "logD", "Human.fup", "Human.Rb2p", "Rat.fup", "Rat.Rb2p"]
    output = ["cas_number", "name", "polarity", "log_d", "fraction_unbound_plasma_human", "blood_to_plasma_ratio_human",
              "fraction_unbound_plasma_rat", "blood_to_plasma_ratio_rat"]
    # the header row is written with `name` first
    path = output_directory / "uchimura2010.csv"
    count = transcribe_sheet(sheet, expected, output, path, "name")
    fix_header(path, ["name"] + [column for column in output if column != "name"])
    print("uchimura2010 rows:", count)

    workbook = openpyxl.load_workbook(paixao_path)
    sheet = workbook["Table 2"]
    expected = ["CAS", "Drug", "Foral", "frenal", "fp", "Clplasma (ml.min−1.kg−1)", "Rb", "Drug Class", "ClH", "Fabs", "Fmet", "Reference"]
    output = ["cas_number", "name", "f_oral", "f_renal", "fp", "cl_plasma_ml_per_min_per_kg", "blood_to_plasma_ratio_rb",
              "drug_class", "cl_h", "f_abs", "f_met", "reference"]
    path = output_directory / "paixao2012_table2.csv"
    count = transcribe_sheet(sheet, expected, output, path, "name")
    fix_header(path, ["name"] + [column for column in output if column != "name"])
    print("paixao2012 Table 2 rows:", count)
    print("paixao2012 reference key:", [[cell for cell in row if cell] for row in workbook["References"].iter_rows(values_only=True)])

    # The file is tab-separated without quoting, one record per line, lines ending in CR LF. Three records
    # (Pcb136, Pcb153, Pcb155) are broken over two lines inside the All.Compound.Names field; such a record is
    # put together again from its pieces, and `source_row` is the number of the line on which the record starts.
    raw_lines = open(httk_path, encoding="utf-8", errors="replace", newline="").read().split("\n")
    header = raw_lines[0].rstrip("\r").split("\t")
    records, pending, pending_start = [], None, None
    for line_number, line in enumerate(raw_lines[1:], start=2):
        piece = line.rstrip("\r")
        if pending is None and piece == "":
            continue
        if pending is None:
            pending, pending_start = piece, line_number
        else:
            pending += piece
        if pending.count("\t") == len(header) - 1:
            records.append((pending_start, line_number, dict(zip(header, pending.split("\t")))))
            pending = None
        elif pending.count("\t") > len(header) - 1:
            raise SystemExit(f"line {line_number}: too many fields")
    if pending is not None:
        raise SystemExit("incomplete last record")
    rows = records
    written = 0
    with open(output_directory / "httk_human_rblood2plasma.csv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["name", "cas_number", "dtxsid", "blood_to_plasma_ratio_human", "reference", "source_row", "note"])
        for first_line, last_line, row in rows:
            value = row["Human.Rblood2plasma"]
            if value in ("", "NA"):
                continue
            reference = row["Human.Rblood2plasma.Reference"]
            notes = []
            if reference in ("", "NA"):
                notes.append("no reference given in the source row")
                reference = ""
            if last_line != first_line:
                notes.append(f"record broken over lines {first_line}-{last_line} of the source file, put together here")
            writer.writerow([row["Compound"], row["CAS"], "" if row["DTXSID"] == "NA" else row["DTXSID"], value, reference, first_line, "; ".join(notes)])
            written += 1
    print("records:", len(rows), "broken records put together:", sum(1 for first, last, _ in rows if first != last))
    print("httk rows with Human.Rblood2plasma:", written, "of", len(rows))
    return 0


def fix_header(path, header) -> None:
    lines = pathlib.Path(path).read_text(encoding="utf-8").split("\n", 1)
    pathlib.Path(path).write_text(",".join(header + ["source_row", "note"]) + "\n" + lines[1], encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
