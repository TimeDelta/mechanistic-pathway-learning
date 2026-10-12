"""Build schulz2020.csv and schulz2020_annotations.csv from the extracted cells and annotations.

Usage: python3 -I schulz_build_csv.py <cells.json> <annotations.json> <out.csv> <out_annotations.csv> <out_summary.json>

Nothing is converted. A number is written to a numeric column only when the printed cell, after its raised
annotation numbers are taken out, has one of the plain forms listed in PARSE_FORMS; the form is named in the
<column>_parse column, and every other cell keeps its text only.
"""
import csv
import json
import re
import sys

SYMBOL_FONT_CHARACTERS = {0xF0B1: "±", 0xF0B3: "≥", 0xF0A3: "≤", 0xF0AE: "→", 0xF06D: "µ", 0xF061: "α"}
VALUE_COLUMNS = ["therapeutic", "toxic", "comatose_fatal"]
# Rows whose name cell prints two entries, one per line, with the values of each entry on the same printed line
# (seen on the rendered pages 5, 6, 20, 29 and 42). Each printed line becomes its own CSV row.
ROWS_WITH_ONE_ENTRY_PER_LINE = {"Amitriptyline", "Aripiprazole", "Clomipramine", "Doxepin", "Imipramine"}

NUMBER = r"\d+(?:\.\d+)?"
PARSE_FORMS = [
    # (name, regular expression, which groups are low and high)
    ("range", rf"^({NUMBER})\s?[-–]\s?({NUMBER})$", (1, 2)),
    ("single", rf"^({NUMBER})$", (1, 1)),
    ("up_to", rf"^-({NUMBER})$", (None, 1)),
    ("less_than", rf"^<\s?({NUMBER})$", (None, 1)),
    ("greater_than", rf"^>\s?({NUMBER})$", (1, None)),
    ("at_least", rf"^≥\s?({NUMBER})$", (1, None)),
    ("approx_range", rf"^appr\.\s?({NUMBER})-({NUMBER})$", (1, 2)),
    ("approx_single", rf"^appr\.\s?({NUMBER})$", (1, 1)),
    ("approx_up_to", rf"^appr\.\s?-({NUMBER})$", (None, 1)),
    ("range_with_bracketed_extension", rf"^(?:\({NUMBER}-\)\s?)?({NUMBER})-({NUMBER})(?:\s?\(-{NUMBER}\))?$", (1, 2)),
    ("single_with_bracketed_extension", rf"^(?:\({NUMBER}-\)\s?)?({NUMBER})(?:\s?\(-{NUMBER}\))?$", (1, 1)),
]
UNIT_PATTERN = re.compile(r"(?<![A-Za-zµμ/])(ng ?/mL|ng/L|ng/g|ng/h/mL|µg/L|μg/L|µg/mL|μg/mL|µmol/L|μmol/L|mmol/L|nmol/L|pmol/L|pg/mL|mg/kg|mg/g|µg/g|μg/g|g/L|IU/L|IU/mL|U/L|U/mL|mU/L|%|‰|ppm|ppb)")
HALF_LIFE_UNIT_PATTERN = re.compile(r"(days?|weeks?|months?|years?|minutes?|\bmin\b|\bh\b)")


def clean(text):
    for private_use_code, character in SYMBOL_FONT_CHARACTERS.items():
        text = text.replace(chr(private_use_code), character)
    return text


def render_lines(lines):
    """(printed, text, annotation numbers, features) of a list of printed lines of one cell.

    printed: raised annotation numbers written as ^{...}; text: the same without them.
    A line that ends in a hyphen with no space after it is joined to the next line without a space.
    """
    printed, text = "", ""
    annotations = []
    features = set()
    kept_lines = [line for line in lines if "".join(run[0] for run in line["runs"]).strip()]
    if len(kept_lines) > 1:
        features.add(f"printed on {len(kept_lines)} lines")
    for line_index, line in enumerate(kept_lines):
        line_printed, line_text = "", ""
        for run_text, kind in line["runs"]:
            run_text = clean(run_text)
            if kind == "sup" and run_text.strip():
                line_printed += "^{" + run_text.strip() + "}" + (" " if run_text.endswith(" ") else "")
                annotations += [int(number) for number in re.findall(r"\d+", run_text)]
                leftover = re.sub(r"[\d, ]", "", run_text)
                if leftover:
                    line_text += leftover + (" " if run_text.endswith(" ") else "")
                    features.add(f"'{leftover}' printed raised with the annotation number")
            else:
                if kind == "sub" and run_text.strip():
                    features.add(f"'{run_text.strip()}' printed as subscript")
                line_printed += run_text
                line_text += run_text
        if line_index > 0:
            previous = printed
            if previous.endswith("-"):
                features.add("line break after a hyphen (joined without a space)")
            else:
                printed, text = printed.rstrip() + " ", text.rstrip() + " "
        printed += line_printed
        text += line_text
    printed = re.sub(r"\s+", " ", printed).strip()
    text = re.sub(r"\s+", " ", text).strip()
    if printed and not text:
        features.add("the cell prints only the annotation number(s), no value")
    return printed, text, annotations, features


def parse_value(text, unit_pattern=UNIT_PATTERN, printed=""):
    """(low, high, form) for a cell text; numbers stay as printed strings."""
    if not text:
        return "", "", ("annotation_only" if printed else "")
    if text == "-":
        return "", "", "dash"
    if re.match(r"^[Ss]ee\b", text):
        return "", "", "cross_reference"
    for name, pattern, (low_group, high_group) in PARSE_FORMS:
        match = re.match(pattern, text)
        if match:
            low = match.group(low_group) if low_group else ""
            high = match.group(high_group) if high_group else ""
            return low, high, name
    if unit_pattern.search(text):
        return "", "", "not_parsed_other_unit"
    if ";" in text:
        return "", "", "not_parsed_several_values"
    return "", "", "not_parsed"


def main():
    cells_path, annotations_path, out_csv, out_annotations_csv, out_summary = sys.argv[1:6]
    bands = json.load(open(cells_path, encoding="utf-8"))
    annotation_data = json.load(open(annotations_path, encoding="utf-8"))
    annotation_text = {a["number"]: a["text"] for a in annotation_data["annotations"]}

    with open(out_annotations_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["annotation", "page", "text"])
        for a in annotation_data["annotations"]:
            writer.writerow([a["number"], a["page"], a["text"]])

    data_bands = [b for b in bands if not any(font[1] == "Arial-BoldMT" for font in b["fonts"])]
    header_bands = len(bands) - len(data_bands)
    rows = []
    unknown_annotations = set()
    for source_row, band in enumerate(data_bands, start=1):
        name_lines = [line for line in band["cells"]["substance"] if "".join(r[0] for r in line["runs"]).strip()]
        first_name = re.sub(r"\s+", " ", "".join(r[0] for r in name_lines[0]["runs"] if r[1] != "sup")).strip() if name_lines else ""
        split = first_name in ROWS_WITH_ONE_ENTRY_PER_LINE and len(name_lines) == 2
        if split:
            groups = []
            for line_index, name_line in enumerate(name_lines):
                group = {"substance": [name_line]}
                for column in VALUE_COLUMNS + ["half_life"]:
                    group[column] = [line for line in band["cells"][column]
                                     if min(range(len(name_lines)), key=lambda i: abs(name_lines[i]["top"] - line["top"])) == line_index]
                    if len([line for line in group[column] if "".join(r[0] for r in line["runs"]).strip()]) > 1:
                        raise SystemExit(f"row {source_row}: more than one printed line of {column} beside one name line")
                # the references cell is one wrapped list for the whole printed row: kept whole, on line 1
                group["references"] = band["cells"]["references"] if line_index == 0 else []
                groups.append(group)
        else:
            groups = [band["cells"]]
        for line_index, group in enumerate(groups, start=1):
            row = {"source_row": source_row, "source_line": line_index, "page": band["page"]}
            notes = []
            all_annotations = []
            name_printed, name_text, name_annotations, name_features = render_lines(group["substance"])
            row["name"] = name_text
            row["name_printed"] = name_printed
            row["name_annotations"] = "; ".join(str(n) for n in name_annotations)
            all_annotations += name_annotations
            notes += [f"name: {feature}" for feature in sorted(name_features)]
            for column in VALUE_COLUMNS:
                printed, text, annotations, features = render_lines(group[column])
                low, high, form = parse_value(text, printed=printed)
                if low and high and float(high) < float(low):
                    notes.append(f"{column}: the printed range has its second number below its first; kept as printed")
                row[f"{column}_printed"] = printed
                row[f"{column}_text"] = text
                row[f"{column}_annotations"] = "; ".join(str(n) for n in annotations)
                row[f"{column}_low_mg_per_l"] = low
                row[f"{column}_high_mg_per_l"] = high
                row[f"{column}_parse"] = form
                all_annotations += annotations
                notes += [f"{column}: {feature}" for feature in sorted(features)]
                units = sorted(set(UNIT_PATTERN.findall(text)))
                if units:
                    notes.append(f"{column}: the cell states its own unit ({', '.join(units)}); text kept, nothing converted")
                annotation_units = sorted(set(unit for n in annotations for unit in UNIT_PATTERN.findall(annotation_text.get(n, ""))))
                if annotation_units:
                    notes.append(f"{column}: annotation(s) {', '.join(str(n) for n in annotations)} mention other units ({', '.join(annotation_units)}); see remarks")
            printed, text, annotations, features = render_lines(group["half_life"])
            half_life_units = sorted(set(HALF_LIFE_UNIT_PATTERN.findall(text)))
            if half_life_units:
                low, high, form = "", "", "not_parsed_other_unit"
                notes.append(f"half_life: the cell states its own unit ({', '.join(half_life_units)}); text kept, not converted to hours")
            else:
                low, high, form = parse_value(text, unit_pattern=re.compile(r"(?!x)x"), printed=printed)
            if low and high and float(high) < float(low):
                notes.append("half_life: the printed range has its second number below its first; kept as printed")
            row["half_life_printed"] = printed
            row["half_life_text"] = text
            row["half_life_annotations"] = "; ".join(str(n) for n in annotations)
            row["half_life_low_h"] = low
            row["half_life_high_h"] = high
            row["half_life_parse"] = form
            all_annotations += annotations
            notes += [f"half_life: {feature}" for feature in sorted(features)]
            printed, text, annotations, features = render_lines(group["references"])
            row["references_printed"] = text
            cross_reference = re.match(r"^[Ss]ee (?:also )?(.+)$", row["therapeutic_text"])
            row["see"] = cross_reference.group(1) if cross_reference else ""
            seen = []
            for number in all_annotations:
                if number not in seen:
                    seen.append(number)
                if number not in annotation_text:
                    unknown_annotations.add(number)
            row["remarks"] = " || ".join(f"[{n}] {annotation_text.get(n, 'ANNOTATION NOT FOUND')}" for n in seen)
            if split:
                notes.append(f"printed row {source_row} holds two entries, one per printed line; this is line {line_index}"
                             + ("" if line_index == 1 else "; cells printed once for the row (on its first line) are in the row of line 1"))
            other_fonts = [font for font in band["fonts"] if font[1] == "SymbolMT"]
            if other_fonts:
                notes.append("a sign in this row is set in the Symbol font (code 0xB3) and written here as ≥")
            row["note"] = "; ".join(notes)
            rows.append(row)

    columns = list(rows[0].keys())
    with open(out_csv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "bands": len(bands), "header_bands": header_bands, "printed_rows": len(data_bands), "csv_rows": len(rows),
        "non_missing": {column: sum(1 for row in rows if str(row[column]) != "") for column in columns},
        "parse_forms": {column: {} for column in VALUE_COLUMNS + ["half_life"]},
        "annotation_numbers_not_found": sorted(unknown_annotations),
    }
    for column in VALUE_COLUMNS + ["half_life"]:
        for row in rows:
            form = row[f"{column}_parse"] or "(empty cell)"
            summary["parse_forms"][column][form] = summary["parse_forms"][column].get(form, 0) + 1
    with open(out_summary, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=1, ensure_ascii=False)
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
