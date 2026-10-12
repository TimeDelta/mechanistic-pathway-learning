"""Extract the table cells of Schulz et al. 2020, Additional file 1 (PDF), with pdfplumber.

Rows are the bands between the horizontal rules that the PDF draws under each table row; columns are the six
printed columns, whose vertical rules sit at the same x positions on every page. Characters printed at 6.5 pt are
raised annotation numbers (superscripts) or, when they sit on the baseline, subscripts (as in "Vitamin B12"); they
are kept apart from the cell text.

Output: JSON list of bands, each with page, band position, fonts and, per column, the printed lines
(top position and runs of [text, kind], kind being "body", "sup" or "sub").
Usage: python3 -I schulz_cells.py <pdf> <first_page> <last_page> <out.json>
"""
import json
import statistics
import sys

import pdfplumber

COLUMN_NAMES = ["substance", "therapeutic", "toxic", "comatose_fatal", "half_life", "references"]
COLUMN_EDGES = [42.9, 233.95, 375.75, 460.75, 574.15, 680.45, 826.1]
SMALL_SIZE = 6.5


def horizontal_rule_positions(page):
    """y positions of the rules under table rows: thin rectangles lying across the substance column."""
    positions = set()
    for rect in page.rects:
        if (rect["bottom"] - rect["top"]) < 1.5 and rect["x0"] < 60 and rect["x1"] > 200:
            positions.add(round(rect["top"], 1))
    merged = []
    for position in sorted(positions):
        if merged and position - merged[-1] < 2.0:
            continue
        merged.append(position)
    return merged


def is_small(char):
    return round(char["size"], 1) == SMALL_SIZE


def lines_of_cell(chars):
    """Printed lines of one cell, top to bottom."""
    body_chars = [c for c in chars if not is_small(c)]
    small_chars = [c for c in chars if is_small(c)]
    line_tops = []
    for char in sorted(body_chars, key=lambda c: c["top"]):
        if not line_tops or char["top"] - line_tops[-1] > 4.0:
            line_tops.append(char["top"])
    if not line_tops and small_chars:
        line_tops = [min(c["top"] for c in small_chars) + 0.7]
    lines = [[] for _ in line_tops]
    for char in chars:
        reference_top = char["top"] + (0.7 if is_small(char) else 0.0)
        index = min(range(len(line_tops)), key=lambda i: abs(line_tops[i] - reference_top))
        lines[index].append(char)
    result = []
    for line_top, line in zip(line_tops, lines):
        line.sort(key=lambda c: c["x0"])
        body_bottoms = [c["bottom"] for c in line if not is_small(c)]
        body_bottom = statistics.median(body_bottoms) if body_bottoms else None
        runs = []
        previous = None
        for char in line:
            if not is_small(char):
                kind = "body"
            elif body_bottom is None or char["bottom"] - body_bottom < -2.0:
                kind = "sup"
            else:
                kind = "sub"
            text = char["text"]
            if previous is not None:
                gap = char["x0"] - previous["x1"]
                if gap > 1.5 and text != " " and not (runs and runs[-1][0].endswith(" ")):
                    runs.append([" ", "body"])  # a gap with no space character printed
            if runs and runs[-1][1] == kind:
                runs[-1][0] += text
            else:
                runs.append([text, kind])
            previous = char
        result.append({"top": round(line_top, 1), "runs": runs})
    return result


def main():
    pdf_path, first_page, last_page, out_path = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    bands = []
    characters_in_table_area = 0
    characters_in_cells = 0
    with pdfplumber.open(pdf_path) as pdf:
        for page_number in range(first_page, last_page + 1):
            page = pdf.pages[page_number - 1]
            rules = horizontal_rule_positions(page)
            if len(rules) < 2:
                print("page", page_number, "has no table rules", file=sys.stderr)
                continue
            in_area = [c for c in page.chars if rules[0] <= (c["top"] + c["bottom"]) / 2 < rules[-1] and c["text"].strip()]
            characters_in_table_area += len(in_area)
            for band_index in range(len(rules) - 1):
                top, bottom = rules[band_index], rules[band_index + 1]
                band_chars = [c for c in page.chars if top <= (c["top"] + c["bottom"]) / 2 < bottom
                              and COLUMN_EDGES[0] <= (c["x0"] + c["x1"]) / 2 < COLUMN_EDGES[-1]]
                cells = {}
                fonts = set()
                for column_index, column_name in enumerate(COLUMN_NAMES):
                    left, right = COLUMN_EDGES[column_index], COLUMN_EDGES[column_index + 1]
                    cell_chars = [c for c in band_chars if left <= (c["x0"] + c["x1"]) / 2 < right]
                    characters_in_cells += sum(1 for c in cell_chars if c["text"].strip())
                    for c in cell_chars:
                        fonts.add((round(c["size"], 1), c["fontname"].split("+")[-1]))
                    cells[column_name] = lines_of_cell(cell_chars)
                bands.append({"page": page_number, "band_top": top, "band_bottom": bottom,
                              "fonts": sorted(fonts), "cells": cells})
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(bands, handle, ensure_ascii=False, indent=0)
    print(len(bands), "bands;", characters_in_table_area, "non-blank characters between the first and last rule;",
          characters_in_cells, "assigned to cells")


if __name__ == "__main__":
    main()
