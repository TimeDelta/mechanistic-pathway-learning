"""Extract the numbered annotations (pages 85-118) of Schulz et al. 2020, Additional file 1 (PDF).

Each annotation starts with its number printed left of x = 114 pt; its text starts at about x = 117 pt and
continuation lines keep that indent. Characters at 8 pt are raised (written ^{...}) or lowered (written inline,
as in Cmax); a "½" printed slightly below the line belongs to the "t" before it.

Usage: python3 -I schulz_annotations.py <pdf> <first_page> <last_page> <out.json>
"""
import json
import re
import statistics
import sys

import pdfplumber

NUMBER_COLUMN_RIGHT = 114.0
# Glyphs set in the Symbol font come out of the PDF as private-use code points (0xF000 + Symbol code).
SYMBOL_FONT_CHARACTERS = {0xF0B1: "±", 0xF0B3: "≥", 0xF0A3: "≤", 0xF0AE: "→", 0xF06D: "µ", 0xF061: "α"}
BODY_SIZE = 12.0


def lines_of_page(page):
    chars = [c for c in page.chars if 60 < c["top"] < 770]
    tops = sorted(set(round(c["top"], 1) for c in chars if round(c["size"], 1) == BODY_SIZE and c["text"].strip() and c["text"] != "½"))
    line_tops = []
    for top in tops:
        if not line_tops or top - line_tops[-1] > 7.0:
            line_tops.append(top)
    lines = [[] for _ in line_tops]
    for char in chars:
        if not line_tops:
            continue
        index = min(range(len(line_tops)), key=lambda i: abs(line_tops[i] - char["top"]))
        if abs(line_tops[index] - char["top"]) > 7.0:
            continue
        lines[index].append(char)
    result = []
    for line_top, line in zip(line_tops, lines):
        line.sort(key=lambda c: c["x0"])
        body_bottoms = [c["bottom"] for c in line if round(c["size"], 1) == BODY_SIZE and c["text"] != "½"]
        body_bottom = statistics.median(body_bottoms)
        text = ""
        in_superscript = False
        for char in line:
            small = round(char["size"], 1) < BODY_SIZE
            superscript = small and (char["bottom"] - body_bottom) < -2.0
            if superscript and not in_superscript:
                text += "^{"
            if not superscript and in_superscript:
                text += "}"
            in_superscript = superscript
            text += char["text"]
        if in_superscript:
            text += "}"
        fonts = sorted(set(c["fontname"].split("+")[-1] for c in line if c["text"].strip()))
        result.append({"top": line_top, "x0": min(c["x0"] for c in line if c["text"].strip()) if any(c["text"].strip() for c in line) else None,
                       "text": text, "fonts": fonts,
                       "number_part": "".join(c["text"] for c in line if c["x0"] < NUMBER_COLUMN_RIGHT and round(c["size"], 1) == BODY_SIZE),
                       "text_part_start": next((i for i, c in enumerate(line) if c["x0"] >= NUMBER_COLUMN_RIGHT), None)})
    return result


def main():
    pdf_path, first_page, last_page, out_path = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    annotations = []
    trailing_text = []
    with pdfplumber.open(pdf_path) as pdf:
        finished = False
        for page_number in range(first_page, last_page + 1):
            for line in lines_of_page(pdf.pages[page_number - 1]):
                if line["x0"] is None:
                    continue
                text = re.sub(r"\s+", " ", line["text"]).strip()
                if text in ("Annotations",):
                    continue
                if text.startswith("Abbreviations:"):
                    finished = True
                if finished:
                    trailing_text.append(text)
                    continue
                number_match = re.match(r"^\s*(\d+)\s+(.*)$", line["text"]) if line["x0"] < NUMBER_COLUMN_RIGHT - 2 and line["number_part"].strip().isdigit() else None
                if number_match:
                    annotations.append({"number": int(number_match.group(1)), "page": page_number,
                                        "lines": [number_match.group(2)]})
                else:
                    if not annotations:
                        print("text before the first annotation:", repr(text), file=sys.stderr)
                        continue
                    annotations[-1]["lines"].append(line["text"])
    for annotation in annotations:
        joined = ""
        for line in annotation["lines"]:
            if joined and joined[-1] in "-–":  # printed line ended at a hyphen or dash with no space after it
                joined += line.lstrip()
            else:
                joined += (" " if joined else "") + line.lstrip()
        for private_use_code, character in SYMBOL_FONT_CHARACTERS.items():
            joined = joined.replace(chr(private_use_code), character)
        annotation["text"] = re.sub(r"\s+", " ", joined).strip()
        del annotation["lines"]
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump({"annotations": annotations, "trailing_text": " ".join(trailing_text)}, handle, ensure_ascii=False, indent=1)
    numbers = [a["number"] for a in annotations]
    print(len(annotations), "annotations; first", numbers[0], "last", numbers[-1],
          "consecutive:", numbers == list(range(1, len(numbers) + 1)))
    private_use = sorted(set(ch for a in annotations for ch in a["text"] if "" <= ch <= ""))
    print("private-use characters left:", [hex(ord(ch)) for ch in private_use])


if __name__ == "__main__":
    main()
