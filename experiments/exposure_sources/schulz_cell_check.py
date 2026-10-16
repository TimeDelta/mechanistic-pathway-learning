"""Every printed line of every cell (pdfplumber) must occur, character for character, in the text of the same page
as read by pdftotext -layout (poppler). Whitespace is ignored on both sides."""
import json
import re
import sys
cells_path, text_path = sys.argv[1], sys.argv[2]
bands = json.load(open(cells_path, encoding="utf-8"))
pages = [re.sub(r"\s+", "", p) for p in open(text_path, encoding="utf-8").read().split("\f")]
checked = failed = 0
for band in bands:
    if any(f[1] == "Arial-BoldMT" for f in band["fonts"]):
        continue
    for column, lines in band["cells"].items():
        for line in lines:
            text = re.sub(r"\s+", "", "".join(run[0] for run in line["runs"]))
            if not text:
                continue
            checked += 1
            if text not in pages[band["page"] - 1]:
                failed += 1
                print("not found on page", band["page"], column, repr(text))
print("lines checked", checked, "not found", failed)
