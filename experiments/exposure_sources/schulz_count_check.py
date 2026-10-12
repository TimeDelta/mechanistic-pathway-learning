"""Compare the name lines found by pdfplumber (cells JSON) with the lines that pdftotext -layout starts in the
first column of the table pages: two independent readings of the same PDF."""
import json, re, sys
cells_path, text_path = sys.argv[1], sys.argv[2]
bands = json.load(open(cells_path, encoding="utf-8"))
pages = open(text_path, encoding="utf-8").read().split("\f")
total_plumber = total_poppler = 0
differences = []
for page_number in range(1, 83):
    plumber_lines = []
    for band in bands:
        if band["page"] != page_number or any(f[1] == "Arial-BoldMT" for f in band["fonts"]):
            continue
        for line in band["cells"]["substance"]:
            text = "".join(run[0] for run in line["runs"])
            if text.strip():
                plumber_lines.append(re.sub(r"\s+", "", text))
    poppler_lines = []
    in_table = False
    for line in pages[page_number - 1].split("\n"):
        if line.startswith(" Substance"):
            in_table = True
            continue
        if not in_table or not line.strip():
            continue
        if line.startswith("Drug concentrations are mentioned") or line.startswith("Note on liability"):
            break
        if re.match(r"^ \S", line):
            name = re.split(r" {2,}", line.strip())[0]
            poppler_lines.append(re.sub(r"\s+", "", name))
    total_plumber += len(plumber_lines)
    total_poppler += len(poppler_lines)
    if plumber_lines != poppler_lines:
        import difflib
        for d in difflib.unified_diff(poppler_lines, plumber_lines, lineterm="", n=0):
            if not d.startswith(("---", "+++", "@@")):
                differences.append((page_number, d))
print("name lines: pdfplumber", total_plumber, "pdftotext", total_poppler)
for page_number, d in differences:
    print(page_number, d)
