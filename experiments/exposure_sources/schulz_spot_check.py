"""Print CSV rows next to the pdftotext -layout lines of the same substance (an independent text extraction)."""
import csv
import re
import sys
csv_path, text_path = sys.argv[1], sys.argv[2]
names = sys.argv[3:]
rows = list(csv.DictReader(open(csv_path, encoding="utf-8")))
pages = open(text_path, encoding="utf-8").read().split("\f")
for name in names:
    matches = [r for r in rows if r["name"].lower().startswith(name.lower())]
    for r in matches:
        print(f"CSV  p{r['page']} name={r['name_printed']!r} ther={r['therapeutic_printed']!r} [{r['therapeutic_low_mg_per_l']}|{r['therapeutic_high_mg_per_l']}|{r['therapeutic_parse']}] tox={r['toxic_printed']!r} [{r['toxic_low_mg_per_l']}|{r['toxic_high_mg_per_l']}] coma={r['comatose_fatal_printed']!r} [{r['comatose_fatal_low_mg_per_l']}|{r['comatose_fatal_high_mg_per_l']}] t½={r['half_life_printed']!r} [{r['half_life_low_h']}|{r['half_life_high_h']}|{r['half_life_parse']}] refs={r['references_printed']!r}")
        page_text = pages[int(r["page"]) - 1].split("\n")
        first_word = r["name"].split(" (")[0][:25]
        for i, line in enumerate(page_text):
            if line.strip().lower().startswith(first_word.lower()):
                print("PDF  " + re.sub(r" {2,}", "  ", line.strip()))
                if i + 1 < len(page_text) and page_text[i + 1].strip():
                    print("PDF+ " + re.sub(r" {2,}", "  ", page_text[i + 1].strip()))
        if r["note"]: print("NOTE", r["note"])
    print()
