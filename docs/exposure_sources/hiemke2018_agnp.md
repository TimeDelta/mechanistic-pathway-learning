# Hiemke et al. 2018 (AGNP consensus guidelines, update 2017): NOT retrieved

> Copied on 11 October 2026 from `data/raw/exposure/hiemke2018_agnp/sources.md`, the note written beside the raw files (which are not in git). Paths such as `raw/...` are relative to that directory. The transcription scripts are in `experiments/exposure_sources/`.

Tried on 12 October 2026 (UTC; the evening of 11 October 2026 local time). **No table was transcribed and no CSV
was written.** `hiemke2018_therapeutic_reference_ranges.csv` and `hiemke2018_dose_related_reference_ranges.csv` do
not exist. Licence: **not verified** (no open licence was found on any route that could be read).

## Citation

Hiemke C, Bergemann N, Clement HW, Conca A, Deckert J, Domschke K, et al. Consensus Guidelines for Therapeutic Drug
Monitoring in Neuropsychopharmacology: Update 2017. Pharmacopsychiatry 2018;51(1-02):9-62.
doi:10.1055/s-0043-116492. PMID 28910830. No PMCID. Erratum: Pharmacopsychiatry 2018;51(1-02):e1,
doi:10.1055/s-0037-1600991, PMID 29390205.

## What is open

Only the abstract (PubMed, Europe PMC). The article is not in PubMed Central and Europe PMC holds no full text.
The publisher's site (Thieme) may let a person read the PDF (one repository record calls the article "bronze"
open access, which means free to read at the publisher without an open licence), but the site's robots file turns
this reader away, so the PDF was not requested.

## Routes tried, in order

| # | Address | How | Result |
|---|---|---|---|
| 1 | `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:"10.1055/s-0043-116492"&resultType=core&format=json` | curl (Europe PMC service) | 200. `isOpenAccess: N`, `inEPMC: N`, `inPMC: N`, `hasPDF: N`, `license: null`, no `pmcid`; the only full-text link is the DOI, "Subscription required". Saved `raw/europepmc_search_hiemke2018.json`. |
| 2 | `https://thieme-connect.de/products/all/doi/10.1055/s-0043-116492` (a search result) | WebFetch | Read before the robots file had been looked at. It returned no article text, no links and no licence, only "© Georg Thieme Verlag KG Stuttgart · New York". Nothing from it was used. |
| 3 | `https://www.thieme-connect.com/products/ejournals/abstract/10.1055/s-0043-116492` | shell script `experiments/read_public_page.py` | **Refused before the request**: "the site's robots file disallows this address for automated readers". |
| 4 | `https://www.thieme-connect.com/products/ejournals/pdf/10.1055/s-0043-116492.pdf`, `https://www.thieme-connect.com/products/ejournals/html/10.1055/s-0043-116492`, `https://www.thieme-connect.de/products/ejournals/pdf/10.1055/s-0043-116492.pdf` | shell script | **Refused before the request**, same reason. The robots file (`raw/thieme_connect_robots.txt`, read with curl) holds `User-agent: ClaudeBot` / `Disallow: /` (also GPTBot, ChatGPT-User, bingbot and others); its `User-agent: *` section does not disallow `/products/ejournals/`. The script checks the rules for ClaudeBot, Claude-User and anthropic-ai as well as `*`, so the site counts as refusing. No further Thieme address was requested by any tool. |
| 5 | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=28910830&retmode=xml` | curl (NCBI E-utilities) | 200. `<CopyrightInformation>© Georg Thieme Verlag KG Stuttgart · New York.</CopyrightInformation>`; article ids are PubMed and DOI only. Saved `raw/pubmed_28910830.xml`. |
| 6 | `https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/?ids=10.1055/s-0043-116492&format=json` | curl (PMC id converter, no e-mail parameter) | 200. "Identifier not found in PMC". Saved `raw/pmc_idconv_hiemke2018.json`. |
| 7 | PubMed connector, `get_copyright_status` for PMID 28910830 | connector | copyright "© Georg Thieme Verlag KG Stuttgart · New York.", `license.type: null`, `is_open_access: false`, source "pubmed". |
| 8 | `https://api.crossref.org/works/10.1055/s-0043-116492` | curl (Crossref, no `mailto`) | 200. No `license` field. One link, to the Thieme PDF, "intended-application: similarity-checking". Saved `raw/crossref_work_hiemke2018.json`. |
| 9 | `https://api.openalex.org/works/doi:10.1055/s-0043-116492` | curl (OpenAlex, no e-mail) | 200. `oa_status: green`; the publisher location is `is_oa: false`, licence none. The only location marked open is the Göttingen repository (`https://resolver.sub.uni-goettingen.de/purl?gro-2/74377`), which OpenAlex labels `license: cc-by`, `version: submittedVersion`. That label was **not verified** (route 10). Saved `raw/openalex_work_hiemke2018.json`. |
| 10 | `https://publications.goettingen-research-online.de/handle/2/74377` and `https://publications.goettingen-research-online.de/bitstream/2/74377/1/document.pdf` | shell script | **Refused by the site: HTTP 403** (3,688 bytes of HTML) for both, although its robots file was served and does not disallow them. Not retried with another tool or user agent. |
| 11 | `https://sonar.ch/global/documents/263484` (a search result) | WebFetch | A record only: its full-text link points to the Thieme PDF; "Open access status: bronze"; no licence stated. |
| 12 | `https://iris.unil.ch/handle/iris/98447`, `https://serval.unil.ch/notice/serval:BIB_4965E3E471A7` (search results) | WebFetch | Records only ("Open Access: Oui/Yes"), no file, no licence. |
| 13 | `https://boris-portal.unibe.ch/entities/publication/7ce8d357-2072-453a-a89a-996b5c34a303/full` (a search result) | WebFetch | Empty page (the site is a script application). |
| 14 | `https://boris-portal.unibe.ch/server/api/core/items/7ce8d357-2072-453a-a89a-996b5c34a303` | shell script | **Refused by the site: a check that the visitor is human** ("Making sure you're not a bot! [...] Protected by Anubis"). The script did not recognise the page as a check and saved it; the saved page was deleted and nothing from it was used. |
| 15 | `https://pred.uni-regensburg.de/id/eprint/15284/` | shell script | 200. "Full text not available from this repository." Saved in `raw/regensburg/`. |
| 16 | `https://epub.ub.uni-muenchen.de/63429/` (a search result) | WebFetch | "Full text not available from 'Open Access LMU'." |
| 17 | `https://freidok.uni-freiburg.de/data/175646` | shell script | 200, a record page built by script; no file link and no licence in the page source. Saved in `raw/other_repositories/`. |
| 18 | `https://mediatum.ub.tum.de/1589862` | shell script | 200, a record only (full text = the DOI). Saved in `raw/other_repositories/`. |
| 19 | `https://publications.rwth-aachen.de/record/757772` | shell script | **Refused by the site: a script challenge page** (248 bytes, `/fast-challenge/`). The script saved it; the saved page was deleted. |
| 20 | `https://hdl.handle.net/21.11116/0000-0002-7002-D` (Max Planck repository) | shell script | **HTTP 403** (94 bytes). |
| 21 | `https://www.ppt-online.de/_Resources/Persistent/6/c/c/5/6cc58517f68e550140086cc560738af382b10035/ppt201803_artikel1410.pdf` (a search result) | WebFetch (summary only; the file was not saved) | This is the German translation: Hefner G, et al. Konsensus-Leitlinien für therapeutisches Drug-Monitoring in der Neuropsychopharmakologie: Update 2017. Psychopharmakotherapie 2018;25(3):92-140 ("Diese Arbeit ist eine leicht gekürzte Übersetzung/adaptierte Fassung von Hiemke C. et al."). It is freely downloadable but carries "Kein Nachdruck, keine Veröffentlichung im Internet oder Intranet ohne Zustimmung des Verlags!" and "© Wissenschaftliche Verlagsgesellschaft Stuttgart"; no open licence. **Not used**: it is not the article asked for and it has no licence. |
| 22 | `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=(TITLE:"therapeutic drug monitoring" AND (TITLE:neuropsychopharmacology OR TITLE:AGNP))&resultType=lite&format=json&pageSize=50` | curl (Europe PMC service) | 200, 16 records. None of the AGNP guideline papers (2004, 2011, 2017 update, its erratum, the 2018 and 2019 summaries) is open access or in PMC. It also lists "Consensus Guidelines for Therapeutic Drug Monitoring in Neuropsychopharmacology: Update 2026", Pharmacopsychiatry, doi:10.1055/a-2860-7861, PMID 42392224, `isOpenAccess: N`. Saved `raw/europepmc_search_agnp_guidelines.json`. |

No address waited on a permission request. Unpaywall was not asked, because it requires an e-mail address in the
request.

## Refusals (respected: no cache, mirror, solved challenge or other user agent was used)

- Thieme (`thieme-connect.com`, `thieme-connect.de`): robots rule against ClaudeBot.
- Göttingen research repository (`publications.goettingen-research-online.de`): HTTP 403.
- Bern repository (`boris-portal.unibe.ch`): check that the visitor is human (Anubis).
- RWTH Aachen publications (`publications.rwth-aachen.de`): script challenge.
- `hdl.handle.net` (Max Planck record): HTTP 403.

## Licence: not verified

Routes tried for a licence statement: the Europe PMC record (none), the PubMed record (copyright line only, quoted
above), the PubMed connector (none), Crossref (none), OpenAlex (none for the publisher's copy; "cc-by" for a
"submittedVersion" in the Göttingen repository, which refused the request), the repository records above (none).
The article itself could not be opened, so whatever licence line it prints was not read.

Addresses for the owner to open:

- `https://www.thieme-connect.com/products/ejournals/pdf/10.1055/s-0043-116492.pdf` (the article; look for the
  licence line on its first page) and `https://doi.org/10.1055/s-0043-116492`.
- `https://resolver.sub.uni-goettingen.de/purl?gro-2/74377` (the repository copy that OpenAlex labels CC BY).
- The German translation at route 21, should the owner decide that it may be used.

If the owner saves the article PDF into `raw/`, the two tables can be transcribed from it in the same way as the
Schulz table (`experiments/exposure_sources/`).

## What the retrieved Schulz table already holds from this guideline

`../schulz2020/schulz2020.csv`: 143 rows cite reference [4] (this guideline) and 169 toxic cells carry annotation
311, "“laboratory alert level”; mostly according to AGNP Consensus Guidelines for therapeutic drug monitoring in
psychiatry: update 2017 [4] [...]". That is Schulz's use of the guideline, in mg/L, and not a transcription of the
guideline's tables; it holds no level of recommendation, no conversion factor and nothing of the dose-related
table (CL/F, F, dosing interval, tmax, DRC factors).

## Text addressed to an AI assistant

None seen in the saved pages or in the fetch summaries.

## SHA-256 of the raw files

| File | SHA-256 |
|---|---|
| `raw/europepmc_search_hiemke2018.json` | `be5d055a1f810fcb43943d819abf0e8a5e398f0f4b2e2c661f08f14df7db2677` |
| `raw/europepmc_search_agnp_guidelines.json` | `ba4f39a868c71296317fb419aa470fb297ac935b19cae6c06587e67747098302` |
| `raw/pubmed_28910830.xml` | `125f76766418cd8471b7dd83cd07909761c7a03f2bec1e36c1d8769b56d3a712` |
| `raw/pmc_idconv_hiemke2018.json` | `97dc4107b0d267ad193480565a166917d72ae72c3a7bf0b888dd38135d202317` |
| `raw/crossref_work_hiemke2018.json` | `40f80d86d908931e2c5430f0f80a1956436e4123ad89a3a6ef2fa4e17266d34b` |
| `raw/openalex_work_hiemke2018.json` | `c6ed145367259b1577bd9ea5ab159c2a89576ff60d95381cb7267689df38026b` |
| `raw/thieme_connect_robots.txt` | `975e3b5171fb02dbdb0268dd2e27c5e9ea3828ac3800c3eeebf35fff47d0b039` |
| `raw/regensburg/2160dad1c6b761ab.html` | `2160dad1c6b761ab191ee46ebe7cbe97caf2bb3e81e12c1bf865fa24ef4121ba` |
| `raw/other_repositories/bdb99f5950aa5d43.html` (Freiburg) | `bdb99f5950aa5d432de966f12663397f394d5a66c15a036767b8062d5fa302d5` |
| `raw/other_repositories/8d2be4c63bc4f34d.html` (TU Munich) | `8d2be4c63bc4f34d45d34dafad7d10cc59b2ea71eb2d662b18953c56aff7e5ac` |

Each page saved by the shell script has a `.json` beside it with its address, time and SHA-256.

## A defect seen in `experiments/read_public_page.py` (not changed)

Its list of challenge markers misses two pages met here: the Anubis page ("Making sure you're not a bot") and the
RWTH page (`/fast-challenge/`). Both came back with HTTP 200 and were saved as if they were content.
