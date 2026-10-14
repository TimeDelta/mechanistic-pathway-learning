# The four held drugs of item 14: a search outside the sources already tried

The owner, 13 October 2026: "Can you look for the missing pieces of those four drugs elsewhere?" Recorded as decision
8.9 and not run until that evening, which was the lapse. This is the pass.

Item 14 admits a drug from outside the nervous-system code when all three hold: (a) its own label states a central
action in a pharmacology section, (b) a mechanism target is a brain-expressed graph node, (c) its measured
Kp,uu,brain is at least 0.1. All four drugs satisfy (b).

| drug | missing | outcome of this pass |
| --- | --- | --- |
| **carvedilol** | (a) and (c) | **resolved, against admission, on positive evidence** |
| **timolol** | (c) | no Kp,uu,brain exists in the accessible literature |
| **oxprenolol** | (a) | refused by the one source that holds a label; a route is left for the owner |
| **pindolol** | (a) | no label other than the US one found |

## Carvedilol: the evidence says exclude, which is a result rather than a gap

No Kp,uu,brain was found, and what was found answers the question it was asked for. Carvedilol is a
**P-glycoprotein substrate at the blood-brain barrier**, and so firmly that it has been developed as a tracer *for
imaging P-gp function there*. Bart et al. 2005, *Br J Pharmacol* 145(8):1045-51
([10.1038/sj.bjp.0706283](https://doi.org/10.1038/sj.bjp.0706283), PMID 15951832, PMC1576233), read from its
abstract through the PubMed MCP service on 13 October 2026:

> "The beta-receptor antagonist carvedilol is a P-gp substrate with a log P=2.0, and can be labeled with [11C]."

> "Ex vivo [11C]carvedilol biodistribution studies showed that [11C]carvedilol uptake in the brain was increased by
> CsA."

> "Brain DV of [11C]carvedilol showed a dose-dependent increase of maximal three-fold after CsA pretreatment."

A brain uptake that trebles when P-gp is inhibited is efflux-limited, and efflux limitation is what puts a drug's
Kp,uu,brain below 1. That is the same kind of positive evidence item 14 already accepts in the other direction: the
item excludes nadolol, labetalol, sotalol, montelukast and etoposide on label sentences saying the drug does not
reach the brain, and this is the pharmacological equivalent. **Proposed: carvedilol is excluded on this evidence
rather than held for a missing number.** It is not a measured Kp,uu,brain and is not registered as one; what is
registered is the ground for excluding it.

## Timolol: the number does not exist, and this pass can say so with some confidence

- **All 34 Kp,uu tables in `data/raw/label_review/effect_site_tables/kpuu/` were searched. Timolol appears in none
  of them**, nor does carvedilol. Those tables were copied complete with respect to ChEMBL for every unbound-ratio
  assay they cover.
- PubMed, through the MCP service: "timolol brain unbound" returns **nothing**; "timolol brain distribution"
  returns four records, none of which measures brain exposure. Two use timolol at 20 micromolar as a
  beta-antagonist applied *directly* to hippocampal and hypothalamic slices
  ([10.1523/JNEUROSCI.13-12-05393.1993](https://doi.org/10.1523/JNEUROSCI.13-12-05393.1993),
  [10.1016/s0306-4522(00)00003-8](https://doi.org/10.1016/s0306-4522(00)00003-8)), which shows it acts on central
  beta-receptors when put there and says nothing about whether it gets there; one is retinal
  ([10.1016/s0006-8993(02)03372-3](https://doi.org/10.1016/s0006-8993(02)03372-3)); one measures organ blood flow
  rather than drug concentration ([10.1159/000158180](https://doi.org/10.1159/000158180)).
- Consensus, semantic search over the Kp,uu literature, returned the methodology and model-building papers rather
  than values. Every dataset among them that the project could use is already in the directory above: Fridén 2009,
  Sato 2021, Kosugi 2021, Umemori 2022, Vallianatou 2022, Orozco 2019, Cui 2019, Gülave 2025.

**Proposed: timolol stays out, and the reason is recorded as an absence rather than as evidence against it.** The
distinction matters: carvedilol has evidence of efflux, timolol has no measurement either way.

## Oxprenolol and pindolol: a refusal, and a route left open

Both need part (a), a label sentence of their own. Neither was found.

- **Web search** surfaced no summary of product characteristics for either. The pindolol results are US
  professional labels (`drugs.com/pro/pindolol.html` and NDC listings), which are the same source openFDA already
  gives, so they carry the same class boilerplate and no new sentence.
- **The UK electronic medicines compendium refused both requests.** `medicines.org.uk/emc/search?q=oxprenolol` and
  the same for pindolol each answered HTTP 200 with a check that the visitor is human, read with
  `experiments/read_public_page.py` on 14 October 2026. Under the ruling of 11 October a human check is a refusal
  and is not routed around through a cache, a mirror or a solved challenge. **Recorded as a refusal, not as "not
  retrieved".**
- **A route for the owner, because the label does exist.** The NHS dictionary of medicines and devices lists
  Trasicor 40 mg tablets (Amdipharm Plc) and Slow-Trasicor 160 mg tablets (Sigma Pharmaceuticals Plc) as actual
  medicinal products, so oxprenolol is a current UK product and has a current summary of product characteristics.
  Its section 5.1, "Pharmacodynamic properties", is the part item 14(a) would read. The owner can open
  `medicines.org.uk` and search either name.

## Routing around the refusal, 14 October 2026

The owner: "Can you please try to route around the refusal instead?" The ruling of 11 October names three routes it
rules out - a cache, a mirror and a solved challenge - and one it rules in: "Record the refusal and use the
provider's other publications." The second was taken and it answered.

- **The Irish regulator, HPRA, holds the record and it reads normally.** `hpra.ie`
  `/find-a-medicine/for-human-use/authorised-medicines/details/-1779`, read 14 October 2026 with
  `experiments/read_public_page.py` (sha256 7a55618dfc93cfdc...). Verbatim from the page: **TRASICOR**, active
  substance **OXPRENOLOL HYDROCHLORIDE**, tablets, licence **PA0028/003/004**, holder **Ciba Geigy Plc**, licence
  issued **29/08/1975**, **"Withdrawn: 31/03/1995"**. The withdrawn-medicines record (sha256 3567ffc166a3478f...)
  carries the same. **No summary of product characteristics is attached to either**, which is what a licence
  withdrawn in 1995 would be expected to look like: it predates online SPC publication.
- **No eMC product page for Trasicor, Slow-Trasicor or pindolol surfaces in search either.** The eMC carries what
  marketing-authorisation holders submit, and a legacy product often is not there at all.

**So the refusal is probably not what is blocking this.** The NHS dictionary of medicines and devices still lists
Trasicor 40 mg and Slow-Trasicor 160 mg as actual medicinal products, so the drug is still dispensed in the United
Kingdom; but an authorisation from 1975 may have no published SPC anywhere, in which case part (a) cannot be
satisfied for oxprenolol by any regulator and the answer is to leave it out rather than to keep searching. That is
a better outcome than a blocked request, because it is a fact about the drug rather than about a web server.

**What was not done, and why it is the owner's to decide.** The three routes the ruling rules out - a cache, a
mirror, a solved human check - were not taken. That prohibition is the owner's own standing ruling, dated 11
October and quoted in `CLAUDE.md`, and it holds "for every session and for every research pass a session starts".
A preregistered study whose source rules carry an unrecorded exception is weaker for it, and the exception would
apply to every future pass and not only to this drug. **If the ruling is to change, it should be amended and dated
like the others**, and this session will make that amendment on the word. What is lost by not taking those routes,
in this instance, is probably nothing: the HPRA record suggests there is no SPC to reach.

## One thing this pass raises that it cannot settle

**Item 14(a) says "the drug's own label", and the screen reads openFDA, which is the United States only.** A UK
summary of product characteristics is a different regulator's label for the same drug, with its pharmacology in
section 5.1 rather than in a "clinical pharmacology" section. By item 14(a)'s wording - "a pharmacology section -
mechanism of action, clinical pharmacology, pharmacodynamics or description" - an SPC's section 5.1 qualifies, and
oxprenolol has **no** US label, so without a non-US label it can never satisfy (a) whatever its pharmacology. But
admitting a second regulator's labels is a change of scope for the rule and not a detail: it widens the evidence
base for every drug the rule governs, not only this one, and the two regulators do not write the same sentences.
**Whether item 14(a) reads a non-US label is the owner's and is recorded as decision 3.1(b).**
