# Plasma substrate binders in the graph

Built by experiments/build_plasma_binder_variant.py from `data/processed/graph_full_neuronal` into `data/processed/graph_full_neuronal_binders`. 27 `binds` edges, sign 0, each from the extracellular copy of a cargo metabolite to its binder, carried on gene nodes, which is the direction and sign the graph's own small-molecule edges use.

The same run built `data/processed/graph_full_neuronal_split_binders` from the matching source graphs. On a gene/protein split graph the carriage targets the protein node, because binding is the protein's property and the gene node there holds only expression, so the counts above hold but the edge endpoints differ.

## Cargo connected, by binder

| binder | cargo | edges |
|---|---|---|
| AFP | Cu(2+) | 1 |
| ALB | (4Z,15Z)-bilirubin IXalpha, Ca(2+), Cu cation, Zn(2+), arachidonate, linoleate, oleate, palmitate, stearate, thyroxine | 10 |
| APOA1 | cholesterol | 1 |
| APOB | cholesterol, cholesteryl ester | 2 |
| ORM1 | progesterone | 1 |
| ORM2 | progesterone | 1 |
| RBP4 | retinol | 1 |
| SERPINA6 | cortisol, progesterone | 2 |
| SERPINA7 | thyroxine, triiodothyronine | 2 |
| SHBG | 5alpha-dihydrotestosterone, estradiol-17beta, testosterone | 3 |
| TTR | L-thyroxine, retinol, triiodothyronine | 3 |

By source: {'pubmed': 17, 'uniprot_binding_site': 8, 'uniprot_function': 2}. Rows read: 30 (22 curated, 8 from UniProt binding-site features).

## Binder genes absent from the graph

None.

## Cargo with no Human-GEM metabolite

None.

## Cargo with no extracellular copy in this graph

None.

## Edges the graph already held

- MAM01615e -> GENE:SERPINA6
- MAM02998e -> GENE:SERPINA7
- MAM02998e -> GENE:TTR

