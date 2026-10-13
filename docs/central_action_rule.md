# Which drugs act centrally, by the measured rule

Written by `experiments/apply_central_action_rule.py` over `data/processed/label_review/evidence_any_atc_cap3` against `data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes`. The rule, approved by the owner on 12 October 2026 for the beta blockers: a drug is admitted when a mechanism target of it is a brain-expressed node of the graph and the drug has a measured Kp,uu,brain of at least 0.1. Nothing is applied by this script; `configs/centrally_acting_drugs.csv` carries the evidence per drug.

924 candidate drugs, 222 of them already in the table in use.

| verdict | all candidates | not already admitted |
| --- | --- | --- |
| admitted | 110 | 40 |
| below the Kp,uu,brain threshold | 28 | 24 |
| target only, no measured Kp,uu,brain | 619 | 492 |
| no brain-expressed target | 167 | 146 |

## Drugs already in the table that the rule would turn away

152 of the 222 admitted drugs. The rule is for admitting drugs the nervous-system filter left out; it is not proposed as a filter on the drugs already in, and this list is here so the owner can see what applying it both ways would cost.

| drug |
| --- |
| Lumateperone |
| acamprosate |
| agomelatine |
| alfentanil |
| almotriptan |
| amifampridine |
| amisulpride |
| amitriptyline |
| amphetamine |
| articaine |
| asenapine |
| aspirin |
| atogepant |
| atomoxetine |
| benperidol |
| benztropine |
| betahistine |
| brexpiprazole |
| brivaracetam |
| bromazepam |
| bromocriptine |
| bupivacaine |
| buprenorphine |
| butorphanol |
| cabergoline |
| cannabidiol |
| cariprazine |
| cenobamate |
| cevimeline |
| chlordiazepoxide |
| chloroprocaine |
| citalopram |
| clobazam |
| clomipramine |
| clonazepam |
| daridorexant |
| desipramine |
| desvenlafaxine |
| diflunisal |
| dihydrocodeine |
| dihydroergotamine |
| disulfiram |
| dothiepin |
| droperidol |
| duloxetine |
| dyclonine |
| eletriptan |
| eptinezumab |
| ergotamine |
| esketamine |
| eslicarbazepine |
| estazolam |
| ethotoin |
| etomidate |
| felbamate |
| fenfluramine |
| fentanyl |
| flunitrazepam |
| fluoxetine |
| flurazepam |
| fluvoxamine |
| fosphenytoin |
| frovatriptan |
| ganaxolone |
| hydromorphone |
| hydroxybutyrate |
| iloperidone |
| imipramine |
| istradefylline |
| ketamine |
| lacosamide |
| lasmiditan |
| lemborexant |
| levodopa |
| lisdexamfetamine |
| lofexidine |
| lorazepam |
| lurasidone |
| maprotiline |
| meperidine |
| mephenytoin |
| mepivacaine |
| methadone |
| methamphetamine |
| methohexital |
| methsuximide |
| methysergide |
| milnacipran |
| moclobemide |
| molindone |
| nalbuphine |
| nalmefene |
| naratriptan |
| nicotine |
| nortriptyline |
| opicapone |
| oxazepam |
| oxcarbazepine |
| paroxetine |
| pentazocine |
| pentobarbital |
| perampanel |
| pilocarpine |
| pimavanserin |
| pimozide |
| pitolisant |
| pramipexole |
| prilocaine |
| procaine |
| prochlorperazine |
| procyclidine |
| promazine |
| propofol |
| protriptyline |
| pyridostigmine |
| quazepam |
| ramelteon |
| rasagiline |
| reboxetine |
| remifentanil |
| remimazolam |
| retigabine |
| riluzole |
| rivastigmine |
| rizatriptan |
| ropinirole |
| ropivacaine |
| rotigotine |
| rufinamide |
| scopolamine |
| secobarbital |
| sertraline |
| sufentanil |
| sulpiride |
| sumatriptan |
| suvorexant |
| tapentadol |
| tasimelteon |
| tetrabenazine |
| tetracaine |
| tranylcypromine |
| triazolam |
| trihexyphenidyl |
| trimethadione |
| ubrogepant |
| valproate |
| venlafaxine |
| vigabatrin |
| ziconotide |
| zolmitriptan |
| zonisamide |
| zopiclone |

## The drugs the rule would admit

| drug | brain-expressed targets | Kp,uu,brain | species | values | source files |
| --- | --- | --- | --- | --- | --- |
| metformin | GPD2, MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND4L, MT-ND5, MT-ND6, NDUFA1, NDUFA10, NDUFA11, NDUFA12, NDUFA13, NDUFA2, NDUFA3, NDUFA5, NDUFA6, NDUFA7, NDUFA8, NDUFA9, NDUFAB1, NDUFAF2, NDUFAF3, NDUFB1, NDUFB10, NDUFB11, NDUFB2, NDUFB3, NDUFB4, NDUFB5, NDUFB6, NDUFB7, NDUFB8, NDUFB9, NDUFC1, NDUFC2, NDUFS1, NDUFS2, NDUFS3, NDUFS4, NDUFS5, NDUFS6, NDUFS7, NDUFS8, NDUFV1, NDUFV2, NDUFV3 | 9.51 (9.51 to 9.51) | other animal | 1 | friden2011_via_chembl.csv |
| nifedipine | CACNA1C, CACNA1D, CACNA1F | 8.64 (1.03 to 16.2) | other animal | 2 | kern2026.csv |
| hydrocodone | OPRD1, OPRK1, OPRM1 | 1.96 (1.96 to 1.96) | other animal | 1 | ma2024.csv |
| flumazenil | GABRA1, GABRA2, GABRA3, GABRA4, GABRA5, GABRA6, GABRB1, GABRB2, GABRB3, GABRD, GABRE, GABRG2, GABRQ | 1.67 (1.67 to 1.67) | human | 1 | sato2021.csv |
| cyclobenzaprine | HTR2A, HTR2C | 1.62 (1.62 to 1.62) | other animal | 1 | ma2024.csv |
| isradipine | CACNA1C, CACNA1D, CACNA1F | 1.42 (1.42 to 1.42) | other animal | 1 | kern2026.csv |
| mecamylamine | CHRNA3 | 1.05 (1.05 to 1.05) | other animal | 1 | gulave2025.csv |
| codeine | OPRD1, OPRK1, OPRM1 | 1 (0.89 to 1.99) | other animal | 5 | friden2009_via_chembl.csv;friden2011_via_chembl.csv;gulave2025.csv;morales2024.csv;sato2021.csv |
| propranolol | ADRB1, ADRB2 | 0.947 (0.947 to 0.947) | non-human primate | 2 | morales2024.csv;sato2021.csv |
| dextromethorphan | GRIN2A, SIGMAR1 | 0.825 (0.825 to 0.825) | other animal | 1 | hu2025.csv |
| loratadine | HRH1 | 0.74 (0.74 to 0.74) | other animal | 1 | morales2024.csv |
| felodipine | CACNA1C, CACNA1D, CACNA1F, NR3C2 | 0.701 (0.701 to 0.701) | other animal | 1 | kern2026.csv |
| metoclopramide | DRD2, HTR3A, HTR4 | 0.669 (0.52 to 0.669) | other animal | 3 | gulave2025.csv;ma2024.csv;morales2024.csv |
| metoprolol | ADRB1 | 0.64 (0.64 to 2.53) | other animal | 7 | friden2009_via_chembl.csv;friden2011_via_chembl.csv;gulave2025.csv;hu2025.csv;kosugi2021.csv;morales2024.csv;sato2021.csv |
| diltiazem | CACNA1C, CACNA1D, CACNA1F | 0.623 (0.345 to 0.9) | other animal | 2 | gulave2025.csv;kern2026.csv |
| omeprazole | ATP4A | 0.622 (0.15 to 1.09) | other animal | 2 | gulave2025.csv;hu2025.csv |
| naproxen | PTGS1, PTGS2 | 0.601 (0.6 to 0.602) | other animal | 2 | morales2024.csv;vallianatou2022.csv |
| pioglitazone | PPARG | 0.561 (0.561 to 0.561) | other animal | 1 | morales2024.csv |
| ibuprofen | PTGS1, PTGS2 | 0.558 (0.556 to 0.56) | other animal | 2 | morales2024.csv;vallianatou2022.csv |
| theophylline | ADORA1, ADORA2A, ADORA2B, ADORA3, PDE3A, PDE3B, PDE4A, PDE4B, PDE4C, PDE4D | 0.52 (0.05 to 0.99) | other animal | 2 | gulave2025.csv;vallianatou2022.csv |
| pindolol | ADRB1, ADRB2 | 0.5 (0.34 to 2.89) | other animal | 7 | friden2009_via_chembl.csv;friden2011_via_chembl.csv;gulave2025.csv;kosugi2021.csv;morales2024.csv;sato2021.csv;vallianatou2022.csv |
| celecoxib | PTGS2 | 0.416 (0.416 to 0.416) | other animal | 1 | morales2024.csv |
| etoricoxib | PTGS2 | 0.388 (0.388 to 0.388) | other animal | 1 | morales2024.csv |
| nimodipine | CACNA1C, CACNA1D, CACNA1F, NR3C2 | 0.385 (0.319 to 0.451) | other animal | 2 | hu2025.csv;kern2026.csv |
| diclofenac | PTGS1, PTGS2 | 0.351 (0.0927 to 0.61) | other animal | 2 | friden2011_via_chembl.csv;morales2024.csv |
| nisoldipine | CACNA1C, CACNA1D, CACNA1F | 0.322 (0.322 to 0.322) | other animal | 1 | kern2026.csv |
| sulfasalazine | ALOX5, PTGS1, PTGS2 | 0.27 (0.27 to 0.27) | other animal | 1 | friden2011_via_chembl.csv |
| prazosin | ADRA1A, ADRA1B, ADRA1D | 0.253 (0.253 to 0.253) | non-human primate | 1 | sato2021.csv |
| pitavastatin | HMGCR | 0.242 (0.242 to 0.242) | non-human primate | 1 | sato2021.csv |
| verapamil | CACNA1C, CACNA1D, CACNA1F | 0.24 (0.24 to 0.24) | human | 1 | sato2021.csv |
| ondansetron | HTR3A | 0.236 (0.19 to 0.283) | other animal | 2 | gulave2025.csv;morales2024.csv |
| valdecoxib | PTGS2 | 0.222 (0.222 to 0.222) | other animal | 1 | morales2024.csv |
| nicardipine | CACNA1C, CACNA1D, CACNA1F | 0.219 (0.219 to 0.219) | other animal | 1 | kern2026.csv |
| flurbiprofen | PTGS1, PTGS2 | 0.206 (0.206 to 0.206) | other animal | 1 | morales2024.csv |
| oxprenolol | ADRB1, ADRB2 | 0.2 (0.2 to 3.25) | other animal | 5 | friden2009_via_chembl.csv;friden2011_via_chembl.csv;morales2024.csv;sato2021.csv;vallianatou2022.csv |
| warfarin | VKORC1 | 0.19 (0.19 to 0.19) | other animal | 1 | ma2024.csv |
| etoposide | TOP2B | 0.164 (0.164 to 0.164) | non-human primate | 1 | sato2021.csv |
| granisetron | HTR3A | 0.13 (0.13 to 0.13) | other animal | 1 | gulave2025.csv |
| ketoprofen | PTGS1, PTGS2 | 0.12 (0.12 to 0.12) | other animal | 1 | gulave2025.csv |
| indomethacin | PTGS1, PTGS2 | 0.11 (0.1 to 0.6) | other animal | 8 | friden2009_via_chembl.csv;friden2011_via_chembl.csv;gulave2025.csv;kosugi2021.csv;langthaler2024.csv;morales2024.csv;sato2021.csv;vallianatou2022.csv |

## Measured and below the threshold

24 drugs, at a threshold of 0.1.

| drug | brain-expressed targets | Kp,uu,brain |
| --- | --- | --- |
| dexamethasone | NR3C1 | 0.0933 |
| ketorolac | PTGS1, PTGS2 | 0.0834 |
| erlotinib | EGFR | 0.076 |
| quinidine | SCN1A, SCN2A, SCN3A, SCN7A, SCN8A, SCN9A | 0.0744 |
| fexofenadine | HRH1 | 0.072 |
| pemetrexed | DHFR, GART, TYMS | 0.059 |
| vinblastine | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 | 0.0564 |
| dantrolene | RYR1 | 0.053 |
| methotrexate | DHFR | 0.041 |
| loperamide | OPRM1 | 0.04 |
| baclofen | GABBR1, GABBR2 | 0.037 |
| nadolol | ADRB1, ADRB2 | 0.037 |
| atenolol | ADRB1 | 0.03 |
| amlodipine | CACNA1C, CACNA1D, CACNA1F | 0.022 |
| colchicine | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 | 0.0205 |
| cimetidine | HRH2 | 0.02 |
| desloratadine | HRH1 | 0.018 |
| ritonavir | CYP3A5 | 0.0104 |
| cetirizine | HRH1 | 0.01 |
| digoxin | ATP1A1, ATP1A2, ATP1A3, ATP1B1, ATP1B2, ATP1B3 | 0.01 |
| salbutamol | ADRB2 | 0.01 |
| lacidipine | CACNA1C, CACNA1D, CACNA1F | 0.008 |
| paclitaxel | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 | 0.0073 |
| lercanidipine | CACNA1C, CACNA1D, CACNA1F | 0.002 |

## A brain-expressed target and no measured Kp,uu,brain

492 drugs. The second half of the rule is unmet for want of a measurement, which is a gap in the exposure tables rather than a finding about the drug, so these are listed and not admitted.

| drug | brain-expressed targets |
| --- | --- |
| 4-AP | KCNA1, KCNA2, KCNA3, KCNA4, KCNA5, KCNB1, KCNB2, KCNC1, KCNC3, KCNC4, KCND1, KCND2, KCND3, KCNG2, KCNH1, KCNH2, KCNH8, KCNQ1, KCNQ2, KCNQ3, KCNQ4, KCNQ5 |
| 4-methylpyrazole | ADH1B |
| 5-ASA | ALOX5, PPARG, PTGS1, PTGS2 |
| 5-FU | TYMS |
| 5-aza-2'-deoxycytidine | DNMT1, DNMT3A |
| A77 | DHODH |
| AMD3100 | CXCR4 |
| BCNU | GSR |
| Bivalirudin | F2 |
| CAS | ACE |
| CAS | ACE |
| Cantril | CHRM1, CHRM3 |
| Carfilzomib | ADRM1, PSMA1, PSMA2, PSMA3, PSMA4, PSMA5, PSMA6, PSMA7, PSMB1, PSMB10, PSMB2, PSMB3, PSMB4, PSMB5, PSMB6, PSMB7, PSMB8, PSMB9, PSMC1, PSMC2, PSMC3, PSMC4, PSMC5, PSMC6, PSMD1, PSMD11, PSMD12, PSMD13, PSMD14, PSMD2, PSMD3, PSMD4, PSMD6, PSMD7, PSMD8, SEM1 |
| DFMO | ODC1 |
| EACA | PLAT |
| FAMP | POLA1, POLA2, POLD1, POLD2, POLD3, POLD4, POLE, POLE2, POLE3, PRIM1, PRIM2, RRM1, RRM2B |
| FTY720 | S1PR1, S1PR2, S1PR3, S1PR5 |
| FdUrd | TYMS |
| ICI | ESR1 |
| L-threo-DOPS | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, ADRB1, ADRB2 |
| Lisinopril | ACE |
| Lutalyse | PTGFR |
| MDV3100 | AR |
| MK-462 | HTR1B, HTR1D |
| N-carbamylglutamate | CPS1 |
| NTBC | HPD |
| Octreotide | SSTR1, SSTR2, SSTR3, SSTR4 |
| Oestrogen | ESR1 |
| PCI-32765 | BTK |
| PGE2 | PTGER3, PTGER4 |
| Romidepsin | HDAC1, HDAC10, HDAC11, HDAC2, HDAC3, HDAC4, HDAC5, HDAC6, HDAC7, HDAC8, HDAC9 |
| SU5416 | FLT1, FLT4, KDR, KIT |
| Teprotumumab | IGF1R |
| abemaciclib | CDK4, CDK6 |
| abiraterone | AR, SRD5A1, SRD5A3 |
| abrocitinib | JAK1 |
| acalabrutinib | BTK |
| acarbose | AMY2A |
| acebutolol | ADRB1 |
| acitretin | RARA, RARB, RARG, RXRA, RXRB, RXRG |
| acrivastine | HRH1 |
| adenosine | ADORA1, ADORA2A, ADORA2B, ADORA3 |
| afamelanotide | MC1R |
| afatinib | EGFR, ERBB2, ERBB4 |
| aldesleukin | IL2RG |
| alectinib | ALK, EML4, RET |
| alemtuzumab | CD52 |
| alendronate | FDPS |
| alfuzosin | ADRA1A, ADRA1B, ADRA1D |
| alosetron | HTR3A |
| alpelisib | PIK3CA |
| alpha-methyl-p-tyrosine | TH |
| ambrisentan | EDNRA, EDNRB |
| amiloride | SCNN1B, SCNN1G |
| aminoglutethimide | CYP11A1 |
| aminophylline | ADORA1, ADORA2A, ADORA2B, ADORA3, PDE3A, PDE3B, PDE4A, PDE4B, PDE4C, PDE4D |
| amiodarone | KCNH2 |
| amivantamab | EGFR, MET |
| amlodipine | CACNA1C, CACNA1D, CACNA1F |
| amsacrine | TOP2B |
| anagrelide | PDE3A, PDE3B |
| apalutamide | AR |
| apraclonidine | ADRA2A, ADRA2B, ADRA2C |
| apremilast | PDE4A, PDE4B, PDE4C, PDE4D |
| aprepitant | TACR1 |
| arsenic trioxide | TXNRD1 |
| asciminib | ABL1, BCR |
| atorvastatin | HMGCR |
| atosiban | OXTR |
| atropine | CHRM1, CHRM2, CHRM3 |
| avacopan | C5AR1 |
| avanafil | PDE5A |
| avapritinib | KIT, PDGFRA |
| avatrombopag | MPL |
| avelumab | CD274 |
| axitinib | FLT1, FLT4, KDR |
| azelastine | HRH1 |
| balsalazide | ALOX5, PPARG, PTGS1, PTGS2 |
| baricitinib | JAK1, JAK2 |
| bazedoxifene | ESR1 |
| belantamab mafodotin | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| belatacept | CD86 |
| belimumab | TNFSF13B |
| belinostat | HDAC1, HDAC10, HDAC11, HDAC2, HDAC3, HDAC4, HDAC5, HDAC6, HDAC7, HDAC8, HDAC9 |
| belumosudil | ROCK1, ROCK2 |
| belzutifan | EPAS1 |
| benazepril | ACE |
| benzphetamine | SLC18A2 |
| benzydamine | PTGS1, PTGS2 |
| bepridil | CACNA1A, CACNA1B, CACNA1C, CACNA1D, CACNA1E, CACNA1F, CACNA1G, CACNA1H, CACNA1I, CACNA2D1, CACNA2D3, CACNB1, CACNB2, CACNB3, CACNB4, CACNG1, CACNG2, CACNG3, CACNG7, CACNG8 |
| betaxolol | ADRB1 |
| bevacizumab | VEGFA |
| bexarotene | RXRA, RXRB, RXRG |
| bezafibrate | PPARA, PPARD, PPARG |
| bicalutamide | AR |
| bimatoprost | PTGFR |
| binimetinib | MAP2K1, MAP2K2 |
| bisoprolol | ADRB1 |
| blinatumomab | CD3D, CD3E |
| bortezomib | ADRM1, PSMA1, PSMA2, PSMA3, PSMA4, PSMA5, PSMA6, PSMA7, PSMB1, PSMB10, PSMB2, PSMB3, PSMB4, PSMB5, PSMB6, PSMB7, PSMB8, PSMB9, PSMC1, PSMC2, PSMC3, PSMC4, PSMC5, PSMC6, PSMD1, PSMD11, PSMD12, PSMD13, PSMD14, PSMD2, PSMD3, PSMD4, PSMD6, PSMD7, PSMD8, SEM1 |
| bosentan | EDNRA, EDNRB |
| brentuximab vedotin | TNFRSF8, TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| brigatinib | ALK, EGFR |
| brimonidine | ADRA2A, ADRA2B, ADRA2C |
| brinzolamide | CA2 |
| brodalumab | IL17RA |
| budesonide | NR3C1 |
| butabarbital | GABRA1, GABRA2, GABRA3, GABRA4, GABRA5, GABRA6, GABRB1, GABRB2, GABRB3, GABRD, GABRE, GABRG2, GABRQ |
| cabazitaxel | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| cabozantinib | KDR, MET |
| capecitabine | TYMS |
| caplacizumab | VWF |
| capmatinib | MET |
| captopril | ACE |
| carbidopa | DDC |
| carbinoxamine | HRH1 |
| carboprost | PTGFR |
| carteolol | ADRB1, ADRB2 |
| carvedilol | ADRA1A, ADRA1B, ADRA1D, ADRB1, ADRB2 |
| celiprolol | ADRA2A, ADRA2B, ADRA2C, ADRB1, ADRB2 |
| ceritinib | ALK, EML4, NPM1 |
| cerivastatin | HMGCR |
| chlorpheniramine | HRH1 |
| chlorzoxazone | KCNN4 |
| ciclesonide | NR3C1 |
| cilostazol | PDE3A |
| clemastine | HRH1 |
| clobetasol | NR3C1 |
| clofarabine | POLA1, POLA2, POLD1, POLD2, POLD3, POLD4, POLE, POLE2, POLE3, PRIM1, PRIM2, RRM1, RRM2B |
| clomiphene | ESR1 |
| clopidogrel | P2RY12 |
| clorazepate | GABRA1, GABRA2, GABRA3, GABRA4, GABRA5, GABRA6, GABRB1, GABRB2, GABRB3, GABRD, GABRE, GABRG2, GABRQ |
| cobicistat | CYP3A5 |
| cobimetinib | MAP2K1, MAP2K2 |
| conivaptan | AVPR2 |
| cortisol | NR3C1 |
| cyclizine | HRH1 |
| cyclopentolate | CHRM3 |
| cyclosporine | PPIA |
| cyproheptadine | HRH1 |
| cyproterone | AR, NR3C1 |
| cytarabine | POLA1, POLA2, POLD1, POLD2, POLD3, POLD4, POLE, POLE2, POLE3, PRIM1, PRIM2 |
| dabigatran | F2 |
| dabrafenib | BRAF |
| daclizumab | IL2RG |
| dacomitinib | EGFR, ERBB2, ERBB4 |
| danazol | AR |
| daratumumab | CD38 |
| darifenacin | CHRM2, CHRM3 |
| darolutamide | AR |
| deflazacort | NR3C1 |
| delta | CNR1 |
| demethyl | NR3C1 |
| desmopressin | AVPR1B, AVPR2 |
| desonide | NR3C1 |
| dexamethasone | NR3C1 |
| dexrazoxane | TOP2B |
| dextrorphan | OPRM1 |
| diazoxide | KCNJ11, KCNJ8 |
| dicyclomine | CHRM1, CHRM3 |
| diethylpropion | SLC6A3 |
| dimenhydrinate | HRH1 |
| dimethyl fumarate | KEAP1 |
| dipyridamole | PDE10A, PDE1A, PDE1B, PDE1C, PDE2A, PDE3A, PDE3B, PDE4A, PDE4B, PDE4C, PDE4D, PDE5A, PDE6B, PDE6C, PDE6D, PDE6H, PDE7A, PDE7B, PDE8A, PDE8B, PDE9A, SLC29A1 |
| disopyramide | SCN1A, SCN2A, SCN3A, SCN7A, SCN8A, SCN9A |
| dobutamine | ADRB1, ADRB2 |
| docetaxel | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| dofetilide | KCNH2 |
| domperidone | DRD2 |
| dopamine | ADRB1 |
| dorzolamide | CA2 |
| doxapram | KCNK3, KCNK9 |
| doxazosin | ADRA1A, ADRA1B, ADRA1D |
| drospirenone | NR3C2 |
| dulaglutide | GLP1R |
| durvalumab | CD274 |
| dutasteride | SRD5A1, SRD5A3 |
| duvelisib | PIK3CD |
| dyphylline | PDE3A, PDE4A, PDE4B, PDE4C, PDE4D |
| eculizumab | C5 |
| edrophonium | ACHE |
| eliglustat | UGCG |
| eltrombopag | MPL |
| eluxadoline | OPRM1 |
| emedastine | HRH1 |
| encorafenib | BRAF |
| entrectinib | ALK, NTRK1, NTRK2, NTRK3, ROS1 |
| ephedrine | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, ADRB1, ADRB2 |
| epinephrine | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, ADRB1, ADRB2 |
| eplerenone | NR3C2 |
| eribulin | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| eslicarbazepine | SCN1A, SCN2A, SCN3A, SCN7A, SCN8A, SCN9A |
| esmolol | ADRB1, ADRB2 |
| esomeprazole | ATP4A |
| estradiol | ESR1 |
| estradiol | ESR1 |
| estradiol | ESR1 |
| estrone | ESR1 |
| ethinyl estradiol | ESR1 |
| etodolac | PTGS2 |
| everolimus | FKBP1A |
| exenatide | GLP1R |
| famotidine | HRH2 |
| fedratinib | FLT3, JAK2 |
| fenofibrate | PPARA |
| fenofibric | PPARA |
| fenoprofen | PTGS1, PTGS2 |
| fesoterodine | CHRM2, CHRM3 |
| filgrastim | CSF3R |
| flavoxate | PDE4A, PDE4B, PDE4C, PDE4D, PDE7A, PDE7B, PDE8A, PDE8B |
| fludrocortisone | NR3C2 |
| flunisolide | NR3C1 |
| fluocinolone | NR3C1 |
| fluoxymesterone | AR |
| fluphenazine | DRD2 |
| flutamide | AR |
| fluticasone | NR3C1 |
| fluticasone | NR3C1 |
| formoterol | ADRB2 |
| fosaprepitant | TACR1 |
| fosinopril | ACE |
| gefitinib | EGFR |
| gemcitabine | POLA1, POLA2, POLD1, POLD2, POLD3, POLD4, POLE, POLE2, POLE3, PRIM1, PRIM2, RRM1, RRM2B |
| gemfibrozil | PPARA |
| gemtuzumab ozogamicin | CD33 |
| gilteritinib | AXL, FLT3 |
| givosiran | ALAS1 |
| glasdegib | SMO |
| glibenclamide | ABCC8, KCNJ11 |
| gliclazide | ABCC8, KCNJ11 |
| glimepiride | ABCC8, KCNJ11 |
| glipizide | ABCC8, KCNJ11 |
| glycopyrrolate | CHRM1, CHRM2, CHRM3, CHRM4, CHRM5 |
| guanfacine | ADRA2A, ADRA2B, ADRA2C |
| guanidinium | KCNA1, KCNA2, KCNA3, KCNA4, KCNA5, KCNB1, KCNB2, KCNC1, KCNC3, KCNC4, KCND1, KCND2, KCND3, KCNG2, KCNH1, KCNH2, KCNH8, KCNQ1, KCNQ2, KCNQ3, KCNQ4, KCNQ5 |
| histamine | HRH1 |
| hydroxychloroquine | TLR7, TLR9 |
| hydroxyurea | RRM1, RRM2B |
| ibandronate | FDPS |
| imiquimod | TLR7 |
| inotuzumab ozogamicin | CD22 |
| insulin | INSR |
| insulin glargine | INSR |
| insulin lispro | INSR |
| interferon alfa-2b | IFNAR1, IFNAR2 |
| interferon beta-1a | IFNAR1, IFNAR2 |
| interferon beta-1b | IFNAR1, IFNAR2 |
| ipratropium | CHRM3 |
| irinotecan | TOP1 |
| isatuximab | CD38 |
| isoproterenol | ADRB1, ADRB2 |
| isosorbide | GUCY1A1, GUCY1A2, GUCY1B1 |
| isosorbide-5-mononitrate | GUCY1A1, GUCY1A2, GUCY1B1 |
| ivabradine | HCN4 |
| ivosidenib | IDH1 |
| ixabepilone | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| ketorolac | PTGS1, PTGS2 |
| ketotifen | HRH1 |
| labetalol | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, ADRB1, ADRB2 |
| lanreotide | SSTR2 |
| lansoprazole | ATP4A |
| lapatinib | EGFR, ERBB2 |
| larotrectinib | NTRK1, NTRK2, NTRK3 |
| lasofoxifene | ESR1 |
| latanoprost | PTGFR |
| leflunomide | DHODH |
| lenalidomide | CRBN, CUL4A, DDB1, RBX1 |
| lenvatinib | FLT1, FLT4, KDR |
| levamlodipine | CACNA1C, CACNA1D, CACNA1F |
| levobunolol | ADRB1, ADRB2 |
| levocetirizine | HRH1 |
| levosimendan | PDE3A, PDE3B, TNNC1, TNNI3, TNNT2 |
| liraglutide | GLP1R |
| lithium | GSK3A, GSK3B, IMPA1 |
| lixisenatide | GLP1R |
| lomitapide | P4HB |
| lonafarnib | FNTA, FNTB |
| lorcaserin | HTR2C |
| lorlatinib | ALK, EML4 |
| lovastatin | HMGCR |
| lubiprostone | CLCN2 |
| luspatercept | TGFB1, TGFB2, TGFB3 |
| mefenamic | PTGS1, PTGS2 |
| meloxicam | PTGS2 |
| metaproterenol | ADRB2 |
| methyldopa | ADRA2A, ADRA2B, ADRA2C |
| methylene | GUCY1A1, GUCY1A2, GUCY1B1 |
| methylergometrine | DRD1 |
| methylnaltrexone | OPRM1 |
| methylprednisolone | NR3C1 |
| methylprednisolone | NR3C1 |
| methyltestosterone | AR |
| metipranolol | ADRB1, ADRB2 |
| metreleptin | LEPR |
| mexiletine | SCN1A, SCN2A, SCN3A, SCN7A, SCN8A, SCN9A |
| midodrine | ADRA1A, ADRA1B, ADRA1D |
| mifepristone | NR3C1 |
| migalastat | GLA |
| miglustat | UGCG |
| minoxidil | ABCC9, KCNJ11 |
| misoprostol | PTGER3 |
| mitapivat | PKLR |
| mitotane | CYP11A1 |
| moexipril | ACE |
| mometasone | NR3C1 |
| monamine | HRH1 |
| moxetumomab pasudotox | CD22, EEF2 |
| moxonidine | NISCH |
| mycophenolate | IMPDH1, IMPDH2 |
| mycophenolic | IMPDH1, IMPDH2 |
| nabilone | CNR1 |
| nabumetone | PTGS2 |
| naloxone | OPRD1, OPRK1, OPRM1 |
| nandrolone | AR |
| naphazoline | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C |
| natalizumab | ITGA4, ITGB1 |
| nateglinide | ABCC8, KCNJ11 |
| nebivolol | ADRB1, ADRB2 |
| neratinib | EGFR, ERBB2, ERBB4 |
| nilotinib | ABL1, BCR |
| nilutamide | AR |
| nintedanib | FGFR1, FGFR2, FGFR3, FGFR4, FLT1, FLT4, KDR, PDGFRA, PDGFRB |
| niraparib | PARP1, PARP2 |
| nitroglycerin | GUCY1A1, GUCY1A2, GUCY1B1 |
| nizatidine | HRH2 |
| norepinephrine | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, ADRB1, ADRB2 |
| ogen | ESR1 |
| olaparib | PARP1, PARP2, PARP3 |
| olopatadine | HRH1 |
| olsalazine | ALOX5, PPARG, PTGS1, PTGS2 |
| orlistat | LIPF, PNLIP |
| osimertinib | EGFR |
| oxandrolone | AR |
| oxaprozin | PTGS1, PTGS2 |
| oxybutynin | CHRM2, CHRM3 |
| oxymetholone | AR |
| palbociclib | CCND1, CDK4, CDK6 |
| palonosetron | HTR3A |
| pamidronate | FDPS |
| panitumumab | EGFR |
| panobinostat | HDAC1, HDAC10, HDAC11, HDAC2, HDAC3, HDAC4, HDAC5, HDAC6, HDAC7, HDAC8, HDAC9 |
| pantoprazole | ATP4A |
| parathyroid hormone | PTH1R |
| parecoxib | PTGS2 |
| peginterferon alfa-2b | IFNAR1, IFNAR2 |
| peginterferon beta-1a | IFNAR1, IFNAR2 |
| pemigatinib | FGFR1, FGFR2, FGFR3 |
| penbutolol | ADRB1, ADRB2 |
| pentostatin | ADA |
| pentoxifylline | ADORA2A, ADORA2B, PDE10A, PDE1A, PDE1B, PDE1C, PDE2A, PDE3A, PDE3B, PDE4A, PDE4B, PDE4C, PDE4D, PDE5A, PDE6B, PDE6C, PDE6D, PDE6H, PDE7A, PDE7B, PDE8A, PDE8B, PDE9A |
| perhexiline | CPT1B, CPT2, MTOR |
| perindopril | ACE |
| perindopril | ACE |
| pertuzumab | ERBB2 |
| phendimetrazine | SLC6A3 |
| phenoxybenzamine | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C |
| phentolamine | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C |
| phenylephrine | ADRA1A, ADRA1B, ADRA1D |
| phenylpropanolamine | ADRA1A |
| pirbuterol | ADRB2 |
| piroxicam | PTGS1, PTGS2 |
| polatuzumab vedotin | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| pomalidomide | CRBN, CUL4A, DDB1, RBX1 |
| ponatinib | ABL1, BCR |
| ponesimod | S1PR1 |
| pralatrexate | DHFR |
| pralidoxime | ACHE |
| pralsetinib | CCDC6, KIF5B, RET |
| prasugrel | P2RY12 |
| pravastatin | HMGCR |
| prednisolone | NR3C1 |
| prednisolone | NR3C1 |
| prednisone | NR3C1 |
| procainamide | SCN1A, SCN2A, SCN3A, SCN7A, SCN8A, SCN9A |
| promethazine | HRH1 |
| propafenone | ADRB1, ADRB2, SCN1A, SCN2A, SCN3A, SCN7A, SCN8A, SCN9A |
| propantheline | CHRM1, CHRM2, CHRM3 |
| prucalopride | HTR4 |
| quinagolide | DRD2 |
| quinapril | ACE |
| rabeprazole | ATP4A |
| raloxifene | ESR1 |
| ramipril | ACE |
| ramucirumab | KDR |
| ranitidine | HRH2 |
| ravulizumab | C5 |
| regadenoson | ADORA2A |
| rescinnamine | SLC18A2 |
| reserpine | SLC18A2 |
| retinoic | RARA, RARB, RARG, RXRA, RXRB, RXRG |
| ribavirin | IMPDH1 |
| ribociclib | CDK4, CDK6 |
| rimonabant | CNR1 |
| ripretinib | KIT, PDGFRA |
| risedronate | FDPS |
| rivaroxaban | F10 |
| rofecoxib | PTGS2 |
| roflumilast | PDE4A, PDE4B, PDE4C, PDE4D |
| rolapitant | TACR1 |
| romiplostim | MPL |
| ropeginterferon alfa-2b | IFNAR1, IFNAR2 |
| rosiglitazone | PPARG |
| rosuvastatin | HMGCR |
| rucaparib | PARP1, PARP2, PARP3 |
| salmeterol | ADRB2 |
| sargramostim | CSF2RA, CSF2RB |
| sarilumab | IL6R |
| satralizumab | IL6R, IL6ST |
| selinexor | XPO1 |
| selumetinib | MAP2K1, MAP2K2 |
| semaglutide | GLP1R |
| setmelanotide | MC4R |
| sibutramine | SLC6A3 |
| sildenafil | PDE5A |
| silodosin | ADRA1A |
| simvastatin | HMGCR |
| siponimod | S1PR1 |
| sirolimus | FKBP1A |
| sitaxsentan | EDNRA, EDNRB |
| solifenacin | CHRM2, CHRM3 |
| somapacitan | GHR |
| somatropin | GHR |
| sonidegib | SMO |
| sotalol | ADRB1, ADRB2, KCNH2 |
| sotorasib | KRAS |
| spironolactone | NR3C2 |
| sulindac | PTGS1, PTGS2 |
| tacrolimus | FKBP1A |
| tadalafil | PDE5A |
| tagraxofusp | EEF2, IL3RA |
| talazoparib | PARP1, PARP2 |
| tamoxifen | ESR1 |
| tamsulosin | ADRA1A, ADRA1B, ADRA1D |
| tazemetostat | EZH2 |
| teduglutide | GLP2R |
| tegaserod | HTR4 |
| teniposide | TOP2B |
| tepotinib | MET |
| terazosin | ADRA1A, ADRA1B, ADRA1D |
| terbutaline | ADRB2 |
| teriparatide | PTH1R |
| testosterone | AR |
| testosterone | AR |
| testosterone | AR |
| testosterone | AR |
| tetrahydrozoline | ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C |
| thalidomide | CRBN, CUL4A, DDB1, RBX1 |
| thyroxine | THRA, THRB |
| ticagrelor | P2RY12 |
| ticlopidine | P2RY12 |
| timolol | ADRB1, ADRB2 |
| tiotropium | CHRM3 |
| tirzepatide | GIPR, GLP1R |
| tivozanib | FLT1, FLT4, KDR, KIT, PDGFRB |
| tizanidine | ADRA2A, ADRA2B, ADRA2C |
| tocilizumab | IL6R |
| tofacitinib | JAK1, JAK2, JAK3, TYK2 |
| tolazamide | ABCC8, KCNJ11 |
| tolmetin | PTGS1, PTGS2 |
| tolterodine | CHRM2, CHRM3 |
| tolvaptan | AVPR2 |
| topotecan | TOP1 |
| toremifene | ESR1 |
| trametinib | MAP2K1, MAP2K2 |
| trandolapril | ACE |
| trastuzumab | ERBB2 |
| travoprost | PTGFR |
| triamcinolone | NR3C1 |
| triamcinolone | NR3C1 |
| triamterene | SCNN1B, SCNN1G |
| triiodothyronine | THRA, THRB |
| trilaciclib | CDK4, CDK6 |
| trimeprazine | HRH1 |
| troglitazone | PPARG |
| tropicamide | CHRM3 |
| trospium | CHRM2, CHRM3 |
| tucatinib | ERBB2 |
| ustekinumab | IL12A |
| vardenafil | PDE5A |
| vedolizumab | ITGA4 |
| vemurafenib | BRAF |
| venetoclax | BCL2 |
| vincristine | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| vinflunine | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| vinorelbine | TUBA1A, TUBA1B, TUBA1C, TUBA4A, TUBB, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6 |
| vismodegib | SMO |
| vorapaxar | F2R |
| voretigene neparvovec | RPE65 |
| vosoritide | NPR2 |
| zanubrutinib | BTK |
| zileuton | ALOX5 |
| zoledronic | FDPS |
