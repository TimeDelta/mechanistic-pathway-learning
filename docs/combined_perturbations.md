# Combination rows of the evidence table

Written by `experiments/build_combined_perturbations.py` from `configs/comparative_checks.csv` on the labels of `data/processed/evidence_full_v4` and the nodes of `data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes`.

106 combinations over 135 statements, 317 positive pairs. 79 statements wrote nothing.

| kind | combinations | positive pairs |
| --- | --- | --- |
| carrier_level | 2 | 29 |
| displacement | 5 | 40 |
| enzyme_activity | 17 | 140 |
| second_drug | 82 | 108 |

## The combinations

| perturbation | kind | symptoms | statement |
| --- | --- | --- | --- |
| alprazolam + loxapine | second_drug | 1 | `pair_positive_alprazolam_with_loxapine_cognitive_impairment` |
| alprazolam + loxapine | second_drug | 1 | `pair_positive_alprazolam_with_loxapine_somnolence_or_hypersomnia` |
| amisulpride + loxapine | second_drug | 1 | `pair_positive_amisulpride_with_loxapine_somnolence_or_hypersomnia` |
| amitriptyline + duloxetine | second_drug | 1 | `pair_positive_amitriptyline_with_duloxetine_cognitive_impairment` |
| amitriptyline + duloxetine | second_drug | 1 | `pair_positive_amitriptyline_with_duloxetine_somnolence_or_hypersomnia` |
| amitriptyline + fluvoxamine | second_drug | 1 | `pair_positive_amitriptyline_with_fluvoxamine_somnolence_or_hypersomnia` |
| amitriptyline + loxapine | second_drug | 1 | `pair_positive_amitriptyline_with_loxapine_cognitive_impairment` |
| amitriptyline + loxapine | second_drug | 1 | `pair_positive_amitriptyline_with_loxapine_somnolence_or_hypersomnia` |
| aripiprazole + CYP2D6 lowered | enzyme_activity | 18 | `enzyme_aripiprazole_CYP2D6` |
| aripiprazole + loxapine | second_drug | 1 | `pair_positive_aripiprazole_with_loxapine_cognitive_impairment` |
| aripiprazole + loxapine | second_drug | 1 | `pair_positive_aripiprazole_with_loxapine_somnolence_or_hypersomnia` |
| atomoxetine + CYP2D6 lowered | enzyme_activity | 11 | `enzyme_atomoxetine_CYP2D6` |
| carbamazepine + valproate + ALB lowered | displacement | 12 | `displacement_carbamazepine_by_valproate` |
| cevimeline + CYP2D6 lowered | enzyme_activity | 5 | `enzyme_cevimeline_CYP2D6` |
| chlordiazepoxide + loxapine | second_drug | 1 | `pair_positive_chlordiazepoxide_with_loxapine_somnolence_or_hypersomnia` |
| chlorpromazine + loxapine | second_drug | 1 | `pair_positive_chlorpromazine_with_loxapine_somnolence_or_hypersomnia` |
| citalopram + CYP2C19 lowered | enzyme_activity | 16 | `enzyme_citalopram_CYP2C19` |
| clobazam + CYP2C19 lowered | enzyme_activity | 10 | `enzyme_clobazam_CYP2C19` |
| clobazam + loxapine | second_drug | 1 | `pair_positive_clobazam_with_loxapine_somnolence_or_hypersomnia` |
| clomipramine + duloxetine | second_drug | 1 | `pair_positive_clomipramine_with_duloxetine_cognitive_impairment` |
| clomipramine + duloxetine | second_drug | 1 | `pair_positive_clomipramine_with_duloxetine_somnolence_or_hypersomnia` |
| clomipramine + fluvoxamine | second_drug | 1 | `pair_positive_clomipramine_with_fluvoxamine_somnolence_or_hypersomnia` |
| clomipramine + loxapine | second_drug | 1 | `pair_positive_clomipramine_with_loxapine_cognitive_impairment` |
| clomipramine + loxapine | second_drug | 1 | `pair_positive_clomipramine_with_loxapine_somnolence_or_hypersomnia` |
| clonazepam + loxapine | second_drug | 1 | `pair_positive_clonazepam_with_loxapine_cognitive_impairment` |
| clonazepam + loxapine | second_drug | 1 | `pair_positive_clonazepam_with_loxapine_somnolence_or_hypersomnia` |
| clozapine + CYP2D6 lowered | enzyme_activity | 9 | `enzyme_clozapine_CYP2D6` |
| clozapine + loxapine | second_drug | 1 | `pair_positive_clozapine_with_loxapine_cognitive_impairment` |
| clozapine + loxapine | second_drug | 1 | `pair_positive_clozapine_with_loxapine_somnolence_or_hypersomnia` |
| dexmedetomidine + loxapine | second_drug | 1 | `pair_positive_dexmedetomidine_with_loxapine_somnolence_or_hypersomnia` |
| diazepam + loxapine | second_drug | 1 | `pair_positive_diazepam_with_loxapine_cognitive_impairment` |
| diazepam + loxapine | second_drug | 1 | `pair_positive_diazepam_with_loxapine_somnolence_or_hypersomnia` |
| diazepam + valproate + ALB lowered | displacement | 8 | `displacement_diazepam_by_valproate` |
| dihydrocodeine + loxapine | second_drug | 1 | `pair_positive_dihydrocodeine_with_loxapine_cognitive_impairment` |
| dihydrocodeine + loxapine | second_drug | 1 | `pair_positive_dihydrocodeine_with_loxapine_somnolence_or_hypersomnia` |
| dothiepin + paroxetine | second_drug | 1 | `pair_positive_dothiepin_with_paroxetine_somnolence_or_hypersomnia` |
| doxepin + CYP2C19 lowered | enzyme_activity | 7 | `enzyme_doxepin_CYP2C19` |
| doxepin + CYP2D6 lowered | enzyme_activity | 7 | `enzyme_doxepin_CYP2D6` |
| doxepin + loxapine | second_drug | 1 | `pair_positive_doxepin_with_loxapine_somnolence_or_hypersomnia` |
| duloxetine + amitriptyline | second_drug | 1 | `pair_positive_duloxetine_with_amitriptyline_cognitive_impairment` |
| duloxetine + amitriptyline | second_drug | 1 | `pair_positive_duloxetine_with_amitriptyline_somnolence_or_hypersomnia` |
| duloxetine + clomipramine | second_drug | 1 | `pair_positive_duloxetine_with_clomipramine_cognitive_impairment` |
| duloxetine + clomipramine | second_drug | 1 | `pair_positive_duloxetine_with_clomipramine_somnolence_or_hypersomnia` |
| duloxetine + imipramine | second_drug | 1 | `pair_positive_duloxetine_with_imipramine_cognitive_impairment` |
| duloxetine + imipramine | second_drug | 1 | `pair_positive_duloxetine_with_imipramine_somnolence_or_hypersomnia` |
| fluphenazine + loxapine | second_drug | 1 | `pair_positive_fluphenazine_with_loxapine_somnolence_or_hypersomnia` |
| flurazepam + loxapine | second_drug | 1 | `pair_positive_flurazepam_with_loxapine_cognitive_impairment` |
| flurazepam + loxapine | second_drug | 1 | `pair_positive_flurazepam_with_loxapine_somnolence_or_hypersomnia` |
| fluvoxamine + CYP2D6 lowered | enzyme_activity | 13 | `enzyme_fluvoxamine_CYP2D6` |
| fluvoxamine + amitriptyline | second_drug | 1 | `pair_positive_fluvoxamine_with_amitriptyline_somnolence_or_hypersomnia` |
| fluvoxamine + clomipramine | second_drug | 1 | `pair_positive_fluvoxamine_with_clomipramine_somnolence_or_hypersomnia` |
| fluvoxamine + imipramine | second_drug | 1 | `pair_positive_fluvoxamine_with_imipramine_somnolence_or_hypersomnia` |
| fosphenytoin + CYP2C9 lowered | enzyme_activity | 5 | `enzyme_fosphenytoin_CYP2C9` |
| gabapentin + morphine | second_drug | 1 | `pair_positive_gabapentin_with_morphine_somnolence_or_hypersomnia` |
| galantamine + CYP2D6 lowered | enzyme_activity | 10 | `enzyme_galantamine_CYP2D6` |
| haloperidol + loxapine | second_drug | 1 | `pair_positive_haloperidol_with_loxapine_cognitive_impairment` |
| haloperidol + loxapine | second_drug | 1 | `pair_positive_haloperidol_with_loxapine_somnolence_or_hypersomnia` |
| haloperidol + olanzapine | second_drug | 1 | `pair_positive_haloperidol_with_olanzapine_parkinsonism` |
| hydromorphone + loxapine | second_drug | 1 | `pair_positive_hydromorphone_with_loxapine_cognitive_impairment` |
| hydromorphone + loxapine | second_drug | 1 | `pair_positive_hydromorphone_with_loxapine_somnolence_or_hypersomnia` |
| hydroxyzine + loxapine | second_drug | 1 | `pair_positive_hydroxyzine_with_loxapine_somnolence_or_hypersomnia` |
| imipramine + duloxetine | second_drug | 1 | `pair_positive_imipramine_with_duloxetine_cognitive_impairment` |
| imipramine + duloxetine | second_drug | 1 | `pair_positive_imipramine_with_duloxetine_somnolence_or_hypersomnia` |
| imipramine + fluvoxamine | second_drug | 1 | `pair_positive_imipramine_with_fluvoxamine_somnolence_or_hypersomnia` |
| imipramine + loxapine | second_drug | 1 | `pair_positive_imipramine_with_loxapine_cognitive_impairment` |
| imipramine + loxapine | second_drug | 1 | `pair_positive_imipramine_with_loxapine_somnolence_or_hypersomnia` |
| lorazepam + loxapine | second_drug | 1 | `pair_positive_lorazepam_with_loxapine_cognitive_impairment` |
| lorazepam + loxapine | second_drug | 1 | `pair_positive_lorazepam_with_loxapine_somnolence_or_hypersomnia` |
| loxapine + alprazolam | second_drug | 1 | `pair_positive_loxapine_with_alprazolam_somnolence_or_hypersomnia` |
| loxapine + amisulpride | second_drug | 1 | `pair_positive_loxapine_with_amisulpride_somnolence_or_hypersomnia` |
| loxapine + amitriptyline | second_drug | 1 | `pair_positive_loxapine_with_amitriptyline_somnolence_or_hypersomnia` |
| loxapine + aripiprazole | second_drug | 1 | `pair_positive_loxapine_with_aripiprazole_somnolence_or_hypersomnia` |
| loxapine + chlordiazepoxide | second_drug | 1 | `pair_positive_loxapine_with_chlordiazepoxide_somnolence_or_hypersomnia` |
| loxapine + chlorpromazine | second_drug | 1 | `pair_positive_loxapine_with_chlorpromazine_somnolence_or_hypersomnia` |
| loxapine + clobazam | second_drug | 1 | `pair_positive_loxapine_with_clobazam_somnolence_or_hypersomnia` |
| loxapine + clomipramine | second_drug | 1 | `pair_positive_loxapine_with_clomipramine_somnolence_or_hypersomnia` |
| loxapine + clonazepam | second_drug | 1 | `pair_positive_loxapine_with_clonazepam_somnolence_or_hypersomnia` |
| loxapine + clozapine | second_drug | 1 | `pair_positive_loxapine_with_clozapine_somnolence_or_hypersomnia` |
| loxapine + dexmedetomidine | second_drug | 1 | `pair_positive_loxapine_with_dexmedetomidine_somnolence_or_hypersomnia` |
| loxapine + diazepam | second_drug | 1 | `pair_positive_loxapine_with_diazepam_somnolence_or_hypersomnia` |
| loxapine + dihydrocodeine | second_drug | 1 | `pair_positive_loxapine_with_dihydrocodeine_somnolence_or_hypersomnia` |
| loxapine + doxepin | second_drug | 1 | `pair_positive_loxapine_with_doxepin_somnolence_or_hypersomnia` |
| loxapine + droperidol | second_drug | 1 | `pair_positive_loxapine_with_droperidol_somnolence_or_hypersomnia` |
| loxapine + fluphenazine | second_drug | 1 | `pair_positive_loxapine_with_fluphenazine_somnolence_or_hypersomnia` |
| loxapine + flurazepam | second_drug | 1 | `pair_positive_loxapine_with_flurazepam_somnolence_or_hypersomnia` |
| loxapine + haloperidol | second_drug | 1 | `pair_positive_loxapine_with_haloperidol_somnolence_or_hypersomnia` |
| loxapine + hydromorphone | second_drug | 1 | `pair_positive_loxapine_with_hydromorphone_somnolence_or_hypersomnia` |
| loxapine + hydroxyzine | second_drug | 1 | `pair_positive_loxapine_with_hydroxyzine_somnolence_or_hypersomnia` |
| loxapine + imipramine | second_drug | 1 | `pair_positive_loxapine_with_imipramine_somnolence_or_hypersomnia` |
| loxapine + lorazepam | second_drug | 1 | `pair_positive_loxapine_with_lorazepam_somnolence_or_hypersomnia` |
| loxapine + midazolam | second_drug | 1 | `pair_positive_loxapine_with_midazolam_somnolence_or_hypersomnia` |
| loxapine + morphine | second_drug | 1 | `pair_positive_loxapine_with_morphine_somnolence_or_hypersomnia` |
| loxapine + olanzapine | second_drug | 1 | `pair_positive_loxapine_with_olanzapine_somnolence_or_hypersomnia` |
| loxapine + oxazepam | second_drug | 1 | `pair_positive_loxapine_with_oxazepam_somnolence_or_hypersomnia` |
| loxapine + oxycodone | second_drug | 1 | `pair_positive_loxapine_with_oxycodone_somnolence_or_hypersomnia` |
| loxapine + perampanel | second_drug | 1 | `pair_positive_loxapine_with_perampanel_somnolence_or_hypersomnia` |
| loxapine + pimozide | second_drug | 1 | `pair_positive_loxapine_with_pimozide_somnolence_or_hypersomnia` |
| loxapine + primidone | second_drug | 1 | `pair_positive_loxapine_with_primidone_somnolence_or_hypersomnia` |
| loxapine + quetiapine | second_drug | 1 | `pair_positive_loxapine_with_quetiapine_somnolence_or_hypersomnia` |
| loxapine + sulpiride | second_drug | 1 | `pair_positive_loxapine_with_sulpiride_somnolence_or_hypersomnia` |
| loxapine + zolpidem | second_drug | 1 | `pair_positive_loxapine_with_zolpidem_somnolence_or_hypersomnia` |
| loxapine + zopiclone | second_drug | 1 | `pair_positive_loxapine_with_zopiclone_somnolence_or_hypersomnia` |
| midazolam + loxapine | second_drug | 1 | `pair_positive_midazolam_with_loxapine_cognitive_impairment` |
| midazolam + loxapine | second_drug | 1 | `pair_positive_midazolam_with_loxapine_somnolence_or_hypersomnia` |
| morphine + gabapentin | second_drug | 1 | `pair_positive_morphine_with_gabapentin_somnolence_or_hypersomnia` |
| morphine + loxapine | second_drug | 1 | `pair_positive_morphine_with_loxapine_cognitive_impairment` |
| morphine + loxapine | second_drug | 1 | `pair_positive_morphine_with_loxapine_somnolence_or_hypersomnia` |
| olanzapine + haloperidol | second_drug | 1 | `pair_positive_olanzapine_with_haloperidol_parkinsonism` |
| olanzapine + loxapine | second_drug | 1 | `pair_positive_olanzapine_with_loxapine_cognitive_impairment` |
| olanzapine + loxapine | second_drug | 1 | `pair_positive_olanzapine_with_loxapine_somnolence_or_hypersomnia` |
| oxazepam + loxapine | second_drug | 1 | `pair_positive_oxazepam_with_loxapine_cognitive_impairment` |
| oxycodone + loxapine | second_drug | 1 | `pair_positive_oxycodone_with_loxapine_cognitive_impairment` |
| oxycodone + loxapine | second_drug | 1 | `pair_positive_oxycodone_with_loxapine_somnolence_or_hypersomnia` |
| paroxetine + dothiepin | second_drug | 1 | `pair_positive_paroxetine_with_dothiepin_somnolence_or_hypersomnia` |
| perphenazine + CYP2D6 lowered | enzyme_activity | 4 | `enzyme_perphenazine_CYP2D6` |
| phenytoin + ALB lowered | carrier_level | 6 | `carrier_phenytoin_albumin_label` |
| phenytoin + ALB lowered | carrier_level | 6 | `carrier_phenytoin_albumin_neurotoxicity_madan` |
| phenytoin + CYP2C9 lowered | enzyme_activity | 6 | `enzyme_phenytoin_CYP2C9` |
| phenytoin + aspirin + ALB lowered | displacement | 6 | `displacement_phenytoin_by_aspirin` |
| phenytoin + valproate + ALB lowered | displacement | 6 | `displacement_phenytoin_by_valproate` |
| pimozide + CYP2D6 lowered | enzyme_activity | 7 | `enzyme_pimozide_CYP2D6` |
| pimozide + loxapine | second_drug | 1 | `pair_positive_pimozide_with_loxapine_somnolence_or_hypersomnia` |
| protriptyline + CYP2D6 lowered | enzyme_activity | 1 | `enzyme_protriptyline_CYP2D6` |
| quetiapine + loxapine | second_drug | 1 | `pair_positive_quetiapine_with_loxapine_cognitive_impairment` |
| quetiapine + loxapine | second_drug | 1 | `pair_positive_quetiapine_with_loxapine_somnolence_or_hypersomnia` |
| tetrabenazine + CYP2D6 lowered | enzyme_activity | 8 | `enzyme_tetrabenazine_CYP2D6` |
| thioridazine + CYP2D6 lowered | enzyme_activity | 3 | `enzyme_thioridazine_CYP2D6` |
| valproate + ALB lowered | carrier_level | 8 | `carrier_valproate_albumin_neurologic_dore` |
| valproate + ALB lowered | carrier_level | 8 | `carrier_valproate_albumin_neurotoxicity_madan` |
| valproate + ALB lowered | carrier_level | 1 | `carrier_valproate_albumin_somnolence` |
| valproate + aspirin + ALB lowered | displacement | 8 | `displacement_valproate_by_aspirin` |
| zolpidem + loxapine | second_drug | 1 | `pair_positive_zolpidem_with_loxapine_cognitive_impairment` |
| zolpidem + loxapine | second_drug | 1 | `pair_positive_zolpidem_with_loxapine_somnolence_or_hypersomnia` |
| zopiclone + loxapine | second_drug | 1 | `pair_positive_zopiclone_with_loxapine_cognitive_impairment` |
| zopiclone + loxapine | second_drug | 1 | `pair_positive_zopiclone_with_loxapine_somnolence_or_hypersomnia` |

## Statements that wrote nothing

| statement | kind | reason |
| --- | --- | --- |
| `enzyme_amifampridine_NAT2` | enzyme_activity | RXCUI:2106338 has no kept positive symptom |
| `enzyme_brexpiprazole_CYP2D6` | enzyme_activity | RXCUI:1658314 has no kept positive symptom |
| `enzyme_brivaracetam_CYP2C19` | enzyme_activity | RXCUI:1739745 has no kept positive symptom |
| `enzyme_iloperidone_CYP2D6` | enzyme_activity | CID100071360 has no kept positive symptom |
| `enzyme_lofexidine_CYP2D6` | enzyme_activity | RXCUI:28863 has no kept positive symptom |
| `enzyme_pitolisant_CYP2D6` | enzyme_activity | RXCUI:2197878 has no kept positive symptom |
| `pair_positive_loxapine_with_midazolam_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_amitriptyline_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_diazepam_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_oxazepam_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_oxazepam_with_loxapine_somnolence_or_hypersomnia` | second_drug | CID100004616 has no kept positive symptom at somnolence_or_hypersomnia |
| `pair_positive_doxepin_with_loxapine_cognitive_impairment` | second_drug | CID100003158 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_doxepin_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_droperidol_with_loxapine_cognitive_impairment` | second_drug | CID100003168 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_droperidol_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_droperidol_with_loxapine_somnolence_or_hypersomnia` | second_drug | CID100003168 has no kept positive symptom at somnolence_or_hypersomnia |
| `pair_positive_loxapine_with_zolpidem_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_pimozide_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_pimozide_with_loxapine_cognitive_impairment` | second_drug | CID100016362 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_primidone_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_primidone_with_loxapine_cognitive_impairment` | second_drug | CID100004909 has no kept positive symptom at cognitive_impairment |
| `pair_positive_primidone_with_loxapine_somnolence_or_hypersomnia` | second_drug | CID100004909 has no kept positive symptom at somnolence_or_hypersomnia |
| `pair_positive_fluphenazine_with_loxapine_cognitive_impairment` | second_drug | CID100003372 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_fluphenazine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_flurazepam_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_aripiprazole_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_haloperidol_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_quetiapine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_hydroxyzine_with_loxapine_cognitive_impairment` | second_drug | CID100003658 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_hydroxyzine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_imipramine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_alprazolam_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_olanzapine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_lorazepam_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_morphine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_oxycodone_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_hydromorphone_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_dihydrocodeine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_zopiclone_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_clobazam_with_loxapine_cognitive_impairment` | second_drug | CID100002789 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_clobazam_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_amisulpride_with_loxapine_cognitive_impairment` | second_drug | CID100002159 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_amisulpride_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_dexmedetomidine_with_loxapine_cognitive_impairment` | second_drug | CID100060612 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_dexmedetomidine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_sulpiride_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_sulpiride_with_loxapine_cognitive_impairment` | second_drug | RXCUI:10239 has no kept positive symptom at cognitive_impairment |
| `pair_positive_sulpiride_with_loxapine_somnolence_or_hypersomnia` | second_drug | RXCUI:10239 has no kept positive symptom at somnolence_or_hypersomnia |
| `pair_positive_loxapine_with_perampanel_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_perampanel_with_loxapine_cognitive_impairment` | second_drug | CID109924495 has no kept positive symptom at cognitive_impairment |
| `pair_positive_perampanel_with_loxapine_somnolence_or_hypersomnia` | second_drug | CID109924495 has no kept positive symptom at somnolence_or_hypersomnia |
| `pair_positive_chlorpromazine_with_loxapine_cognitive_impairment` | second_drug | CID100002726 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_chlorpromazine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_clomipramine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_clonazepam_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_clozapine_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_positive_chlordiazepoxide_with_loxapine_cognitive_impairment` | second_drug | CID100002712 has no kept positive symptom at cognitive_impairment |
| `pair_positive_loxapine_with_chlordiazepoxide_cognitive_impairment` | second_drug | CID100003964 has no kept positive symptom at cognitive_impairment |
| `pair_negative_sumatriptan_with_zolpidem_irritability_or_aggression` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_eletriptan_with_trihexyphenidyl_cognitive_impairment` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_diazepam_with_eletriptan_parkinsonism` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_acetaminophen_with_flurazepam_psychomotor_agitation` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_duloxetine_with_flurazepam_fatigue` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_lorazepam_with_rivastigmine_parkinsonism` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_sumatriptan_with_zopiclone_somnolence_or_hypersomnia` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_citalopram_with_clobazam_fatigue` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_clomipramine_with_gabapentin_psychomotor_agitation` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_imipramine_with_zolpidem_fatigue` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_citalopram_with_primidone_cognitive_impairment` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_sumatriptan_with_zolpidem_cognitive_impairment` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_loxapine_with_valproate_insomnia` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_methylphenidate_with_rizatriptan_insomnia` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_cabergoline_with_venlafaxine_cognitive_impairment` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_pyridostigmine_with_ropinirole_cognitive_impairment` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_flurazepam_with_sertraline_cognitive_impairment` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_flurazepam_with_sertraline_insomnia` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_clonazepam_with_rotigotine_fatigue` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_paliperidone_with_sumatriptan_somnolence_or_hypersomnia` | second_drug | expectation 'none': a negative control is not written as a label |
| `pair_negative_acetaminophen_with_alfentanil_irritability_or_aggression` | second_drug | expectation 'none': a negative control is not written as a label |
