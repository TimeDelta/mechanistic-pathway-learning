# What a drug's own label says about acting centrally

Read from openFDA by `experiments/read_label_central_action.py`, pharmacology sections only (mechanism of action, clinical pharmacology, pharmacodynamics, description). The adverse-reaction, warning and indication sections are deliberately not read: a drug admitted for the psychiatric events on its label would select the drug set on the outcome the study predicts. A positive verdict is a shortlist to read, not a decision: class boilerplate makes every beta blocker positive, and a combination-product label can put a drug's only match on a sentence about the other ingredient.

56 drugs read.

| verdict | drugs |
| --- | --- |
| no such sentence | 24 |
| states a central action | 17 |
| no label section read | 7 |
| states that it does not reach the brain | 6 |
| only an adverse-event sentence | 2 |

## states a central action

| drug | section | quote | openFDA set |
| --- | --- | --- | --- |
| bisoprolol | clinical_pharmacology | Factors which may be involved include: 1) Decreased cardiac output, 2) Inhibition of renin release by the kidneys, 3) Diminution of tonic sympathetic outflow from the vasomotor centers in the brain. | 084308e4-3b8b-40d6-acbf-d775772a8a73 |
| codeine | clinical_pharmacology | The analgesic properties of codeine have been speculated to come from its conversion to morphine, although the exact mechanism of analgesic action remains unknown. 12.2 Pharmacodynamics Effects on the Central Nervous System: Codeine produces respiratory depression by direct action on brain stem respiratory centers. | 010905f9-3bcb-4b50-9fe8-a3ad0010f14c |
| cyclobenzaprine | clinical_pharmacology | It is ineffective in muscle spasm due to central nervous system disease. | 00b20a8c-debd-e930-e063-6394a90a7ed0 |
| felodipine | clinical_pharmacology | Animal studies have demonstrated that felodipine crosses the blood-brain barrier and the placenta. | 05feb0a1-8a0c-4d8a-9b86-5e700ea6fd64 |
| granisetron | clinical_pharmacology | Serotonin receptors of the 5-HT 3 type are located peripherally on vagal nerve terminals and centrally in the chemoreceptor trigger zone of the area postrema. | 2758eb4d-2fe0-4f74-b3f5-b4903009e97c |
| hydrocodone | mechanism_of_action | However, specific CNS opioid receptors for endogenous compounds with opioid-like activity have been identified throughout the brain and spinal cord and are thought to play a role in the analgesic effects of this drug. | 025eb0f7-d002-4ce0-87b8-71461ab18356 |
| indomethacin | clinical_pharmacology | Indomethacin has been found to cross the blood-brain barrier and the placenta, and appears in breast milk. | 009097a5-2c1e-4f5d-8054-896cf896cb3d |
| metformin | clinical_pharmacology | In mice, no glipizide or metabolites were detectable autoradiographically in the brain or spinal cord of males or females, nor in the fetuses of pregnant females. | 03d9ce02-0559-494c-a53d-de9b9a64b606 |
| metoprolol | clinical_pharmacology | Metoprolol is also known to cross the blood brain barrier following oral administration and CSF concentrations close to that observed in plasma have been reported. | 00940cc5-d2eb-4841-9138-de97d7b1c674 |
| nebivolol | mechanism_of_action | Possible factors that may be involved include: (1) decreased heart rate, (2) decreased myocardial contractility, (3) diminution of tonic sympathetic outflow to the periphery from cerebral vasomotor centers,(4) suppression of renin activity and (5) vasodilation and decreased peripheral vascular resistance. | 04a64dad-6922-4737-ab46-98c185782b38 |
| nimodipine | clinical_pharmacology | In animal experiments, nimodipine had a greater effect on cerebral arteries than on arteries elsewhere in the body perhaps because it is highly lipophilic, allowing it to cross the blood-brain barrier; concentrations of nimodipine as high as 12.5 ng/mL have been detected in the cerebrospinal fluid of nimodipine-treated subarachnoid hemorrhage (SAH) patients. | a323ebd8-ef1b-439b-93d1-ba73a5b0127d |
| omeprazole | clinical_pharmacology | Other Effects Systemic effects of omeprazole in the CNS, cardiovascular and respiratory systems have not been found to date. | 0344268f-d8e6-48d5-b1e1-f750ed1d8a52 |
| ondansetron | mechanism_of_action | Serotonin receptors of the 5-HT 3 type are present both peripherally on vagal nerve terminals and centrally in the chemoreceptor trigger zone of the area postrema. | 00327696-c496-4c83-a63e-9e29fd6246d4 |
| palonosetron | mechanism_of_action | Cancer chemotherapy may be associated with a high incidence of nausea and vomiting, particularly when certain agents, such as cisplatin, are used. 5-HT 3 receptors are located on the nerve terminals of the vagus in the periphery and centrally in the chemoreceptor trigger zone of the area postrema. | 0a781796-bb74-4c5f-8d26-9ce520b818f7 |
| pindolol | clinical_pharmacology | The mechanism of the antihypertensive effects of beta-blocking agents has not been established, but several mechanisms have been postulated: 1) an effect on the central nervous system resulting in a reduced sympathetic outflow to the periphery, 2) competitive antagonism of catecholamines at peripheral (especially cardiac) adrenergic receptor sites, leading to decreased cardiac output, 3) an inhibition of renin release. | 0daed7b3-6017-4611-8c3a-2f390470f8bc |
| propranolol | clinical_pharmacology | Factors that may contribute to the antihypertensive action include: (1) decreased cardiac output, (2) inhibition of renin release by the kidneys, and (3) diminution of tonic sympathetic nerve outflow from vasomotor centers in the brain. | 013a485a-4227-46ae-b19f-c550302fcb96 |
| timolol | clinical_pharmacology | Possible mechanisms of action include reduction in cardiac output, reduction in plasma renin activity, and a central nervous system sympatholytic action. | 0bc40b2c-65eb-4095-87e6-4752d5b19a3a |

## states that it does not reach the brain

| drug | section | quote | openFDA set |
| --- | --- | --- | --- |
| etoposide | clinical_pharmacology | Although it is detectable in CSF and intracerebral tumors, the concentrations are lower than in extracerebral tumors and in plasma. | 0dddbe44-b4f4-431f-9011-8a54abbd4cbf |
| flumazenil | clinical_pharmacology | Flumazenil does not antagonize the central nervous system effects of drugs affecting GABA-ergic neurons by means other than the benzodiazepine receptor (including ethanol, barbiturates, or general anesthetics) and does not reverse the effects of opioids. | 1855bb03-392f-4a0d-8361-24fc2a0e635e |
| labetalol | clinical_pharmacology | Only negligible amounts of the drug crossed the blood-brain barrier in animal studies. | 04f9dcdf-4023-48ee-9e7d-ca9366d9db72 |
| montelukast | clinical_pharmacology | Studies in rats with radiolabeled montelukast indicate minimal distribution across the blood-brain barrier. | 07d5fbea-b183-41c9-a20f-30da1f44ab1b |
| nadolol | clinical_pharmacology | In dogs, only minimal amounts of nadolol were detected in the brain relative to amounts in blood and other organs and tissues. | 0a2acdc0-aab6-462e-a377-67a0f3d55ec6 |
| sotalol | clinical_pharmacology | Sotalol crosses the blood brain barrier poorly. | 0fd4afc7-dbf8-46be-9657-aa83bbb2c026 |

## only an adverse-event sentence

ibuprofen, theophylline


## no such sentence

acebutolol, atenolol, betaxolol, carvedilol, celecoxib, diclofenac, diltiazem, flurbiprofen, isotretinoin, isradipine, ketoprofen, mecamylamine, metoclopramide, naproxen, nicardipine, nifedipine, nisoldipine, pioglitazone, pitavastatin, prazosin, ribavirin, sulfasalazine, verapamil, warfarin


## no label section read

celiprolol, dextromethorphan, etoricoxib, loratadine, oxprenolol, penbutolol, valdecoxib

