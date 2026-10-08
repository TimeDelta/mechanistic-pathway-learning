# Node descriptor columns

Version 1, generated 2026-10-08 by experiments/write_descriptor_column_doc.py from the code at commit f60f5ea. Rerun the script after any descriptor table changes; it reads the column names and their placement from the tables.

Every node gets the same descriptor columns. Each node type fills its own block and is 0 in the others, so the encoder's input layer acts as one linear map per node type and nothing is learned per node (mechanistic_pathway_learning/graph/node_descriptors.py). A column is listed under a node type below when at least one node of that type has a non-zero value in it.

## Columns per node type in the confirmatory table (`data/processed/node_descriptors/full_neuronal_descriptors_brain_expression.parquet`)

### gene (92 columns)

- `protein_rrr_1`: component 1 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_2`: component 2 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_3`: component 3 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_4`: component 4 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_5`: component 5 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_6`: component 6 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_7`: component 7 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_8`: component 8 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_9`: component 9 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_10`: component 10 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_11`: component 11 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_12`: component 12 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_13`: component 13 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_14`: component 14 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_15`: component 15 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_16`: component 16 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_17`: component 17 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_18`: component 18 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_19`: component 19 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_20`: component 20 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_21`: component 21 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_22`: component 22 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_23`: component 23 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_24`: component 24 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_25`: component 25 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_26`: component 26 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_27`: component 27 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_28`: component 28 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_29`: component 29 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_30`: component 30 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_31`: component 31 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_32`: component 32 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_33`: component 33 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_34`: component 34 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_35`: component 35 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_36`: component 36 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_37`: component 37 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_38`: component 38 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_39`: component 39 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_40`: component 40 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_41`: component 41 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_42`: component 42 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_43`: component 43 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_44`: component 44 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_45`: component 45 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_46`: component 46 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_47`: component 47 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_48`: component 48 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_49`: component 49 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_50`: component 50 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_51`: component 51 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_52`: component 52 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_53`: component 53 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_54`: component 54 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_55`: component 55 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_56`: component 56 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_57`: component 57 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_58`: component 58 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_59`: component 59 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_60`: component 60 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_61`: component 61 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_62`: component 62 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_63`: component 63 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_rrr_64`: component 64 of the reduced-rank regression of ESM-2 (esm2_t12_35M_UR50D) embeddings onto EC, GO function, GO component, UniProt location and Pfam annotations; standardised
- `protein_has_protein_descriptors`: 1 when the node's protein has an ESM-2 embedding; 0 leaves the components at 0
- `gene_brain_gtex_brain_max`: largest GTEx v10 median TPM over the 13 brain tissues, log1p, standardised
- `gene_brain_gtex_other_max`: largest GTEx v10 median TPM over every other tissue, log1p, standardised
- `gene_brain_region_amygdala`: Human Protein Atlas consensus nTPM in the amygdala, log1p, standardised
- `gene_brain_region_basal_ganglia`: Human Protein Atlas consensus nTPM in the basal ganglia, log1p, standardised
- `gene_brain_region_cerebellum`: Human Protein Atlas consensus nTPM in the cerebellum, log1p, standardised
- `gene_brain_region_cerebral_cortex`: Human Protein Atlas consensus nTPM in the cerebral cortex, log1p, standardised
- `gene_brain_region_choroid_plexus`: Human Protein Atlas consensus nTPM in the choroid plexus, log1p, standardised
- `gene_brain_region_hippocampal_formation`: Human Protein Atlas consensus nTPM in the hippocampal formation, log1p, standardised
- `gene_brain_region_hypothalamus`: Human Protein Atlas consensus nTPM in the hypothalamus, log1p, standardised
- `gene_brain_region_medulla_oblongata`: Human Protein Atlas consensus nTPM in the medulla oblongata, log1p, standardised
- `gene_brain_region_midbrain`: Human Protein Atlas consensus nTPM in the midbrain, log1p, standardised
- `gene_brain_region_pons`: Human Protein Atlas consensus nTPM in the pons, log1p, standardised
- `gene_brain_region_spinal_cord`: Human Protein Atlas consensus nTPM in the spinal cord, log1p, standardised
- `gene_brain_region_thalamus`: Human Protein Atlas consensus nTPM in the thalamus, log1p, standardised
- `gene_brain_region_white_matter`: Human Protein Atlas consensus nTPM in the white matter, log1p, standardised
- `gene_brain_class_astrocyte`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class astrocyte, log1p, standardised
- `gene_brain_class_cortical_interneuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class cortical interneuron, log1p, standardised
- `gene_brain_class_ependymal_choroid`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class ependymal choroid, log1p, standardised
- `gene_brain_class_excitatory_neuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class excitatory neuron, log1p, standardised
- `gene_brain_class_immune`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class immune, log1p, standardised
- `gene_brain_class_medium_spiny_neuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class medium spiny neuron, log1p, standardised
- `gene_brain_class_oligodendrocyte_lineage`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class oligodendrocyte lineage, log1p, standardised
- `gene_brain_class_other_inhibitory_neuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class other inhibitory neuron, log1p, standardised
- `gene_brain_class_other_neuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class other neuron, log1p, standardised
- `gene_brain_class_vascular`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class vascular, log1p, standardised
- `gene_brain_class_dopaminergic_neuron`: dopaminergic cluster 395 of Siletti et al. 2023 (CELLxGENE counts put on the HPA nCPM scale), log1p, standardised
- `gene_brain_has_expression`: 1 when the gene (for a reaction: a gene in its rule) is in the Human Protein Atlas region table

### membrane_potential (0 columns)

No descriptor column; all 0.

### metabolite (10 columns)

- `metabolite_log_molecular_weight`: log molecular weight (RDKit MolWt), standardised
- `metabolite_logp`: Crippen logP, standardised
- `metabolite_net_charge`: net charge in Human-GEM (SBML fbc:charge, near pH 7.3), standardised; known for every metabolite
- `metabolite_log_polar_surface_area`: log(1 + topological polar surface area), standardised
- `metabolite_log_hydrogen_bond_donors`: log(1 + hydrogen-bond donors), standardised
- `metabolite_log_hydrogen_bond_acceptors`: log(1 + hydrogen-bond acceptors), standardised
- `metabolite_log_rotatable_bonds`: log(1 + rotatable bonds), standardised
- `metabolite_log_rings`: log(1 + rings), standardised
- `metabolite_has_structure`: 1 when a complete SMILES gave the properties above; 0 leaves them at the mean
- `metabolite_partial_structure`: 1 when the SMILES has R-group atoms (*), e.g. the acyl chain of an acyl-CoA

### protein_entity (0 columns)

No descriptor column; all 0.

### reaction (35 columns)

- `reaction_ec_class_1`: 1 when an annotated EC number of the reaction starts with 1 (oxidoreductases)
- `reaction_ec_class_2`: 1 when an annotated EC number of the reaction starts with 2 (transferases)
- `reaction_ec_class_3`: 1 when an annotated EC number of the reaction starts with 3 (hydrolases)
- `reaction_ec_class_4`: 1 when an annotated EC number of the reaction starts with 4 (lyases)
- `reaction_ec_class_5`: 1 when an annotated EC number of the reaction starts with 5 (isomerases)
- `reaction_ec_class_6`: 1 when an annotated EC number of the reaction starts with 6 (ligases)
- `reaction_ec_class_7`: 1 when an annotated EC number of the reaction starts with 7 (translocases)
- `reaction_has_ec`: 1 when the reaction has an annotated EC number
- `reaction_brain_gtex_brain_max`: largest GTEx v10 median TPM over the 13 brain tissues, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_gtex_other_max`: largest GTEx v10 median TPM over every other tissue, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_amygdala`: Human Protein Atlas consensus nTPM in the amygdala, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_basal_ganglia`: Human Protein Atlas consensus nTPM in the basal ganglia, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_cerebellum`: Human Protein Atlas consensus nTPM in the cerebellum, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_cerebral_cortex`: Human Protein Atlas consensus nTPM in the cerebral cortex, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_choroid_plexus`: Human Protein Atlas consensus nTPM in the choroid plexus, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_hippocampal_formation`: Human Protein Atlas consensus nTPM in the hippocampal formation, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_hypothalamus`: Human Protein Atlas consensus nTPM in the hypothalamus, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_medulla_oblongata`: Human Protein Atlas consensus nTPM in the medulla oblongata, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_midbrain`: Human Protein Atlas consensus nTPM in the midbrain, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_pons`: Human Protein Atlas consensus nTPM in the pons, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_spinal_cord`: Human Protein Atlas consensus nTPM in the spinal cord, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_thalamus`: Human Protein Atlas consensus nTPM in the thalamus, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_region_white_matter`: Human Protein Atlas consensus nTPM in the white matter, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_astrocyte`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class astrocyte, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_cortical_interneuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class cortical interneuron, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_ependymal_choroid`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class ependymal choroid, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_excitatory_neuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class excitatory neuron, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_immune`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class immune, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_medium_spiny_neuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class medium spiny neuron, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_oligodendrocyte_lineage`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class oligodendrocyte lineage, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_other_inhibitory_neuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class other inhibitory neuron, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_other_neuron`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class other neuron, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_vascular`: Human Protein Atlas single-nucleus nCPM, largest over the cluster types of the class vascular, log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_class_dopaminergic_neuron`: dopaminergic cluster 395 of Siletti et al. 2023 (CELLxGENE counts put on the HPA nCPM scale), log1p, standardised; from the gene rule (minimum over and, maximum over or)
- `reaction_brain_has_expression`: 1 when the gene (for a reaction: a gene in its rule) is in the Human Protein Atlas region table; from the gene rule (minimum over and, maximum over or)

## Where the split graphs differ (`data/processed/node_descriptors/full_neuronal_split_descriptors_brain_expression.parquet`)

With genes and proteins split (docs/gene_protein_split.md), the protein block sits on the protein nodes and the gene nodes keep the brain-expression block. Node types whose columns differ from the table above:

- gene: 27 columns, block `gene_brain`
- protein: 65 columns, block `protein`

## Structural features before the descriptors

The encoders read the descriptors after the structural features of ExperimentData.structural_node_features() (experiments/run_main_model.py, node_feature_matrix). The code builds these by position, without names; the names below are given here for reading, in order, for `data/processed/graph_full_neuronal`:

1. `type_gene`
2. `type_membrane_potential`
3. `type_metabolite`
4. `type_protein_entity`
5. `type_reaction`
6. `compartment_c`
7. `compartment_e`
8. `compartment_g`
9. `compartment_i`
10. `compartment_l`
11. `compartment_m`
12. `compartment_n`
13. `compartment_r`
14. `compartment_v`
15. `compartment_x`
16. `log1p_degree`
17. `is_currency`
18. `is_transport`
19. `is_reversible`
20. `log1p_gtex_brain_median_tpm_max`
21. `brain_expressed`

## Blocks

experiments/run_main_model.py --drop-descriptor-blocks leaves out whole blocks, by column prefix (DESCRIPTOR_BLOCK_PREFIXES in node_descriptors.py):

- `metabolite`: 10 columns
- `reaction_ec`: 8 columns
- `reaction_brain`: 27 columns
- `protein`: 65 columns
- `gene_brain`: 27 columns

## Placement in every descriptor table

Cells give how many of a block's columns are non-zero on at least one node of the type (blank: none).

### `data/processed/node_descriptors/full_neuronal_descriptors_brain_expression.parquet`

Used by: the confirmatory configurations (FULL_GRAPH_NODE_PROPERTIES in experiments/run_main_model_batch.py). Graph: `data/processed/graph_full_neuronal`. 137 columns. SHA-256 `14033cb35ad97098`.

| node type | nodes | metabolite | reaction_ec | reaction_brain | protein | gene_brain |
|---|---|---|---|---|---|---|
| gene | 12810 |  |  |  | all | all |
| membrane_potential | 1 |  |  |  |  |  |
| metabolite | 8580 | all |  |  |  |  |
| protein_entity | 1681 |  |  |  |  |  |
| reaction | 13793 |  | all | all |  |  |

### `data/processed/node_descriptors/full_neuronal_descriptors_brain_expression_resolved.parquet`

Used by: merged full graph with the protein block resolved per entry (twin of the full split). Graph: `data/processed/graph_full_neuronal`. 137 columns. SHA-256 `9fee213fb2c0716e`.

| node type | nodes | metabolite | reaction_ec | reaction_brain | protein | gene_brain |
|---|---|---|---|---|---|---|
| gene | 12810 |  |  |  | all | all |
| membrane_potential | 1 |  |  |  |  |  |
| metabolite | 8580 | all |  |  |  |  |
| protein_entity | 1681 |  |  |  |  |  |
| reaction | 13793 |  | all | all |  |  |

### `data/processed/node_descriptors/full_neuronal_split_descriptors_brain_expression.parquet`

Used by: full graph with genes and proteins split. Graph: `data/processed/graph_full_neuronal_split`. 137 columns. SHA-256 `27ed5e89ccb369d9`.

| node type | nodes | metabolite | reaction_ec | reaction_brain | protein | gene_brain |
|---|---|---|---|---|---|---|
| gene | 12810 |  |  |  |  | all |
| membrane_potential | 1 |  |  |  |  |  |
| metabolite | 8580 | all |  |  |  |  |
| protein | 12709 |  |  |  | all |  |
| protein_entity | 1681 |  |  |  |  |  |
| reaction | 13793 |  | all | all |  |  |

### `data/processed/node_descriptors/slice_descriptors_brain_expression.parquet`

Used by: slice. Graph: `data/processed/graph`. 137 columns. SHA-256 `87b2f34ad6bca3f3`.

| node type | nodes | metabolite | reaction_ec | reaction_brain | protein | gene_brain |
|---|---|---|---|---|---|---|
| gene | 2848 |  |  |  | all | all |
| metabolite | 8460 | all |  |  |  |  |
| reaction | 12877 |  | all | all |  |  |

### `data/processed/node_descriptors/slice_descriptors_brain_expression_resolved.parquet`

Used by: slice with the protein block resolved per entry (twin of the split slice). Graph: `data/processed/graph`. 137 columns. SHA-256 `e35af6f2b98b5489`.

| node type | nodes | metabolite | reaction_ec | reaction_brain | protein | gene_brain |
|---|---|---|---|---|---|---|
| gene | 2848 |  |  |  | all | all |
| metabolite | 8460 | all |  |  |  |  |
| reaction | 12877 |  | all | all |  |  |

### `data/processed/node_descriptors/slice_split_descriptors_brain_expression.parquet`

Used by: slice with genes and proteins split. Graph: `data/processed/graph_split`. 137 columns. SHA-256 `fef2a1b250760a87`.

| node type | nodes | metabolite | reaction_ec | reaction_brain | protein | gene_brain |
|---|---|---|---|---|---|---|
| gene | 2848 |  |  |  |  | all |
| metabolite | 8460 | all |  |  |  |  |
| protein | 2837 |  |  |  | all |  |
| reaction | 12877 |  | all | all |  |  |

### `data/processed/graph/node_descriptors.parquet`

Used by: slice configurations without brain expression. Graph: `data/processed/graph`. 83 columns. SHA-256 `5390d78f662d1387`.

| node type | nodes | metabolite | reaction_ec | protein |
|---|---|---|---|---|
| gene | 2848 |  |  | all |
| metabolite | 8460 | all |  |  |
| reaction | 12877 |  | all |  |

### `data/processed/graph_full_neuronal/node_descriptors.parquet`

Used by: full graph without brain expression. Graph: `data/processed/graph_full_neuronal`. 83 columns. SHA-256 `60afb642110187a7`.

| node type | nodes | metabolite | reaction_ec | protein |
|---|---|---|---|---|
| gene | 12810 |  |  | all |
| membrane_potential | 1 |  |  |  |
| metabolite | 8580 | all |  |  |
| protein_entity | 1681 |  |  |  |
| reaction | 13793 |  | all |  |

### `data/processed/graph_split/node_descriptors.parquet`

Used by: split slice without brain expression. Graph: `data/processed/graph_split`. 83 columns. SHA-256 `1339881da3cb767e`.

| node type | nodes | metabolite | reaction_ec | protein |
|---|---|---|---|---|
| gene | 2848 |  |  |  |
| metabolite | 8460 | all |  |  |
| protein | 2837 |  |  | all |
| reaction | 12877 |  | all |  |

### `data/processed/graph_full_neuronal_split/node_descriptors.parquet`

Used by: full split without brain expression. Graph: `data/processed/graph_full_neuronal_split`. 83 columns. SHA-256 `fd34a8c4d7ba38d5`.

| node type | nodes | metabolite | reaction_ec | protein |
|---|---|---|---|---|
| gene | 12810 |  |  |  |
| membrane_potential | 1 |  |  |  |
| metabolite | 8580 | all |  |  |
| protein | 12709 |  |  | all |
| protein_entity | 1681 |  |  |  |
| reaction | 13793 |  | all |  |

## Changes

| version | date | change |
|---|---|---|
| 1 | 2026-10-08 | First version: the 137 descriptor columns of the brain-expression tables and the 83 of the graph directories' own tables. |
