# EvoResist
**Evolution-Guided Prioritization of Drug Resistance Mutations Enhances Molecular Prediction of Tuberculosis Drug Susceptibility**

---

## Abstract

### Background
Drug susceptibility testing (DST) is crucial in guiding the selection of effective treatment regimens for tuberculosis. Since culture-based drug susceptibility testing typically requires 4-8 weeks, molecular diagnostics using mutational markers provide an alternative approach for rapidly predicting phenotypic drug resistance. However, prioritizing mutations to enhance phenotypic DST (pDST) prediction remains challenging.

### Methods
We propose leveraging convergent evolution to prioritize mutations in known drug-resistance genes. We analyzed whole-genome sequencing (WGS) data from 106,069 Mycobacterium tuberculosis (Mtb) isolates to identify mutations that have been more frequently targeted in the bacterial population. We applied the following criteria to prioritize drug-resistant mutations: 1) Mutations with >= 4 convergent events in drug-resistance genes; 2) Frameshift indels and premature stop mutations within genes known to confer drug resistance via loss-of-function; and 3) Exclusion of mutations that occurred in ancestral nodes (phylogenetic mutations). We compared the performance of the prioritized mutations to the Group 1 and Group 2 mutations from the WHO Catalogue of drug-resistant mutations in three previously published datasets: WHO dataset (36,270 isolates), a Chinese dataset (9,051 isolates), and a Spanish dataset (701 isolates).

### Findings
We identified 1,827 convergent mutations, along with 4,332 indels or premature stop codons, for predicting phenotypic drug resistance (median number of convergent events: 6; interquartile range [IQR]: 4-12). Applying these mutations to predict drug susceptibility in the WHO dataset demonstrated significant improvements in sensitivity (1.2%-12.9%) for predicting phenotypic resistance in first-line drugs, and 0.5%-34.4% for second-line or new drugs compared to the WHO Catalogue Group 1 mutations, at a minimum (0.3%-7.2%) reduction in specificity. Significant improvement was also observed in two separate datasets from China and Spain. Additionally, we showed that both biological and technical factors can affect the specificity of certain mutations. For example, lineage background effects specificity of certain drug-resistance mutations (e.g., *inhA* -779G>T only exhibited low specificity in lineage 4); the low specificity of some mutations (e.g., *gyrA* Ala90Val) was due to quality of certain data source, cross-referenced by the similarly poor performance of well-established resistance markers katG S315T and rpoB S450L in those same datasets. 

### Interpretation
We present an evolution-based, phenotype-independent framework for prioritizing drug-resistance mutations and demonstrate that convergent mutations in resistance genes are promising markers for predicting phenotypic drug resistance. These findings underscore the value of incorporating evolutionary signatures derived from population genomics to guide the development of more targeted and effective molecular diagnostics.

---

## Repository Contents

This repository contains scripts associated with the study and Supplementary Table 5:

- **`Scripts/`** – Analysis and processing scripts  
  - `EvoResist_pipeline.sh`: Pipeline of EvoResist.
  - `Infer_ancestral_mutations.pl`: Infer the ancestral mutations.
  - `Simulation_JC96.py`: Run 100 simulations of mutational events across the whole genome using the Jukes–Cantor 1996 (JC96) model.
  - `Permutation_test.py`: Perform a permutation test to evaluate whether the difference between the observed event count and 100 simulated event counts is statistically significant.
  - `Threshold_selection.R`: Select the optimized mutational event threshold.
  - `R_S_count.py`: For each variant, count the number of isolates carrying the variant that have pDST results classified as Resistant (R) or Sensitive (S).
  - `Sens_Spec_change.py`: For each variant, compute the change in prediction sensitivity and specificity compared to models using mutations from WHO G1 or WHO G1 + G2.
  - `Lineage_specificity_difference_test.R`: Test whether prediction specificity differs significantly among Lineages 1, 2, 3, and 4.
  - `Variant_difference_test.R`: Test if the prevance of variants within primary drug-resistant genes differ statistically significant among Lineage 1, 2, 3, and 4.

- **`TableS5_Mutations_and_their_event_number_across_genome.tsv`** - Supplementary table 5.
