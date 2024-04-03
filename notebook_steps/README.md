# Methods

This section explains the methods for the paper. <br>
We used CPTAC's data on ccrcc, luad, hnscc, lscc, and ucec tissues.<br>
Proteomic data was obtained from "umich" database.<br>
Transcriptomic data was obained from "washu" database.<br>

## Table of Contents

- [Δ_corr](#Δ_corr)
- [Transmutation Effects](#Transmutation-Effects)
- [Proteomic and Transcriptomic Tables](#Proteomic-and-Transcriptomic-Differential-Expression)
- [Cancer Stages and Correlations](#Cancer-Stages-and-Correlations)
- [Set Enrichment](#Set-Enrichment)

## Δ_corr


We calculated Spearman correlations between rna and protein abundances for all genes found in the tissues for both healthy and diseased samples. Δ_corr is the the difference between the healthy and diseased correlations. A permutation label-swap test was done for determining the significance of each gene's Δ_corr. A table showing all the Δ_corr and p-values can be found [here](./data/delta_correlation_df_with_significance.csv)<br>
The correlation scripts and tables used to create that table can be found in this folder [here](./data/Scripts_to_Make_Cancer_Delta_Corr_and_P_value_Dataframe). <br>
That table is used to produce the figures 2, 4, 5, S2, and S4. <br>

## Transmutation Effects

The transmutation effects scripts and tables can be found in the data folder [here](./data/Scripts_to_make_transmutation_effects_dataframes). <br>
You can find a table with the transmutation data in this [csv](./data/transmutation_df.csv)
This was used to create figure 3.<br>

## Proteomic and Transcriptomic Differential Expression

Proteomic and Transcriptomic data was obtained from CPTAC and the script can be found [here](./data/Make_Cancer_stages_correlations_df.ipynb).<br>
Proteomic differential expression can be found in this [csv](./data/Proteomics_differential_expression_df.csv)<br>
Transcriptomic differential expression can be found in this [csv](./data/Transcriptomics_differential_expression_df.csv)<br>
This was used to create figure 2, S2.<br>

## Cancer Stages and Correlations

Cancer stages were obtained from CPTAC and the script can be found [here](./data/Scripts_to_make_transmutation_effects_dataframes). <br>
We colapsed the different stage name into 1-4 stages. (Ex: ["IA", "Stage 1", "Stage 1B", "Stage IB", "Stage IA" or "Stage IA3"] into "Stage I")<br>
We then calculated the mRNA/Correlation in each of the stages and saved it as a [csv](./data/Cancer_stages_correlations.csv)<br>
This was used to create figure 1.<br>

A complete tumor-normal Correlation script can be found [here](./data/Make_Tumor-Normal_Correlation_Dataframe.ipynb)<br>
The table is saved in this [csv](./data/tumor_normal_correlation_df.csv)<br>
This was used to create figure 4.

## Set Enrichment

The Script calculate the Set Enrichment can be found [here](./data/Scripts_to_make_transmutation_effects_dataframes). <br>
We used the KEGG source for the enrichment. You can use the script to find tissues with altered pathways found using significant Δ_corr genes.<br>
This script's code is used to generate figure 5.