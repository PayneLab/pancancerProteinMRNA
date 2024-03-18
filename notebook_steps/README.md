# Methods

This section explains the methods for the paper. <br>
We used CPTAC's data on ccrcc, luad, hnscc, lscc, and ucec tissues.

## Table of Contents

- [Δ_corr](#Δ_corr)
- [Transmutation Effects](#transmutation effects)
- [Proteomic and Transcriptomic Tables](#proteomic and transcriptomic tables)


## Δ_corr


We calculated Spearman correlations between rna and protein abundances for all genes found in the tissues for both healthy and diseased samples. Δ_corr is the the difference between the healthy and diseased correlations. A permutation label-swap test was done for determining the significance of each gene's Δ_corr. A table showing all the Δ_corr and p-values can be found [here](./data/delta_correlation_df_with_significance.csv)<br>
The correlation scripts and tables used to create that table can be found in this folder [here](./data/Scripts_to_Make_Cancer_Delta_Corr_and_P_value_Dataframe). <br>
That table is used to produce the figures 2, 4, 5, S2, and S4. <br>

## Transmutation Effects

The transmutation effects scripts and tables can be found in the data folder [here](./data/Scripts_to_make_transmutation_effects_dataframes). <br>


## Proteomic and Transcriptomic Tables


