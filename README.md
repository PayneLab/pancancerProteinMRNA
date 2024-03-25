# Activating cancer hallmarks through changes in mRNA/protein regulation

As a diverse family of diseases, cancer is unified by a set of common dysfunctions, such as limitless growth potential and an insensitivity to anti-growth signals. These shared overarching biological processes have been termed the hallmarks of cancer. To better understand the root cause of cellular dysregulation, intense molecular characterization of tumors have utilized DNA, RNA and protein measurement techniques to produce proteogenomic data. In large cancer cohort studies, genomic and proteogenomic data have frequently identified many cancer hallmarks, including cell cycle and cell signaling. However, altered metabolism - a known cancer hallmark - is not as clearly identified in mutation screens of simple differential expression analyses. Here we introduce a new computational method to identify changes in cellular regulation by focusing on the mRNA/protein relationship. We create a metric, Δ_corr, to capture when the mRNA/protein correlation changes significantly between tumor and normal tissues, and show that it is distinct from differential expression and also not associated with DNA mutation profiles. Our method clearly highlights altered metabolic pathways across multiple tumor types. Δ_corr gives researchers a new perspective on dysfunction of tumor cells and introduces a novel method for proteogenomic data integration.

## Table of Contents

- [Requirements](#requirements)
- [Δ_corr](#Δ_corr)
- [Figures](#figures)


## Requirements

- pip install -r requirements.txt
- pip install -e ./pcprutils

## Δ_corr

Δ_corr is the difference between the healthy and diseased correlations of RNA and protein abundances. The correlation scripts and tables can be found in the data folder [here](./notebook_steps/data/Scripts_to_Make_Cancer_Delta_Corr_and_P_value_Dataframe). <br>
The csv files contain Δ_corr tables calculated with the scripts in the folder. <br>
The scripts were run in BYU's supercomputer. To run this localy: <br>
Run download_pancan.py and pass as argument a token for downloading non published datasets.<br>
Use the command in the found in the cancer bash files to calculate the tables. However, note that the it may take some time to finish the run.

## Transmutation Effects

The transmutation effects scripts and tables can be found in the data folder [here](./notebook_steps/data/Scripts_to_make_transmutation_effects_dataframes). <br>
The csv files contain the transmutation effect tables calculated with the scripts in the folder. <br>
The scripts were run in BYU's supercomputer. To run this localy: <br>
Run download_pancan.py and pass as argument a token for downloading non published datasets.<br> 
Use the command in the cancer bash files. Again, note that the it may take some time to finish running.

## Figures

All figures are produced by the scripts found [here](./notebook_step/Figures). <br>
All the scripts that produce the images can be found [here](./notebook_steps). <br>

## pcprutils
Collection of useful functions to load cancers, fetch proteomic and transcriptomic data, and calculate correlations and permutation scores.<br>


