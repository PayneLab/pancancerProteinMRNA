import pandas as pd
import cptac
import statistics as st
import cptac.utils as ut
import numpy as np
from scipy import stats
import scipy
import warnings
import csv
import sys
import copy


warnings.filterwarnings('ignore')

input_cancer_type = sys.argv[1]
mutated_gene = sys.argv[2]
num_permutations = int(sys.argv[3])

if input_cancer_type == "ccrcc":
    #cptac.download('ccrcc')
    cancer = cptac.Ccrcc()
elif input_cancer_type == "en":
    cancer = cptac.Endometrial()
elif input_cancer_type == "luad":
    cancer = cptac.Luad()
elif input_cancer_type == "hnscc":
    cancer  = cptac.Hnscc()
elif input_cancer_type == "lscc":
    cancer = cptac.Lscc()
    
def permute(df,original_correlation,Type, label_1, label_2, column_one, column_two, permutation_times):
    permutation_list = []
    permu_df = copy.deepcopy(df)

    for i in range(permutation_times):
        permu_df[Type] = np.random.permutation(permu_df[Type])
        permu_is_label_1 = permu_df[Type] == label_1
        permu_is_label_2 = permu_df[Type] == label_2
        label_1_correlation,label_1_pval = scipy.stats.pearsonr(permu_df[permu_is_label_1][column_one], permu_df[permu_is_label_1][column_two])
        label_2_correlation,label_2_pval = scipy.stats.pearsonr(permu_df[permu_is_label_2][column_one], permu_df[permu_is_label_2][column_two])
        delta = label_1_correlation - label_2_correlation
        permutation_list.append(delta)

    z_score = (original_correlation - np.mean(permutation_list)) / np.std(permutation_list)
    p_val = scipy.stats.norm.sf(abs(z_score))*2
    return p_val

def get_omics_df(cancer):
    transcriptomics_df = cancer.get_transcriptomics(tissue_type='tumor')
    proteomics_df = cancer.get_proteomics(tissue_type='tumor')
    if isinstance(proteomics_df.columns, pd.MultiIndex):
        proteomics_df = proteomics_df.droplevel('Database_ID', axis = 1)
    if isinstance(transcriptomics_df.columns, pd.MultiIndex):
        transcriptomics_df = transcriptomics_df.droplevel('Database_ID', axis = 1)
    proteomics_df['patient_ID'] = proteomics_df.index
    transcriptomics_df['patient_ID'] = transcriptomics_df.index
    transcriptomics_df = transcriptomics_df.melt(id_vars='patient_ID', var_name = 'gene', value_name='transcriptomics')
    proteomics_df = proteomics_df.melt(id_vars='patient_ID', var_name = 'gene', value_name='proteomics')
    mutation_df = cancer.get_somatic_mutation()
    print(pd.unique(mutation_df.Mutation))
    mutation_df = mutation_df[mutation_df.Gene == mutated_gene]
    print(mutation_df)
    omics_df = pd.merge(transcriptomics_df, proteomics_df, how = 'inner')
    omics_df['mutation_status'] = omics_df.patient_ID.isin(mutation_df.index)
    omics_df = omics_df.dropna()
    print(omics_df)
    return omics_df

def get_corr_df(omics_df):
    mut_corrs = []
    mut_p_vals = []
    non_mut_corrs = []
    non_mut_p_vals = []
    corr_diffs = []
    corr_diff_pvals = []
    genes = []
    for gene in pd.unique(omics_df.gene):
        df = omics_df[omics_df.gene == gene]
        mut_df = df[df.mutation_status == True]
        non_mut_df = df[df.mutation_status == False]
        if len(mut_df) < 4 or len(non_mut_df) < 4:
            continue
        mut_r, mut_p = stats.pearsonr(mut_df.transcriptomics, mut_df.proteomics)
        non_mut_r, non_mut_p = stats.pearsonr(non_mut_df.transcriptomics, non_mut_df.proteomics)
        mut_corrs.append(mut_r)
        mut_p_vals.append(mut_p)
        non_mut_corrs.append(non_mut_r)
        non_mut_p_vals.append(non_mut_p)
        corr_diff = mut_r - non_mut_r
        corr_diffs.append(corr_diff)
        diff_p_val = permute(df, corr_diff,'mutation_status', True, False, 'transcriptomics', 'proteomics', num_permutations)
 
        corr_diff_pvals.append(diff_p_val)
        genes.append(gene)
    correlation_df = pd.DataFrame({'gene': genes, 'mutated_correlation': mut_corrs, 'non_mutated_correlation': non_mut_corrs,
                                   'non_mutated_p_vals': non_mut_p_vals, 'mutated_p_vals': mut_p_vals,
                                   'delta_correlation': corr_diffs, 'delta_correlation_pval': corr_diff_pvals})
    return correlation_df

omics_df = get_omics_df(cancer)
omics_df['cancer'] = [input_cancer_type] * len(omics_df)
corr_df = get_corr_df(omics_df)
corr_df['cancer'] = [input_cancer_type] * len(corr_df)
corr_df['sig_mut_pvals'] = corr_df.mutated_p_vals < 0.05 / len(corr_df)
corr_df['sig_wt_pvals'] = corr_df.non_mutated_p_vals < 0.05 / len(corr_df)
corr_df['sig_delta_corr_pval'] = corr_df.delta_correlation_pval < 0.05 / len(corr_df)
df_name = input_cancer_type + "_" + mutated_gene + "_mutation_effects_permutation.csv"
corr_df.to_csv(df_name)