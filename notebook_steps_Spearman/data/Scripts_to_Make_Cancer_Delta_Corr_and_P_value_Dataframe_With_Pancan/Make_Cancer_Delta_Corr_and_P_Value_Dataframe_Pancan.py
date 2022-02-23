import cptac.pancan as pc
import scipy.stats as stats
import numpy as np
import copy
import cptac.utils as ut
import statsmodels.stats.multitest as ssm
import pandas as pd
import warnings
import sys, os

def get_prot_trans_df(cancer):

    prot_normal_df = cancer.get_proteomics("pdc", 'normal')
    if isinstance(prot_normal_df.columns, pd.MultiIndex):
        prot_normal_df = ut.reduce_multiindex(df= prot_normal_df, levels_to_drop = 'Database_ID', quiet=True)

    trans_normal_df = cancer.get_transcriptomics("washu", 'normal')
    if isinstance(trans_normal_df.columns, pd.MultiIndex):
        trans_normal_df = ut.reduce_multiindex(df = trans_normal_df, levels_to_drop='Database_ID', quiet=True)

    prot_normal_df['Name'] = prot_normal_df.index
    trans_normal_df['Name'] = trans_normal_df.index

    prot_normal_df = prot_normal_df.melt(id_vars = 'Name', var_name = 'Gene', value_name = 'Proteomics')
    trans_normal_df = trans_normal_df.melt(id_vars = 'Name', var_name = 'Gene', value_name = 'Transcriptomics')

    prot_tumor_df = cancer.get_proteomics("pdc", 'tumor')
    if isinstance(prot_tumor_df.columns, pd.MultiIndex):
        prot_tumor_df = ut.reduce_multiindex(df= prot_tumor_df, levels_to_drop = 'Database_ID', quiet=True)

    trans_tumor_df = cancer.get_transcriptomics("washu", 'tumor')
    if isinstance(trans_tumor_df.columns, pd.MultiIndex):
        trans_tumor_df = ut.reduce_multiindex(df = trans_tumor_df, levels_to_drop='Database_ID', quiet=True)

    prot_tumor_df['Name'] = prot_tumor_df.index
    trans_tumor_df['Name'] = trans_tumor_df.index

    prot_tumor_df = prot_tumor_df.melt(id_vars = 'Name', var_name = 'Gene', value_name = 'Proteomics')
    trans_tumor_df = trans_tumor_df.melt(id_vars = 'Name', var_name = 'Gene', value_name = 'Transcriptomics')

    prot_tumor_df['Tissue'] = ['Tumor'] * len(prot_tumor_df)
    prot_normal_df['Tissue'] = ['Normal'] * len(prot_normal_df)
    trans_tumor_df['Tissue'] = ['Tumor'] * len(trans_tumor_df)
    trans_normal_df['Tissue'] = ['Normal'] * len(trans_normal_df)

    prot_df = pd.concat([prot_tumor_df, prot_normal_df])
    trans_df = pd.concat([trans_tumor_df, trans_normal_df])

    return(pd.merge(prot_df, trans_df).dropna())

def permutate(df, column = 'Tissue', label1 = 'Tumor', label2 ='Normal', cutoff = 15, num_permutations = 10000, return_perm_list = False):
    delta_corr = delta_correlation(df, column = column, label1 = label1, label2 = label2, cutoff= cutoff)
    perm_delta_corrs = []
    for i in range(0, num_permutations):
        df[column] = np.random.permutation(df[column])
        perm_delta_corr = delta_correlation(df, column = column, label1 = label1, label2 = label2, cutoff = cutoff)
        perm_delta_corrs.append(perm_delta_corr)
    if return_perm_list:
        return perm_delta_corrs
    z_score = (delta_corr - np.mean(perm_delta_corrs)) / np.std(perm_delta_corrs)
    p_val = stats.norm.sf(abs(z_score)) * 2  
    return(delta_corr, p_val)

def delta_correlation(df, column = 'Tissue', label1 = 'Tumor', label2 ='Normal', cutoff = 15):
    normal_corr = df[df[column] == label2].corr(method = 'spearman',min_periods = cutoff ).iloc[0][1]
    tumor_corr = df[df[column] == label1].corr(method = 'spearman',min_periods = cutoff).iloc[0][1]
    delta_corr = tumor_corr - normal_corr
    return delta_corr

warnings.filterwarnings('ignore')
input_cancer_type = sys.argv[1]
input_permutation_number = int(sys.argv[2])
token = sys.argv[3]
cutoff = 15

if input_cancer_type == "CCRCC":
    cancer = pc.PancanCcrcc()
elif input_cancer_type == "Endometrial":
    cancer = pc.PancanUcec()
    cutoff = 10
elif input_cancer_type == "LUAD":
    cancer = pc.PancanLuad()
elif input_cancer_type == "HNSCC":
    cancer = pc.PancanHnscc()
elif input_cancer_type == "LSCC":
    cancer = pc.PancanLscc()


prot_trans_df = get_prot_trans_df(cancer)
perm_results = prot_trans_df.groupby('Gene').apply(lambda x: permutate(df = x, cutoff = cutoff, num_permutations = input_permutation_number))
df = pd.DataFrame.from_records(perm_results, columns = ['Delta_Correlation', 'P_Value'])
df.index = perm_results.index
df = df.dropna()
df.reset_index(inplace = True)
df['FDR'] = ssm.fdrcorrection(df.P_Value)[1]
df['Cancer'] = [input_cancer_type] * len(df)
file_name = input_cancer_type + '_delta_corr_pval_df.csv'
df.to_csv(file_name, index = False)
  