import cptac
import scipy.stats as stats
import numpy as np
import pandas as pd
import copy
import cptac.utils as ut

# 
def get_prot_trans_df(cancer):
    prot_normal_df = cancer.get_proteomics('normal')
    if isinstance(prot_normal_df.columns, pd.MultiIndex):
        prot_normal_df = ut.reduce_multiindex(df= prot_normal_df, levels_to_drop = 'Database_ID')
    trans_normal_df = cancer.get_transcriptomics('normal')
    if isinstance(trans_normal_df.columns, pd.MultiIndex):
        trans_normal_df = ut.reduce_multiindex(df = trans_normal_df, levels_to_drop='Database_ID')
    prot_normal_df['Patient_ID'] = prot_normal_df.index
    trans_normal_df['Patient_ID'] = trans_normal_df.index
    prot_normal_df = prot_normal_df.melt(id_vars = 'Patient_ID', var_name = 'Gene', value_name = 'Proteomics')
    trans_normal_df = trans_normal_df.melt(id_vars = 'Patient_ID', var_name = 'Gene', value_name = 'Transcriptomics')
    prot_tumor_df = cancer.get_proteomics('tumor')
    if isinstance(prot_tumor_df.columns, pd.MultiIndex):
        prot_tumor_df = ut.reduce_multiindex(df= prot_tumor_df, levels_to_drop = 'Database_ID')
    trans_tumor_df = cancer.get_transcriptomics('tumor')
    if isinstance(trans_tumor_df.columns, pd.MultiIndex):
        trans_tumor_df = ut.reduce_multiindex(df = trans_tumor_df, levels_to_drop='Database_ID')
    prot_tumor_df['Patient_ID'] = prot_tumor_df.index
    trans_tumor_df['Patient_ID'] = trans_tumor_df.index
    prot_tumor_df = prot_tumor_df.melt(id_vars = 'Patient_ID', var_name = 'Gene', value_name = 'Proteomics')
    trans_tumor_df = trans_tumor_df.melt(id_vars = 'Patient_ID', var_name = 'Gene', value_name = 'Transcriptomics')
    prot_tumor_df['Tissue'] = ['tumor'] * len(prot_tumor_df)
    prot_normal_df['Tissue'] = ['normal'] * len(prot_normal_df)
    trans_tumor_df['Tissue'] = ['tumor'] * len(trans_tumor_df)
    trans_normal_df['Tissue'] = ['normal'] * len(trans_normal_df)
    prot_df = pd.concat([prot_tumor_df, prot_normal_df])
    trans_df = pd.concat([trans_tumor_df, trans_normal_df])
    return(pd.merge(prot_df, trans_df).dropna())

def permutate(df, column = 'Tissue', label1 = 'tumor', label2 ='normal', cutoff = 15, num_permutations = 10000):
    delta_corr = delta_correlation(df, column = column, label1 = label1, label2 = label2)
    perm_delta_corrs = []
    for i in range(0, num_permutations):
        df.Tissue = np.random.permutation(df.Tissue)
        perm_delta_corr = delta_correlation(df, column = column, label1 = label1, label2 = label2)
        perm_delta_corrs.append(perm_delta_corr)
    z_score = (delta_corr - np.mean(perm_delta_corrs)) / np.std(perm_delta_corrs)
    p_val = stats.norm.sf(abs(z_score)) * 2
    
    return(delta_corr, p_val)

def delta_correlation(df, column = 'Tissue', label1 = 'tumor', label2 ='normal', cutoff = 15):
    normal_corr = df[df[column] == label2].corr(method = 'pearson',min_periods = cutoff ).iloc[0][1]
    tumor_corr = df[df[column] == label1].corr(method = 'pearson',min_periods = cutoff).iloc[0][1]
    delta_corr = tumor_corr - normal_corr
    return delta_corr