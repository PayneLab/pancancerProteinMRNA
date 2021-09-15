import cptac
import scipy.stats as stats
import numpy as np
import pandas as pd
import copy
import cptac.utils as ut


def load_cancers(include_pdac = False):
    ccrcc = cptac.Ccrcc()
    en = cptac.Endometrial()
    luad = cptac.Luad()
    hnscc  = cptac.Hnscc()
    lscc = cptac.Lscc()
    cancers = [ccrcc, en, luad, hnscc, lscc]
    cancer_names = ['CCRCC', 'Endometrial', 'LUAD', 'HNSCC', 'LSCC']
    if include_pdac:
        pdac = cptac.Pdac()
        cancers.append(pdac)
        cancer_names.append('PDAC')
    return cancers, cancer_names

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
