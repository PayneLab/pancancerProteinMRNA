import cptac
import scipy.stats as stats
import numpy as np
import pandas as pd
import copy
import cptac.utils as ut
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as robj

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
    delta_corr = delta_correlation(df, column = column, label1 = label1, label2 = label2, cutoff= cutoff)
    perm_delta_corrs = []
    for i in range(0, num_permutations):
        df[column] = np.random.permutation(df[column])
        perm_delta_corr = delta_correlation(df, column = column, label1 = label1, label2 = label2, cutoff = cutoff)
        perm_delta_corrs.append(perm_delta_corr)
   # print(perm_delta_corrs)
    z_score = (delta_corr - np.mean(perm_delta_corrs)) / np.std(perm_delta_corrs)
   # print(z_score)
    p_val = stats.norm.sf(abs(z_score)) * 2
    
    return(delta_corr, p_val)

def delta_correlation(df, column = 'Tissue', label1 = 'tumor', label2 ='normal', cutoff = 15):
    normal_corr = df[df[column] == label2].corr(method = 'pearson',min_periods = cutoff ).iloc[0][1]
    tumor_corr = df[df[column] == label1].corr(method = 'pearson',min_periods = cutoff).iloc[0][1]
    delta_corr = tumor_corr - normal_corr
    return delta_corr

def linear_model(data, Input, Output, Condition):
    try:
        stats = importr('stats')
        base = importr('base')
        pandas2ri.activate()
        r_df = pandas2ri.py2rpy(data)
        pandas2ri.deactivate()
        formula = '{y}~{x}*{condition}'.format(y = Output, x = Input, condition = Condition)
        lm = stats.lm(formula, r_df)
        summary = (base.summary(lm))
        results = summary.rx2('coefficients')
        results_df = base.as_data_frame_matrix(results)
        py_results_df = pd.DataFrame(results_df).transpose()
        py_results_df.columns = results_df.colnames
        py_results_df.index = results_df.rownames
        return(py_results_df)
    except:
        return(pd.DataFrame({}))

def regression(df):
    lm_df= linear_model(df, 'Transcriptomics', 'Proteomics', 'Tissue')
    d = dict()
    d['interaction_coeff'] = float('NaN')
    d['condition_coeff'] = float('NaN')
    d['transcript_coeff'] = float('NaN')
    d['intercept'] = float('NaN')
    d['interaction_pval'] = float('NaN')
    d['condition_pval'] = float('NaN')
    d['transcript_pval'] = float('NaN')
    d['intercept_pval'] = float('NaN')
    if len(lm_df) == 4:
        d['interaction_coeff'] = lm_df['Estimate'][3]
        d['condition_coeff'] = lm_df['Estimate'][2]
        d['transcript_coeff'] = lm_df['Estimate'][1]
        d['intercept'] = lm_df['Estimate'][0]
        d['interaction_pval'] = lm_df['Pr(>|t|)'][3]
        d['condition_pval'] = lm_df['Pr(>|t|)'][2]
        d['transcript_pval'] = lm_df['Pr(>|t|)'][1]
        d['intercept_pval'] = lm_df['Pr(>|t|)'][0]
    return d
        