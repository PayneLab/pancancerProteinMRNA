import cptac
import pandas as pd
from scipy import stats
import numpy as np
import statsmodels.stats.multitest as ssm
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as robj
import pandas as pd


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

    
    
def calculate_regression (full_df, cutoff = 10):
    cancer_dfs = []
    for cancer in pd.unique(full_df.Cancer):
        print(cancer)
        rows = []
        for gene in list(pd.unique(full_df.Gene)):
            print(gene)
            d = {}
            df = full_df[full_df.Gene == gene]
            df = df[df.Cancer == cancer]
            df = df[['Tissue', 'Proteomics', 'Transcriptomics']]
            if len(df[df.Tissue == 'normal']) < cutoff or len(df[df.Tissue == 'tumor']) < cutoff:
                print(gene + 'not enough')
                continue
            lm_df= linear_model(df, 'Transcriptomics', 'Proteomics', 'Tissue')
            d['gene'] = gene
            d['cancer'] = cancer
            if len(lm_df) == 4:
                d['interaction_coeff'] = lm_df['Estimate'][3]
                d['condition_coeff'] = lm_df['Estimate'][2]
                d['transcript_coeff'] = lm_df['Estimate'][1]
                d['intercept'] = lm_df['Estimate'][0]
                d['interaction_pval'] = lm_df['Pr(>|t|)'][3]
                d['condition_pval'] = lm_df['Pr(>|t|)'][2]
                d['transcript_pval'] = lm_df['Pr(>|t|)'][1]
                d['intercept_pval'] = lm_df['Pr(>|t|)'][0]
                rows.append(d)
            else:
                print(gene + 'linear model error')
                continue
        cancer_df = pd.DataFrame(rows)        
        for column in cancer_df:
            if 'pval' in column:
                old_pvals = list(cancer_df[column])
                adj_pvals = list(ssm.fdrcorrection(old_pvals)[1])
                cancer_df[column] = adj_pvals    
        cancer_dfs.append(cancer_df)    
    lm_df = pd.concat(cancer_dfs)
    return(lm_df)