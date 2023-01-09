import cptac 
import pandas as pd
import os, sys

currentdir = os.path.dirname(os.path.realpath('make_regression_df.py'))
parentdir = os.path.dirname(currentdir)
parentdir = os.path.dirname(parentdir)
sys.path.append(parentdir) #RPY2_CFFI_MODE=ABI
import Delta_Correlation as dc
import statsmodels.stats.multitest as ssm
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as robj

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

def compute_regression(input_cancer_type):

    if input_cancer_type == "CCRCC":
        cancer = cptac.Ccrcc()
    elif input_cancer_type == "Endometrial":
        cancer = cptac.Endometrial()
    elif input_cancer_type == "LUAD":
        cancer = cptac.Luad()
    elif input_cancer_type == "HNSCC":
        cancer  = cptac.Hnscc()
    elif input_cancer_type == "LSCC":
        cancer = cptac.Lscc()
    elif input_cancer_type == "PDAC":
        cancer = cptac.Pdac()

    df = dc.get_prot_trans_df(cancer)
    results = df.groupby('Gene').apply(regression)
    reg_df = pd.DataFrame(list(results))
    reg_df.index = results.index
    reg_df.reset_index(inplace = True)
    reg_df = reg_df.dropna()
    reg_df['interaction_FDR'] = ssm.fdrcorrection(reg_df['interaction_pval'])[1]
    reg_df['condition_FDR'] = ssm.fdrcorrection(reg_df['condition_pval'])[1]
    reg_df['intercept_FDR'] = ssm.fdrcorrection(reg_df['intercept_pval'])[1]
    reg_df['Cancer'] = [input_cancer_type] * len(reg_df)

    file_name = input_cancer_type + '_regressions.csv'
    reg_df.to_csv(file_name, index = False)
