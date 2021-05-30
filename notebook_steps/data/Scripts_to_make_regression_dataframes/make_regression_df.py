import cptac 
import pandas as pd
import os, sys
currentdir = os.path.dirname(os.path.realpath('make_regression_df.py'))
parentdir = os.path.dirname(currentdir)
currentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import Delta_Correlation as dc
import statsmodels.stats.multitest as ssm

input_cancer_type = sys.argv[1]

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
results = df.groupby('Gene').apply(dc.regression)
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
