import Delta_Correlation as dc
import cptac
import statsmodels.stats.multitest as ssm
import pandas as pd
import warnings
import sys

warnings.filterwarnings('ignore')

input_cancer_type = sys.argv[1]
input_permutation_number = int(sys.argv[2])
cutoff = 15

if input_cancer_type == "CCRCC":
    cancer = cptac.Ccrcc()
elif input_cancer_type == "Endometrial":
    cancer = cptac.Endometrial()
    cutoff = 10
elif input_cancer_type == "LUAD":
    cancer = cptac.Luad()
elif input_cancer_type == "HNSCC":
    cancer  = cptac.Hnscc()
elif input_cancer_type == "LSCC":
    cancer = cptac.Lscc()
    cancer  = cptac.Hnscc()
elif input_cancer_type == "PDAC":
    cancer = cptac.Pdac()

prot_trans_df = dc.get_prot_trans_df(cancer)
perm_results = prot_trans_df.groupby('Gene').apply(lambda x: dc.permutate(df = x, cutoff = cutoff, num_permutations = input_permutation_number))
df = pd.DataFrame.from_records(perm_results, columns = ['Delta_Correlation', 'P_Value'])
df.index = perm_results.index
df = df.dropna()
df.reset_index(inplace = True)
df['FDR'] = ssm.fdrcorrection(df.P_Value)[1]
df['Cancer'] = [input_cancer_type] * len(df)
file_name = input_cancer_type + '_delta_corr_pval_df.csv'
df.to_csv(file_name, index = False)
  