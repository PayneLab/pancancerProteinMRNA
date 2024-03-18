import cptac.pancan as pc
import statsmodels.stats.multitest as ssm
import pandas as pd
import warnings
import sys, os

warnings.filterwarnings('ignore')
import pcprutils as dc

input_cancer_type = sys.argv[1]
mutated_gene = sys.argv[2]
input_permutation_number = int(sys.argv[3])


cutoff = 15

if input_cancer_type == "CCRCC":
    cancer = pc.PancanCcrcc(no_internet = True)
elif input_cancer_type == "Endometrial":
    cancer = pc.PancanUcec(no_internet = True)
    cutoff = 10
elif input_cancer_type == "LUAD":
    cancer = pc.PancanLuad(no_internet = True)
elif input_cancer_type == "HNSCC":
    cancer = pc.PancanHnscc(no_internet = True)
elif input_cancer_type == "LSCC":
    cancer = pc.PancanLscc(no_internet = True)
    
mutation_df = cancer.get_somatic_mutation()
mutation_df = mutation_df[mutation_df.Gene == mutated_gene]
mutation_df = mutation_df[mutation_df.Mutation != 'Silent']
mutation_df = mutation_df[mutation_df.Mutation != 'RNA']
mutation_df = mutation_df[mutation_df.Mutation != 'synonymous SNV']
mutation_df.reset_index(inplace = True)
prot_trans_df = dc.get_prot_trans_df(cancer)
prot_trans_df = prot_trans_df[prot_trans_df.Tissue == 'Tumor']
prot_trans_df['Mutation'] = prot_trans_df.Name.isin(mutation_df.Patient_ID)
prot_trans_df = prot_trans_df.drop(columns = 'Tissue')
perm_results = prot_trans_df.groupby('Gene').apply(lambda x: dc.permutate(df = x, cutoff = cutoff, num_permutations = input_permutation_number, column = 'Mutation', label1 = True, label2 = False))
df = pd.DataFrame.from_records(perm_results, columns = ['Delta_Correlation', 'P_Value'])
df.index = perm_results.index
df = df.dropna()
df.reset_index(inplace = True)
df['FDR'] = ssm.fdrcorrection(df.P_Value)[1]
df['Cancer'] = [input_cancer_type] * len(df)
df['Mutated_Gene'] = [mutated_gene] * len(df)
file_name = input_cancer_type + '_' + mutated_gene + '_transmutation_effects.csv'
df.to_csv(file_name, index = False)
