#have command line parameter that says "breast cancer" and makes sure to save it out to my computer
import cptac
import scipy
import sys #first argument is cancer type, second is number of permutations
from scipy import stats
#import seaborn as sns
#import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
#import statistics
# import parse_correlations_dataframe as get_corr
import copy
import csv
# import get_correlations
import cptac.utils as ut
import warnings
warnings.filterwarnings("ignore")

input_cancer_type = sys.argv[1]
if input_cancer_type == "ccrcc":
    ccrcc = cptac.Ccrcc()
    cancer_list = [ccrcc]
elif input_cancer_type == "en":
    en = cptac.Endometrial()
    cancer_list = [en]
elif input_cancer_type == "luad":
    luad = cptac.Luad()
    cancer_list = [luad]
elif input_cancer_type == "hnscc":
    hnscc  = cptac.Hnscc()
    cancer_list = [hnscc]
elif input_cancer_type == "lscc":
    lscc = cptac.Lscc()
    cancer_list = [lscc]
# brca = cptac.Brca()
# ccrcc = cptac.Ccrcc()
# colon = cptac.Colon()
# en = cptac.Endometrial()
# gbm = cptac.Gbm()
# luad = cptac.Luad()
# ovarian = cptac.Ovarian()
# hnscc  = cptac.Hnscc()
# lscc = cptac.Lscc()

# type_dict = {brca:"brca",ccrcc:"ccrcc",colon:"colon",en:"endometrial",gbm:"gbm",luad:"luad",
#                   ovarian:"ovarian",hnscc:"hnscc",lscc:"lscc"}

# cancer_from_string = {"brca":brca,"ccrcc":ccrcc,"colon":colon,"en":en,"gbm":gbm,"luad":luad,
#                   "ovarian":ovarian,"hnscc":hnscc,"lscc":lscc}

#cancer_list = [ccrcc, en, luad, hnscc, lscc]
with open('gene_list.csv', newline='') as in_file:
    reader = csv.reader(in_file)
    gene_list = list(reader)[0]

# https://link.springer.com/article/10.3758/s13428-012-0289-7
def compare_correlations(r1, r2, n1, n2):
    rp1 = np.arctanh(r1)
    rp2 = np.arctanh(r2)

    if n1 < 4 or n2 < 4:
        return(np.nan)
    Sr12 = math.sqrt((1/(n1-3))+(1/(n2-3)))
    z = (rp1-rp2) / Sr12
    p = scipy.stats.norm.sf(abs(z))*2
    return (p)

def permute(df,original_correlation, label_1, label_2, column_one, column_two, permutation_times):
    permutation_list = []
    permu_df = copy.deepcopy(df)

    for i in range(permutation_times):
        permu_df["Type"] = np.random.permutation(permu_df["Type"])
        permu_is_label_1 = permu_df["Type"] == label_1
        permu_is_label_2 = permu_df["Type"] == label_2
        label_1_correlation,label_1_pval = scipy.stats.pearsonr(permu_df[permu_is_label_1][column_one], permu_df[permu_is_label_1][column_two])
        label_2_correlation,label_2_pval = scipy.stats.pearsonr(permu_df[permu_is_label_2][column_one], permu_df[permu_is_label_2][column_two])
        delta = label_1_correlation - label_2_correlation
        permutation_list.append(delta)

    z_score = (original_correlation - np.mean(permutation_list)) / np.std(permutation_list)
    p_val = scipy.stats.norm.sf(abs(z_score))*2
    return p_val

tot_diff_list = []
tot_pval_list = []
tot_perm_list = []
for cancer in cancer_list:
    cancer_diff_list = [sys.argv[1]]
    cancer_pval_list = [sys.argv[1]]
    cancer_perm_list = [sys.argv[1]]

    tumor_cancer_df = cancer.join_omics_to_omics("transcriptomics","proteomics",tissue_type="tumor",quiet=True)
    if isinstance(tumor_cancer_df.columns, pd.MultiIndex):
        tumor_cancer_df = ut.reduce_multiindex(df = tumor_cancer_df, levels_to_drop="Database_ID",quiet=True)

    normal_cancer_df = cancer.join_omics_to_omics("transcriptomics","proteomics",tissue_type="normal",quiet=True)
    if isinstance(normal_cancer_df.columns, pd.MultiIndex):
        normal_cancer_df = ut.reduce_multiindex(df = normal_cancer_df, levels_to_drop="Database_ID",quiet=True)

    for gene in gene_list:
        gene_trans = gene + "_transcriptomics"
        gene_prot = gene + "_proteomics"
        gene_in_tumor = gene_trans in tumor_cancer_df.columns and gene_prot in tumor_cancer_df.columns
        gene_in_normal = gene_trans in normal_cancer_df.columns and gene_prot in normal_cancer_df.columns

        if not(gene_in_tumor and gene_in_normal):
            cancer_diff_list.append(np.nan)
            cancer_pval_list.append(np.nan)
            cancer_perm_list.append(np.nan)
            continue

        tumor_df = tumor_cancer_df[[gene_trans,gene_prot]]
        #The following line takes care of the problem that arises when reducing the multi-index like we did earlier.
        #There are sometimes multiple columns of proteomics and sometimes multiple of transcriptomics, this
        #takes the first one.
        if isinstance(tumor_df[gene_trans], pd.core.frame.DataFrame) or isinstance(tumor_df[gene_prot], pd.core.frame.DataFrame): #This is to take first column of multi-index
            trans_col = tumor_df[gene_trans]
            if isinstance(tumor_df[gene_trans], pd.core.frame.DataFrame):
                trans_col = trans_col.iloc[:,0]
            prot_col = tumor_df[gene_prot]
            if isinstance(tumor_df[gene_prot], pd.core.frame.DataFrame):
                prot_col = prot_col.iloc[:,0]
            frame = {gene_trans : trans_col, gene_prot : prot_col}
            tumor_df = pd.DataFrame(frame)
        tumor_df = tumor_df.dropna()
        num_tumor = len(tumor_df)
        tumor_corr = tumor_df.corr().iloc[0][1]

        normal_df = normal_cancer_df[[gene_trans,gene_prot]]
        if isinstance(normal_df[gene_trans], pd.core.frame.DataFrame) or isinstance(normal_df[gene_prot], pd.core.frame.DataFrame): #This is to take first column of multi-index
            trans_col = normal_df[gene_trans]
            if isinstance(normal_df[gene_trans], pd.core.frame.DataFrame):
                trans_col = trans_col.iloc[:,0]
            prot_col = normal_df[gene_prot]
            if isinstance(normal_df[gene_prot], pd.core.frame.DataFrame):
                prot_col = prot_col.iloc[:,0]
            frame = {gene_trans : trans_col, gene_prot : prot_col}
            normal_df = pd.DataFrame(frame)
        normal_df = normal_df.dropna()
        num_normal = len(normal_df)
        normal_corr = normal_df.corr().iloc[0][1]

        corr_diff = tumor_corr - normal_corr
        cancer_diff_list.append(corr_diff)

        gene_pval = compare_correlations(tumor_corr, normal_corr, num_tumor, num_normal)
        cancer_pval_list.append(gene_pval)

        #Permutations
        if num_tumor < 4 or num_normal < 4:
            cancer_perm_list.append(np.nan)
            continue
        tumor_label_list = ['tumor'] * len(tumor_df)
        tumor_df["Type"] = tumor_label_list

        normal_label_list = ['normal'] * len(normal_df)
        normal_df["Type"] = normal_label_list

        perm_list = [tumor_df,normal_df]
        perm_df = pd.concat(perm_list)

        column_one = perm_df.columns[0]
        column_two = perm_df.columns[1]
        perm_val = permute(perm_df,corr_diff,"tumor","normal",column_one,column_two,int(sys.argv[2])) # Look at how long for 1000, then calculate
        #Break into units for cancer (command line for bash shell)
        cancer_perm_list.append(perm_val)

    tot_diff_list.append(cancer_diff_list)
    tot_pval_list.append(cancer_pval_list)
    tot_perm_list.append(cancer_perm_list)
labels = ["Cancer"]
labels.extend(gene_list)
df = pd.DataFrame.from_records(tot_diff_list,columns=labels)
df2 = pd.DataFrame.from_records(tot_pval_list,columns=labels)
df3 = pd.DataFrame.from_records(tot_perm_list,columns=labels)

# df.to_csv("corr_diff.csv",index=False)
# df2.to_csv("p_val.csv",index=False)
out_string = sys.argv[1] + "_permutation_" + sys.argv[2] + ".csv"
df3.to_csv(out_string,index=False)
