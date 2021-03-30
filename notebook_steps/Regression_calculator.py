import cptac
import pandas as pd
from scipy import stats
import numpy as np
import statsmodels.stats.multitest as ssm
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as robj
import pandas as pd


#ccrcc = cptac.Ccrcc()
#en = cptac.Endometrial()
#luad = cptac.Luad()
#hnscc  = cptac.Hnscc()
#lscc = cptac.Lscc()


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
        print(data)

    
    
def calculate_regression (cancer_list, gene_list, cancer_names):
    dfs = []
    for gene in gene_list:        
        for cancer_type, cancer_name in zip(cancer_list, cancer_names):
            tumor = cancer_type.multi_join({'proteomics': gene, 'transcriptomics': gene}, tissue_type= 'tumor', flatten = True)
            normal = cancer_type.multi_join({'proteomics': gene, 'transcriptomics': gene}, tissue_type= 'normal', flatten = True)
            tumor = tumor.dropna()
            normal = normal.dropna()
            if len(normal) < 2 or len(tumor) < 2:
                return float("NaN"), float("NaN"), float("NaN")
            tumor.columns = ['proteomics', 'transcriptomics']
            normal.columns = ['proteomics', 'transcriptomics']
            groups = ['tumor'] * len(tumor)
            groups.extend(['normal']*len(normal))
            prot_list = list(tumor['proteomics'])
            prot_list.extend(list(normal['proteomics']))
            trans_list = list(tumor['transcriptomics'])
            trans_list.extend(list(normal['transcriptomics']))
            gene_df = pd.DataFrame({'Type': groups, 'Proteomics': prot_list, 'Transcriptomics': trans_list})
            gene_df['Cancer'] = [cancer_name] * len(gene_df)
            gene_df['Gene'] = [gene] * len(gene_df)
            dfs.append(gene_df)    
    full_df = pd.concat(dfs)
    full_df = full_df.rename(columns ={'Type': 'Tissue'})
   # print(full_df)
    cancer_dfs = []
    for cancer in pd.unique(full_df.Cancer):
      #  print(cancer)
        rows = []
        for gene in list(pd.unique(full_df.Gene)):
          #  print(gene)
            d = {}
            df = full_df[full_df.Gene == gene]
            df = df[df.Cancer == cancer]
            if len(df) < 4:
                continue
            df = df[['Tissue', 'Proteomics', 'Transcriptomics']]
            lm_df= linear_model(df, 'Transcriptomics', 'Proteomics', 'Tissue')
            d['gene'] = gene
            d['cancer'] = cancer
            d['interaction_coeff'] = lm_df['Estimate'][3]
            d['condition_coeff'] = lm_df['Estimate'][2]
            d['transcript_coeff'] = lm_df['Estimate'][1]
            d['intercept'] = lm_df['Estimate'][0]
            d['interaction_pval'] = lm_df['Pr(>|t|)'][3]
            d['condition_pval'] = lm_df['Pr(>|t|)'][2]
            d['transcript_pval'] = lm_df['Pr(>|t|)'][1]
            d['intercept_pval'] = lm_df['Pr(>|t|)'][0]
            rows.append(d)
        cancer_df = pd.DataFrame(rows)
        for column in cancer_df:
            if 'pval' in column:
                old_pvals = list(cancer_df[column])
                adj_pvals = list(ssm.fdrcorrection(old_pvals)[1])
                cancer_df[column] = adj_pvals    
        cancer_dfs.append(cancer_df)
    
    lm_df = pd.concat(cancer_dfs)
    return(lm_df)