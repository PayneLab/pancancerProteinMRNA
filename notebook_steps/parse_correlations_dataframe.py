import pandas as pd
import load_data
def get_data(cancer_type="all", tissue_type="all", genes = "all", get_corr=True, get_r2=False, get_p=False):
    df = load_data.load_correlation_df()
    #df = df.rename(columns = {'Unnamed: 0': "Cancer Type"})
    if cancer_type != "all":
        df = df[df['Cancer_Type'].str.lower().isin(cancer_type)]  
    if tissue_type != "all":
        df = df[df['Tissue_Type'].str.lower().isin(tissue_type)]
    if genes != "all":
        df = df[df['Gene'].isin(genes)]
    if not get_corr:
        df = df.drop(columns = ["Correlation"])
    if not get_r2:
        df = df.drop(columns = ["P-value"])
    if not get_p:
        df = df.drop(columns = ["R-squared"])
    return (df)
        