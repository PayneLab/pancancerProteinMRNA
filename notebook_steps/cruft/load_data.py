import pandas as pd

def load_proteins_29():
    file="data/proteins_29.csv"
    df = pd.read_csv(file)
    df = df.drop(columns=["Protein ID", "Gene ID"])
    df = df.set_index("Gene name")
    df = df.apply(pd.to_numeric)
    return df
def load_correlation_df():
    df = pd.read_csv("data/correlations_dataframe.csv")
    return df
