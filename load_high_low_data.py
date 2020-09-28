import pandas as pd

def load_from_file():
    file="proteins_29.csv"
    df = pd.read_csv(file)
    df = df.drop(columns=["Protein ID", "Gene ID"])
    df = df.set_index("Gene name")
    df = df.apply(pd.to_numeric)
    return df
