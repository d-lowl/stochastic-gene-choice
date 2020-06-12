#!/bin/env python3
#%%
import pandas as pd
from tqdm import tqdm
try:
    if __IPYTHON__:
        familyset_filename = "../HGNC_families.csv"
        dataset_filename = "../lcl_european.csv"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Convert familyset from XLSX to CSV")
    parser.add_argument("familyset_filename")
    parser.add_argument("dataset_filename")
    args = parser.parse_args()
    familyset_filename = args.familyset_filename
    dataset_filename = args.dataset_filename
#%%
df = pd.read_csv(familyset_filename)
dataset_df = pd.read_csv(dataset_filename, index_col="gene_id")

#%%
df = df.loc[:,["family_id", "gene_symbol"]]#.drop_duplicates().dropna()

# %%
overlap = dataset_df.index.intersection(df.gene_symbol)
df = df.set_index("gene_symbol").loc[overlap,:].reset_index().drop_duplicates()
# %%
fam_counts = df.groupby("family_id").count().iloc[:,0]
for i, N in tqdm(fam_counts.iteritems()):
    df.loc[df.family_id == i, "N"] = N

# #%% Pick the smallest family assuming the lowest in the heirarchy
# df = df.drop_duplicates().sort_values("N")
# df = df.drop_duplicates("mgi_symbol",keep="first").sort_values(by="mgi_symbol")

# # %% Recalculate family sizes again
# fam_counts = df.groupby("family_rm_sub").count().iloc[:,0]
# for i, N in fam_counts.iteritems():
#     df.loc[df.family_rm_sub == i, "N"] = N

#%% Drop families smaller then 5
df = df[df.N >= 5]

#%%
df.columns = ["gene_id", "family_id", "N"]
#%%
df.to_csv("familyset.csv", index=False)