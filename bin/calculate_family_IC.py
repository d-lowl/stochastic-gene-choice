#!/bin/env python3
#%%
#%% Code so the same scripts can be used as notebooks and as part of a pipeline
try:
    if __IPYTHON__:
        import sys
        sys.path.insert(0,'./bin')
except:
    pass

from pipeline_utils.ic_utils import get_IC, get_observed_variance, get_pb_variance
import pandas as pd
import numpy as np
from scipy.stats.mstats import gmean
from tqdm import tqdm
tqdm.pandas()

try:
    if __IPYTHON__:
        familyset_filename = "../familyset.csv"
        dataset_filename = "../lcl_european.csv"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Calculate IC for a set of families")
    parser.add_argument("familyset_filename")
    parser.add_argument("dichotomised_dataset_filename")
    args = parser.parse_args()
    familyset_filename = args.familyset_filename
    dichotomised_dataset_filename = args.dichotomised_dataset_filename

# %%
dichotomised_df = pd.read_csv(dichotomised_dataset_filename, index_col="gene_id")

#%%
families_df = pd.read_csv(familyset_filename)

#%%
def calculate_IC(group):
    entry = {}
    entry["family_id"] = group.iloc[0,1]
    gene_ids = group.gene_id
    entry["family_size"] = len(gene_ids)
    genes = dichotomised_df.reindex(gene_ids).dropna(how="all").copy()
    entry['n_genes'] = genes.shape[0]
    entry['obs_var'] = get_observed_variance(genes)
    entry['pb_var'] = get_pb_variance(genes)
    entry['ic'] = get_IC(entry['obs_var'], entry['pb_var'])
    entry['mean_expression'] = genes.sum(axis=0).mean()
    return pd.Series(entry)

#%%
print("Calculate IC and report")
report = families_df.groupby("family_id").progress_apply(calculate_IC)
# 621 ms

#%%
report.to_csv("family_IC.csv", index=False)