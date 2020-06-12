#!/bin/env python3
try:
    if __IPYTHON__:
        import sys
        sys.path.insert(0,'./bin')
except:
    pass

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
    parser = argparse.ArgumentParser(description="Dichotomise genes with family-wise thresholds")
    parser.add_argument("familyset_filename")
    parser.add_argument("dataset_filename")
    parser.add_argument("type")
    args = parser.parse_args()
    familyset_filename = args.familyset_filename
    dataset_filename = args.dataset_filename
    dich_type = args.type

# %%
df = pd.read_csv(dataset_filename, index_col="gene_id")

#%%
families_df = pd.read_csv(familyset_filename)

#%%
def _geomean05(genes):
    _genes = genes.values
    bottom = [np.mean(np.sort(row[row >= 0.5])[:3]) if len(row[row >= 0.5]) < 120 else np.quantile(row[row >= 0.5], 0.025) for row in _genes]
    top = [np.mean(np.sort(row[row >= 0.5])[-3:]) if len(row[row >= 0.5]) < 120 else np.quantile(row[row >= 0.5], 0.975) for row in _genes]
    if(len(bottom) == 0 or len(top) == 0):
        return np.nan
    return gmean([np.nanmin(bottom), np.nanmax(top)])

def _fm05(genes):
    _genes = genes.values
    top = [np.mean(np.sort(row[row >= 0.5])[-3:]) if len(row[row >= 0.5]) < 120 else np.quantile(row[row >= 0.5], 0.975) for row in _genes]
    if(len(top) == 0):
        return np.nan
    return np.nanmax(top) / 10

def get_threshold(genes):
    if dich_type == "geomean05":
        return _geomean05(genes)
    elif dich_type == "fm05":
        return _fm05(genes)
    elif dich_type == "one":
        return 1.0
    else:
        raise Exception("Unknown dichotomisation type")

#%% Get thresholds 
def get_threshold_wrapper(group):
    entry = {}
    gene_ids = group.gene_id
    genes = df.reindex(gene_ids).dropna().copy()
    threshold = get_threshold(genes)
    entry['uncorrected_threshold'] = threshold
    threshold = max(0.5, threshold)
    entry['threshold'] = threshold

    return pd.Series(entry)

print("Calculate thresholds")
thresholds = families_df.groupby("family_id").progress_apply(get_threshold_wrapper)

#%% Merge thresholds
families_df = families_df.merge(thresholds, left_on="family_id", right_index=True)
# 8.45 ms

#%% Dichotomise
def dichotomise():
    _thresholds = families_df.set_index("gene_id").threshold
    transposed_values = df.loc[_thresholds.index].T.values

    transposed_values[transposed_values < _thresholds.values] = 0
    transposed_values[transposed_values > 0] = 1

    return pd.DataFrame(transposed_values.T, index=_thresholds.index, columns=df.columns)

print("Dichotomise genes")
dichotomised_df = dichotomise()

dichotomised_df.to_csv("dichotomised_genes.csv")
thresholds.to_csv("family_thresholds.csv")