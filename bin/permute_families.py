#!/bin/env python3
#%%
#%% Code so the same scripts can be used as notebooks and as part of a pipeline
try:
    if __IPYTHON__:
        import sys
        sys.path.insert(0,'./bin')
except:
    pass

import pandas as pd
import numpy as np

try:
    if __IPYTHON__:
        familyset_filename = "../somatosensory_rpkm_suppl/merged_hgnc_biomart_processed1_corrected_scattered/intermediate/familyset.csv"
        dataset_filename = "../somatosensory_rpkm_suppl.csv"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Calculate IC for a set of families")
    parser.add_argument("familyset_filename")
    parser.add_argument("dataset_filename")
    args = parser.parse_args()
    familyset_filename = args.familyset_filename
    dataset_filename = args.dataset_filename

# %%
df = pd.read_csv(dataset_filename, index_col="gene_id")
families_df = pd.read_csv(familyset_filename)

# %%
families_df.loc[:,"gene_id"] = families_df.gene_id.sample(frac=1.0).reset_index(drop=True)

# %%
families_df.to_csv("shuffled_families.csv",index=False)