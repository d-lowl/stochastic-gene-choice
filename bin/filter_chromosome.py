#!/bin/env python
#%% Code so the same scripts can be used as notebooks and as part of a pipeline
try:
    if __IPYTHON__:
        import sys
        sys.path.insert(0,'./bin')
except:
    pass

#%%
from pipeline_utils.gff_utils import read_gff3, get_genes_from_gff3
import numpy as np

try:
    if __IPYTHON__:
        args = {}
        gff_in = "GRCm38.p6_genomic.gff"
        gff_out = "intermediate/chr18_filtered.csv"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Filter a chromosome from duplicate genes and predicted genes")
    parser.add_argument("gff_in")
    parser.add_argument("gff_out")
    args = parser.parse_args()
    gff_in = args.gff_in
    gff_out = args.gff_out
# %%
chr18_df = get_genes_from_gff3(read_gff3(gff_in)).set_index("Name")
print("All genes:",chr18_df.shape[0])

#%%
chr18_df = chr18_df[chr18_df.seqid == "NC_000084.6"]
print("Chr18 genes (before filtering):",chr18_df.shape[0])

# %% Remove duplicates 
chr18_df = chr18_df.loc[~chr18_df.index.duplicated(keep='first')]

# %% Remove predicted genes
chr18_df = chr18_df.loc[~chr18_df.description.str.startswith("predicted gene").fillna(True).astype(bool)]
print("Chr18 genes (after filtering):",chr18_df.shape[0])

# %%
chr18_df.to_csv(gff_out)