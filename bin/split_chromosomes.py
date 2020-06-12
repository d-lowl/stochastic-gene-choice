#!/bin/env python3
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
import re

try:
    if __IPYTHON__:
        args = {}
        gff_in = "GRCm38.p6_genomic.gff"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Split gff file into chromosomes and filter")
    parser.add_argument("gff_in")
    args = parser.parse_args()
    gff_in = args.gff_in

# %%
df = read_gff3(gff_in)


# %% extract chromosome sequence ID
df_slice = df[df.type == "region"]

# %%
def get_chromosome_name(x):
    try:
        return re.search(r'(?:chromosome=)(.*?)(?:;)', x).groups()[0]
    except:
        return None
df_slice["chromosome"] = df_slice.attributes.apply(get_chromosome_name)

# %%
df_regions = df[
    (df.seqid.str.startswith("NC_")) & 
    (df.type == "region") &
    (df.source == "RefSeq") & 
    (df.attributes.str.contains("chromosome"))
]

# %%
df_regions.loc[:,"chromosome"] = df_slice.loc[:,"attributes"].apply(get_chromosome_name)

# %%
genes = get_genes_from_gff3(df).set_index("Name")

# %%
for (_id, chr_name) in df_regions.loc[:,["seqid", "chromosome"]].values:
    chr_df = genes[genes.seqid == _id]
    print("Chr{} genes (before filtering):".format(chr_name),chr_df.shape[0])

    # %% Remove duplicates 
    chr_df = chr_df.loc[~chr_df.index.duplicated(keep='first')]

    # %% Remove predicted genes
    chr_df = chr_df.loc[~chr_df.description.str.startswith("predicted gene").fillna(True).astype(bool)]
    print("Chr{} genes (after filtering):".format(chr_name),chr_df.shape[0])

    chr_df.to_csv("chr{}_filtered.csv".format(chr_name))

# %%
