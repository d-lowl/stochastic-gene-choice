#!/bin/env python3
#%% Code so the same scripts can be used as notebooks and as part of a pipeline
try:
    if __IPYTHON__:
        import sys
        sys.path.insert(0,'./bin')
except:
    pass

#%%
import pandas as pd
import numpy as np
from tqdm import tqdm
from pipeline_utils.gff_utils import get_genes_from_gff3, read_gff3
from pipeline_utils.ic_utils import get_IC, get_observed_variance, get_pb_variance
tqdm.pandas()
try:
    if __IPYTHON__:
        dichotomised_filename = "intermediate/dichotomised_genes.csv"
        chromosome_filename = "intermediate/filtered_chromosome.csv"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Calculate optimal bimodality thresholds")
    parser.add_argument("dichotomised_filename")
    parser.add_argument("chromosome_filename")
    parser.add_argument("chr_name")
    args = parser.parse_args()
    dichotomised_filename = args.dichotomised_filename
    chromosome_filename = args.chromosome_filename
    chr_name = args.chr_name

# %%
df = pd.read_csv(dichotomised_filename, index_col="gene_id")
chr_df = pd.read_csv(chromosome_filename, index_col="Name")
ic_filename = "stage1_{}_IC.csv".format(chr_name)

# %%
ic_df = pd.DataFrame(columns=["start","end","stretch","n_genes","obs_var","pb_var","ic"])
total_genes = chr_df.shape[0]
for stretch in [7,14,21]:
    for start, end in tqdm(zip(range(0,total_genes-stretch),range(stretch,total_genes))):
        entry = {
            "start": start,
            "end": end,
            "stretch": stretch
        }
        try:
            gene_names = chr_df.index[start:end]
            df_slice = df.loc[gene_names,:].dropna(how="all")
            entry['n_genes'] = df_slice.shape[0]
            entry['obs_var'] = get_observed_variance(df_slice)
            entry['pb_var'] = get_pb_variance(df_slice)
            entry['ic'] = get_IC(entry['obs_var'], entry['pb_var'])
            entry['mean_expression'] = df_slice.sum(axis=0).mean()
        except Exception as e:
            print(e)

        ic_df = ic_df.append(entry,ignore_index=True)

#%%
ic_df.to_csv(ic_filename,index=False)