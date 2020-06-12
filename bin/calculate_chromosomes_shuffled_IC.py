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
from pipeline_utils.ic_utils import get_IC, get_observed_variance, get_pb_variance
from pipeline_utils.nextflow_utils import map_chromosome_filename, split_filenames
tqdm.pandas()
try:
    if __IPYTHON__:
        dichotomised_filename = "intermediate/vrs/filtered_dichotomised_genes.csv"
        chromosome_filenames = "[/home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr10_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr11_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr12_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr13_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr14_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr15_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr16_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr17_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr18_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr19_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr1_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr2_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr3_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr4_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr5_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr6_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr7_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr8_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr9_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chrX_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chrY_filtered.csv]"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Calculate IC for shuffled genome")
    parser.add_argument("dichotomised_filename")
    parser.add_argument("chromosome_filenames")
    args = parser.parse_args()
    dichotomised_filename = args.dichotomised_filename
    chromosome_filenames = args.chromosome_filenames

# %%
chromosome_filenames = {map_chromosome_filename(x): x for x in split_filenames(chromosome_filenames)}
chrs = {k: pd.read_csv(v) for k, v in chromosome_filenames.items()}
dichotomised_df = pd.read_csv(dichotomised_filename, index_col="gene_id")

#%%
_chrs = []
for k in chrs:
    chrs[k].loc[:,"chromosome"] = k
    _chrs += [chrs[k].loc[:,["Name", "chromosome"]]]
chrs_df = pd.concat(_chrs).set_index("Name")
# %%
gene_intersection = chrs_df.index.intersection(dichotomised_df.index)
chrs_df["new_name"] = chrs_df.index
chrs_df.loc[gene_intersection,"new_name"] = chrs_df.loc[gene_intersection].reset_index().Name.sample(frac=1).values
# %%
ics = []
for chr_name in chrs.keys():#list(range(1,20)) + ["X", "Y"]:
    genes = chrs_df[chrs_df.chromosome == str(chr_name)].new_name.values
    total_genes = len(genes)
    for stretch in [7,14,21]:
        for start, end in tqdm(zip(range(0,total_genes-stretch),range(stretch,total_genes))):
            entry = {
                "start": start,
                "end": end,
                "stretch": stretch
            }
            try:
                gene_names = genes[start:end]
                df_slice = dichotomised_df.loc[gene_names,:].dropna(how="all")
                entry['n_genes'] = df_slice.shape[0]
                entry['obs_var'] = get_observed_variance(df_slice)
                entry['pb_var'] = get_pb_variance(df_slice)
                entry['ic'] = get_IC(entry['obs_var'], entry['pb_var'])
                entry['mean_expression'] = df_slice.sum(axis=0).mean()
            except Exception as e:
                print(e)

            ics += [entry]

# %%
ic_df = pd.DataFrame(ics)

# %%
ic_df.to_csv("shuffled_ic.csv", index=False)