#!/bin/env python3
#%%
try:
    if __IPYTHON__:
        import sys
        sys.path.insert(0,'./bin')
except:
    pass

#%%
import pandas as pd
from pipeline_utils.nextflow_utils import split_filenames
try:
    if __IPYTHON__:
        shuffled_chromosomes_IC_filenames = "[/home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr10_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr11_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr12_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr13_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr14_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr15_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr16_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr17_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr18_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr19_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr1_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr2_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr3_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr4_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr5_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr6_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr7_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr8_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chr9_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chrX_filtered.csv, /home/iakovl0000/stochastic/work/4f/0b8306e38994c5bbdca60e86973e10/chrY_filtered.csv]"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Collect and combine the outputs of chromosomes shuffled IC calculations")
    parser.add_argument("shuffled_chromosomes_IC_filenames")
    args = parser.parse_args()
    shuffled_chromosomes_IC_filenames = args.shuffled_chromosomes_IC_filenames
    
#%%
shuffled_dfs = [pd.read_csv(x) for x in split_filenames(shuffled_chromosomes_IC_filenames)]

for i in range(len(shuffled_dfs)):
    shuffled_dfs[i]["iteration"] = i

pd.concat(shuffled_dfs).to_csv("shuffled_IC.csv", index=False)