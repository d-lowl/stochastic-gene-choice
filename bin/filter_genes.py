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
from pipeline_utils.gff_utils import get_genes_from_gff3, read_gff3
from scipy.stats import kurtosis
try:
    if __IPYTHON__:
        args = {}
        all_genes_filename = "somatosensory_rpkm_suppl.csv"
        filtered_all_genes_filename = "intermediate/filtered_all_genes.csv"
        dichotomisation_type = "cvrs"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Filter the dataset of genes by mean expression and kurtoisis bimodality index")
    parser.add_argument("all_genes_filename")
    parser.add_argument("dichotomised_genes_filename")
    parser.add_argument("filtered_dichotomised_genes_filename")
    parser.add_argument("type")
    parser.add_argument("--no-bi", action="store_true")
    args = parser.parse_args()
    all_genes_filename = args.all_genes_filename
    dichotomised_genes_filename = args.dichotomised_genes_filename
    filtered_dichotomised_genes_filename = args.filtered_dichotomised_genes_filename
    dichotomisation_type = args.type
# %%
df = pd.read_csv(all_genes_filename, index_col="gene_id")
dichotomised_df = pd.read_csv(dichotomised_genes_filename, index_col="gene_id")
# chr18 = pd.read_csv(chromosome_filename, index_col="Name")

# %% Filter by mean
# print("Genes removed due to low mean expression value (< 2**-6):",(df.mean(axis=1) < 2**-6).sum())
# filtered_df = df[df.mean(axis=1) >= 2**-6]


def BI_kurtosis(df):
    n = df.shape[1]
    K = df.kurt(axis=1)
    skewness = df.skew(axis=1)
    return (skewness ** 2 + 1) / (K + (3 * (n-1) ** 2) / ((n-2) * (n-3)))
# %%
if args.no_bi:
    filtered_df = df
else:
    BI = BI_kurtosis(df)
    print("Genes removed due to small Kurtosis BI (< 0.55):",(BI < 0.55).sum())
    filtered_df = df[BI >= 0.55]

#%%
# if dichotomisation_type == "cvrs":
#     print("Genes removed due to small skewness (<= 0):",(filtered_df.skew(axis=1) <= 0).sum())
#     filtered_df = filtered_df[filtered_df.skew(axis=1) > 0]

# %%
print("Total number of genes selected:", filtered_df.shape[0])
dichotomised_df.loc[filtered_df.index].to_csv(filtered_dichotomised_genes_filename)
