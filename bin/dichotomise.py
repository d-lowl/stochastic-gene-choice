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
import matplotlib.pyplot as plt
from tqdm import tqdm
from pipeline_utils.vrs_utils import grange, get_vrs, get_cvrs
tqdm.pandas()
import swifter
from scipy.stats.mstats import gmean
try:
    if __IPYTHON__:
        genes_filename = "../dopaminergic.csv"
        bimodal_filename = "../intermediate/_bimodal_genes.csv"
        thresholds_filename = "../intermediate/_optimal_thresholds.csv"
        dichotomisation_type = "vrs"
except:
    import argparse
    parser = argparse.ArgumentParser(description="Calculate optimal bimodality thresholds")
    parser.add_argument("genes_filename")
    parser.add_argument("bimodal_filename")
    parser.add_argument("thresholds_filename")
    parser.add_argument("type")
    args = parser.parse_args()
    genes_filename = args.genes_filename
    bimodal_filename = args.bimodal_filename
    thresholds_filename = args.thresholds_filename
    dichotomisation_type = args.type

# %%
df = pd.read_csv(genes_filename, index_col="gene_id")

# %%
def find_threshold(_gene):
    if(len(_gene[_gene > 0]) == 0):
        return pd.Series({
            "value": np.nan,
            "threshold": 0.5,
            "uncorrected_threshold": np.nan
        })
    if(dichotomisation_type == "3max"):
        if len(_gene[_gene > 0]) < 120:
            threshold = _gene[_gene > 0].sort_values()[-3:].mean() / 10.0
        else:
            threshold = _gene[_gene > 0].quantile(0.975) / 10.0
        return pd.Series({
            "value": np.nan,
            "threshold": np.nanmax((threshold, 0.5)),
            "uncorrected_threshold": threshold
        })
    elif(dichotomisation_type == "geomean"):
        if len(_gene[_gene > 0]) < 120:
            threshold = gmean([
                _gene[_gene > 0].sort_values()[-3:].mean(), # Top three
                _gene[_gene > 0].sort_values()[:3].mean(), # Bottom three
            ])
        else:
            threshold = gmean([
                _gene[_gene > 0].quantile(0.975), # Top quantile
                _gene[_gene > 0].quantile(0.025), # Bottom quantile
            ])
        return pd.Series({
            "value": np.nan,
            "threshold": np.nanmax((threshold, 0.5)),
            "uncorrected_threshold": threshold
        })

    threshold_range = grange(*np.nanquantile(_gene[_gene > 0.0].values,[0.025, 0.975]))
    if(len(threshold_range) == 0):
        value_min = np.nan
        threshold = np.nan
    else:
        if(dichotomisation_type == "vrs"):
            values = [get_vrs(_gene, x) for x in threshold_range]
        elif(dichotomisation_type == "cvrs"):
            values = [get_cvrs(_gene, x) for x in threshold_range]
        else:
            raise Exception("Invalid dichotomisation type")
        if np.isnan(values).all() or (np.nanargmin(values) >= len(threshold_range)):
            value_min = np.nan
            threshold = np.nan
        else:
            value_min = np.nanmin(values)
            threshold = threshold_range[np.nanargmin(values)]
    return pd.Series({
        "value": value_min,
        "threshold": np.nanmax((threshold, 0.5)),
        "uncorrected_threshold": threshold
    }, name = _gene.name)
# %%
if dichotomisation_type == "3max" or dichotomisation_type == "geomean":
    thresholds_df = df.swifter.progress_bar(enable=True).apply(find_threshold, axis=1)
else:
    thresholds_df = df.progress_apply(find_threshold, axis=1)

thresholds_df.to_csv(thresholds_filename)

# %%
transposed_values = df.T.values

transposed_values[transposed_values < thresholds_df.threshold.values] = 0
transposed_values[transposed_values > 0] = 1

# for gene_id, (_, threshold, _) in tqdm(thresholds_df.iterrows(), total=thresholds_df.shape[0]):
#     if np.isnan(threshold):
#         if gene_id in df.index:
#             df.drop(gene_id,axis=0,inplace=True)
#         continue
#     df.loc[gene_id,df.loc[gene_id,:] < threshold] = 0
#     df.loc[gene_id,df.loc[gene_id,:] >= threshold] = 1
df = pd.DataFrame(transposed_values.T, index=df.index, columns=df.columns)

# %%
df.to_csv(bimodal_filename)
