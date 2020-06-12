import numpy as np
from scipy.stats.mstats import gmean

def get_individual_gtme(_gene):
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
    return threshold

def get_family_gtme_05(genes):
    _genes = genes.values
    bottom = [np.mean(np.sort(row[row >= 0.5])[:3]) if len(row[row >= 0.5]) < 120 else np.quantile(row[row >= 0.5], 0.025) for row in _genes]
    top = [np.mean(np.sort(row[row >= 0.5])[-3:]) if len(row[row >= 0.5]) < 120 else np.quantile(row[row >= 0.5], 0.975) for row in _genes]
    if(len(bottom) == 0 or len(top) == 0):
        return np.nan
    return gmean([np.nanmin(bottom), np.nanmax(top)])

def get_family_gtme_00(genes):
    _genes = genes.values
    bottom = [np.mean(np.sort(row[row > 0])[:3]) if len(row[row > 0]) < 120 else np.quantile(row[row > 0], 0.025) for row in _genes]
    top = [np.mean(np.sort(row[row > 0])[-3:]) if len(row[row > 0]) < 120 else np.quantile(row[row > 0], 0.975) for row in _genes]
    if(len(bottom) == 0 or len(top) == 0):
        return np.nan
    return gmean([np.nanmin(bottom), np.nanmax(top)])