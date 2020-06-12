import numpy as np
import pandas as pd

def get_observed_variance(genes):
    return genes.sum(axis=0).var()

def get_pb_variance(genes):
    fx = genes.sum(axis=1) / genes.shape[1]
    return np.sum([(1-p)*p for p in fx])

def get_IC(obs_var, pb_var):
    return obs_var / pb_var

def bootstrap_IC(genes, bootstrap_n = 10000):
    bootstrap_populations = [genes.sample(frac=1.0, replace=True, axis=1) for _ in range(bootstrap_n)]
    pb_vars = np.array([get_pb_variance(x) for x in bootstrap_populations])
    obs_vars = np.array([get_observed_variance(x) for x in bootstrap_populations])
    ics = obs_vars / pb_vars

    return pd.Series({
        "mean": np.nanmean(ics),
        "lower": np.nanquantile(ics, q=0.025),
        "median": np.nanquantile(ics, q=0.5),
        "upper": np.nanquantile(ics, q=0.975),
    })