import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def grange(a,b,step=1.2):
    """
    Returns a list between a and b of geometrically progressing series 
    """
    r = []
    while(a < b):
        r += [a]
        a *= step
    return r

#Treat genes separetely, therefore input is a Series
def get_var(gene):
    return ((gene - gene.mean()) ** 2).sum()
    # return (df.subtract(df.mean(axis=1), axis=0) ** 2).sum(axis=1)

def get_vrs(gene, threshold):
    try:
        var = get_var(gene)
        var1 = get_var(gene[gene < threshold])
        var2 = get_var(gene[gene >= threshold])
        return (var1 + var2) / var
    except:
        return np.nan

def get_cvrs(gene, threshold):
    try:
        var = get_var(gene) / (gene.mean() ** 2)
        var1 = get_var(gene[gene < threshold]) / (gene[gene < threshold].mean() ** 2)
        var2 = get_var(gene[gene >= threshold]) / (gene[gene >= threshold].mean() ** 2)
        return (var1 + var2) / var
    except:
        return np.nan