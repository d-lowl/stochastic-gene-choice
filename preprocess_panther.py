#!/bin/env python3
#%%
import pandas as pd
import re

# %%
df = pd.read_csv("PTHR15.0_mouse.tsv", delimiter="\t", header=0)
df.columns = [
    "gene_id",
    "protein_id",
    "sf_id",
    "family_name",
    "subfamily_name",
    "molecular_function",
    "biological_process",
    "cellular_components",
    "protein_class",
    "pathway"
]

# %%
df.loc[:,"family_id"] = df.sf_id.apply(lambda x: x.split(":")[0])

# %% Extract MGI ids
mgi_id_pattern = re.compile(r'(?:MGI=)(\d+)')
def extract_mgi_id(x):
    match = mgi_id_pattern.search(x)
    if(match is not None):
        return "MGI:{}".format(match.group(1))
df.loc[:,"mgi_gene_id"] = df.gene_id.apply(extract_mgi_id)

#%% Extract Ensembl Ids
ensembl_id_pattern = re.compile(r'(?:Ensembl=)(ENSMUSG\d+)')
def extract_ensembl_id(x):
    match = ensembl_id_pattern.search(x)
    if(match is not None):
        return match.group(1)
df.loc[:,"ensembl_gene_id"] = df.gene_id.apply(extract_ensembl_id)

#%% Extract Uniprot Ids
uniprot_id_pattern = re.compile(r'(?:Gene=)(.+_MOUSE)(?:\|)')
def extract_uniprot_id(x):
    match = uniprot_id_pattern.search(x)
    if(match is not None):
        return match.group(1)

df.loc[:,"uniprot_id"] = df.gene_id.apply(extract_uniprot_id)

#%% Extract EntrezGene Ids
entrez_id_pattern = re.compile(r'(?:GeneID=)(\d+)(?:\|)')
def extract_entrez_id(x):
    match = entrez_id_pattern.search(x)
    if(match is not None):
        return match.group(1)

df.loc[:,"entrez_id"] = df.gene_id.apply(extract_entrez_id)

#%% Extract gene ids, treat as symbols
gene_symbol_pattern = re.compile(r'(?:Gene=)(.+)(?<!_MOUSE)(?:\|)')
def extract_gene_symbol(x):
    match = gene_symbol_pattern.search(x)
    if(match is not None):
        return match.group(1)

df.loc[:,"gene_symbol"] = df.gene_id.apply(extract_gene_symbol)
#%%
df[~df.gene_symbol.isna()].shape
# %%
# df.loc[:,["gene_id","gene_symbol", "entrez_id", "ensembl_id", "mgi_gene_id", "uniprot_id", "family_id", "sf_id", "family_name", "subfamily_name"]]
df.loc[:,["gene_symbol", "entrez_id", "ensembl_gene_id", "mgi_gene_id", "uniprot_id", "family_id", "sf_id", "family_name", "subfamily_name"]].to_csv("preprocessed_panther.csv", index=False)

#%% After gene symbols were converted
df = pd.read_csv("db2db_panther.csv")

# %%
df = df.replace({"-": None})

#%% Drop where no gene symbol is present
df = df.dropna(subset=["gene_symbol"])

#%%
df[df.duplicated("gene_symbol", keep=False)].sort_values("gene_symbol")
#%% Manually fix one duplicate
df.loc[20529,"gene_symbol"] = "A430078G23Rik"

#%% Other three duplicates come from different accession ids from multiple databases, drop
df = df.drop_duplicates("gene_symbol", keep="first")

#%%
df.to_csv("clean_panther.csv", index=False)