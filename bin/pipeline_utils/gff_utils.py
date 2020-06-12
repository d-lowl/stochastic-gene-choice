import pandas as pd
import re

def read_gff3(filename):
    with open(filename, 'r') as f:
        gff3_lines = list(f)
    gff3 = pd.read_csv(
        filename,
        sep="\t",
        skiprows=(lambda x: gff3_lines[x].startswith("#")),
        header=None,
        names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    )
    return gff3

def get_genes_from_gff3(gff3):
    df_slice = gff3[(gff3.type == "gene") & (gff3.source.str.contains("RefSeq") | gff3.source.str.contains("Curated"))].copy()
    def get_attribute(x, attribute_name):
        try:
            return re.search(r'(?:{}=)(.*?)(?:;)'.format(attribute_name), x).groups()[0]
        except:
            return ""
    for attr in ["Name", "description"]:
        df_slice[attr] = df_slice.attributes.apply(lambda x: get_attribute(x, attr))
    return df_slice.drop("attributes", axis=1)