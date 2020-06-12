#!/bin/env python3
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Get a set of family sizes")
parser.add_argument("familyset_filename")
args = parser.parse_args()
familyset_filename = args.familyset_filename

df = pd.read_csv(familyset_filename)
sizes = df.N.unique()

print("\n".join(sizes.astype(int).astype(str)))