#!/usr/bin/env python
import pandas as pd
import glob

dfs = []
for f in glob.glob("./*_topSNP.txt"):
    dfs.append(pd.read_csv(f, sep="\t")[["SNP","gene"]])

strong_pairs = (
    pd.concat(dfs)
    .drop_duplicates()
    .reset_index(drop=True))

strong_pairs.to_csv("./strong_pairs.global.txt", sep="\t", index=False)

