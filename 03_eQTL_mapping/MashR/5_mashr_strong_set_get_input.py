#!/usr/bin/env python
import pandas as pd
import glob
import os

all_results = []
tissues = ["Mu", "BF", "Li", "AF"]

for tissue in tissues:
    tissue_files = glob.glob(f"./association/{tissue}_*_strong.qassoc")
    
    for file in tissue_files:
        basename = os.path.basename(file)
        gene = basename.replace(f"{tissue}_", "").replace("_strong.qassoc", "")
        df = pd.read_csv(file, delim_whitespace=True)
        df = df[["SNP", "BETA", "SE"]].copy()
        df["Z"] = df["BETA"] / df["SE"]
        df = df[["SNP", "Z"]]
        df["gene"] = gene
        df["tissue"] = tissue
        all_results.append(df)

long_df = pd.concat(all_results, ignore_index=True)
wide_df = long_df.pivot_table(index=["SNP", "gene"], columns="tissue", values="Z")
wide_df = wide_df.reindex(columns=tissues)
wide_df.columns = [f"Z_{t}" for t in tissues]
wide_df = wide_df.reset_index()
wide_df.to_csv("./strong_set.z.txt", sep="\t", index=False, na_rep="NA")

