#!/usr/bin/env python
import pandas as pd
import glob
import os

tissue_order = ["Mu", "BF", "Li", "AF"]
records = []
for file in glob.glob("./association/*_random.qassoc"):
    basename = os.path.basename(file)
    tissue, gene, _ = basename.split("_", 2)
    if tissue not in tissue_order:
        continue
    df = pd.read_csv(file, delim_whitespace=True)
    df = df[["SNP", "BETA", "SE"]].copy()
    df["Z"] = df["BETA"] / df["SE"]
    df = df[["SNP", "Z"]]
    df["gene"] = gene
    df["tissue"] = tissue
    records.append(df)

long_df = pd.concat(records, ignore_index=True)
wide_df = long_df.pivot_table(index=["SNP", "gene"],columns="tissue",values="Z")
wide_df = wide_df.reindex(columns=tissue_order)
wide_df.columns = [f"Z_{t}" for t in tissue_order]
wide_df = wide_df.reset_index()
wide_df.to_csv("./mashr_input_random_set.z.txt",sep="\t",index=False,na_rep="NA")

