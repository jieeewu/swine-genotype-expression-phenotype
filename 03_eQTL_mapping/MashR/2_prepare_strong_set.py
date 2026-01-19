#!/usr/bin/env python
import argparse
import pandas as pd
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument("--tissue", required=True)
args = parser.parse_args()

tissue = args.tissue

def extract_top_snps(folder):
    results = []
    for f in glob.glob(os.path.join(folder, "*.qassoc")):
        gene = os.path.basename(f).replace("Observe.", "").replace(".qassoc", "")
        df = pd.read_csv(f, delim_whitespace=True)
        if "P" not in df.columns or df.empty:
            continue
        top = df.loc[df["P"].idxmin()].copy()
        top["gene"] = gene
        results.append(top)
    return pd.DataFrame(results)

input_dir = f"./eQTL_mapping/eQTL_output/{tissue}"
out_dir = "./eQTL_mapping/MASHR"
os.makedirs(out_dir, exist_ok=True)

out = extract_top_snps(input_dir)
out.to_csv(f"{out_dir}/{tissue}_topSNP.txt", sep="\t", index=False)

