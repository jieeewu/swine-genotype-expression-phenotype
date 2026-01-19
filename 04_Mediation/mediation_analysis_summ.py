#!/usr/bin/env python
import argparse
import pandas as pd

def classify_eQTL(input_tss_file, input_gene_snp_file, input_snp_pos_file, output_file):
    tss_df = pd.read_csv(input_tss_file, sep='\t')
    gene_snp_df = pd.read_csv(input_gene_snp_file, sep='\t')
    gene_snp_df = gene_snp_df.set_index('gene').snps.str.split(',', expand=True).stack().reset_index(level=1, drop=True).reset_index(name='snp')
    snp_pos_df = pd.read_csv(input_snp_pos_file, sep='\t', usecols=['snpID', 'chr', 'pos'])
    snp_pos_df = snp_pos_df.rename(columns={'snpID': 'snp'})
    gene_snp_tss_df = pd.merge(gene_snp_df, tss_df, on='gene')
    merged_df = pd.merge(gene_snp_tss_df, snp_pos_df, on='snp')
    print(merged_df.columns)
	merged_df['distance'] = merged_df.apply(lambda row: abs(row['pos'] - row['TSS']), axis=1)
    print(merged_df.columns)
    print(merged_df['chr_x'].unique())
    print(merged_df['chr_y'].unique())
    merged_df['eQTL_type'] = merged_df.apply(lambda row: 'cis-eQTL' if row['chr_x'] == row['chr_y'] and row['distance'] <= 1000000 else 'trans-eQTL', axis=1)
    print(merged_df.columns)
    merged_df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='classify eQTLs as cis or trans based on distance to gene TSS')
    parser.add_argument('tss_file', help='file containing gene TSS information')
    parser.add_argument('gene_snp_file', help='file containing gene-SNP information')
    parser.add_argument('snp_pos_file', help='file containing SNP position information')
    parser.add_argument('output_file', help='output file')
    args = parser.parse_args()
    classify_eQTL(args.tss_file, args.gene_snp_file, args.snp_pos_file, args.output_file)
