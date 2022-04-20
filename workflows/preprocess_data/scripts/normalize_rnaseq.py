# Script purpose
# --------------
#
# Outline
# -------
#

import os
import pandas as pd
import numpy as np
import argparse
import gc

# variables
SAVE_PARAMS = {"index": False, "compression": "gzip", "sep":"\t"}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_adhesion_index'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
counts_file = os.path.join(PREP_DIR,'genexpr_counts','TCGA','KICH.tsv')
counts_file = os.path.join(PREP_DIR,'genexpr_counts','TabulaSapiens','smartseq2','Bone_Marrow.parquet.gz')
gene_info_file = os.path.join(PREP_DIR,"references","gene_info-hg38.tsv.gz")
geneid_col = "symbol"
genelen_col = "length"
method = "TPM"
"""

##### FUNCTIONS #####
def load_data(counts_file, gene_info_file):
    counts = pd.read_table(counts_file, index_col=0)
    
    if gene_info_file is not None:
        gene_info = pd.read_table(gene_info_file)
    else:
        gene_info = None

    gc.collect()

    return counts, gene_info


def normalize(counts, method, gene_info, geneid_col, genelen_col):
    if method == "CPM":
        norm = 1e6 * (counts / counts.sum(axis=0))
    elif method == "TPM":
        gene_info = gene_info[[geneid_col,genelen_col]].dropna()
        gene_lengths = gene_info.set_index(geneid_col)[genelen_col]
        
        common_genes = list(set(gene_lengths.index).intersection(set(counts.index)))
        
        norm = counts.loc[common_genes] / gene_lengths[
            counts.loc[common_genes].index
        ].values.reshape(-1, 1)
        
        norm = norm.replace([np.inf, -np.inf], np.nan)
        norm = 1e6 * (norm / norm.sum(axis=0))

    return norm


def parse_args():
    """
    For scripting.
    """
    # required
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", type=str)
    parser.add_argument("--counts_file", type=str)
    parser.add_argument("--gene_info_file", type=str, default=None)
    parser.add_argument("--geneid_col", type=str, default=None)
    parser.add_argument("--genelen_col", type=str, default=None)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()
    return args


def main():
    # parse arguments
    args = parse_args()
    method = args.method
    counts_file = args.counts_file
    gene_info_file = args.gene_info_file
    geneid_col = args.geneid_col
    genelen_col = args.genelen_col
    output_file = args.output_file

    # read and prepare data
    print("Reading data...")
    counts, gene_info = load_data(counts_file, gene_info_file)

    # normalize
    print("Normalizing...")
    norm = normalize(counts, method, gene_info, geneid_col, genelen_col)

    # save
    print("Saving...")
    norm.reset_index().to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
