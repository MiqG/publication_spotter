#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Count mutations per sample by effect. 
#

import pandas as pd
import argparse

"""
Development
-----------
import os
ROOT = os.path.join('~','projects','publication_splicing_dependency')
RAW_DIR = os.path.join(ROOT,'data','raw')
snv_file = os.path.join(RAW_DIR,'DepMap','achilles_ccle','CCLE_mutations.csv')
"""

##### FUNCTIONS #####
def load_data(snv_file):
    if snv_file.endswith('.csv'):
        snv = pd.read_csv(snv_file, low_memory=False)
    else:
        snv = pd.read_table(snv_file, low_memory=False)
    return snv


def count_mutations_per_gene_and_effect(snv, id_col, groupby_cols):
    """
    mutations per gene and effect
    if a sample has multiple mutations in 
    the same gene it is counted only once
    """
    
    counts = snv[[id_col]+groupby_cols]\
                .drop_duplicates()\
                .groupby(groupby_cols)\
                .size()\
                .reset_index()\
                .rename(columns={0:'n'})
    
    return counts


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--snv_file',type=str)
    parser.add_argument('--id_col',type=str)
    parser.add_argument('--groupby_cols',type=str)
    parser.add_argument('--output_file',type=str)
    
    args = parser.parse_args()
    return args


def main():
    # parse arguments
    args = parse_args()
    
    snv_file = args.snv_file
    id_col = args.id_col
    groupby_cols = args.groupby_cols.split(',')
    output_file = args.output_file
    
    # load
    snv = load_data(snv_file)
    result = count_mutations_per_gene_and_effect(snv, id_col, groupby_cols)
    
    # save
    result.to_csv(output_file, sep='\t', index=False, compression='gzip')


##### SCRIPT #####
if __name__ == '__main__':
    main()
    print('Done!')

