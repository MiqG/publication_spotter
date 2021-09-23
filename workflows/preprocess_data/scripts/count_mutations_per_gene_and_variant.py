#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Count mutations per sample by effect. 
#

import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm

"""
Development
-----------
import os
ROOT = os.path.join('~','projects','publication_splicing_dependency')
RAW_DIR = os.path.join(ROOT,'data','raw')
annotation_file = os.path.join(RAW_DIR,'ENSEMBL','gene_annotation-hg19.tsv.gz')
snv_file = os.path.join(RAW_DIR,'DepMap','achilles_ccle','CCLE_mutations.csv')
id_col = 'DepMap_ID'
gene_col = 'Hugo_Symbol'
effect_col = 'Variant_Classification'
"""

##### FUNCTIONS #####
def load_data(snv_file, annotation_file):
    if snv_file.endswith('.csv'):
        snv = pd.read_csv(snv_file, low_memory=False)
    else:
        snv = pd.read_table(snv_file, low_memory=False)
        
    annot = pd.read_table(annotation_file)
    
    return snv, annot


def count_mutations_per_gene_and_effect(snv, id_col, gene_col, effect_col):
    """
    mutations per gene and effect
    if a sample has multiple mutations in 
    the same gene it is counted only once
    """
    counts = snv[[id_col, gene_col, effect_col]]\
                .drop_duplicates()\
                .groupby([gene_col, effect_col])\
                .size()\
                .reset_index()\
                .rename(columns={0:'n'})
    
    return counts


def compute_mutation_entropy(df_gene_effect):
    positions = df_gene_effect['Start_position']
    start, end = int(df_gene_effect['Start'].values[0]), int(df_gene_effect['End'].values[0])
    
    mut_positions, mut_counts = np.unique(positions, return_counts=True)
    mut_freq = pd.Series(index=np.arange(start, end+1), dtype=np.int)
    mut_freq[mut_positions] = mut_counts
    mut_freq = mut_freq.fillna(0)
    mut_freq = mut_freq.reset_index()
    mut_freq.columns = ['position','count']
    mut_freq['p'] = mut_freq['count'] / len(mut_freq) # probability
    mut_freq['partial_entropy'] = mut_freq['p']*np.log2(mut_freq['p']) # in bits
    entropy = -mut_freq['partial_entropy'].sum()
    
    return pd.Series([entropy, len(mut_freq)], index=['entropy','length'])


def in_interval(pos, start, end):
    return (start<=pos) & (pos<=end)
    

def compute_mutation_entropy_per_gene_and_effect(snv, annot, gene_col, effect_col):
    cols_oi = ['Start','End','length']
    df = pd.merge(snv, annot[['Gene']+cols_oi].rename(columns={'Gene':gene_col}), how='left', on='Hugo_Symbol')
    df = df[[gene_col, effect_col, 'Start_position', 'End_position']+cols_oi].dropna()
    df['in_interval'] = in_interval(df['Start_position'], df['Start'], df['End'])
    df = df.loc[df['in_interval']]
    
    entropies = df.groupby([gene_col, effect_col])\
                  .apply(compute_mutation_entropy)\
                  .reset_index()
    entropies['min_entropy'] = -(1/entropies['length'])*np.log2(1/entropies['length'])
    entropies['max_entropy'] = -np.log2(1/entropies['length']) # -len*(1/len)*log2(1/len)
    entropies['max_rel_entropy'] = entropies['entropy'] / entropies['max_entropy']
    entropies['min_rel_entropy'] = entropies['entropy'] / entropies['min_entropy']
    
    return entropies
        

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--snv_file',type=str)
    parser.add_argument('--annotation_file',type=str)
    parser.add_argument('--id_col',type=str)
    parser.add_argument('--gene_col',type=str)
    parser.add_argument('--effect_col',type=str)
    parser.add_argument('--output_file',type=str)
    
    args = parser.parse_args()
    return args


def main():
    # parse arguments
    args = parse_args()
    
    snv_file = args.snv_file
    annotation_file = args.annotation_file
    id_col = args.id_col
    gene_col = args.gene_col
    effect_col = args.effect_col
    output_file = args.output_file
    
    # load
    print('Loading data...')
    snv, annot = load_data(snv_file, annotation_file)
    print('Counting mutations per gene...')
    mut_freq = count_mutations_per_gene_and_effect(snv, id_col, gene_col, effect_col)
    print('Computing entropy per gene...')
    mut_entropy = compute_mutation_entropy_per_gene_and_effect(snv, annot, gene_col, effect_col)
    result = pd.merge(mut_freq, mut_entropy, how='left', on=[gene_col,effect_col])
    result['mut_freq_per_kb'] = (result['n'] / result['length']) * 1000
    
    # save
    result.to_csv(output_file, sep='\t', index=False, compression='gzip')


##### SCRIPT #####
if __name__ == '__main__':
    main()
    print('Done!')

