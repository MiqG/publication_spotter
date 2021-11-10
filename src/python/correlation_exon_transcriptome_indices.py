#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
# Last Update: 2021-02-09
#
# Script purpose
# --------------
# Compute spearman correlation of exon PSI with their sample's proliferation index.
#
# Outline
# -------
# 1. Load data: exon PSI, proliferation indices
# 2. Compute Spearman correlation, keeping coefficient and p-value, p-adj, and logs.
# 3. Save

from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np
import argparse

# variables
PROPERTIES_OI = ['MKI67', 'mitotic_index', 
                 'stemness',
                 'ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE', 
                 'DEPTH-only_tumor',
                 'DEPTH-tumor_and_normal', 
                 'aneuploidy', 'CA20',
                 'median_promotes_adhesion', 'median_represses_adhesion',
                 'adhesion_index']

"""
Development
-----------

"""

##### FUNCTIONS #####
def load_data(omic_file, sample_properties_file, metadata_file, sample_col):
    psi = pd.read_table(omic_file, index_col=0)
    sample_properties = pd.read_table(sample_properties_file, index_col=sample_col).dropna(axis=1, how='all')
    
    if metadata_file is not None: 
        metadata = pd.read_table(metadata_file, index_col=sample_col)
    else: 
        metadata = None
    
    return psi, sample_properties, metadata


def process_inputs(omic_file, sample_properties_file, metadata_file, sample_col, subset_col, subset_values):
    psi, sample_properties, metadata = load_data(omic_file, sample_properties_file, metadata_file, sample_col)
    
    # subset samples
    if subset_col is not None:
        subset_samples = list(metadata.index[metadata[subset_col].isin(subset_values)])
        psi = psi[subset_samples]        
    
    common_samples = set(psi.columns).intersection(sample_properties.index)
    return psi[common_samples], sample_properties.loc[common_samples]


def get_cor_func(cor_method):
    cor_funcs = {
        'spearman': spearmanr
    }
    return cor_funcs[cor_method]


def cor(x, sample_index, cor_method):
    cor_func = get_cor_func(cor_method)
    try:
        corr, pvalue = cor_func(x, sample_index, nan_policy='omit')
    except:
        corr, pvalue = np.nan, np.nan
    output = {'correlation': corr, 'pvalue': pvalue}    
    return pd.Series(output)


def prepare_output(result, padj_method, cor_method):
    # correct for multiple testing
    result['padj'] = np.nan
    is_missing = result['pvalue'].isnull()
    if sum(~is_missing)>1:
        result.loc[~is_missing,'padj'] = multipletests(result.loc[~is_missing,'pvalue'].values, method=padj_method)[1]
    
    # log transform p-values
    result['log10_pvalue'] = - np.log10(result['pvalue'])
    result['log10_padj'] = - np.log10(result['padj'])
    
    # add more info
    result['test_func'] = cor_method
    result['padj_metod'] = padj_method
    
    return result.reset_index()


def compute_correlation(omic_file, sample_properties_file, metadata_file, sample_col, subset_col, subset_values, cor_method, padj_method):
    # load
    psi, sample_properties = process_inputs(omic_file, sample_properties_file, metadata_file, sample_col, subset_col, subset_values)
    
    # correlate
    results = {}
    properties_oi = list(set(sample_properties.columns).intersection(PROPERTIES_OI))
    for index in properties_oi:
        sample_index = sample_properties[index]
        result = psi.apply(lambda x: cor(x, sample_index, cor_method), axis = 1)
        results[index] = result
    
    # output
    result = pd.concat(
        [prepare_output(result, padj_method, cor_method).assign(index_name=index)                
         for index, result in results.items()]
    )
    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--omic_file',type=str)
    parser.add_argument('--sample_properties_file',type=str)
    parser.add_argument('--metadata_file',type=str,default=None)
    parser.add_argument('--sample_col',type=str)
    parser.add_argument('--subset_col',type=str,default=None)
    parser.add_argument('--subset_values',type=str,default=None)
    parser.add_argument('--output_file',type=str)
    parser.add_argument('--cor_method',type=str)
    parser.add_argument('--padj_method',type=str)
    
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    
    omic_file = args.omic_file
    sample_properties_file = args.sample_properties_file
    metadata_file = args.metadata_file
    sample_col = args.sample_col
    subset_col = args.subset_col
    subset_values = args.subset_values.split(',') if args.subset_values is not None else args.subset_values 
    cor_method = args.cor_method
    padj_method = args.padj_method
    output_file = args.output_file
        
    result = compute_correlation(omic_file, sample_properties_file, metadata_file, 
                                 sample_col, subset_col, subset_values, 
                                 cor_method, padj_method)
    
    # save
    result.to_csv(output_file, sep="\t", compression="gzip", index=False)

    
##### SCRIPT #####
if __name__ == '__main__':
    main()
    print('Done!')
