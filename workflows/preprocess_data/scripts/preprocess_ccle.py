#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Preprocess and clean CCLE data.
#
# Outline
# -------
# 1. exchange to official sample identifiers
# 2. Drop rows with all missing values
# 3. Add cancer types to metadata

import argparse
import pandas as pd
import numpy as np
import string

# variables
SAVE_PARAMS = {'sep':'\t', 'index':False, 'compression':'gzip'}
THRESH_NON_MISSING = 50

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')

# metadata
dataset = 'metadata'
sample_info_file = os.path.join(RAW_DIR,'DepMap','achilles_ccle','sample_info.csv')
ccle_cancer_types_file = os.path.join(RAW_DIR,'articles','Yu2019','ccle_metadata.xls')
sample_annot_file = os.path.join(RAW_DIR,'CCLE','ENA_filereport-PRJNA523380-CCLE.tsv')
mat_file, metadata_file, event_type = '', '', ''

# psi
dataset = 'event_psi'
event_type = 'ALTD'
mat_file = os.path.join(RAW_DIR,'CCLE','vast_out','PSI-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab.gz')
metadata_file = os.path.join(PREP_DIR,'metadata','CCLE.tsv.gz')
sample_info_file, ccle_cancer_types_file, sample_annot_file = '','',''

# genexpr
dataset = 'genexpr_tpm'
mat_file = os.path.join(RAW_DIR,'CCLE','vast_out','TPM-hg38-1019.tab.gz')
metadata_file = os.path.join(PREP_DIR,'metadata','CCLE.tsv.gz')
sample_info_file, ccle_cancer_types_file, sample_annot_file, event_type = '','','',''
"""


##### FUNCTIONS #####
def load_data(
        dataset, 
        sample_info_file, ccle_cancer_types_file, sample_annot_file,
        metadata_file, mat_file, event_type
    ):
    
    if dataset=='metadata':
        data = {
            'sample_info': pd.read_csv(sample_info_file),
            'cancer_types': pd.read_excel(ccle_cancer_types_file),
            'sample_annot': pd.read_table(sample_annot_file),
        }
        
    elif dataset=='event_psi':
        data = {
            'metadata': pd.read_table(metadata_file)
        }
        psi = pd.read_table(mat_file, index_col=0)
        psi = psi.loc[psi.index.str.contains(event_type)].copy()
        data['psi'] = psi
        
    elif dataset=='genexpr_tpm':
        data = {
            'metadata': pd.read_table(metadata_file),
            'genexpr': pd.read_table(mat_file, index_col=[0,1])
        }
        
    else:
        print('Invalid dataset.')

    return data


def preprocess_ccle(data, dataset):
    if dataset=='metadata':
        sample_info = data['sample_info']
        sample_annot = data['sample_annot'].rename(
                columns={"sample_alias": "CCLE_Name"}
            )
        sample_annot = sample_annot.loc[sample_annot['library_strategy']=='RNA-Seq']
        cancer_types = data['cancer_types'].rename(
                columns={"CCLE_name": "CCLE_Name", "disease": "cancer_type"}
            )
        
        # combine
        metadata = pd.merge(sample_info, cancer_types, on='CCLE_Name', how='left')
        metadata = pd.merge(metadata, sample_annot[['run_accession','CCLE_Name']], on='CCLE_Name', how='left')
        
        output = metadata
        
    elif dataset=='event_psi':
        psi = data['psi']
        metadata = data['metadata']
        
        # drop rows missing all values
        is_na = psi.isnull()
        non_missing = is_na.shape[1] - is_na.sum(1)
        to_keep = non_missing >= THRESH_NON_MISSING
        psi = psi.loc[to_keep].copy()
        
        # remove suffix from vast-tools
        psi.columns = [c.replace('_1','') for c in psi.columns]
        
        # subset
        common_samples = set(psi.columns).intersection(metadata['run_accession'])
        psi = psi[common_samples]
        metadata = metadata.loc[metadata['run_accession'].isin(common_samples)]

        # rename
        psi = psi.rename(columns = metadata.set_index('run_accession')['DepMap_ID'].to_dict())
        
        output = psi.reset_index()
    
    elif dataset=='genexpr_tpm':
        genexpr = data['genexpr']
        metadata = data['metadata']
        
        # log-transform
        genexpr = np.log2(genexpr + 1)
        
        # remove suffix from vast-tools
        genexpr.columns = [c.replace('_1','') for c in genexpr.columns]
        
        # subset
        common_samples = set(genexpr.columns).intersection(metadata['run_accession'])
        genexpr = genexpr[common_samples]
        metadata = metadata.loc[metadata['run_accession'].isin(common_samples)]

        # rename
        genexpr = genexpr.rename(columns = metadata.set_index('run_accession')['DepMap_ID'].to_dict())
        
        output = genexpr.reset_index().drop(columns='ID')
        
    return output


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, default="")
    # metadata
    parser.add_argument("--sample_info_file", type=str, default="")
    parser.add_argument("--ccle_cancer_types_file", type=str, default="")
    parser.add_argument("--sample_annot_file", type=str, default="")
    # PSI or genexpr matrices
    parser.add_argument("--metadata_file", type=str, default="")
    parser.add_argument("--event_type", type=str, default="")
    parser.add_argument("--mat_file", type=str, default="")
    # generic
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    dataset = args.dataset
    sample_info_file = args.sample_info_file
    ccle_cancer_types_file = args.ccle_cancer_types_file
    sample_annot_file = args.sample_annot_file
    metadata_file = args.metadata_file
    event_type = args.event_type
    mat_file = args.mat_file
    output_file = args.output_file

    # load
    data = load_data(
        dataset, 
        sample_info_file, ccle_cancer_types_file, sample_annot_file,
        metadata_file, mat_file, event_type
    )
    
    output = preprocess_ccle(data, dataset)
    
    # save
    output.to_csv(output_file, **SAVE_PARAMS)
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
