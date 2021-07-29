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
import string

# variables
SAVE_PARAMS = {'sep':'\t', 'index':False}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
raw_psi_file = os.path.join(RAW_DIR,'inhouse','CCLE','CCLE_ALL_INCLUSION_LEVELS_FULL-Hs2935-minN_1-minSD_0-noVLOW-min_ALT_use25-Tidy.tab')
raw_metadata_file = os.path.join(ROOT,'data','raw','DepMap','achilles_ccle','sample_info.csv')
ccle_cancertypes_file = os.path.join(ROOT,'data','raw','articles','Yu2019','ccle_metadata.xls')
"""

##### FUNCTIONS #####
def load_data(psi_file, metadata_file, ccle_cancertypes_file):
    psi = pd.read_table(psi_file, index_col=0)
    metadata = pd.read_csv(metadata_file)
    cancertypes = pd.read_excel(ccle_cancertypes_file)

    return psi, metadata, cancertypes


def format_cell_names(x):
    x = pd.Series(x)
    x = (
        x.str.replace("^.+?\\.", "")
        .str.replace("\\.+?.", "")
        .str.replace(rf'[{string.punctuation}]', '')
        .str.upper()
    )

    return x.values


def preprocess_ccle(psi, metadata, cancertypes):
    # add cancertypes to metadata
    metadata = pd.merge(
        metadata,
        cancertypes.rename(
            columns={"CCLE_name": "CCLE_Name", "disease": "cancer_type"}
        ),
        on="CCLE_Name",
    )

    # drop rows missing all values
    is_na = psi.isnull()
    is_all_na = is_na.sum(1) == is_na.shape[0]
    psi = psi.loc[~is_all_na]

    # rename psi columns
    psi.columns = format_cell_names(psi.columns)
    
    # subset
    common_samples = set(psi.columns).intersection(metadata['stripped_cell_line_name'])
    psi = psi[common_samples]
    metadata = metadata.loc[metadata['stripped_cell_line_name'].isin(common_samples)]
    
    # rename
    psi = psi.rename(columns = metadata.set_index('stripped_cell_line_name')['DepMap_ID'].to_dict())
    return psi, metadata


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_psi_file", type=str)
    parser.add_argument("--raw_metadata_file", type=str)
    parser.add_argument("--ccle_cancertypes_file", type=str)
    parser.add_argument("--prep_psi_file", type=str)
    parser.add_argument("--prep_metadata_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    raw_psi_file = args.raw_psi_file
    raw_metadata_file = args.raw_metadata_file
    ccle_cancertypes_file = args.ccle_cancertypes_file
    prep_psi_file = args.prep_psi_file
    prep_metadata_file = args.prep_metadata_file

    # load
    psi, metadata, cancertypes = load_data(raw_psi_file, raw_metadata_file, ccle_cancertypes_file)
    
    psi, metadata = preprocess_ccle(psi, metadata, cancertypes)
    
    # save
    psi.reset_index().to_csv(prep_psi_file, compression='gzip', **SAVE_PARAMS)
    metadata.to_csv(prep_metadata_file, **SAVE_PARAMS)
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
