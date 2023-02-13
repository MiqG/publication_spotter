#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Preprocess and clean CCLE data.

import argparse
import pandas as pd

# variables
SAVE_PARAMS = {'sep':'\t', 'index':False, 'compression':'gzip'}
THRESH_NON_MISSING = 50


"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
"""

##### FUNCTIONS #####
def read_file(input_file, **kws):
    if 'csv' in input_file:
        data = pd.read_csv(input_file, **kws)
    elif 'tsv' in input_file:
        data = pd.read_table(input_file, **kws)
    else:
        print('Wrong input file format.')
        data = None
    return data


def load_data(depmap_file, metadata_file):
    depmap = read_file(depmap_file, index_col=0)
    metadata = read_file(metadata_file)
    
    return depmap, metadata


def preprocess_depmap(depmap, metadata):
    # strip gene names
    depmap.index = [symbol for symbol, entrez in depmap.index.str.split(" ")]
    
    # rename samples
    depmap = depmap.rename(columns=metadata.set_index("CCLE_Name")["DepMap_ID"].to_dict())
    
    # drop rows missing many values
    is_na = depmap.isnull()
    non_missing = is_na.shape[1] - is_na.sum(1)
    to_keep = non_missing >= THRESH_NON_MISSING
    depmap = depmap.loc[to_keep].copy()
    
    return depmap


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_depmap_file", type=str)
    parser.add_argument("--raw_metadata_file", type=str)
    parser.add_argument("--prep_depmap_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    raw_depmap_file = args.raw_depmap_file
    raw_metadata_file = args.raw_metadata_file
    prep_depmap_file = args.prep_depmap_file

    # load
    print('Loading data...')
    depmap, metadata = load_data(raw_depmap_file, raw_metadata_file)
    print('Preprocessing data...')
    depmap = preprocess_depmap(depmap, metadata)
    
    # save
    print('Saving data...')
    depmap.reset_index().to_csv(prep_depmap_file, **SAVE_PARAMS)
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
