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


def load_data(achilles_file, metadata_file):
    achilles = read_file(achilles_file, index_col=0)
    metadata = read_file(metadata_file)
    
    return achilles, metadata


def preprocess_achilles(achilles, metadata):
    # strip gene names
    achilles.index = [symbol for symbol, entrez in achilles.index.str.split(" ")]
    
    # rename samples
    achilles = achilles.rename(columns=metadata.set_index("CCLE_Name")["DepMap_ID"].to_dict())
    
    # drop rows missing many values
    is_na = achilles.isnull()
    non_missing = is_na.shape[1] - is_na.sum(1)
    to_keep = non_missing >= THRESH_NON_MISSING
    achilles = achilles.loc[to_keep].copy()
    
    return achilles


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_achilles_file", type=str)
    parser.add_argument("--raw_metadata_file", type=str)
    parser.add_argument("--prep_achilles_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    raw_achilles_file = args.raw_achilles_file
    raw_metadata_file = args.raw_metadata_file
    prep_achilles_file = args.prep_achilles_file

    # load
    print('Loading data...')
    achilles, metadata = load_data(raw_achilles_file, raw_metadata_file)
    print('Preprocessing data...')
    achilles = preprocess_achilles(achilles, metadata)
    
    # save
    print('Saving data...')
    achilles.reset_index().to_csv(prep_achilles_file, **SAVE_PARAMS)
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
