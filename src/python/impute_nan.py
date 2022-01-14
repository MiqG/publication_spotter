# impute missing values with various methods:
# - knn
# - pmm

import argparse
import json
import pandas as pd
import numpy as np

# imputation methods
from sklearn.impute import KNNImputer

# variables
SAVE_PARAMS = {'sep':'\t', 'compression':'gzip', 'index':False}

"""
Development
-----------
import os
from sso_targets import config
import json
cancer = 'LUAD'
prep_clean_tcga_dir = os.path.join(config.ROOT,'data','prep','clean','TCGA')
prep_imputed_tcga_dir = os.path.join(config.ROOT,'data','prep','imputed','TCGA')

input_file = os.path.join(prep_clean_tcga_dir,'exon_psi',cancer+'.tsv')
output_file = os.path.join(prep_imputed_tcga_dir,'exon_psi',cancer+'.tsv')
method = 'knn'
method_kws = json.loads('{"n_neighbors":5}')
features_as_rows = True
"""

##### FUNCTIONS #####
def load_data(input_file, rows_oi_file=None):
    if 'csv' in input_file:
        data = pd.read_csv(input_file, index_col=0)
    elif 'tsv' in input_file:
        data = pd.read_table(input_file, index_col=0)
    else:
        print('Wrong input file format.')
        data = None
        
    # subset
    if rows_oi_file is not None:
        rows_oi = pd.read_table(rows_oi_file, header=None)[0].values

        print("Subsetting %s rows..." % len(rows_oi))
        data = data.loc[data.index.isin(rows_oi)].copy()
    
    return data


def get_imputation_method(method, method_kws):
    if method == 'knn':
        imputation_method = KNNImputer(**method_kws)
    return imputation_method


def impute_nan(data, method, method_kws, features_as_rows=True):
    method = get_imputation_method(method, method_kws)
    output = np.full(data.shape, np.nan)

    if features_as_rows:
        # do not consider features missing all values
        all_missing = data.isnull().all(axis=1)
        imputed = method.fit_transform(data.loc[~all_missing].T).T
        output[~all_missing,:] = imputed
    else:
        # do not consider features missing all values
        all_missing = data.isnull().all(axis=0)
        imputed = method.fit_transform(data[~all_missing])
        output[:,~all_missing] = imputed
    
    output = pd.DataFrame(output, index = data.index, columns = data.columns)
    return output


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--output_file', type=str)
    parser.add_argument('--rows_oi_file', type=str, default=None)
    parser.add_argument('--method', type=str)
    parser.add_argument('--method_kws', type=str,
                        help='Example: \'{"n_neighbors":10, "weights":"distance"}\'')
    parser.add_argument('--features_as_rows', type=bool, default=True)
    args = parser.parse_args()
    return args
    
    
def main():
    # parse arguments
    args = parse_args()
    input_file = args.input_file
    output_file = args.output_file
    rows_oi_file = args.rows_oi_file
    method = args.method
    method_kws = json.loads(args.method_kws)
    features_as_rows = args.features_as_rows
        
    # run
    print('Loading data...')
    data = load_data(input_file, rows_oi_file)
    
    print('Imputing data...')
    result = impute_nan(data, method, method_kws, features_as_rows)
    
    # save
    result.reset_index().to_csv(output_file, **SAVE_PARAMS)
    
    
##### SCRIPT #####
if __name__ == '__main__':
    main()
    print('Done!')
