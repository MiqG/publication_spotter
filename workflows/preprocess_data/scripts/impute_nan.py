# impute missing values with various methods:
# - knn
# - pmm

import argparse
import json
import pandas as pd

# imputation methods
from sklearn.impute import KNNImputer
from autoimpute.imputations import MultipleImputer

# variables
SAVE_PARAMS = {'sep':'\t', 'compression':'gzip', 'index':False}

##### FUNCTIONS #####
def load_data(input_file):
    data = pd.read_table(input_file, index_col=0)
    return data


def get_imputation_method(method, method_kws):
    if method == 'knn':
        imputation_method = KNNImputer(**method_kws)
    elif method == 'pmm':
        imputation_method = MultipleImputer(strategy=method,**method_kws)
    return imputation_method


def impute_nan(data, method, method_kws, features_as_rows=True):
    method = get_imputation_method(method, method_kws)
    if features_as_rows:
        imputed = method.fit_transform(data.T).T
    else:
        imputed = method.fit_transform(data)
    imputed = pd.DataFrame(imputed, index = data.index, columns = data.columns)
    return imputed


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_file', type=str)
    parser.add_argument('--output_file', type=str)
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
    method = args.method
    method_kws = json.loads(args.method_kws)
    features_as_rows = args.features_as_rows
        
    # run
    print('Loading data...')
    data = load_data(input_file)
    print('Imputing data...')
    result = impute_nan(data, method, method_kws, features_as_rows)
    
    # save
    result.reset_index().to_csv(output_file, **SAVE_PARAMS)
    
    
##### SCRIPT #####
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
method_kws = json.loads('{"n_neighbors":10, "weights":"distance"}')
features_as_rows = True
"""
if __name__ == '__main__':
    main()
    print('Done!')
