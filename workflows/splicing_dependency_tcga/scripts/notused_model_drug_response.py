#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
# Last Update: 2021-02-09
#
# Script purpose
# --------------
#
#

import argparse
import gc
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
from tqdm import tqdm
from limix.qtl import scan


"""
Development
-----------
import os
ROOT = '~/projects/publication_splicing_dependency'
PREP_DIR = os.path.join(ROOT,'data','prep')
MODELS_DIR = os.path.join(ROOT,'results','model_splicing_dependency')
RESULTS_DIR = os.path.join(ROOT,'results','splicing_dependency_tcga')
spldep_file = os.path.join(RESULTS_DIR,'files','LGG','imputed_splicing_dependency_mean-EX.tsv.gz')
psi_file = os.path.join(PREP_DIR,'event_psi_imputed','LGG.tsv.gz')
metadata_file = os.path.join(PREP_DIR,'Moiso2021','drug_response.tsv.gz')
selected_models_file = os.path.join(MODELS_DIR,'files','selected_models-EX.txt')
method = 'RandomForestClassifier'
method_kws = {'n_jobs':1}
"""

##### FUNCTIONS #####
def load_data(spldep_file, psi_file, metadata_file, selected_models_file):
    spldep = pd.read_table(spldep_file, index_col=0)
    psi = pd.read_table(psi_file, index_col=0)
    metadata = pd.read_table(metadata_file)
    selected_events = list(pd.read_table(selected_models_file, header=None)[0].values)

    # drop undetected & uninformative events
    spldep = spldep.loc[spldep.std(axis=1) > 0]
    psi = psi.loc[psi.std(axis=1) > 0]

    # subset
    ## events
    common_events = set(spldep.index).intersection(psi.index)

    ## samples
    idx = metadata["sample_type"] == "Primary Tumor"
    metadata = metadata.loc[idx].dropna()
    common_samples = (
        set(spldep.columns).intersection(psi.columns).intersection(metadata["sampleID"])
    )

    spldep = spldep.loc[common_events, common_samples].copy()
    psi = psi.loc[common_events, common_samples].copy()
    metadata = (
        metadata.loc[
            metadata["sampleID"].isin(common_samples),
            ["sampleID", "treatment", "response"],
        ]
        .drop_duplicates()
        .set_index("sampleID")
    )
    gc.collect()

    return spldep, psi, metadata, selected_events


def fit_model(spldep, psi, metadata, selected_events, method, method_kws):
    y = pd.get_dummies(metadata.loc[spldep.columns,'response'])['NON RESPONDER']
    avail_events = set(selected_events).intersection(spldep.index).intersection(psi.index)
    print('Total selected events found: %s' % len(avail_events))
    
    # fit estimator with imputed Splicing Dependencies
    ## using selected variables
    X = spldep.loc[avail_events].T
    estimator = eval(method)(**method_kws)
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33)
    estimator.fit(X_train, y_train)
    estimator.score(X_test, y_test)
    
    # training
    pred_proba = estimator.predict_proba(X_train)[:,1]
    y_true = y_train
    fpr, tpr, thresholds = metrics.roc_curve(y_true, pred_proba)
    sns.lineplot(x=fpr, y=tpr)
    metrics.roc_auc_score(y_true, pred_proba)

    # test
    pred_proba = estimator.predict_proba(X_test)[:,1]
    y_true = y_test
    fpr, tpr, thresholds = metrics.roc_curve(y_true, pred_proba)
    sns.lineplot(x=fpr, y=tpr)
    metrics.roc_auc_score(y_true, pred_proba)
    
    ## using permuted variables
    sel_events = np.random.choice(list(set(spldep.index) - set(avail_events)), len(avail_events), replace=False)
    X = spldep.loc[sel_events].T
    estimator = eval(method)(**method_kws)
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33)
    estimator.fit(X_train, y_train)
    estimator.score(X_test, y_test)
    
    # training
    pred_proba = estimator.predict_proba(X_train)[:,1]
    y_true = y_train
    fpr, tpr, thresholds = metrics.roc_curve(y_true, pred_proba)
    sns.lineplot(x=fpr, y=tpr)
    metrics.roc_auc_score(y_true, pred_proba)

    # test
    pred_proba = estimator.predict_proba(X_test)[:,1]
    y_true = y_test
    fpr, tpr, thresholds = metrics.roc_curve(y_true, pred_proba)
    sns.lineplot(x=fpr, y=tpr)
    metrics.roc_auc_score(y_true, pred_proba)
    
    # fit estimator with imputed PSI
    ## using selected variables
    
    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spldep_file", type=str)
    parser.add_argument("--psi_file", type=str)
    parser.add_argument("--metadata_file", type=str, default=None)
    parser.add_argument("--selected_models_file", type=str, default=None)
    parser.add_argument("--method", type=str, default=None)
    parser.add_argument("--method_kws", type=str, default=None)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    spldep_file = args.spldep_file
    psi_file = args.psi_file
    selected_models_file = args.selected_models_file
    metadata_file = args.metadata_file
    method = args.method
    method_kws = args.method_kws
    output_file = args.output_file

    print("Loading data...")
    spldep, psi, metadata, selected_events = load_data(
        spldep_file, psi_file, metadata_file, selected_models_file
    )

    print("Fitting model...")
    result = fit_model(spldep, psi, metadata, selected_events, method, method_kws)

    print("Saving results...")
    result.to_csv(output_file, sep="\t", compression="gzip", index=False)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
