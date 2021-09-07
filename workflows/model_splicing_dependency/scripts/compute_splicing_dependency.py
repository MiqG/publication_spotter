#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Compute splicing dependency using the parameters from the fitted linear
# models to predict gene dependency from splicing PSI and gene expression
# log2(TPM+1)
#
# gene_dependency = intercept + psi + genexpr + psi*genexpr
# splicing_dependency = intercept + psi + psi*genexprs

import pandas as pd
import argparse

"""
Development
-----------
import os
ROOT = '~/projects/publication_splicing_dependency'
PREP_DIR = os.path.join(ROOT,'data','prep')
RESULTS_DIR = os.path.join(ROOT,'results','model_splicing_dependency')
psi_file = os.path.join(PREP_DIR,'event_psi','CCLE-EX.tsv.gz')
genexpr_file = os.path.join(PREP_DIR,'genexpr_tpm','CCLE.tsv.gz')
models_file = os.path.join(RESULTS_DIR,'files','models_gene_dependency-EX.tsv.gz')
"""

##### FUNCTIONS #####
def load_data(models_file, psi_file, genexpr_file):
    models = pd.read_table(models_file)
    psi = pd.read_table(psi_file, index_col=0)
    genexpr = pd.read_table(genexpr_file, index_col=0)

    # common elements
    common_samples = set(psi.columns).intersection(genexpr.columns)
    psi = psi.loc[:, common_samples].copy()
    genexpr = genexpr.loc[:, common_samples].copy()

    return models, psi, genexpr


def compute_splicing_dependency(models, psi, genexpr):
    # prep
    event_mean = models["event_mean"].values.reshape(-1, 1)
    event_std = models["event_std"].values.reshape(-1, 1)
    gene_mean = models["gene_mean"].values.reshape(-1, 1)
    gene_std = models["gene_std"].values.reshape(-1, 1)

    # subset
    E = psi.loc[models["EVENT"]].values
    G = genexpr.loc[models["GENE"]].values

    # standardize
    E = (E - event_mean) / event_std
    G = (G - gene_mean) / gene_mean
    EG = E * G

    # compute
    ## coefficients
    intercept = models["intercept_coefficient"].values.reshape(-1, 1)
    e = models["event_coefficient"].values.reshape(-1, 1)
    g = models["gene_coefficient"].values.reshape(-1, 1)
    eg = models["interaction_coefficient"].values.reshape(-1, 1)

    ## computation
    splicing_dependency = intercept + e * E + g * G + eg * EG
    splicing_dependency = pd.DataFrame(
        splicing_dependency, index=models["EVENT"], columns=psi.columns
    )

    return splicing_dependency


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psi_file", type=str)
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--models_file", type=str)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    psi_file = args.psi_file
    genexpr_file = args.genexpr_file
    models_file = args.models_file
    output_file = args.output_file

    print("Loading data...")
    models, psi, genexpr = load_data(models_file, psi_file, genexpr_file)

    print("Computing splicing dependencies...")
    result = compute_splicing_dependency(models, psi, genexpr)

    print("Saving results...")
    result.reset_index().to_csv(output_file, sep="\t", compression="gzip", index=False)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
