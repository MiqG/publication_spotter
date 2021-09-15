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
import numpy as np
import pymc3 as pm
from tqdm import tqdm
import argparse

SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}

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

    # drop events that cannot be inferred
    idx_todrop = (
        models[
            [
                "event_coefficient",
                "event_stderr",
                "gene_coefficient",
                "gene_stderr",
                "intercept_coefficient",
                "intercept_stderr",
            ]
        ]
        .isin([np.nan, -np.inf, np.inf])
        .any(axis=1)
    )
    models = models.loc[~idx_todrop]

    # subset
    ## common event - genes
    common_genes = set(models["ENSEMBL"]).intersection(genexpr.index)
    common_events = set(
        models.loc[models["ENSEMBL"].isin(common_genes), "EVENT"]
    ).intersection(psi.index)

    ## common elements
    common_samples = set(psi.columns).intersection(genexpr.columns)

    models = models.loc[
        models["EVENT"].isin(common_events) & models["ENSEMBL"].isin(common_genes)
    ].copy()
    psi = psi.loc[common_events, common_samples].copy()
    genexpr = genexpr.loc[common_genes, common_samples].copy()

    return models, psi, genexpr


def sample_splicing_dependency(summary, x_psi, x_genexpr, size=1000):
    # prep input data
    ## standardize event splicing and gene expression
    E = x_psi.values
    G = x_genexpr.values
    E = (E - summary["event_mean"]) / summary["event_std"]
    G = (G - summary["gene_mean"]) / summary["gene_mean"]
    EG = E * G
    ## save with the right shape
    event_vals = E.reshape(-1, 1)
    interaction_vals = EG.reshape(-1, 1)

    # coefficients
    coefs_event = (
        pm.Normal.dist(mu=summary["event_coefficient"], sigma=summary["event_stderr"])
        .random(size=size)
        .reshape(1, -1)
    )
    coefs_interaction = (
        pm.Normal.dist(
            mu=summary["interaction_coefficient"], sigma=summary["interaction_stderr"]
        )
        .random(size=size)
        .reshape(1, -1)
    )
    coefs_intercept = (
        pm.Normal.dist(
            mu=summary["intercept_coefficient"], sigma=summary["intercept_stderr"]
        )
        .random(size=size)
        .reshape(1, -1)
    )

    # predict y
    y_pred = (
        coefs_intercept
        + coefs_event * event_vals
        + coefs_interaction * interaction_vals
    )

    # summarize predictions
    mean = pd.Series(y_pred.mean(axis=1), name=x_psi.name, index=x_psi.index)
    std = pd.Series(y_pred.std(axis=1), name=x_psi.name, index=x_psi.index)

    return mean, std


def compute_splicing_dependency(models, psi, genexpr):
    models = models.set_index(["EVENT", "ENSEMBL"])

    spldep_mean = []
    spldep_std = []
    for event, ensembl in tqdm(models.index):
        tmp_mean, tmp_std = sample_splicing_dependency(
            models.loc[(event, ensembl)], psi.loc[event], genexpr.loc[ensembl]
        )
        spldep_mean.append(tmp_mean)
        spldep_std.append(tmp_std)

        del tmp_mean, tmp_std

    spldep_mean = pd.DataFrame(spldep_mean)
    spldep_std = pd.DataFrame(spldep_std)

    return spldep_mean, spldep_std


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psi_file", type=str)
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--models_file", type=str)
    parser.add_argument("--spldep_mean_file", type=str)
    parser.add_argument("--spldep_std_file", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    psi_file = args.psi_file
    genexpr_file = args.genexpr_file
    models_file = args.models_file
    spldep_mean = args.spldep_mean_file
    spldep_std = args.spldep_std_file

    print("Loading data...")
    models, psi, genexpr = load_data(models_file, psi_file, genexpr_file)

    print("Computing splicing dependencies...")
    spldep_mean, spldep_std = compute_splicing_dependency(models, psi, genexpr)

    print("Saving results...")
    spldep_mean.reset_index().to_csv(spldep_mean_file, **SAVE_PARAMS)
    spldep_std.reset_index().to_csv(
        spldep_std_file, sep="\t", compression="gzip", index=False
    )


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
