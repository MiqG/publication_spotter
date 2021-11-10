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

import os
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from joblib import Parallel, delayed

SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}

"""
Development
-----------
ROOT = '~/projects/publication_splicing_dependency'
PREP_DIR = os.path.join(ROOT,'data','prep')
RESULTS_DIR = os.path.join(ROOT,'results','model_splicing_dependency')
psi_file = os.path.join(PREP_DIR,'event_psi','LUAD.tsv')
genexpr_file = os.path.join(PREP_DIR,'gene_counts','LUAD.tsv')
models_file = os.path.join(RESULTS_DIR,'files','models_gene_dependency-EX','model_summaries.tsv.gz')
models_coefs_dir = os.path.join(RESULTS_DIR,'files','models_gene_dependency-EX')
gene_lengths_file = os.path.join(ROOT,'data','raw','VastDB','assemblies','Hs2','EXPRESSION','Hs2_mRNA-50-SS.eff')
normalize_counts = True
"""

##### FUNCTIONS #####
def compute_tpm(gene_counts, gene_lengths):
    X = gene_counts / gene_lengths.loc[gene_counts.index].values
    tpm = 1e6 * X / X.sum(axis=0)

    return np.log2(tpm + 1)


def load_data(
    models_file,
    models_coefs_dir,
    psi_file,
    genexpr_file,
    normalize_counts,
    gene_lengths_file,
):
    # read
    models = pd.read_table(models_file).set_index(["EVENT", "ENSEMBL"])
    coefs = {
        "event": pd.read_pickle(
            os.path.join(models_coefs_dir, "coefs_event.pickle.gz")
        ),
        "gene": pd.read_pickle(os.path.join(models_coefs_dir, "coefs_gene.pickle.gz")),
        "interaction": pd.read_pickle(
            os.path.join(models_coefs_dir, "coefs_interaction.pickle.gz")
        ),
        "intercept": pd.read_pickle(
            os.path.join(models_coefs_dir, "coefs_intercept.pickle.gz")
        ),
    }
    psi = pd.read_table(psi_file, index_col=0)
    genexpr = pd.read_table(genexpr_file, index_col=0)

    # subset
    ## common event - genes
    event_gene = coefs["event"][["EVENT", "GENE", "ENSEMBL"]]
    events_avail = set(event_gene["EVENT"])
    genes_avail = set(event_gene["ENSEMBL"])

    common_events = events_avail.intersection(psi.index).intersection(
        event_gene.loc[event_gene["ENSEMBL"].isin(genexpr.index), "EVENT"]
    )
    common_genes = genes_avail.intersection(genexpr.index).intersection(
        event_gene.loc[event_gene["EVENT"].isin(psi.index), "ENSEMBL"]
    )

    ## common elements
    common_samples = set(psi.columns).intersection(genexpr.columns)

    coefs = {
        k: coef.loc[
            coef["EVENT"].isin(common_events) & coef["ENSEMBL"].isin(common_genes)
        ].copy()
        for k, coef in coefs.items()
    }
    psi = psi.loc[common_events, common_samples].copy()
    genexpr = genexpr.loc[common_genes, common_samples].copy()

    if normalize_counts:
        print("Normalizing counts...")
        gene_lengths = pd.read_table(gene_lengths_file, header=None, index_col=0)
        gene_lengths.columns = ["length"]
        genexpr = compute_tpm(genexpr, gene_lengths)

    # standardize PSI and TPMs
    event_mean = models.loc[(psi.index, slice(None)), "event_mean"].values.reshape(
        -1, 1
    )
    event_std = models.loc[(psi.index, slice(None)), "event_std"].values.reshape(-1, 1)
    psi = (psi - event_mean) / event_std

    gene_mean = (
        models.loc[(slice(None), genexpr.index), "gene_mean"]
        .reset_index(["ENSEMBL"])
        .drop_duplicates()["gene_mean"]
        .values.reshape(-1, 1)
    )
    gene_std = (
        models.loc[(slice(None), genexpr.index), "gene_std"]
        .reset_index(["ENSEMBL"])
        .drop_duplicates()["gene_std"]
        .values.reshape(-1, 1)
    )
    genexpr = (genexpr - gene_mean) / gene_std

    return coefs, psi, genexpr


def compute_single_splicing_dependency(
    b_event, b_gene, b_interaction, b_intercept, x_psi, x_genexpr
):

    samples = x_psi.index
    event = x_psi.name

    PSI = x_psi.values.reshape(1, -1)
    TPM = x_genexpr.values.reshape(1, -1)
    PROD = PSI * TPM

    # compute
    y = b_intercept + b_event * PSI + b_gene * TPM + b_interaction * PROD

    # summarize
    mean = pd.Series(np.mean(y, axis=0), index=samples, name=event)
    median = pd.Series(np.median(y, axis=0), index=samples, name=event)
    std = pd.Series(np.std(y, axis=0), index=samples, name=event)

    summary = {"mean": mean, "median": median, "std": std}

    return summary


def compute_splicing_dependency(coefs, psi, genexpr, n_jobs):
    # unpack coefficients
    coefs_event = coefs["event"].drop(columns=["GENE"]).set_index(["EVENT", "ENSEMBL"])
    coefs_gene = coefs["gene"].drop(columns=["GENE"]).set_index(["EVENT", "ENSEMBL"])
    coefs_interaction = (
        coefs["interaction"].drop(columns=["GENE"]).set_index(["EVENT", "ENSEMBL"])
    )
    coefs_intercept = (
        coefs["intercept"].drop(columns=["GENE"]).set_index(["EVENT", "ENSEMBL"])
    )

    # predict splicing dependency for each combination of parameters
    event_gene = coefs_event.index.to_frame()

    result = Parallel(n_jobs=n_jobs)(
        delayed(compute_single_splicing_dependency)(
            b_event=coefs_event.loc[(event, gene)].values.reshape(-1, 1),
            b_gene=coefs_gene.loc[(event, gene)].values.reshape(-1, 1),
            b_interaction=coefs_interaction.loc[(event, gene)].values.reshape(-1, 1),
            b_intercept=coefs_intercept.loc[(event, gene)].values.reshape(-1, 1),
            x_psi=psi.loc[event],
            x_genexpr=genexpr.loc[gene],
        )
        for event, gene in tqdm(event_gene[["EVENT", "ENSEMBL"]].values)
    )
    spldep_mean = pd.DataFrame([r["mean"] for r in result])
    spldep_median = pd.DataFrame([r["median"] for r in result])
    spldep_std = pd.DataFrame([r["std"] for r in result])

    return spldep_mean, spldep_median, spldep_std


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psi_file", type=str)
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--models_file", type=str)
    parser.add_argument("--models_coefs_dir", type=str)
    parser.add_argument("--normalize_counts", type=bool, default=False)
    parser.add_argument("--spldep_mean_file", type=str)
    parser.add_argument("--spldep_median_file", type=str)
    parser.add_argument("--spldep_std_file", type=str)
    parser.add_argument("--gene_lengths_file", type=str)
    parser.add_argument("--n_jobs", type=int)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    psi_file = args.psi_file
    genexpr_file = args.genexpr_file
    models_file = args.models_file
    models_coefs_dir = args.models_coefs_dir
    normalize_counts = args.normalize_counts
    gene_lengths_file = args.gene_lengths_file
    spldep_mean_file = args.spldep_mean_file
    spldep_median_file = args.spldep_median_file
    spldep_std_file = args.spldep_std_file
    n_jobs = args.n_jobs

    print("Loading data...")
    coefs, psi, genexpr = load_data(
        models_file,
        models_coefs_dir,
        psi_file,
        genexpr_file,
        normalize_counts,
        gene_lengths_file,
    )

    print("Computing splicing dependencies...")
    spldep_mean, spldep_median, spldep_std = compute_splicing_dependency(
        coefs, psi, genexpr, n_jobs
    )

    print("Saving results...")
    spldep_mean.reset_index().to_csv(spldep_mean_file, **SAVE_PARAMS)
    spldep_median.reset_index().to_csv(spldep_median_file, **SAVE_PARAMS)
    spldep_std.reset_index().to_csv(spldep_std_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
