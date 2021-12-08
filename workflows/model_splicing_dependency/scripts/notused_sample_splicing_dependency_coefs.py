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
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse

"""
Development
-----------
ROOT = '~/projects/publication_splicing_dependency'
PREP_DIR = os.path.join(ROOT,'data','prep')
RESULTS_DIR = os.path.join(ROOT,'results','model_splicing_dependency')
models_file = os.path.join(RESULTS_DIR,'files','models_gene_dependency-EX.tsv.gz')
size = 1000
random_seed = 12345
n_jobs = 5
"""

##### FUNCTIONS #####
def load_data(models_file):
    models = pd.read_table(models_file)

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

    return models


def sample_single_splicing_dependency_coefs(summary, size, random_seed):

    np.random.seed(random_seed)

    coefs_event = np.random.normal(
        loc=summary["event_coefficient"], scale=summary["event_stderr"], size=size
    )
    coefs_gene = np.random.normal(
        loc=summary["gene_coefficient"], scale=summary["gene_stderr"], size=size
    )
    coefs_interaction = np.random.normal(
        loc=summary["interaction_coefficient"],
        scale=summary["interaction_stderr"],
        size=size,
    )
    coefs_intercept = np.random.normal(
        loc=summary["intercept_coefficient"],
        scale=summary["intercept_stderr"],
        size=size,
    )

    coefs = {
        "EVENT": summary["EVENT"],
        "GENE": summary["GENE"],
        "ENSEMBL": summary["ENSEMBL"],
        "event": coefs_event,
        "gene": coefs_gene,
        "interaction": coefs_interaction,
        "intercept": coefs_intercept,
    }

    return coefs


def get_coefs(result, coef_oi, size):
    index = ["EVENT", "GENE", "ENSEMBL"] + list(range(size))
    coefs = pd.DataFrame(
        [
            pd.Series(
                [res["EVENT"], res["GENE"], res["ENSEMBL"]] + list(res[coef_oi]),
                index=index,
            )
            for res in result
        ]
    )
    return coefs


def sample_splicing_dependency_coefs(models, size, random_seed, n_jobs):

    result = Parallel(n_jobs=n_jobs)(
        delayed(sample_single_splicing_dependency_coefs)(row, size, random_seed)
        for index, row in tqdm(models.iterrows())
    )

    event = get_coefs(result, "event", size)
    gene = get_coefs(result, "gene", size)
    interaction = get_coefs(result, "interaction", size)
    intercept = get_coefs(result, "intercept", size)

    return event, gene, interaction, intercept


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models_file", type=str)
    parser.add_argument("--output_dir", type=str)
    parser.add_argument("--size", type=int, default=1000)
    parser.add_argument("--random_seed", type=int, default=None)
    parser.add_argument("--n_jobs", type=int, default=None)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    models_file = args.models_file
    output_dir = args.output_dir
    size = args.size
    random_seed = args.random_seed
    n_jobs = args.n_jobs

    os.mkdir(output_dir)

    print("Loading data...")
    models = load_data(models_file)

    print("Sampling coefficients...")
    event, gene, interaction, intercept = sample_splicing_dependency_coefs(
        models, size, random_seed, n_jobs
    )

    print("Saving results...")
    event.to_pickle(os.path.join(output_dir, "event.pickle.gz"))
    gene.to_pickle(os.path.join(output_dir, "gene.pickle.gz"))
    interaction.to_pickle(os.path.join(output_dir, "interaction.pickle.gz"))
    intercept.to_pickle(os.path.join(output_dir, "intercept.pickle.gz"))


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
