#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
# Last Update: 2021-02-09
#
# Script purpose
# --------------
# Create a linear model to compute the association between exon PSI and
# cell fitness independently of gene expression.
#

import argparse
import gc
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
from tqdm import tqdm


"""
Development
-----------
import os
ROOT = '~/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
psi_file = os.path.join(PREP_DIR,'event_psi','CCLE-EX.tsv.gz')
genexpr_file = os.path.join(PREP_DIR,'genexpr_tpm','CCLE.tsv.gz')
rnai_file = os.path.join(PREP_DIR,'demeter2','CCLE.tsv.gz')
annotation_file = os.path.join(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
"""

##### FUNCTIONS #####
def load_data(psi_file, genexpr_file, rnai_file, annotation_file):
    psi = pd.read_table(psi_file, index_col=0)
    genexpr = pd.read_table(genexpr_file, index_col=0)
    annotation = pd.read_table(annotation_file)
    rnai = pd.read_table(rnai_file, index_col=0)

    gene_annot = annotation[["ENSEMBL", "GENE"]].drop_duplicates().dropna()

    # drop undetected & uninformative events
    psi = psi.dropna(thresh=2)
    psi = psi.loc[psi.std(axis=1) != 0]

    # subset
    common_samples = (
        set(rnai.columns).intersection(psi.columns).intersection(genexpr.columns)
    )

    common_genes = set(
        gene_annot.loc[gene_annot["GENE"].isin(rnai.index), "ENSEMBL"]
    ).intersection(genexpr.index)

    common_events = set(psi.index).intersection(
        annotation.loc[annotation["ENSEMBL"].isin(common_genes), "EVENT"]
    )

    psi = psi.loc[common_events, common_samples]
    genexpr = genexpr.loc[common_genes, common_samples]
    rnai = rnai.loc[
        set(gene_annot.set_index("ENSEMBL").loc[common_genes, "GENE"]), common_samples
    ]
    annotation = annotation.loc[annotation["EVENT"].isin(common_events)]

    gc.collect()

    return psi, genexpr, rnai, annotation


def fit_olsmodel(y, X):
    event, gene = X.columns[:2]

    # split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1)

    # fit linear model to training data
    model = sm.OLS(y_train, X_train).fit()

    # log-likelihood test
    model_null = sm.OLS(y_train, X_train[[gene, "intercept"]]).fit()
    lr_stat, lr_pvalue, _ = model.compare_lr_test(model_null)

    # score using test data
    prediction = model.predict(X_test)
    pearson_coef, pearson_pvalue = stats.pearsonr(prediction, y_test)
    spearman_coef, spearman_pvalue = stats.spearmanr(prediction, y_test)

    # prepare output
    summary = pd.DataFrame(
        {
            "coefficient": model.params,
            "stderr": model.bse,
            "zscore": model.params / model.bse,
            "pvalue": model.pvalues,
            "n_obs": model.nobs,
            "rsquared": model.rsquared,
            "pearson_correlation": pearson_coef,
            "pearson_pvalue": pearson_pvalue,
            "spearman_correlation": spearman_coef,
            "spearman_pvalue": spearman_pvalue,
            "lr_stat": lr_stat,
            "lr_pvalue": lr_pvalue,
        }
    )

    return summary


def fit_single_model(y, X, method):
    methods = {
        "OLS": fit_olsmodel,
    }
    summary = methods[method](y, X)
    return summary


def fit_model(x_psi, x_genexpr, y_rnai, method):

    X = pd.DataFrame([x_psi, x_genexpr]).T
    y = y_rnai

    # dropna
    is_nan = X.isnull().any(1) | y.isnull()
    X = X.loc[~is_nan].copy()
    y = y[~is_nan].copy()

    try:
        # standardize features
        X.values[:, :] = StandardScaler().fit_transform(X)
        X["interaction"] = X[x_psi.name] * X[x_genexpr.name]
        X["intercept"] = 1.0

        summary = fit_single_model(y, X, method)

    except:
        X["interaction"] = np.nan
        X["intercept"] = np.nan

        summary = pd.DataFrame(
            {
                "coefficient": [np.nan] * X.shape[1],
                "stderr": np.nan,
                "zscore": np.nan,
                "pvalue": np.nan,
                "n_obs": np.nan,
                "rsquared": np.nan,
                "pearson_correlation": np.nan,
                "pearson_pvalue": np.nan,
                "spearman_correlation": np.nan,
                "spearman_pvalue": np.nan,
                "lr_stat": np.nan,
                "lr_pvalue": np.nan,
            },
            index=X.columns,
        )

    summary = pd.Series(
        {
            "EVENT": x_psi.name,
            "ENSEMBL": x_genexpr.name,
            "GENE": y_rnai.name,
            "event_coefficient": summary.loc[x_psi.name, "coefficient"],
            "event_stderr": summary.loc[x_psi.name, "stderr"],
            "event_zscore": summary.loc[x_psi.name, "zscore"],
            "event_pvalue": summary.loc[x_psi.name, "pvalue"],
            "gene_coefficient": summary.loc[x_genexpr.name, "coefficient"],
            "gene_stderr": summary.loc[x_genexpr.name, "stderr"],
            "gene_zscore": summary.loc[x_genexpr.name, "zscore"],
            "gene_pvalue": summary.loc[x_genexpr.name, "pvalue"],
            "interaction_coefficient": summary.loc["interaction", "coefficient"],
            "interaction_stderr": summary.loc["interaction", "stderr"],
            "interaction_zscore": summary.loc["interaction", "zscore"],
            "interaction_pvalue": summary.loc["interaction", "pvalue"],
            "intercept_coefficient": summary.loc["intercept", "coefficient"],
            "intercept_stderr": summary.loc["intercept", "stderr"],
            "intercept_zscore": summary.loc["intercept", "zscore"],
            "intercept_pvalue": summary.loc["intercept", "pvalue"],
            "rsquared": summary.loc[x_genexpr.name, "rsquared"],
            "pearson_correlation": summary.loc[x_genexpr.name, "pearson_correlation"],
            "pearson_pvalue": summary.loc[x_genexpr.name, "pearson_pvalue"],
            "spearman_correlation": summary.loc[x_genexpr.name, "spearman_correlation"],
            "spearman_pvalue": summary.loc[x_genexpr.name, "spearman_pvalue"],
            "lr_stat": summary.loc[x_genexpr.name, "lr_stat"],
            "lr_pvalue": summary.loc[x_genexpr.name, "lr_pvalue"],
            "n_obs": summary.loc[x_genexpr.name, "n_obs"],
            "event_mean": x_psi.mean(),
            "event_std": x_psi.std(),
            "gene_mean": x_genexpr.mean(),
            "gene_std": x_genexpr.std(),
        }
    )

    return summary


def fit_models(psi, genexpr, rnai, annotation, n_jobs):
    results = Parallel(n_jobs=n_jobs)(
        delayed(fit_model)(
            psi.loc[event], genexpr.loc[ensembl], rnai.loc[gene], method="OLS"
        )
        for event, ensembl, gene in tqdm(annotation.values)
    )
    results = pd.DataFrame(results)
    results["lr_padj"] = np.nan
    idx = ~results["lr_pvalue"].isnull()
    results.loc[idx, "lr_padj"] = sm.stats.multipletests(
        results.loc[idx, "lr_pvalue"], method="fdr_bh"
    )[1]

    return results


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psi_file", type=str)
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--annotation_file", type=str)
    parser.add_argument("--rnai_file", type=str)
    parser.add_argument("--n_jobs", type=int)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    psi_file = args.psi_file
    genexpr_file = args.genexpr_file
    annotation_file = args.annotation_file
    rnai_file = args.rnai_file
    n_jobs = args.n_jobs
    output_file = args.output_file

    print("Loading data...")
    psi, genexpr, rnai, annotation = load_data(
        psi_file, genexpr_file, rnai_file, annotation_file
    )

    print("Fitting models...")
    result = fit_models(psi, genexpr, rnai, annotation, n_jobs)

    print("Saving results...")
    result.to_csv(output_file, sep="\t", compression="gzip", index=False)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
