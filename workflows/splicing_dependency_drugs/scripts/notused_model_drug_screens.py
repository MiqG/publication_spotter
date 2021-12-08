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
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
from tqdm import tqdm
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs


"""
Development
-----------
import os
ROOT = '~/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
drug_file = os.path.join(RAW_DIR,'DepMap','gdsc','sanger-dose-response.csv')
spldep_file = os.path.join(ROOT,'results','model_splicing_dependency','files','splicing_dependency-EX','mean.tsv.gz')
selected_models_file = os.path.join(ROOT,'results','model_splicing_dependency','files','selected_models-EX.txt')
annotation_file = os.path.join(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
"""

##### FUNCTIONS #####
def load_data(spldep_file, drug_file, annotation_file, selected_models_file=None):
    spldep = pd.read_table(spldep_file, index_col=0)
    drug = pd.read_table(drug_file)
    annotation = pd.read_table(annotation_file)

    # drop undetected & uninformative features
    spldep = spldep.dropna(thresh=2)
    spldep = spldep.loc[spldep.std(axis=1) != 0]

    # subset
    ## events
    if selected_models_file is not None:
        selected_models = list(
            pd.read_table(selected_models_file, header=None)[0].values
        )
        idx = spldep.index.isin(selected_models)
        spldep = spldep.loc[idx].copy()

    ## samples
    common_samples = set(spldep.columns).intersection(drug["ARXSPAN_ID"])
    spldep = spldep.loc[:, common_samples].copy()
    drug = drug.loc[drug["ARXSPAN_ID"].isin(common_samples)].copy()
    annotation = annotation.loc[annotation["EVENT"].isin(spldep.index)].copy()

    gc.collect()

    return spldep, drug, annotation


def get_drug_pcs(drug):
    drugmat = drug.pivot_table(
        index="DRUG_ID",
        columns="ARXSPAN_ID",
        values="IC50_PUBLISHED",
        aggfunc=np.median,
    )
    drugmat = drugmat.apply(lambda x: x.fillna(np.median(x.dropna())), axis=1)
    pca = PCA(1)
    pca.fit(drugmat)
    pcs = pd.DataFrame(
        pca.components_.T,
        index=drugmat.columns,
        columns=["PC%s" % (n + 1) for n in range(pca.n_components_)],
    )
    return pcs


def get_model_lr_info(model):
    llf = model.lml()
    rank = np.linalg.matrix_rank(model.X)
    df_resid = model.nsamples - rank
    return llf, df_resid


def compare_lr_test(model_null, model_alt):
    llf_null, df_null = get_model_lr_info(model_null)
    llf_alt, df_alt = get_model_lr_info(model_alt)

    lrdf = df_null - df_alt
    lrstat = -2 * (llf_null - llf_alt)
    lr_pvalue = stats.chi2.sf(lrstat, lrdf)

    return lrstat, lr_pvalue, lrdf


def fit_limixmodel(y, X, sigma):
    event = X.columns[0]

    # fit model
    if sigma is not None:
        QS = economic_qs(sigma.loc[y.index, y.index])
    else:
        QS = None
    model = LMM(y, X, QS)
    model.fit(verbose=False)

    # get rsquared
    pred = np.dot(X, model.beta)
    rsquared = metrics.r2_score(y, pred)

    # likelihood ratio test
    model_null = LMM(y, X[["PC1", "intercept"]], QS)
    model_null.fit(verbose=False)
    lr_stat, lr_pvalue, lr_df = compare_lr_test(model_null, model)

    # score using test data
    pearson_coef, pearson_pvalue = stats.pearsonr(pred, y)
    spearman_coef, spearman_pvalue = stats.spearmanr(pred, y)

    # prepare output
    ## get info
    params = pd.Series(model.beta, index=X.columns)
    scanner = model.get_fast_scanner()
    bse = pd.Series(scanner.null_beta_se, index=X.columns)
    zscores = pd.Series(params / bse, index=X.columns)
    pvalues = pd.Series(stats.norm.sf(np.abs(zscores)) * 2, index=X.columns)

    # make summary
    summary = {
        "DRUG_ID": np.nan,
        "EVENT": np.nan,
        "ENSEMBL": np.nan,
        "GENE": np.nan,
        "spldep_coefficient": params[event],
        "spldep_stderr": bse[event],
        "spldep_zscore": zscores[event],
        "spldep_pvalue": pvalues[event],
        "growth_coefficient": params["PC1"],
        "growth_stderr": bse["PC1"],
        "growth_zscore": zscores["PC1"],
        "growth_pvalue": pvalues["PC1"],
        "intercept_coefficient": params["intercept"],
        "intercept_stderr": bse["intercept"],
        "intercept_zscore": zscores["intercept"],
        "intercept_pvalue": pvalues["intercept"],
        "n_obs": model.nsamples,
        "rsquared": rsquared,
        "pearson_correlation": pearson_coef,
        "pearson_pvalue": pearson_pvalue,
        "spearman_correlation": spearman_coef,
        "spearman_pvalue": spearman_pvalue,
        "lr_stat": lr_stat,
        "lr_pvalue": lr_pvalue,
        "lr_df": lr_df,
    }
    summary = pd.Series(summary)

    return summary


def fit_single_model(y, X, sigma, method):
    methods = {"limix": fit_limixmodel}
    summary = methods[method](y, X, sigma)
    return summary


def fit_model(y_drug, x_spldep, x_pcs, sigma, ensembl, gene, method):

    X = pd.concat([x_spldep, x_pcs["PC1"]], axis=1)
    y = y_drug

    # dropna
    is_nan = X.isnull().any(1) | y.isnull()
    X = X.loc[~is_nan].copy()
    y = y[~is_nan].copy()

    try:
        # standardize features
        X.values[:, :] = StandardScaler().fit_transform(X)
        X["intercept"] = 1.0

        summary = fit_single_model(y, X, sigma, method)

    except:
        X["intercept"] = np.nan

        # create empty summary
        summary = pd.Series(
            np.nan,
            index=[
                "DRUG_ID",
                "EVENT",
                "ENSEMBL",
                "GENE",
                "spldep_coefficient",
                "spldep_stderr",
                "spldep_zscore",
                "spldep_pvalue",
                "growth_coefficient",  # PC1
                "growth_stderr",
                "growth_zscore",
                "growth_pvalue",
                "intercept_coefficient",
                "intercept_stderr",
                "intercept_zscore",
                "intercept_pvalue",
                "n_obs",
                "rsquared",
                "pearson_correlation",
                "pearson_pvalue",
                "spearman_correlation",
                "spearman_pvalue",
                "lr_stat",
                "lr_pvalue",
                "lr_df",
            ],
        )
    # update
    summary["DRUG_ID"] = y_drug.name
    summary["EVENT"] = x_spldep.name
    summary["ENSEMBL"] = ensembl
    summary["GENE"] = gene
    # add
    summary["spldep_mean"] = x_spldep.mean()
    summary["spldep_std"] = x_spldep.std()
    summary["growth_mean"] = x_pcs.mean().values[0]  # PC1
    summary["growth_std"] = x_pcs.std().values[0]  # PC1

    return summary


def fit_models(spldep, drug, pcs, annotation, n_jobs):
    sigma = spldep.cov()
    drugs_oi = drug["DRUG_ID"].unique()
    results = []
    for drug_oi in drugs_oi:
        print(drug_oi)

        # prepare drug target variable
        y_drug = drug.loc[drug["DRUG_ID"] == drug_oi]
        y_drug = pd.Series(
            np.log(y_drug["IC50_PUBLISHED"].values),  # log-normalize
            index=y_drug["ARXSPAN_ID"].values,
            name=drug_oi,
        )

        # run against all events
        res = Parallel(n_jobs=n_jobs)(
            delayed(fit_model)(
                y_drug,
                spldep.loc[event, y_drug.index],
                pcs.loc[y_drug.index],
                sigma,
                ensembl,
                gene,
                method="limix",
            )
            for event, ensembl, gene in tqdm(annotation.values)
        )
        res = pd.DataFrame(res)

        # compute adjusted p-values
        res["lr_padj"] = np.nan
        idx = ~res["lr_pvalue"].isnull()
        res.loc[idx, "lr_padj"] = sm.stats.multipletests(
            res.loc[idx, "lr_pvalue"], method="fdr_bh"
        )[1]
        # save
        results.append(res)

    results = pd.concat(results)
    results["DRUG_ID"] = results["DRUG_ID"].astype("int") # formatting

    return results


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spldep_file", type=str)
    parser.add_argument("--drug_file", type=str)
    parser.add_argument("--selected_models_file", type=str, default=None)
    parser.add_argument("--annotation_file", type=str)
    parser.add_argument("--n_jobs", type=int)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    spldep_file = args.spldep_file
    drug_file = args.drug_file
    selected_models_file = args.selected_models_file
    annotation_file = args.annotation_file
    n_jobs = args.n_jobs
    output_file = args.output_file

    print("Loading data...")
    spldep, drug, annotation = load_data(
        spldep_file, drug_file, annotation_file, selected_models_file
    )

    print("Computing drug profiles PC1...")
    pcs = get_drug_pcs(drug)

    print("Fitting models...")
    result = fit_models(spldep, drug, pcs, annotation, n_jobs)

    print("Saving results...")
    result.to_csv(output_file, sep="\t", compression="gzip", index=False)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
