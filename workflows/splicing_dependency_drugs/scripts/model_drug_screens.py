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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
from tqdm import tqdm
from limix.qtl import scan


"""
Development
-----------
import os
ROOT = '~/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
drug_file = os.path.join(RAW_DIR,'DepMap','gdsc','sanger-dose-response.csv')
spldep_file = os.path.join(ROOT,'results','model_splicing_dependency','files','splicing_dependency_mean-EX.tsv.gz')
selected_models_file = os.path.join(ROOT,'results','model_splicing_dependency','files','selected_models-EX.txt')
annotation_file = os.path.join(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
"""

##### FUNCTIONS #####
def load_data(spldep_file, drug_file, annotation_file, selected_models_file=None):
    spldep = pd.read_table(spldep_file, index_col=0)
    drug = pd.read_csv(drug_file)
    annotation = pd.read_table(annotation_file)

    # drop undetected & uninformative events
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
        index="DRUG_NAME",
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
    
    
def fit_limixmodel(y, X, sigma):
    event = X.columns[0]

    # Linear Mixed Model
    m = X[["intercept","PC1"]]  # covariates
    K = sigma.loc[y.index, y.index]  # kinship
    lmm = scan(X[[event]], y, K=K, M=m, lik="normal", verbose=False)
    model = lmm.effsizes["h2"].set_index("effect_name")
    lr_pvalue = lmm.stats["pv20"].values[0]

    # prepare output
    summary = pd.DataFrame(
        {
            "coefficient": model["effsize"],
            "stderr": model["effsize_se"],
            "zscore": model["effsize"] / model["effsize_se"],
            "n_obs": len(y),
            "lr_pvalue": lr_pvalue,
        }
    )

    return summary


def fit_single_model(y, X, sigma, method):
    methods = {"limix": fit_limixmodel}
    summary = methods[method](y, X, sigma)
    return summary


def fit_model(y_drug, x_spldep, x_pcs, sigma, ensembl, gene, method):

    X = pd.concat([x_spldep, x_pcs['PC1']], axis=1)
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
        X["interaction"] = np.nan
        X["intercept"] = np.nan

        summary = pd.DataFrame(
            {
                "coefficient": [np.nan] * X.shape[1],
                "stderr": np.nan,
                "zscore": np.nan,
                "n_obs": np.nan,
                "lr_pvalue": np.nan,
            },
            index=X.columns,
        )

    summary = pd.Series(
        {
            "drug_name": y_drug.name,
            "EVENT": x_spldep.name,
            "ENSEMBL": ensembl,
            "GENE": gene,
            "spldep_coefficient": summary.loc[x_spldep.name, "coefficient"],
            "spldep_stderr": summary.loc[x_spldep.name, "stderr"],
            "spldep_zscore": summary.loc[x_spldep.name, "zscore"],
            "intercept_coefficient": summary.loc["intercept", "coefficient"],
            "intercept_stderr": summary.loc["intercept", "stderr"],
            "intercept_zscore": summary.loc["intercept", "zscore"],
            "lr_pvalue": summary.loc[x_spldep.name, "lr_pvalue"],
            "n_obs": summary.loc[x_spldep.name, "n_obs"],
            "spldep_mean": x_spldep.mean(),
            "spldep_std": x_spldep.std(),
        }
    )

    return summary


def fit_models(spldep, drug, pcs, annotation, n_jobs):
    sigma = spldep.cov()
    drugs_oi = drug["DRUG_NAME"].unique()
    results = []
    for drug_oi in drugs_oi:
        print(drug_oi)

        # prepare drug target variable
        y_drug = drug.loc[drug["DRUG_NAME"] == drug_oi]
        y_drug = pd.Series(
            np.log(y_drug["IC50_PUBLISHED"].values), # log-normalize
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
