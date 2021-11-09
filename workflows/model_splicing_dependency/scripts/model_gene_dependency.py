#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
# Last Update: 2021-02-09
#
# Script purpose
# --------------
# Create a linear model to compute the association between exon splicing and
# cell fitness independently of gene expression.
#

import os
import argparse
import gc
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
from tqdm import tqdm
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs


# variables
METHOD = "OLS"
RANDOM_SEED = 1234
TEST_SIZE = 0.15
SAVE_PARAMS = {"sep":"\t", "compression":"gzip", "index":False}
N_ITERATIONS = 100

"""
Development
-----------
ROOT = '~/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
splicing_file = os.path.join(PREP_DIR,'event_psi','CCLE-EX.tsv.gz')
genexpr_file = os.path.join(PREP_DIR,'genexpr_tpm','CCLE.tsv.gz')
gene_dependency_file = os.path.join(PREP_DIR,'demeter2','CCLE.tsv.gz')
mapping_file = os.path.join(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
n_jobs=10
n_iterations=2
"""

##### FUNCTIONS #####
def load_data(splicing_file, genexpr_file, gene_dependency_file, mapping_file):
    splicing = pd.read_table(splicing_file, index_col=0)
    genexpr = pd.read_table(genexpr_file, index_col=0)
    mapping = pd.read_table(mapping_file)
    gene_dependency = pd.read_table(gene_dependency_file, index_col=0)

    gene_annot = mapping[["ENSEMBL", "GENE"]].drop_duplicates().dropna()

    # drop undetected & uninformative events
    splicing = splicing.dropna(thresh=2)
    splicing = splicing.loc[splicing.std(axis=1) != 0]

    # subset
    common_samples = (
        set(gene_dependency.columns).intersection(splicing.columns).intersection(genexpr.columns)
    )

    common_genes = set(
        gene_annot.loc[gene_annot["GENE"].isin(gene_dependency.index), "ENSEMBL"]
    ).intersection(genexpr.index)

    common_events = set(splicing.index).intersection(
        mapping.loc[mapping["ENSEMBL"].isin(common_genes), "EVENT"]
    )

    splicing = splicing.loc[common_events, common_samples]
    genexpr = genexpr.loc[common_genes, common_samples]
    gene_dependency = gene_dependency.loc[
        set(gene_annot.set_index("ENSEMBL").loc[common_genes, "GENE"]), common_samples
    ]
    mapping = mapping.loc[mapping["EVENT"].isin(common_events)]
    
    gc.collect()
    
    # standardize splicing and gene expression
    splicing_mean = splicing.mean(axis=1).values.reshape(-1, 1)
    splicing_std = splicing.std(axis=1).values.reshape(-1, 1)
    splicing = (splicing - splicing_mean) / splicing_std
    
    genexpr_mean = genexpr.mean(axis=1).values.reshape(-1, 1)
    genexpr_std = genexpr.std(axis=1).values.reshape(-1, 1)
    genexpr = (genexpr - genexpr_mean) / genexpr_std
    
    gc.collect()
    
    return splicing, genexpr, gene_dependency, mapping


def get_summary_stats(df, col_oi):
    summary_stats = {
        col_oi + "_mean": np.mean(df[col_oi]),
        col_oi + "_median": np.median(df[col_oi]),
        col_oi + "_std": np.std(df[col_oi]),
        col_oi + "_q25": np.quantile(df[col_oi], 0.25),
        col_oi + "_q75": np.quantile(df[col_oi], 0.75),
    }
    return summary_stats


def get_model_lr_info(model):
    llf = model.lml()
    rank = np.linalg.matrix_rank(model.X)
    df_resid = model.nsamples - rank
    return llf, df_resid
    
    
def compare_lr_test(model_null, model_alt):
    llf_null, df_null = get_model_lr_info(model_null)
    llf_alt, df_alt = get_model_lr_info(model_alt)
    
    lrdf = (df_null - df_alt)
    lrstat = -2*(llf_null - llf_alt)
    lr_pvalue = stats.chi2.sf(lrstat, lrdf)
    
    return lrstat, lr_pvalue, lrdf
    
    
def fit_limixmodel(y, X, sigma, n_iterations=N_ITERATIONS):
    event, gene = X.columns[:2]

    summaries = []
    for i in range(n_iterations):
        # split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE)
        
        # fit linear models with population structure
        if sigma is not None:
            QS = economic_qs(sigma.loc[y_train.index, y_train.index])
        else:
            QS = None
        model = LMM(y_train, X_train, QS)
        model.fit(verbose=False)       
        
        # get rsquared
        pred_train = np.dot(X_train, model.beta)
        rsquared = metrics.r2_score(y_train, pred_train)
        
        # likelihood ratio test
        model_null = LMM(y_train, X_train[[gene,"intercept"]], QS)
        model_null.fit(verbose=False)
        lr_stat, lr_pvalue, lr_df = compare_lr_test(model_null, model)
        
        # score using test data
        pred_test = np.dot(X_test, model.beta)
        pearson_coef, pearson_pvalue = stats.pearsonr(pred_test, y_test)
        spearman_coef, spearman_pvalue = stats.spearmanr(pred_test, y_test)
        
        # prepare output
        ## get info
        params = pd.Series(model.beta, index=X_train.columns)
        scanner = model.get_fast_scanner()
        bse = pd.Series(scanner.null_beta_se, index=X_train.columns)
        zscores = pd.Series(params / bse, index=X_train.columns)
        pvalues = pd.Series(stats.norm.sf(np.abs(zscores)) * 2, index=X_train.columns)
        
        ## make summary
        summary_it = {
            "iteration": i,
            "event_coefficient": params[event],
            "event_stderr": bse[event],
            "event_zscore": zscores[event],
            "event_pvalue": pvalues[event],
            "gene_coefficient": params[gene],
            "gene_stderr": bse[gene],
            "gene_zscore": zscores[gene],
            "gene_pvalue": pvalues[gene],
            "interaction_coefficient": params["interaction"],
            "interaction_stderr": bse["interaction"],
            "interaction_zscore": zscores["interaction"],
            "interaction_pvalue": pvalues["interaction"],
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
        summaries.append(summary_it)
        
    summaries = pd.DataFrame(summaries)

    # compute average likelihood-ratio test
    avg_lr_stat = np.mean(summaries["lr_stat"])
    avg_lr_df = np.round(summaries["lr_df"].mean())
    lr_pvalue = stats.chi2.sf(avg_lr_stat, avg_lr_df)

    # prepare output
    ## summary
    summary = {"EVENT": event, "ENSEMBL": gene, "GENE": y.name, "n_obs": model.nsamples}
    summary.update(get_summary_stats(summaries, "event_coefficient"))
    summary.update(get_summary_stats(summaries, "gene_coefficient"))
    summary.update(get_summary_stats(summaries, "interaction_coefficient"))
    summary.update(get_summary_stats(summaries, "intercept_coefficient"))
    summary.update(get_summary_stats(summaries, "rsquared"))
    summary.update(get_summary_stats(summaries, "pearson_correlation"))
    summary.update(get_summary_stats(summaries, "spearman_correlation"))
    summary.update(get_summary_stats(summaries, "lr_stat"))
    summary.update(
        {"lr_df": lr_df, "lr_pvalue": lr_pvalue,}
    )
    summary = pd.Series(summary)
    ## empirical distributions of coefficients
    coefs = {
        "EVENT": summary["EVENT"],
        "GENE": summary["GENE"],
        "ENSEMBL": summary["ENSEMBL"],
        "event": summaries["event_coefficient"].values,
        "gene": summaries["gene_coefficient"].values,
        "interaction": summaries["interaction_coefficient"].values,
        "intercept": summaries["intercept_coefficient"].values,
    }

    return summary, coefs


def fit_olsmodel(y, X, sigma=None, n_iterations=N_ITERATIONS):
    # np.random.seed(RANDOM_SEED) Doesn't work
    
    event, gene = X.columns[:2]

    summaries = []
    for i in range(n_iterations):
        # split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=TEST_SIZE)

        # fit linear model to training data
        model = sm.OLS(y_train, X_train).fit()

        # log-likelihood test
        model_null = sm.OLS(y_train, X_train[[gene, "intercept"]]).fit()
        lr_stat, lr_pvalue, lr_df = model.compare_lr_test(model_null)

        # score using test data
        prediction = model.predict(X_test)
        pearson_coef, pearson_pvalue = stats.pearsonr(prediction, y_test)
        spearman_coef, spearman_pvalue = stats.spearmanr(prediction, y_test)

        # prepare output
        summary_it = {
            "iteration": i,
            "event_coefficient": model.params[event],
            "event_stderr": model.bse[event],
            "event_zscore": model.params[event] / model.bse[event],
            "event_pvalue": model.pvalues[event],
            "gene_coefficient": model.params[gene],
            "gene_stderr": model.bse[gene],
            "gene_zscore": model.params[gene] / model.bse[gene],
            "gene_pvalue": model.pvalues[gene],
            "interaction_coefficient": model.params["interaction"],
            "interaction_stderr": model.bse["interaction"],
            "interaction_zscore": model.params["interaction"]
            / model.bse["interaction"],
            "interaction_pvalue": model.pvalues["interaction"],
            "intercept_coefficient": model.params["intercept"],
            "intercept_stderr": model.bse["intercept"],
            "intercept_zscore": model.params["intercept"] / model.bse["intercept"],
            "intercept_pvalue": model.pvalues["intercept"],
            "n_obs": model.nobs,
            "rsquared": model.rsquared,
            "pearson_correlation": pearson_coef,
            "pearson_pvalue": pearson_pvalue,
            "spearman_correlation": spearman_coef,
            "spearman_pvalue": spearman_pvalue,
            "lr_stat": lr_stat,
            "lr_pvalue": lr_pvalue,
            "lr_df": lr_df,
        }
        summaries.append(summary_it)
        
    summaries = pd.DataFrame(summaries)

    # compute average likelihood-ratio test
    avg_lr_stat = np.mean(summaries["lr_stat"])
    avg_lr_df = np.round(summaries["lr_df"].mean())
    lr_pvalue = stats.chi2.sf(avg_lr_stat, avg_lr_df)

    # prepare output
    ## summary
    summary = {"EVENT": event, "ENSEMBL": gene, "GENE": y.name, "n_obs": model.nobs}
    summary.update(get_summary_stats(summaries, "event_coefficient"))
    summary.update(get_summary_stats(summaries, "gene_coefficient"))
    summary.update(get_summary_stats(summaries, "interaction_coefficient"))
    summary.update(get_summary_stats(summaries, "intercept_coefficient"))
    summary.update(get_summary_stats(summaries, "rsquared"))
    summary.update(get_summary_stats(summaries, "pearson_correlation"))
    summary.update(get_summary_stats(summaries, "spearman_correlation"))
    summary.update(get_summary_stats(summaries, "lr_stat"))
    summary.update(
        {"lr_df": lr_df, "lr_pvalue": lr_pvalue,}
    )
    summary = pd.Series(summary)
    ## empirical distributions of coefficients
    coefs = {
        "EVENT": summary["EVENT"],
        "GENE": summary["GENE"],
        "ENSEMBL": summary["ENSEMBL"],
        "event": summaries["event_coefficient"].values,
        "gene": summaries["gene_coefficient"].values,
        "interaction": summaries["interaction_coefficient"].values,
        "intercept": summaries["intercept_coefficient"].values,
    }

    return summary, coefs


def fit_single_model(y, X, sigma=None, n_iterations=N_ITERATIONS, method=METHOD):
    methods = {
        "OLS": fit_olsmodel,
        "limix": fit_limixmodel
    }
    summary, coefs = methods[method](y, X, sigma, n_iterations)
    return summary, coefs


def fit_model(x_splicing, x_genexpr, y_gene_dependency, sigma=None, n_iterations=N_ITERATIONS, method=METHOD):

    X = pd.DataFrame([x_splicing, x_genexpr]).T
    y = y_gene_dependency

    # dropna
    is_nan = X.isnull().any(1) | y.isnull()
    X = X.loc[~is_nan].copy()
    y = y[~is_nan].copy()

    try:
        # add interaction and intercept
        X["interaction"] = X[x_splicing.name] * X[x_genexpr.name]
        X["intercept"] = 1.0

        summary, coefs = fit_single_model(y, X, sigma, n_iterations, method)

    except:
        X["interaction"] = np.nan
        X["intercept"] = np.nan

        # create empy summary
        summary = pd.Series(
            np.nan,
            index=[
                "EVENT",
                "ENSEMBL",
                "GENE",
                "n_obs",
                "event_coefficient_mean",
                "event_coefficient_median",
                "event_coefficient_std",
                "event_coefficient_q25",
                "event_coefficient_q75",
                "gene_coefficient_mean",
                "gene_coefficient_median",
                "gene_coefficient_std",
                "gene_coefficient_q25",
                "gene_coefficient_q75",
                "interaction_coefficient_mean",
                "interaction_coefficient_median",
                "interaction_coefficient_std",
                "interaction_coefficient_q25",
                "interaction_coefficient_q75",
                "intercept_coefficient_mean",
                "intercept_coefficient_median",
                "intercept_coefficient_std",
                "intercept_coefficient_q25",
                "intercept_coefficient_q75",
                "rsquared_mean",
                "rsquared_median",
                "rsquared_std",
                "rsquared_q25",
                "rsquared_q75",
                "pearson_correlation_mean",
                "pearson_correlation_median",
                "pearson_correlation_std",
                "pearson_correlation_q25",
                "pearson_correlation_q75",
                "spearman_correlation_mean",
                "spearman_correlation_median",
                "spearman_correlation_std",
                "spearman_correlation_q25",
                "spearman_correlation_q75",
                "lr_stat_mean",
                "lr_stat_median",
                "lr_stat_std",
                "lr_stat_q25",
                "lr_stat_q75",
                "lr_df",
                "lr_pvalue",
            ],
        )
        summary["EVENT"] = x_splicing.name
        summary["ENSEMBL"] = x_genexpr.name
        summary["GENE"] = y_gene_dependency.name

        # create empty empirical coefficients
        coefs = {
            "EVENT": summary["EVENT"],
            "GENE": summary["GENE"],
            "ENSEMBL": summary["ENSEMBL"],
            "event": np.full(n_iterations, np.nan),
            "gene": np.full(n_iterations, np.nan),
            "interaction": np.full(n_iterations, np.nan),
            "intercept": np.full(n_iterations, np.nan),
        }

    # add some more info to the summary
    summary["event_mean"] = x_splicing.mean()
    summary["event_std"] = x_splicing.std()
    summary["gene_mean"] = x_genexpr.mean()
    summary["gene_std"] = x_genexpr.std()

    return summary, coefs


def get_coefs(res, coef_oi, size):
    index = ["EVENT", "GENE", "ENSEMBL"] + list(range(size))
    coefs = pd.Series(
        [res["EVENT"], res["GENE"], res["ENSEMBL"]] + list(res[coef_oi]), index=index,
    )
    return coefs


def fit_models(splicing, genexpr, gene_dependency, mapping, n_iterations=N_ITERATIONS, method=METHOD, n_jobs=None):
    # compute covariance matrix using splicing and TPMs
    print('Computing covariance matrix first...')
    # sigma = pd.concat([splicing,genexpr]).cov()
    sigma = None
    
    # fit models
    results = Parallel(n_jobs=n_jobs)(
        delayed(fit_model)(
            splicing.loc[event],
            genexpr.loc[ensembl],
            gene_dependency.loc[gene],
            sigma = sigma,
            n_iterations = n_iterations,
            method = method
        )
        for event, ensembl, gene in tqdm(mapping.values[:100])
    )
    
    # split results
    summaries = []
    coefs_event = []
    coefs_gene = []
    coefs_interaction = []
    coefs_intercept = []
    for summary, coefs in results:
        summaries.append(summary)
        coefs_event.append(get_coefs(coefs, "event", n_iterations))
        coefs_gene.append(get_coefs(coefs, "gene", n_iterations))
        coefs_interaction.append(get_coefs(coefs, "interaction", n_iterations))
        coefs_intercept.append(get_coefs(coefs, "intercept", n_iterations))

    summaries = pd.DataFrame(summaries)
    coefs_event = pd.DataFrame(coefs_event)
    coefs_gene = pd.DataFrame(coefs_gene)
    coefs_interaction = pd.DataFrame(coefs_interaction)
    coefs_intercept = pd.DataFrame(coefs_intercept)

    # add FDR correction to model summaries
    summaries["lr_padj"] = np.nan
    idx = ~summaries["lr_pvalue"].isnull()
    summaries.loc[idx, "lr_padj"] = sm.stats.multipletests(
        summaries.loc[idx, "lr_pvalue"], method="fdr_bh"
    )[1]
    
    # add more info
    summaries["n_iterations"] = n_iterations
    summaries["method"] = method
    
    return summaries, coefs_event, coefs_gene, coefs_interaction, coefs_intercept


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--splicing_file", type=str)
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--mapping_file", type=str)
    parser.add_argument("--gene_dependency_file", type=str)
    parser.add_argument("--n_iterations", type=int)
    parser.add_argument("--n_jobs", type=int)
    parser.add_argument("--output_dir", type=str)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    splicing_file = args.splicing_file
    genexpr_file = args.genexpr_file
    mapping_file = args.mapping_file
    gene_dependency_file = args.gene_dependency_file
    n_iterations = args.n_iterations
    n_jobs = args.n_jobs
    output_dir = args.output_dir
    
    os.mkdir(output_dir)
    
    print("Loading data...")
    splicing, genexpr, gene_dependency, mapping = load_data(
        splicing_file, genexpr_file, gene_dependency_file, mapping_file
    )

    print("Fitting models...")
    summaries, coefs_event, coefs_gene, coefs_interaction, coefs_intercept = fit_models(
        splicing, genexpr, gene_dependency, mapping, n_iterations, n_jobs
    )

    print("Saving results...")
    summaries.to_csv(os.path.join(output_dir,'model_summaries.tsv.gz'), **SAVE_PARAMS)
    coefs_event.to_pickle(os.path.join(output_dir, "coefs_event.pickle.gz"))
    coefs_gene.to_pickle(os.path.join(output_dir, "coefs_gene.pickle.gz"))
    coefs_interaction.to_pickle(os.path.join(output_dir, "coefs_interaction.pickle.gz"))
    coefs_intercept.to_pickle(os.path.join(output_dir, "coefs_intercept.pickle.gz"))
    
##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
