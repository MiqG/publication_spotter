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

import statsmodels.regression.linear_model as lm
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
from scipy import stats
import argparse
import gc
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

"""
Development
-----------
import os
ROOT = '~/projects/publication_splicing_dependency'
psi_file = os.path.join(ROOT,'data','prep','exon_psi','CCLE.tsv.gz')
genexpr_file = os.path.join(ROOT,'data','raw','DepMap','achilles_ccle','CCLE_expression_transposed.tsv.gz')
annotation_file = os.path.join(ROOT,'data','raw','VastDB','EVENT_INFO-hg38_noseqs.tsv')
rnai_file =  os.path.join(ROOT,'data','raw','DepMap','demeter2','D2_combined_gene_dep_scores.csv')
metadata_file = os.path.join(ROOT,'data','prep','metadata','CCLE.tsv')
"""

##### FUNCTIONS #####
def load_data(psi_file, genexpr_file, rnai_file, metadata_file, annotation_file):
    psi = pd.read_table(psi_file, index_col=0)
    genexpr = pd.read_table(genexpr_file, index_col=0)
    annotation = pd.read_table(annotation_file)
    rnai = pd.read_csv(rnai_file, index_col=0)
    metadata = pd.read_table(metadata_file)
    
    # drop undetected & uninformative events
    psi = psi.dropna(thresh=2)
    psi = psi.loc[psi.std(axis=1)!=0]
    
    # strip gene names
    genexpr.index = [symbol for symbol, entrez in genexpr.index.str.split(" ")]
    rnai.index = [symbol for symbol, entrez in rnai.index.str.split(" ")]

    # rename samples
    rnai = rnai.rename(columns=metadata.set_index("CCLE_Name")["DepMap_ID"].to_dict())

    # subset
    common_samples = (
        set(rnai.columns).intersection(psi.columns).intersection(genexpr.columns)
    )
    common_genes = set(rnai.index).intersection(genexpr.index)
    common_events = set(psi.index).intersection(
        annotation.loc[annotation["GENE"].isin(common_genes), "EVENT"]
    )

    psi = psi.loc[common_events, common_samples]
    genexpr = genexpr.loc[common_genes, common_samples]
    rnai = rnai.loc[common_genes, common_samples]
    annotation = annotation.loc[annotation["EVENT"].isin(common_events)]

    gc.collect()

    return psi, genexpr, rnai, annotation


def fit_model(x_psi, x_genexpr, y_rnai):

    X = pd.DataFrame([x_psi, x_genexpr]).T
    y = y_rnai
    
    # dropna
    is_nan = X.isnull().any(1) | y.isnull()
    X = X.loc[~is_nan]
    y = y[~is_nan]
    
    try:
        # standardize features
        X.values[:, :] = StandardScaler().fit_transform(X)
        X['interaction'] = X[x_psi.name]*X[x_genexpr.name]
        X['intercept'] = 1.0
        
        # fit linear model
        model = lm.OLS(y, X).fit()
        
        # score
        pearson_coef, pearson_pvalue = stats.pearsonr(model.predict(X),y)
        spearman_coef, spearman_pvalue = stats.spearmanr(model.predict(X),y)
        
        # prepare output
        summary = pd.DataFrame({
            "coefficient": model.params, 
            "stderr": model.bse,
            "zscore": model.params / model.bse,
            "pvalue": model.pvalues,
            "n_obs": model.nobs,
            "rsquared": model.rsquared,
            "pearson_coef": pearson_coef,
            "pearson_pvalue": pearson_pvalue,
        })
        
    except:
        X['interaction'] = np.nan
        X['intercept'] = np.nan
        
        pearson_coef, pearson_pvalue = np.nan, np.nan
        spearman_coef, spearman_pvalue = np.nan, np.nan
        
        summary = pd.DataFrame({
            "coefficient": [np.nan]*X.shape[1], 
            "stderr": np.nan,
            "zscore": np.nan,
            "pvalue": np.nan,
            "n_obs": np.nan,
            "rsquared": np.nan,
            "pearson_coef": np.nan,
            "pearson_pvalue": np.nan,
            },
            index=X.columns
        )
        
    summary = pd.Series({
        'EVENT': x_psi.name,
        'GENE': x_genexpr.name,
        'event_coefficient': summary.loc[x_psi.name,'coefficient'],
        'event_stderr': summary.loc[x_psi.name,'stderr'],
        'event_zscore': summary.loc[x_psi.name,'zscore'],
        'event_pvalue': summary.loc[x_psi.name,'pvalue'],
        'gene_coefficient': summary.loc[x_genexpr.name,'coefficient'],
        'gene_stderr': summary.loc[x_genexpr.name,'stderr'],
        'gene_zscore': summary.loc[x_genexpr.name,'zscore'],
        'gene_pvalue': summary.loc[x_genexpr.name,'pvalue'],
        'interaction_coefficient': summary.loc['interaction','coefficient'],
        'interaction_stderr': summary.loc['interaction','stderr'],
        'interaction_zscore': summary.loc['interaction','zscore'],
        'interaction_pvalue': summary.loc['interaction','pvalue'],
        'intercept_coefficient': summary.loc['intercept','coefficient'],
        'intercept_stderr': summary.loc['intercept','stderr'],
        'intercept_zscore': summary.loc['intercept','zscore'],
        'intercept_pvalue': summary.loc['intercept','pvalue'],
        'n_obs': summary.loc[x_genexpr.name,'n_obs'],
        'rsquared': summary.loc[x_genexpr.name,'rsquared'],
        'pearson_coefficient': pearson_coef,
        'pearson_pvalue': pearson_pvalue,
        'spearman_coefficient': spearman_coef,
        'spearman_pvalue': spearman_pvalue,
        'event_mean': x_psi.mean(),
        'event_std': x_psi.std(),
        'gene_mean': x_genexpr.mean(),
        'gene_std': x_genexpr.std()
    })
    
    return summary


def fit_models(psi, genexpr, rnai, annotation, n_jobs):
    # compute sample covariance
    results = Parallel(n_jobs=n_jobs)(
        delayed(fit_model)(psi.loc[event], genexpr.loc[gene], rnai.loc[gene])
        for gene, event in tqdm(annotation[["GENE", "EVENT"]].values)
    )
    results = pd.DataFrame(results)
    return results


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--psi_file", type=str)
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--annotation_file", type=str)
    parser.add_argument("--rnai_file", type=str)
    parser.add_argument("--metadata_file", type=str)
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
    metadata_file = args.metadata_file
    n_jobs = args.n_jobs
    output_file = args.output_file
    
    print('Loading data...')
    psi, genexpr, rnai, annotation = load_data(
        psi_file, genexpr_file, rnai_file, metadata_file, annotation_file
    )
    
    print('Fitting models...')
    result = fit_models(psi, genexpr, rnai, annotation, n_jobs)
    
    print('Saving results...')
    result.to_csv(output_file, sep='\t', compression='gzip', index=False)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
