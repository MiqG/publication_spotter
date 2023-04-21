#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# For each drug, measure the shortest path length from each of its targets
# to each of the genes significantly associated to the drug sensitivity profile

import argparse
import pandas as pd
import numpy as np
import networkx as nx
from joblib import Parallel, delayed
from tqdm import tqdm

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
FDR_THRESH = 0.1
N_JOBS = 1
RANDOM_SEED = 1234

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
RESULTS_DIR = os.path.join(ROOT,'results','splicing_dependency_drugs')
ppi_file = os.path.join(PREP_DIR,'ppi','STRINGDB.tsv.gz')
drug_targets_file = os.path.join(PREP_DIR,'drug_screens','drug_targets.tsv.gz')
drug_associations_file = os.path.join(RESULTS_DIR,'files','model_summaries_drug_response-EX.tsv.gz')
thresh_fdr = 0.1
n_jobs=10
n_random_sources = 100
"""

##### FUNCTIONS #####
def load_data(ppi_file, drug_targets_file, drug_associations_file):
    # load
    ppi = pd.read_table(ppi_file)
    drug_targets = pd.read_table(drug_targets_file)

    if drug_associations_file is not None:
        drug_associations = pd.read_table(drug_associations_file)
    else:
        drug_associations = None

    return ppi, drug_targets, drug_associations


def prepare_data(ppi, drug_targets, drug_associations, thresh_fdr):
    # make network
    ppi = nx.from_pandas_edgelist(ppi)
    avail_genes = list(ppi.nodes())
    print("Total nodes in PPI network: %s" % len(avail_genes))

    # drug - target pairs
    ## drop genes not in the network
    drug_targets = drug_targets.loc[drug_targets["TARGET"].isin(avail_genes)]

    # significant drug - event interactions at the gene level
    if drug_associations is not None:
        drug_associations = drug_associations.loc[
            drug_associations["lr_padj"] < thresh_fdr
        ].copy()
        drug_associations = drug_associations[["DRUG_ID", "GENE"]].drop_duplicates()
        ## drop genes not in the network
        drug_associations = drug_associations.loc[
            drug_associations["GENE"].isin(avail_genes)
        ]

    return ppi, drug_targets, drug_associations


def single_shortest_path(G, source, target, weight=None):
    try:
        l = nx.shortest_path_length(G, source=source, target=target, weight=weight)
    except nx.NetworkXNoPath:
        print("Unreachable from %s to %s." % (source, target))
        l = np.nan

    return {"source": source, "target": target, "shortest_path_length": l}


def compute_shortest_paths_real(ppi, drug_targets, drug_associations, n_jobs=None):
    """
    Shortest path length from drug targets to significant associations or
    to random genes.
    """
    results = []
    for drug_oi in drug_targets["DRUG_ID"].unique():
        print(drug_oi)
        targets = drug_targets.loc[drug_targets["DRUG_ID"] == drug_oi, "TARGET"].values
        sources = drug_associations.loc[
            drug_associations["DRUG_ID"] == drug_oi, "GENE"
        ].values
        pairs = [(source, target) for target in targets for source in sources]

        result = Parallel(n_jobs=n_jobs)(
            delayed(single_shortest_path)(ppi, s, t, weight="distance")
            for s, t in tqdm(pairs)
        )
        result = pd.DataFrame(result)
        result["type"] = "real"
        result["DRUG_ID"] = drug_oi
        results.append(result)

    results = pd.concat(results)

    return results


def compute_shortest_paths_random(
    ppi, drug_targets, n_random_sources, n_jobs=None, random_seed=RANDOM_SEED
):
    """
    Shortest path length from drug targets to random genes.
    """
    np.random.seed(random_seed)

    print("Random paths...")
    targets = drug_targets["TARGET"].unique()
    pairs = []
    for target in targets:
        sources = np.random.choice(ppi.nodes(), size=n_random_sources, replace=False)
        pairs.append([(source, target) for source in sources])
    pairs = sum(pairs, [])

    result = Parallel(n_jobs=n_jobs)(
        delayed(single_shortest_path)(ppi, s, t, weight="distance")
        for s, t in tqdm(pairs)
    )
    result = pd.DataFrame(result)
    result["type"] = "random"
    result["n_random_sources"] = n_random_sources

    return result


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ppi_file", type=str)
    parser.add_argument("--drug_targets_file", type=str)
    parser.add_argument("--drug_associations_file", type=str, default=None)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--thresh_fdr", type=float, default=FDR_THRESH)
    parser.add_argument("--n_jobs", type=int, default=N_JOBS)
    parser.add_argument("--n_random_sources", type=int)
    parser.add_argument("--mode", type=str, default=None)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    ppi_file = args.ppi_file
    drug_targets_file = args.drug_targets_file
    drug_associations_file = args.drug_associations_file
    thresh_fdr = args.thresh_fdr
    n_jobs = args.n_jobs
    n_random_sources = args.n_random_sources
    mode = args.mode
    output_file = args.output_file

    print(args)

    assert mode is not None

    # load
    print("Loading data...")
    ppi, drug_targets, drug_associations = load_data(
        ppi_file, drug_targets_file, drug_associations_file
    )

    # prepare
    ppi, drug_targets, drug_associations = prepare_data(
        ppi, drug_targets, drug_associations, thresh_fdr
    )

    # compute shortest paths
    if mode == "real":
        result = compute_shortest_paths_real(
            ppi, drug_targets, drug_associations, n_jobs
        )
    elif mode == "random":
        result = compute_shortest_paths_random(
            ppi, drug_targets, n_random_sources, n_jobs
        )

    # save
    print("Saving data...")
    result.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
