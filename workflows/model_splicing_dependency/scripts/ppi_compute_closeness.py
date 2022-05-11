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
import itertools

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
N_JOBS = 1
RANDOM_SEED = 1234
N_RANDOM_IT = 100

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
RESULTS_DIR = os.path.join(ROOT,'results','model_splicing_dependency')

ppi_file = os.path.join(PREP_DIR,'ppi','STRINGDB.tsv.gz')
seed_nodes_file = os.path.join(RESULTS_DIR,'files','selected_models-EX-genes.txt')
seed_nodes_file = os.path.join(RESULTS_DIR,'files','random_model_selection-EX-1000its','genes','386.txt')
sink_nodes_file = os.path.join(PREP_DIR,'gene_sets','cancer_gene_census.txt')
n_jobs = 10
"""

##### FUNCTIONS #####
def load_data(ppi_file, seed_nodes_file, sink_nodes_file):
    # load
    ppi = pd.read_table(ppi_file)
    seed_nodes = pd.read_table(seed_nodes_file, header=None)[0].to_list()
    sink_nodes = pd.read_table(sink_nodes_file, header=None)[0].to_list()

    return ppi, seed_nodes, sink_nodes


def prepare_data(ppi, seed_nodes, sink_nodes):
    # make network
    ppi = nx.from_pandas_edgelist(ppi)
    avail_genes = list(ppi.nodes())
    print("Total nodes in PPI network: %s" % len(avail_genes))

    # subset gene sets of interest
    seed_nodes_clean = set(seed_nodes).intersection(avail_genes)
    print(
        "%s seed nodes out of %s were found in the network"
        % (len(seed_nodes_clean), len(seed_nodes))
    )
    sink_nodes_clean = set(sink_nodes).intersection(avail_genes)
    print(
        "%s seed nodes out of %s were found in the network"
        % (len(sink_nodes_clean), len(sink_nodes))
    )

    return ppi, seed_nodes_clean, sink_nodes_clean


def single_shortest_path(G, source, target, weight=None):
    try:
        l = nx.shortest_path_length(G, source=source, target=target, weight=weight)
    except nx.NetworkXNoPath:
        print("Unreachable from %s to %s." % (source, target))
        l = np.nan

    return {"seed": source, "sink": target, "shortest_path_length": l}


def compute_shortest_paths(ppi, seed_nodes, sink_nodes, n_jobs=None):
    """
    Shortest path length from drug targets to significant associations or
    to random genes.
    """
    # real
    pairs = [(seed, sink) for sink in sink_nodes for seed in seed_nodes]
    result = Parallel(n_jobs=n_jobs)(
        delayed(single_shortest_path)(ppi, seed, sink, weight="distance")
        for seed, sink in tqdm(pairs)
    )
    result = pd.DataFrame(result)

    # keep only shortest path for each sink
    idx = result.dropna().groupby("sink")["shortest_path_length"].idxmin()
    result = result.loc[idx]

    return result


def compute_shortest_paths_random(
    ppi,
    seed_nodes,
    sink_nodes,
    n_jobs=None,
    n_random_it=N_RANDOM_IT,
    random_seed=RANDOM_SEED,
):
    """
    Shortest path length from drug targets to significant associations or
    to random genes.
    """
    np.random.seed(random_seed)
    random_seed_nodes = {
        "it"
        + str(i): np.random.choice(
            list(ppi.nodes()), size=len(seed_nodes), replace=False
        )
        for i in range(n_random_it)
    }

    results = []
    for it, seed_nodes in tqdm(random_seed_nodes.items()):
        print(it)

        pairs = [(seed, sink) for sink in sink_nodes for seed in seed_nodes]
        result = Parallel(n_jobs=n_jobs)(
            delayed(single_shortest_path)(ppi, seed, sink, weight="distance")
            for seed, sink in tqdm(pairs)
        )
        result = pd.DataFrame(result)

        # keep only shortest path for each sink
        idx = result.groupby("sink")["shortest_path_length"].idxmin()
        result = result.loc[idx]
        result["type"] = "random"
        result["iteration"] = it

        results.append(result)

    results = pd.concat(results)

    return results


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ppi_file", type=str)
    parser.add_argument("--seed_nodes_file", type=str)
    parser.add_argument("--sink_nodes_file", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--n_jobs", type=int, default=N_JOBS)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    ppi_file = args.ppi_file
    seed_nodes_file = args.seed_nodes_file
    sink_nodes_file = args.sink_nodes_file
    n_jobs = args.n_jobs
    output_file = args.output_file

    print(args)

    # load
    print("Loading data...")
    ppi, seed_nodes, sink_nodes = load_data(ppi_file, seed_nodes_file, sink_nodes_file)

    # prepare
    ppi, seed_nodes, sink_nodes = prepare_data(ppi, seed_nodes, sink_nodes)
    
    print("Computing shortest paths...")
    # compute shortest paths
    result = compute_shortest_paths(ppi, seed_nodes, sink_nodes, n_jobs)

    # save
    print("Saving data...")
    result.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
