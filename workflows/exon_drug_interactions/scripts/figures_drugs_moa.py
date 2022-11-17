#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Make figures of pathways.

import argparse
import pandas as pd
import numpy as np
import networkx as nx
from joblib import Parallel, delayed
from tqdm import tqdm
import matplotlib.pyplot as plt

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
FDR_THRESH = 0.1
N_JOBS = 1
RANDOM_SEED = 1234

DRUGS_OI = ["AZD4547", "NUTLIN-3A (-)"]
SOURCES = {"AZD4547": ["FGFR2", "FGFR1", "FGFR3"], "NUTLIN-3A (-)": ["MDM2"]}
TARGETS = {
    "AZD4547": [
        "FGFR2",
        "NFATC2",
        "IFI27",
        "GPNMB",
        "BLNK",
        "PARVG",
        "MYB",
        "CBX2",
        "STAT5A",
        "EZH2",
    ],
    "NUTLIN-3A (-)": ["KIF21B", "MDM2", "PARVG", "TJP3", "ANKRD23", "MDM4", "GNA11"],
}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
RESULTS_DIR = os.path.join(ROOT,'results','splicing_dependency_drugs')
ppi_file = os.path.join(PREP_DIR,'ppi','STRINGDB.tsv.gz')
figs_dir = os.path.join(RESULTS_DIR,"figures","ppi_moa")
n_jobs = 10
"""

##### FUNCTIONS #####
def load_data(ppi_file):
    # load
    ppi = pd.read_table(ppi_file)

    return ppi


def make_network(ppi):
    # make network
    ppi = nx.from_pandas_edgelist(ppi)
    avail_genes = list(ppi.nodes())
    print("Total nodes in PPI network: %s" % len(avail_genes))

    return ppi


def single_shortest_path_length(G, source, target, weight=None):
    try:
        l = nx.shortest_path_length(G, source=source, target=target, weight=weight)
    except:
        print("Unreachable from %s to %s." % (source, target))
        l = np.nan

    return {"source": source, "target": target, "shortest_path_length": l}


def single_shortest_path(G, source, target, weight=None):
    try:
        paths = nx.shortest_path(G, source=source, target=target, weight=weight)
        # paths = nx.all_shortest_paths(G, source=source, target=target, weight=weight)
        # paths = list(set(sum(paths,[])))
    except nx.NetworkXNoPath:
        print("Unreachable from %s to %s." % (source, target))
        paths = []

    return paths


def get_subnetworks(ppi, sources, targets, n_jobs=None):
    """
    Shortest path length from drug targets to significant associations or
    to random genes.
    """
    pairs = [(source, target) for target in targets for source in sources]

    path_lengths = Parallel(n_jobs=n_jobs)(
        delayed(single_shortest_path_length)(ppi, s, t) for s, t in tqdm(pairs)
    )
    path_lengths = pd.DataFrame(path_lengths).dropna()
    shortest = (
        path_lengths.groupby("target")["shortest_path_length"].min().reset_index()
    )
    shortest = pd.merge(
        shortest, path_lengths, how="left", on=["target", "shortest_path_length"]
    )
    shortest_pairs = shortest[["source", "target"]]

    subnetworks = {}
    for source in shortest_pairs["source"].unique():
        nodes_oi = Parallel(n_jobs=n_jobs)(
            delayed(single_shortest_path)(ppi, s, t)
            for s, t in tqdm(
                shortest_pairs.loc[shortest_pairs["source"] == source].values
            )
        )
        nodes_oi = list(set(sum(nodes_oi, [])))

        subnetwork = nx.subgraph(ppi, nodes_oi)
        subnetworks[source] = subnetwork

    return subnetworks


def plot_network(
    G, sources, targets, seed=0, k=2.5, width=5, height=5, output_file="tmp.pdf"
):
    # make properties of nodes for plotting
    info_nodes = pd.DataFrame({"node": G.nodes()})
    info_nodes["is_source"] = info_nodes["node"].isin(sources)
    info_nodes["is_target"] = info_nodes["node"].isin(targets)
    info_nodes["is_mediator"] = (~info_nodes["is_source"]) & (~info_nodes["is_target"])
    info_nodes["size"] = 25
    info_nodes.loc[info_nodes["is_source"], "size"] = 50
    info_nodes["node_shape"] = "o"
    info_nodes.loc[info_nodes["is_source"], "node_shape"] = "H"
    info_nodes["color"] = "grey"
    info_nodes.loc[info_nodes["is_target"], "color"] = "#6e8377ff"
    info_nodes.loc[
        info_nodes["is_source"] & info_nodes["is_target"], "color"
    ] = "#6AC2BF"
    info_nodes["font_weight"] = "normal"
    info_nodes.loc[info_nodes["is_source"], "font_weight"] = "bold"

    cm = 1 / 2.54
    fig = plt.figure(figsize=(width * cm, height * cm), frameon=False)  # width x height

    pos = nx.spring_layout(G, seed=seed, k=k)

    # all
    nx.draw_networkx_edges(G, pos, width=0.3, edge_color="gray")  # edge
    nx.draw_networkx_labels(
        G, pos, font_size=6, font_weight="normal", font_family="Arial", verticalalignment="bottom"
    )
    # sources
    idx = info_nodes["is_source"]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=info_nodes.loc[idx, "node"].to_list(),
        node_color=info_nodes.loc[idx, "color"],
        node_size=info_nodes.loc[idx, "size"],
        node_shape=info_nodes.loc[idx, "node_shape"].unique()[0],
        edgecolors="black",
        linewidths=0.25,
        alpha=0.9,
   )
    # targets
    idx = info_nodes["is_target"] & ~info_nodes["is_source"]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=info_nodes.loc[idx, "node"].to_list(),
        node_color=info_nodes.loc[idx, "color"],
        node_size=info_nodes.loc[idx, "size"],
        node_shape=info_nodes.loc[idx, "node_shape"].unique()[0],
        linewidths=0,
        alpha=0.9,
    )
    # mediators
    idx = info_nodes["is_mediator"]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=info_nodes.loc[idx, "node"].to_list(),
        node_color=info_nodes.loc[idx, "color"],
        node_size=info_nodes.loc[idx, "size"],
        node_shape=info_nodes.loc[idx, "node_shape"].unique()[0],
        linewidths=0,
        alpha=0.5,
    )
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_file, pad_inches=0, bbox_inches="tight")


def make_plots(subnets, figs_dir):

    # Nutlin
    drug_oi = "NUTLIN-3A (-)"
    plot_network(
        subnets[drug_oi]["MDM2"],
        SOURCES[drug_oi],
        TARGETS[drug_oi],
        seed=0,
        k=2,
        width=4,
        height=4,
        output_file=os.path.join(figs_dir, "NUTLIN-MDM2.pdf"),
    )

    # AZD4547
    drug_oi = "AZD4547"
#     plot_network(
#         subnets[drug_oi]["FGFR1"],
#         SOURCES[drug_oi],
#         TARGETS[drug_oi],
#         seed=0,
#         k=1,
#         width=4,
#         height=4,
#         output_file=os.path.join(figs_dir, "AZD4547-FGFR1.pdf"),
#     )
    plot_network(
        subnets[drug_oi]["FGFR2"],
        SOURCES[drug_oi],
        TARGETS[drug_oi],
        seed=1,
        k=1.5,
        width=4,
        height=4,
        output_file=os.path.join(figs_dir, "AZD4547-FGFR2.pdf"),
    )
#     plot_network(
#         subnets[drug_oi]["FGFR3"],
#         SOURCES[drug_oi],
#         TARGETS[drug_oi],
#         seed=2,
#         k=2,
#         width=4,
#         height=4,
#         output_file=os.path.join(figs_dir, "AZD4547-FGFR3.pdf"),
#     )


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

    os.makedirs(figs_dir, exist_ok=True)

    # load
    print("Loading data...")
    ppi = load_data(ppi_file)

    # prepare
    ppi = make_network(ppi)

    # compute shortest paths to get subnetwork with relevant nodes
    subnets = {
        drug: get_subnetworks(ppi, SOURCES[drug][:1], TARGETS[drug], n_jobs)
        for drug in DRUGS_OI
    }

    # visualize
    make_plots(subnets, figs_dir)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
