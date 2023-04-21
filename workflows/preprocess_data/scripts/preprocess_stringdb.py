#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Preprocess and clean CCLE data based on Goncalves 2020 (https://github.com/EmanuelGoncalves/dtrace/blob/a81e379e6ca511cc013aa25b0c7309f8fd0a5f16/dtrace/DataImporter.py#L481)

import argparse
import pandas as pd
import networkx as nx

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
SCORE_THRESH = 900

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
raw_ppi_file = os.path.join(RAW_DIR,'STRINGDB','9606.protein.links.full.v11.5.txt.gz')
raw_aliases_file = os.path.join(RAW_DIR,'STRINGDB','9606.protein.aliases.v11.5.txt.gz')
"""

##### FUNCTIONS #####
def load_data(raw_ppi_file, raw_aliases_file):
    ppi = pd.read_table(raw_ppi_file, sep=" ")
    aliases = pd.read_table(raw_aliases_file)

    return ppi, aliases


def preprocess_stringdb(ppi, aliases):
    # map ENSP to gene symbol
    gmap = (
        aliases.loc[aliases["source"] == "BioMart_HUGO"]
        .set_index("#string_protein_id")["alias"]
        .to_dict()
    )

    # Filter by moderate confidence
    ppi = ppi[ppi["combined_score"] > SCORE_THRESH].copy()

    # Filter and map to gene symbol
    avail = ppi["protein1"].isin(gmap.keys()) & ppi["protein2"].isin(gmap.keys())
    ppi = ppi.loc[avail].copy()

    ppi["protein1"] = [gmap[p1] for p1 in ppi["protein1"]]
    ppi["protein2"] = [gmap[p2] for p2 in ppi["protein2"]]

    # if duplicated, keep only edges with maximum weight
    ppi = ppi.groupby(["protein1", "protein2"])["combined_score"].max().reset_index()

    # remove self loops
    is_selfloop = ppi["protein1"] == ppi["protein2"]
    ppi = ppi.loc[~is_selfloop].copy()

    # standardize names
    ppi = ppi.rename(
        columns={
            "protein1": "source",
            "protein2": "target",
            "combined_score": "stringdb_score",
        }
    )

    return ppi


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw_ppi_file", type=str)
    parser.add_argument("--raw_aliases_file", type=str)
    parser.add_argument("--prep_ppi_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    raw_ppi_file = args.raw_ppi_file
    raw_aliases_file = args.raw_aliases_file
    prep_ppi_file = args.prep_ppi_file

    # load
    print("Loading data...")
    ppi, aliases = load_data(raw_ppi_file, raw_aliases_file)
    print("Preprocessing data...")
    ppi = preprocess_stringdb(ppi, aliases)

    # save
    print("Saving data...")
    ppi.to_csv(prep_ppi_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
