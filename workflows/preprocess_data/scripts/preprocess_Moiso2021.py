#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Preprocess and clean CCLE data.

import argparse
import pandas as pd

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
metadata_file = os.path.join(PREP_DIR,'metadata','PANCAN.tsv.gz')
response_file = os.path.join(RAW_DIR,'articles','Moiso2021','drug_response_filtered.tsv')
"""

##### FUNCTIONS #####
def load_data(metadata_file, response_file):
    metadata = pd.read_table(metadata_file, low_memory=False)
    response = pd.read_table(response_file)

    return metadata, response


def preprocess_response(metadata, response):
    metadata["patient_ID"] = metadata["sampleID"].str[0:12]
    response = pd.merge(
        response,
        metadata[["patient_ID", "sampleID", "sample_type"]],
        how="left",
        on="patient_ID",
    )

    return response


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata_file", type=str)
    parser.add_argument("--response_file", type=str)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    metadata_file = args.metadata_file
    response_file = args.response_file
    output_file = args.output_file

    # load
    print("Loading data...")
    metadata, response = load_data(metadata_file, response_file)
    print("Preprocessing data...")
    response = preprocess_response(metadata, response)

    # save
    print("Saving data...")
    response.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
