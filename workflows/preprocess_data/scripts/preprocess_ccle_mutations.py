# Script purpose
# --------------
# - lift mutations coordinates from hg19 to hg38
# - small changes
#
# Notes
# -----
# There were 81 coordinates that could not be lifted.

import argparse
import numpy as np
import pandas as pd
from liftover import get_lifter
import gc
from tqdm import tqdm

"""
Development
-----------
import os 
ROOT = '~/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
input_file = os.path.join(RAW_DIR,'DepMap','achilles_ccle','CCLE_mutations.csv')
"""

##### FUNCTIONS #####
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    # unpack arguments
    args = parse_args()
    input_file = args.input_file
    output_file = args.output_file

    # load
    print("Loading data...")
    mutations = pd.read_csv(input_file, low_memory=False)
    gc.collect()

    # prepreocess
    print("Preprocessing...")
    mutations["Chromosome"] = "chr" + mutations["Chromosome"]

    ## hg38 to hg19
    converter = get_lifter("hg19", "hg38")
    new_starts = []
    new_ends = []
    for idx, (chrom, start, end) in tqdm(mutations[
        ["Chromosome", "Start_position", "End_position"]
    ].iterrows(), total=len(mutations)):

        new_start = converter[chrom][int(start)]
        new_end = converter[chrom][int(end)]

        if len(new_start) == 0:
            new_starts.append(np.nan)
        else:
            new_starts.append(int(new_start[0][1]))

        if len(new_end) == 0:
            new_ends.append(np.nan)
        else:
            new_ends.append(int(new_end[0][1]))

            
    mutations["lifted_start_hg38"] = pd.Series(new_starts).astype("Int64")
    mutations["lifted_end_hg38"] = pd.Series(new_ends).astype("Int64")
    
    # save
    print("Saving...")
    mutations.to_csv(output_file, sep="\t", compression="gzip", index=False)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
