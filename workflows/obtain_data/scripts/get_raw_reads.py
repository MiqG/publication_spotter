#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#

import argparse
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
N_JOBS = 1

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_spotter'
RAW_DIR = os.path.join(ROOT,'data','raw')
inclusion_levels_file = os.path.join(RAW_DIR,"TCGA","GBM","vast_out","INCLUSION_LEVELS_FULL-hg38-175.tab.gz")
inclusion_levels_file = os.path.join(RAW_DIR,"articles",'Thomas2020','vast_out','INCLUSION_LEVELS_FULL-hg38-2.tab.gz')
n_jobs = 10
chunksize = 1000
"""

##### FUNCTIONS #####
def get_total_rows(filename):
    return len(pd.read_table(filename, usecols=[0]))
    
def get_inclusion_reads(x):
    # score3 contains INC1=INC2=EXC or INC=EXC (reads mapping inclusion and exclusion junctions)
    ## some scores are not in the correct order
    scores = x.str.replace(",=",",").str.split(",", expand=True)
    score3 = scores[2]
    is_alt_or_int = scores.index.str.contains("ALT|INT")
    if any(is_alt_or_int):
        score3[is_alt_or_int] = scores.loc[is_alt_or_int,3]
        
    n_ss = score3.str.count("=")
    idx_two_ss = n_ss==2
    idx_one_ss = n_ss==1
    
    inclusion_reads = score3.copy()
    if any(idx_two_ss):
        inclusion_reads.loc[idx_two_ss] = inclusion_reads.loc[idx_two_ss].str.split("=", expand=True)[[0,1]].replace("",None).astype("int64").sum(axis=1)
    if any(idx_one_ss):        
        inclusion_reads.loc[idx_one_ss] = inclusion_reads.loc[idx_one_ss].str.split("=", expand=True)[0].replace("",None).astype("int64")
    
    return inclusion_reads


def get_exclusion_reads(x):
    # score3 contains INC1=INC2=EXC or INC=EXC (reads mapping inclusion and exclusion junctions)
    ## some scores are not in the correct order
    scores = x.str.replace(",=",",").str.split(",", expand=True)
    score3 = scores[2]
    is_alt_or_int = scores.index.str.contains("ALT|INT")
    if any(is_alt_or_int):
        score3[is_alt_or_int] = scores.loc[is_alt_or_int,3]
        
    n_ss = score3.str.count("=")
    idx_two_ss = n_ss==2
    idx_one_ss = n_ss==1
    
    exclusion_reads = score3.copy()
    if any(idx_two_ss):
        exclusion_reads.loc[idx_two_ss] = exclusion_reads.loc[idx_two_ss].str.split("=", expand=True)[2].replace("",None).astype("int64")
    if any(idx_one_ss):        
        exclusion_reads.loc[idx_one_ss] = exclusion_reads.loc[idx_one_ss].str.split("=", expand=True)[1].replace("",None).astype("int64")
    
    return exclusion_reads

    
def get_raw_reads(inclusion_levels):

    # subset event and quality columns
    cols_oi = ["EVENT"] + [c for c in inclusion_levels.columns if c.endswith("-Q")]
    inclusion_levels = inclusion_levels[cols_oi].set_index("EVENT")

    # get read counts
    inclusion_reads = inclusion_levels.apply(get_inclusion_reads, axis=0)
    exclusion_reads = inclusion_levels.apply(get_exclusion_reads, axis=0)
    total_reads = inclusion_reads + exclusion_reads
    
    return inclusion_reads, exclusion_reads, total_reads

    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inclusion_levels_file", type=str)
    parser.add_argument("--n_jobs", type=int, default=N_JOBS)
    parser.add_argument("--chunksize", type=int, default=1000)
    parser.add_argument("--inclusion_reads_file", type=str)
    parser.add_argument("--exclusion_reads_file", type=str)
    parser.add_argument("--total_reads_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    inclusion_levels_file = args.inclusion_levels_file
    n_jobs = args.n_jobs
    chunksize = args.chunksize
    inclusion_reads_file = args.inclusion_reads_file
    exclusion_reads_file = args.exclusion_reads_file
    total_reads_file = args.total_reads_file

    print(args)

    # load
    print("Loading data...")
    inclusion_levels_chunks = pd.read_table(inclusion_levels_file, chunksize=chunksize)
    
    # preprocess
    print("Preprocessing...")
    results = Parallel(n_jobs=n_jobs)(
        delayed(get_raw_reads)(
            inclusion_levels_chunk
        )
        for inclusion_levels_chunk in tqdm(inclusion_levels_chunks)
    )
    
    # prepare outputs
    inclusion_reads = []
    exclusion_reads = []
    total_reads = []
    for result in results:
        inclusion_reads_chunk, exclusion_reads_chunk, total_reads_chunk = result
        
        inclusion_reads.append(inclusion_reads_chunk)
        exclusion_reads.append(exclusion_reads_chunk)
        total_reads.append(total_reads_chunk)
        
    inclusion_reads = pd.concat(inclusion_reads)
    exclusion_reads = pd.concat(exclusion_reads)
    total_reads = pd.concat(total_reads)
    
    # save
    print("Saving data...")
    inclusion_reads.reset_index().to_csv(inclusion_reads_file, **SAVE_PARAMS)
    exclusion_reads.reset_index().to_csv(exclusion_reads_file, **SAVE_PARAMS)
    total_reads.reset_index().to_csv(total_reads_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
