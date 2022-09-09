#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Given a VastDB exon, get all isoforms including it or not.

import argparse
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm
import gget
from Bio.Seq import Seq
from difflib import SequenceMatcher

# variables

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
event_oi = "HsaEX0050345"
annotation_file = os.path.join(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
event_info_file = os.path.join(RAW_DIR,"VastDB","EVENT_INFO-hg38.tab.gz")
"""


##### FUNCTIONS #####
def fasta_to_pandas(fasta_list):
    # identifiers and sequences are interleaved
    ids = [fasta_list[i] for i in range(0, len(fasta_list), 2)]
    seqs = [fasta_list[i] for i in range(1, len(fasta_list), 2)]

    # make df
    df = pd.DataFrame({"id": ids, "aa": seqs})

    return df


def longestSubstring(str1, str2):
    # initialize SequenceMatcher object with
    # input string
    seqMatch = SequenceMatcher(None,str1,str2)

    # find match of longest sub-string
    # output will be like Match(a=0, b=0, size=5)
    match = seqMatch.find_longest_match(0, len(str1), 0, len(str2))

    # print longest substring
    if (match.size!=0):
        match_seq = str1[match.a: match.a + match.size]
    else:
        print ('No longest common sub-string found')
        match_seq = None
        
    return match_seq
            
            
def get_proteoforms(event_oi, annot, event_info):
    # get ensembl id
    gene_oi = annot.loc[annot["EVENT"] == event_oi, "ENSEMBL"].to_list()

    # get aa sequences of the gene
    canonical = gget.seq(gene_oi, translate=True, isoforms=False)
    isoforms = gget.seq(gene_oi, translate=True, isoforms=True)

    # transform to a dataframe
    df = fasta_to_pandas(isoforms)

    # drop duplicate sequences as multiple transcripts may
    # lead to the same AA sequence
    df = df.loc[~df["aa"].duplicated()].copy()

    # indicate the canonical sequence
    df["is_canonical"] = df["aa"] == canonical[1]

    # indicate if it contains the event of interest
    idx = event_info["EVENT"] == event_oi
    event_dna = event_info.loc[idx, "Seq_A"].values[0]
    event_strand = event_info.loc[idx, "REF_CO"].str[-1].values[0]
    
    ## the ORF can have different start points
    possible_aa = [str(Seq(event_dna[i:]).transcribe().translate()) for i in range(0,3)]
    match_seqs = []
    for pos_aa in possible_aa:
        match_seq = pd.DataFrame(df["aa"].apply(lambda aa: longestSubstring(aa, pos_aa)))
        match_seq["len"] = match_seq["aa"].str.len()
        match_seq["full_aa"] = pos_aa
        match_seqs.append(match_seq)
    match_seqs = pd.concat(match_seqs).reset_index()
    idx = match_seqs["len"].idxmax()
    event_aa = match_seqs.loc[idx,"aa"] # get longest match with
    event_aa_full = match_seqs.loc[idx,"full_aa"]
    
    ## indicate
    df["event_included"] = df["aa"].str.contains(event_aa, regex=True)
    df["event_aa"] = event_aa
    df["event_aa_full"] = event_aa_full

    # add info to id
    df["id_new"] = [
        "%s | is_canonical=%s | event_included=%s"
        % (row["id"], row["is_canonical"], row["event_included"])
        for idx, row in df.iterrows()
    ]

    # delete event
    df["aa_noevent"] = df["aa"].str.replace(event_aa, "")

    # prepare fasta of the canonical with and without the event sequence
    #     if any(df["is_canonical"] & df["event_included"]):
    #         sel = df.loc[df["is_canonical"]]
    #         fasta = [
    #             sel["id_new"].values[0] + " | event_included=True",
    #             sel["aa"].values[0],
    #             sel["id_new"].values[0] + " | event_included=False",
    #             sel["aa_noevent"].values[0],
    #         ]
        
    if any(df["event_included"]):
        sel = df.loc[df["event_included"]]
        
        fasta = []
        for idx, row in sel.iterrows():
            fasta.append([
                row["id_new"] + " | event_included=True",
                row["aa"],
                row["id_new"] + " | event_included=False",
                row["aa_noevent"],
            ])
        fasta = sum(fasta,[])
        
    else:
        print("The canonical proteoform does not contain the event!")

    return df, fasta


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--event_oi", type=str)
    parser.add_argument("--annotation_file", type=str)
    parser.add_argument("--event_info_file", type=str)
    parser.add_argument("--output_df_file", type=str)
    parser.add_argument("--output_fasta_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    event_oi = args.event_oi
    annotation_file = args.annotation_file
    event_info_file = args.event_info_file
    output_df_file = args.output_df_file
    output_fasta_file = args.output_fasta_file

    # load
    print("Loading data...")
    annot = pd.read_table(annotation_file)
    event_info = pd.read_table(event_info_file)

    # get the proteoforms of interest
    print("Getting proteoforms...")
    df, fasta = get_proteoforms(event_oi, annot, event_info)

    # save
    print("Saving proteoforms...")
    df.to_csv(output_df_file, sep="\t", compression="gzip", index=False)
    with open(output_fasta_file, "w") as f:
        f.writelines("\n".join(fasta) + "\n")


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
