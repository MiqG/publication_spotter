#
# Author: Miquel Anglada Girotto
# Contact: miquel[dot]anglada[at]crg[dot]eu
# Last Update: 2021-02-26
#
# Script purpose
# --------------
# Count number of mutations that mapped to every event from a column of a table.


from sso_targets import config, util
import pandas as pd
import numpy as np
import argparse

"""
Development
-----------
import os
ROOT = '~/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
annotation_file = os.path.join(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
snv_event_file = os.path.join(PREP_DIR,'event_snv','CCLE-EX.tsv.gz')
snv_gene_file = os.path.join(RAW_DIR,'DepMap','achilles_ccle','CCLE_mutations.csv')
id_col = 'DepMap_ID'
gene_col = 'Hugo_Symbol'
event_col = 'EVENT'
effect_col = 'Variant_Classification'
"""


##### FUNCTIONS #####
def read_snv(snv_file):
    if snv_file.endswith(".csv"):
        snv = pd.read_csv(snv_file, low_memory=False)
    else:
        snv = pd.read_table(snv_file, low_memory=False)

    return snv


def load_data(snv_event_file, snv_gene_file, annotation_file):
    snv_event = read_snv(snv_event_file)
    snv_gene = read_snv(snv_gene_file)
    annot = pd.read_table(annotation_file)

    return snv_event, snv_gene, annot


def count_mutations_per_item_and_effect(snv, id_col, item_col, effect_col):
    """
    mutations per event and effect
    if a sample has multiple mutations in 
    the same gene it is counted only once
    """
    counts = (
        snv[[id_col, item_col, effect_col]]
        .drop_duplicates()
        .groupby([item_col, effect_col])
        .size()
        .reset_index()
        .rename(columns={0: "n"})
    )

    return counts


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--snv_event_file", type=str)
    parser.add_argument("--snv_gene_file", type=str)
    parser.add_argument("--annotation_file", type=str)
    parser.add_argument("--id_col", type=str)
    parser.add_argument("--gene_col", type=str)
    parser.add_argument("--event_col", type=str)
    parser.add_argument("--effect_col", type=str)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()
    return args


def main():
    # parse arguments
    args = parse_args()
    snv_event_file = args.snv_event_file
    snv_gene_file = args.snv_gene_file
    annotation_file = args.annotation_file
    output_file = args.output_file
    id_col = args.id_col
    gene_col = args.gene_col
    event_col = args.event_col
    effect_col = args.effect_col

    # load data
    snv_event, snv_gene, annot = load_data(
        snv_event_file, snv_gene_file, annotation_file
    )

    # count mapped mutations
    mut_freq_event = count_mutations_per_item_and_effect(
        snv_event, id_col, event_col, effect_col
    )
    mut_freq_gene = count_mutations_per_item_and_effect(
        snv_gene, id_col, gene_col, effect_col
    )

    # for each variant, how many mutations in the gene mapped on the exon?
    mut_freq_event = pd.merge(
        mut_freq_event, annot[["EVENT", "GENE"]], on="EVENT", how="left"
    )
    mut_freq_gene = mut_freq_gene.rename(columns={gene_col: "GENE", "n": "total"})

    mut_freq_event = pd.merge(
        mut_freq_event, mut_freq_gene, on=[effect_col, "GENE"], how="left"
    )
    mut_freq_event["ratio"] = mut_freq_event["n"] / mut_freq_event["total"]

    # save
    mut_freq_event.to_csv(output_file, sep="\t", index=False, compression="gzip")


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
