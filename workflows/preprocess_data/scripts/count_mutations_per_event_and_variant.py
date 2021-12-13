#
# Author: Miquel Anglada Girotto
# Contact: miquel[dot]anglada[at]crg[dot]eu
# Last Update: 2021-02-26
#
# Script purpose
# --------------
# Count number of mutations that mapped to every event from a column of a table.
# Given its length and the total number of mutations that hit the gene,
# is an exon more/less mutated than expected?


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
annot_events_file = os.path.join(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
annot_genes_file = os.path.join(RAW_DIR,'ENSEMBL','gene_annotation-hg19.tsv.gz')
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


def load_data(snv_event_file, snv_gene_file, annot_events_file, annot_genes_file):
    snv_event = read_snv(snv_event_file)
    snv_gene = read_snv(snv_gene_file)
    annot_events = pd.read_table(annot_events_file).rename(
        columns={"LE_o": "length_event"}
    )
    annot_genes = pd.read_table(annot_genes_file).rename(
        columns={"Gene": "GENE", "length": "length_gene"}
    )

    return snv_event, snv_gene, annot_events, annot_genes


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
    parser.add_argument("--annot_events_file", type=str)
    parser.add_argument("--annot_genes_file", type=str)
    parser.add_argument("--id_col", type=str)
    parser.add_argument("--event_col", type=str)
    parser.add_argument("--gene_col", type=str)
    parser.add_argument("--effect_col", type=str)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()
    return args


def main():
    # parse arguments
    args = parse_args()
    snv_event_file = args.snv_event_file
    snv_gene_file = args.snv_gene_file
    annot_events_file = args.annot_events_file
    annot_genes_file = args.annot_genes_file
    output_file = args.output_file
    id_col = args.id_col
    event_col = args.event_col
    gene_col = args.gene_col
    effect_col = args.effect_col

    # load data
    snv_event, snv_gene, annot_events, annot_genes = load_data(
        snv_event_file, snv_gene_file, annot_events_file, annot_genes_file
    )

    # count mapped mutations
    mut_freq_event = count_mutations_per_item_and_effect(
        snv_event, id_col, event_col, effect_col
    )
    mut_freq_gene = count_mutations_per_item_and_effect(
        snv_gene, id_col, gene_col, effect_col
    )

    # for each variant, what is the mutation frequency per kilobase?
    ## events
    mut_freq_event = pd.merge(
        mut_freq_event,
        annot_events[["EVENT", "GENE", "length_event"]],
        on="EVENT",
        how="left",
    )
    mut_freq_event["event_mut_freq_per_kb"] = (
        mut_freq_event["n"] / mut_freq_event["length_event"]
    ) * 1000
    ## genes
    mut_freq_gene = pd.merge(
        mut_freq_gene.rename(columns={gene_col: "GENE", "n": "total"}),
        annot_genes[["GENE", "length_gene"]],
        on="GENE",
        how="left",
    )
    mut_freq_gene["gene_mut_freq_per_kb"] = (
        mut_freq_gene["total"] / mut_freq_gene["length_gene"]
    ) * 1000

    # for each variant, how many mutations would I expect given the event length?
    mut_freq_event["event_kb"] = mut_freq_event["length_event"] / 1000  # from bp to kb
    mut_freq_event = pd.merge(
        mut_freq_event, mut_freq_gene, how="left", on=["Variant_Classification", "GENE"]
    )
    mut_freq_event["expected_mutations"] = mut_freq_event["event_kb"] * mut_freq_event["gene_mut_freq_per_kb"]
    mut_freq_event["fc_mutations"] = np.log2(mut_freq_event["n"] / mut_freq_event["expected_mutations"])
    
    # save
    mut_freq_event.to_csv(output_file, sep="\t", index=False, compression="gzip")


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
