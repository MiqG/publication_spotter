# Script purpose
# --------------
# Map mutations to events. Results in a mutations file that has new columns starting
# with "EVENT*". Note that mutations could map to multiple events or to none.
# We used a MARGIN of 2 to consider mutations on splice donor or acceptor sites,
# but used the original coordinates for the result.

import argparse
import pandas as pd
import pyranges as pr
import gc


"""
Development
-----------
import os 
ROOT = '~/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
mutations_file = os.path.join(PREP_DIR,'mutations','CCLE.tsv.gz')
annotation_file = os.path.join(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
margin = 100
event_type = 'EX'
"""


##### FUNCTIONS #####
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations_file", type=str)
    parser.add_argument("--annotation_file", type=str)
    parser.add_argument("--event_type", type=str)
    parser.add_argument("--margin", type=int)
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    # unpack arguments
    args = parse_args()
    mutations_file = args.mutations_file
    annotation_file = args.annotation_file
    event_type = args.event_type
    margin = args.margin
    output_file = args.output_file

    # load
    print("Loading data...")
    mutations = pd.read_table(mutations_file, low_memory=False)
    events = pd.read_table(annotation_file)
    gc.collect()

    # prepreocess
    print("Preprocessing...")
    ## only event type of interest
    events = events.loc[events["EVENT"].str.contains(event_type)].copy()
    ## process event coordinates
    events["Chromosome"] = (
        events["COORD_o"].str.split(":").str[0]  #.str.replace("chr", "")
    )
    events["Start"] = events["COORD_o"].str.split(":").str[1].str.split("-").str[0]
    events["End"] = events["COORD_o"].str.split(":").str[1].str.split("-").str[1]
    ## process mutations
    mutations = mutations.rename(
        columns={"lifted_start_hg38": "Start", "lifted_end_hg38": "End"}
    )
    mutations = mutations.loc[~mutations[["Start", "End"]].isnull().any(axis=1)]

    # intersect events with variants
    X = pr.PyRanges(events[["EVENT", "Chromosome", "Start", "End"]], int64=True)
    Y = pr.PyRanges(
        mutations[
            [
                "Hugo_Symbol",
                "Entrez_Gene_Id",
                "Chromosome",
                "Start",
                "End",
                "Reference_Allele",
                "Tumor_Seq_Allele1",
                "Variant_Type",
                "Variant_Classification",
                "DepMap_ID"
            ]
        ],
        int64=True,
    )
    mapping = X.join(Y, slack=margin, how="left", report_overlap=True)
    
    # prepare outputs
    result = mapping.df.rename(columns={"Start_b":"Start_mutation", "End_b":"End_mutation"})
    result = result.loc[result["Reference_Allele"] != "-1"]
    
    # save
    result.to_csv(output_file, sep="\t", compression="gzip", index=False)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
