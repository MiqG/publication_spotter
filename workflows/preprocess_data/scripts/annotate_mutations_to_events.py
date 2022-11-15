# Script purpose
# --------------
# Map mutations to events. Results in a mutations file that has new columns starting
# with "EVENT*". Note that mutations could map to multiple events or to none.
# We used a MARGIN of 2 to consider mutations on splice donor or acceptor sites,
# but used the original coordinates for the result.
# Consider the whole neighboring intron as margin

import argparse
import numpy as np
import pandas as pd
import pyranges as pr
import gc


"""
Development
-----------
import os 
ROOT = '~/projects/publication_spotter'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
mutations_file = os.path.join(PREP_DIR,'mutations','CCLE.tsv.gz')
annotation_file = os.path.join(RAW_DIR,'VastDB','EVENT_INFO-hg38_noseqs.tsv')
margin = 100
event_type = 'EX'
"""


##### FUNCTIONS #####
def get_abs_min(X):
    cols = np.argmin(np.abs(X).values, axis=1)
    rows = np.arange(len(cols))
    d = pd.Series(X.values[rows, cols], index=X.index)

    return d


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations_file", type=str)
    parser.add_argument("--annotation_file", type=str)
    parser.add_argument("--event_type", type=str)
    parser.add_argument("--margin", type=int, default=None)
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
        events["COORD_o"].str.split(":").str[0]  # .str.replace("chr", "")
    )
    events["EVENT_start"] = (
        events["COORD_o"].str.split(":").str[1].str.split("-").str[0]
    ).astype("int")
    events["EVENT_end"] = (
        events["COORD_o"].str.split(":").str[1].str.split("-").str[1].astype("int")
    )

    if margin is not None:
        events["Start"] = events["EVENT_start"]
        events["End"] = events["EVENT_end"]
    else:
        # start and end will be the next splice site junctions
        events["Start"] = events["CO_C1"].str.split(":").str[1].str.split("-").str[1]
        events["End"] = events["CO_C2"].str.split(":").str[1].str.split("-").str[0]
        margin = 0

    ## process mutations
    mutations = mutations.rename(
        columns={"lifted_start_hg38": "Start", "lifted_end_hg38": "End"}
    )
    mutations = mutations.loc[~mutations[["Start", "End"]].isnull().any(axis=1)]

    # intersect events with variants
    X = pr.PyRanges(
        events[["EVENT", "EVENT_start", "EVENT_end", "Chromosome", "Start", "End"]],
        int64=True,
    )
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
                "DepMap_ID",
            ]
        ],
        int64=True,
    )
    mapping = X.join(Y, slack=margin, how="left", report_overlap=True)

    # prepare outputs
    result = mapping.df.rename(
        columns={"Start_b": "Start_mutation", "End_b": "End_mutation"}
    )
    result = result.loc[result["Reference_Allele"] != "-1"]
    result["distance_to_3ss"] = get_abs_min(
        pd.concat(
            [
                result["EVENT_start"] - result["Start_mutation"],
                result["EVENT_start"] - result["End_mutation"],
            ],
            axis=1,
        )
    )
    result["distance_to_5ss"] = get_abs_min(
        pd.concat(
            [
                result["EVENT_end"] - result["Start_mutation"],
                result["EVENT_end"] - result["End_mutation"],
            ],
            axis=1,
        )
    )
    result["distance_to_closest_ss"] = get_abs_min(
        result[["distance_to_3ss", "distance_to_5ss"]]
    )

    # save
    result.to_csv(output_file, sep="\t", compression="gzip", index=False)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
