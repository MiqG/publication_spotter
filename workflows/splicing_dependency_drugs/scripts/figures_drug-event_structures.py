#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# For each drug, measure the shortest path length from each of its targets
# to each of the genes significantly associated to the drug sensitivity profile

import os
import argparse
import copy
import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.Data import IUPACData
from Bio.Seq import Seq
from yasara_macro import yasara_macro as ys

# variables
SAVE_PARAMS = "X=1024,Y=1024,Zoom=1.0,Atoms=Balls,LabelShadow=No,SecAlpha=100,Display=Yes,Outline=0.0,Background=No" # fog, lighning, 
ORIENTATION = {
    "HsaEX0038400": {
        "inc": "X=-8.257, Y=-26.379, Z=71.680, Alpha=-36.979, Beta=-139.042, Gamma=10.970",
        "exc": "X=0.626, Y=-22.848, Z=64.930, Alpha=42.282, Beta=-126.465, Gamma=-7.557"
    }
}


"""
Development
-----------
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,'results','splicing_dependency_drugs')
pdb_inc_file = os.path.join(SUPPORT_DIR,"MDM2_HsaEX0038400-iso0-included_unrelaxed_rank_1_model_3.pdb")
pdb_exc_file = os.path.join(SUPPORT_DIR,"MDM2_HsaEX0038400-iso0-excluded_unrelaxed_rank_1_model_3.pdb")
metadata_file = os.path.join(RESULTS_DIR,'files','structure_inference','proteoforms_merged.tsv.gz')
gene_oi = "MDM2"
seq_id = "MDM2_HsaEX0038400-iso0"
margin = 10
yasara_path = "~/yasara/yasara"
figs_dir = "toremove"
"""

##### FUNCTIONS #####
def load_data(metadata_file, pdb_inc_file, pdb_exc_file):
    metadata = pd.read_table(metadata_file)
    df_inc = PandasPdb().read_pdb(pdb_inc_file)
    df_exc = PandasPdb().read_pdb(pdb_exc_file)
    pdb_inc = open(pdb_inc_file, "r").read()
    pdb_exc = open(pdb_exc_file, "r").read()

    return metadata, df_inc, df_exc, pdb_inc, pdb_exc


def prepare_data(df_inc, df_exc):
    # make proper dataframes
    ## included
    df_inc = df_inc.df["ATOM"][
        ["residue_name", "chain_id", "residue_number"]
    ].drop_duplicates()
    df_inc["residue_title"] = df_inc["residue_name"].str.title()
    df_inc["residue_symbol"] = df_inc["residue_title"].map(
        IUPACData.protein_letters_3to1
    )
    seq_inc = "".join(df_inc.loc[df_inc["chain_id"] == "A", "residue_symbol"].to_list())
    res_inc = df_inc.loc[df_inc["chain_id"] == "A", "residue_number"].to_list()

    ## excluded
    df_exc = df_exc.df["ATOM"][
        ["residue_name", "chain_id", "residue_number"]
    ].drop_duplicates()
    df_exc["residue_title"] = df_exc["residue_name"].str.title()
    df_exc["residue_symbol"] = df_exc["residue_title"].map(
        IUPACData.protein_letters_3to1
    )
    seq_exc = "".join(df_exc.loc[df_exc["chain_id"] == "A", "residue_symbol"].to_list())
    res_exc = df_exc.loc[df_exc["chain_id"] == "A", "residue_number"].to_list()

    print(seq_inc, df_inc.groupby(["chain_id"]).size(), df_inc.head())
    print(seq_exc, df_exc.groupby(["chain_id"]).size(), df_exc.head())

    return df_inc, seq_inc, res_inc, df_exc, seq_exc, res_exc


def get_positions_aa(metadata, seq_inc):
    # exon sequences
    seqs_oi = (
        metadata.loc[metadata["GENE"] == gene_oi, ["EVENT", "event_aa"]]
        .drop_duplicates()
        .set_index("EVENT")
        .to_dict()["event_aa"]
    )

    # get locations of events
    ## find
    locs_oi_inc = {event: seq_inc.find(aa) for event, aa in seqs_oi.items()}
    ## remove not found
    locs_oi_inc = {
        event: loc_start
        for event, loc_start in locs_oi_inc.items()
        if loc_start != (-1)
    }
    ## add end
    locs_oi_inc = {
        event: {"start": loc_start, "end": loc_start + len(seqs_oi[event])}
        for event, loc_start in locs_oi_inc.items()
    }
    ## add margins
    locs_oi_inc = {
        event: {
            "start": loc["start"],
            "end": loc["end"],
            "upstream_start": loc["start"] - margin - 1,
            "upstream_end": loc["start"] - 1,
            "downstream_start": loc["end"] + 1,
            "downstream_end": loc["end"] + margin + 1,
        }
        for event, loc in locs_oi_inc.items()
    }

    ## locations of excluded isoform
    locs_oi_exc = copy.deepcopy(locs_oi_inc)
    for event in locs_oi_exc.keys():
        locs_oi_exc[event]["downstream_start"] = (
            locs_oi_exc[event]["downstream_start"] - len(seqs_oi[event]) - 1
        )
        locs_oi_exc[event]["downstream_end"] = (
            locs_oi_exc[event]["downstream_end"] - len(seqs_oi[event]) - 1
        )

    print(seqs_oi)
    print(locs_oi_inc)
    print(locs_oi_exc)

    return seqs_oi, locs_oi_inc, locs_oi_exc


def run_yasara(yasara_path, pdb_inc_file, pdb_exc_file, locs_oi_inc, locs_oi_exc):
    events = list(locs_oi_inc.keys())

    # included
    cmd_init = [
        # load
        "LoadPDB {pdb_file}".format(pdb_file=pdb_inc_file),
        # style
        "HideObj 1",
        "ShowSecStrObj all,Style=Tube",
        # color
        "ColorObj all,Gray",  # baseline
    ]

    cmd_coloring = [
        [
            "ColorRes {start}-{end},Blue".format(
                start=locs_oi_inc[event]["upstream_start"],
                end=locs_oi_inc[event]["upstream_end"],
            ),  # upstream
            "ColorRes {start}-{end},154".format(
                start=locs_oi_inc[event]["start"], end=locs_oi_inc[event]["end"]
            ),  # exon
            "ColorRes {start}-{end},Green".format(
                start=locs_oi_inc[event]["downstream_start"],
                end=locs_oi_inc[event]["downstream_end"],
            ),  # downstream
            "PosOriObj 01, {orientation}".format(orientation=ORIENTATION[event]["inc"]) # orientation
        ]
        for event in events
    ]
    cmd_coloring = sum(cmd_coloring,[])
    cmd_save = [
        # save
        ## scene
        "SaveSce {output_file}".format(
            output_file=os.path.join(figs_dir, "structure-inc.sce")
        ),
        ## image
        "RayTrace {output_file},{save_params}".format(
            output_file=os.path.join(figs_dir, "structure-inc.png"),
            save_params=SAVE_PARAMS,
        ),
    ]

    cmd = cmd_init + cmd_coloring + cmd_save

    macro = ys.interface(yasara_path)
    macro.write(cmd)
    exit = macro.submit()
    assert exit == 0
    
    # excluded
    cmd_init = [
        # load
        "LoadPDB {pdb_file}".format(pdb_file=pdb_exc_file),
        # style
        "HideObj 1",
        "ShowSecStrObj all,Style=Tube",
        # color
        "ColorObj all,Gray",  # baseline
    ]
    cmd_coloring = [
        [
            "ColorRes {start}-{end},Blue".format(
                start=locs_oi_exc[event]["upstream_start"],
                end=locs_oi_exc[event]["upstream_end"],
            ),  # upstream
            "ColorRes {start}-{end},Green".format(
                start=locs_oi_exc[event]["downstream_start"],
                end=locs_oi_exc[event]["downstream_end"],
            ),  # downstream
            "PosOriObj 01, {orientation}".format(orientation=ORIENTATION[event]["exc"]) # orientation

        ]
        for event in events
    ]
    cmd_coloring = sum(cmd_coloring,[])
    cmd_save = [
        # save
        ## scene
        "SaveSce {output_file}".format(
            output_file=os.path.join(figs_dir, "structure-exc.sce")
        ),
        ## image
        "RayTrace {output_file},{save_params}".format(
            output_file=os.path.join(figs_dir, "structure-exc.png"),
            save_params=SAVE_PARAMS,
        ),
    ]
    cmd = cmd_init + cmd_coloring + cmd_save

    macro = ys.interface(yasara_path)
    macro.write(cmd)
    exit = macro.submit()
    assert exit == 0

    return 1


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_inc_file", type=str)
    parser.add_argument("--pdb_exc_file", type=str)
    parser.add_argument("--metadata_file", type=str)
    parser.add_argument("--gene_oi", type=str)
    parser.add_argument("--seq_id", type=str)
    parser.add_argument("--margin", type=int)
    parser.add_argument("--yasara_path", type=str)
    parser.add_argument("--figs_dir", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    pdb_inc_file = args.pdb_inc_file
    pdb_exc_file = args.pdb_exc_file
    metadata_file = args.metadata_file
    gene_oi = args.gene_oi
    seq_id = args.seq_id
    margin = args.margin
    yasara_path = args.yasara_path
    figs_dir = args.figs_dir

    # load
    print("Loading data...")
    metadata, df_inc, df_exc, pdb_inc, pdb_exc = load_data(
        metadata_file, pdb_inc_file, pdb_exc_file
    )

    # prepare
    print("Preparing...")
    df_inc, seq_inc, res_inc, df_exc, seq_exc, res_exc = prepare_data(df_inc, df_exc)
    seqs_oi, locs_oi_inc, locs_oi_exc = get_positions_aa(metadata, seq_inc)

    # visualize
    print("Saving data...")
    os.makedirs(figs_dir, exist_ok=True)
    exit = run_yasara(yasara_path, pdb_inc_file, pdb_exc_file, locs_oi_inc, locs_oi_exc)
    assert exit==1
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
