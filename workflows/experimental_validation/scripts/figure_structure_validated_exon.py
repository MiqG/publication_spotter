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
import pandas as pd
from yasara_macro import yasara_macro as ys

# variables
SAVE_PARAMS = "X=4024,Y=4024,Zoom=1.0,Atoms=Balls,LabelShadow=No,SecAlpha=100,Display=Yes,Outline=0.0,Background=No" # fog, lighning, 
ORIENTATION = {
    # obtained in Yasara app with 'PosOriObj 1'
    "FNBP1": "X=13.888, Y=-7.311, Z=177.957, Alpha=49.549, Beta=-88.327, Gamma=-9.494",
    "PPP1R12A": "X=-5.706, Y=7.584, Z=188.305, Alpha=-86.538, Beta=102.258, Gamma=35.111",
    "KRAS": "X=-2.970, Y=1.420, Z=88.013, Alpha=-56.989, Beta=-153.165, Gamma=-72.566",
    "PRPF18": "X=-7.272, Y=-4.613, Z=125.544, Alpha=170.086, Beta=77.046, Gamma=44.936",
    "YAP1": "X=-11.976, Y=0.804, Z=149.225, Alpha=-52.342, Beta=115.258, Gamma=-83.962",
    "VLDLR": "X=12.456, Y=-3.589, Z=147.055, Alpha=-59.008, Beta=-94.040, Gamma=60.905",
    "RCC1": "X=7.608, Y=2.071, Z=115.215, Alpha=-97.592, Beta=-51.287, Gamma=12.338",
    "NUP85": "X=-11.376, Y=0.327, Z=154.464, Alpha=-149.661, Beta=-106.654, Gamma=-45.563"
}


"""
Development
-----------
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
SUPPORT_DIR = os.path.join(ROOT,"support")
RESULTS_DIR = os.path.join(ROOT,'results','splicing_dependency_drugs')
pdb_file = os.path.join(SUPPORT_DIR,"structures","ESMfold","FNBP1includedENSP00000413625.pdb")
exon_info_file = os.path.join(ROOT,"support","exon_mapping_validated_exons.tsv")
gene_oi = "ENSP00000413625"
yasara_path = "~/yasara/yasara"
figs_dir = "toremove"
"""

##### FUNCTIONS #####
def run_yasara(yasara_path, pdb_file, locs_oi_inc, figs_dir):
    events = list(locs_oi_inc.keys())
    pdb_name = os.path.basename(pdb_file).replace(".pdb","")

    # included
    cmd_init = [
        # load
        "LoadPDB {pdb_file}".format(pdb_file=pdb_file),
        # style
        "HideObj 1",
        "ShowSecStrObj all,Style=Ribbon",
        # color
        "ColorObj all,Gray",  # baseline
        # light
        "LightSource Alpha=000,Gamma=000,Ambience=000,Ambience2=000,Shadow=000",
        # fog
        "Fog Density=000"
    ]

    cmd_coloring = [
        [
            "ColorRes {start}-{end},a9a9a9".format(
                start=locs_oi_inc[event]["upstream_start"],
                end=locs_oi_inc[event]["upstream_end"],
            ),  # upstream
            "ColorRes {start}-{end},ffa500".format(
                start=locs_oi_inc[event]["start"], end=locs_oi_inc[event]["end"]
            ),  # exon
            "ColorRes {start}-{end},a9a9a9".format(
                start=locs_oi_inc[event]["downstream_start"],
                end=locs_oi_inc[event]["downstream_end"],
            ),  # downstream
            "PosOriObj 01, {orientation}".format(orientation=ORIENTATION[event]) # orientation
        ]
        for event in events
    ]
    cmd_coloring = sum(cmd_coloring,[])
    cmd_save = [
        # save
        ## scene
        "SaveSce {output_file}".format(
            output_file=os.path.join(figs_dir, pdb_name+".sce")
        ),
        ## image
        "RayTrace {output_file},{save_params}".format(
            output_file=os.path.join(figs_dir, pdb_name+".png"),
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
    parser.add_argument("--pdb_file", type=str)
    parser.add_argument("--pdb_exc_file", type=str)
    parser.add_argument("--exon_info_file", type=str)
    parser.add_argument("--yasara_path", type=str)
    parser.add_argument("--gene_oi", type=str)
    parser.add_argument("--figs_dir", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    pdb_file = args.pdb_file
    exon_info_file = args.exon_info_file
    gene_oi = args.gene_oi
    yasara_path = args.yasara_path
    figs_dir = args.figs_dir

    # load
    print("Loading data...")
    exon_info = pd.read_table(exon_info_file)
    exon_info = exon_info.loc[exon_info["gene"]==gene_oi]
    locs_oi_inc = {gene_oi:
        {
            "upstream_start": 1,
            "upstream_end": max(1, exon_info["exon_proteoform_start"].values[0] - 1),
            "start": exon_info["exon_proteoform_start"].values[0],
            "end": exon_info["exon_proteoform_end"].values[0],
            "downstream_start": min(
                exon_info["main_inclusion_proteoform_length"].values[0], 
                exon_info["exon_proteoform_end"].values[0] + 1
            ),
            "downstream_end": exon_info["main_inclusion_proteoform_length"].values[0]
        }              
    }
    
    # visualize
    print("Saving data...")
    os.makedirs(figs_dir, exist_ok=True)
    exit = run_yasara(yasara_path, pdb_file, locs_oi_inc, figs_dir)
    assert exit==1
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")