# Script purpose
# --------------
# Compute several transcriptome indices:
#  - expression of KI67
#  - Mitotic Index: average expression over 9 genes (CDKN3, ILF2, KDELR2, RFC4,
#    TOP2A, MCM3, KPNA2, CKS2, and CDC2), from Yang et al. 2016 (doi: https://doi.org/10.1186/s13059-016-1064-3)
#  - Stemness score: computed as in TCGAbiolinks::TCGAanalyze_Stemness(),
#    a function that generates the mRNAsi  using the stemness signature (12'955 genes)
#    trained on the PCBC's dataset as described in Malta et al. 2018 (doi: 10.1016/j.cell.2018.03.034)
#  - Tumor purity: obtained from TCGAbiolinks::Tumor.purity object.
#  - DEPTH: tumor heterogeneity
#
# Outline
# -------
# 1. Loads
#    - genexpr_file
#    - annotation_file
#    - translate_from_col
#    - translate_to_col

import os
import pandas as pd
import numpy as np
import argparse

from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri

base = importr("base")
TCGAbiolinks = importr("TCGAbiolinks")
DEPTH = importr("DEPTH")

# stemness signature as series
STEMSIG = r["PCBC_stemSig"]  # stemness signature
STEMSIG = pd.Series(STEMSIG, index=STEMSIG.names)

# tumor purity
pandas2ri.activate()
tp = r["Tumor.purity"]
char_cols = ["Sample.ID", "Cancer.type"]
for key in set(tp.columns) - set(char_cols):
    tp[key] = [x.replace(",", ".") for x in tp[key].values]
    tp[key] = tp[key].astype(float)

TUMORPURITY = tp.set_index("Sample.ID").copy()
TUMORPURITY.drop(columns=["Cancer.type"], inplace=True)

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}
MITOTIC_INDEX_SYMBOLS = [
    "CDKN3",
    "ILF2",
    "KDELR2",
    "RFC4",
    "TOP2A",
    "MCM3",
    "KPNA2",
    "CKS2",
    "CDK1",
]  # CDC2 was the previous name for CDK1
ADHESION_SIGNATURE = pd.read_table(
    os.path.join(
        "~",
        "projects",
        "sso_targets",
        "results",
        "adhesion_index",
        "files",
        "adhesion_signature.tsv",
    )
)

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
PREP_DIR = os.path.join(ROOT,'data','prep')
genexpr_file = os.path.join(PREP_DIR,'genexpr_tpm','CCLE.tsv.gz')
annotation_file = os.path.join(RAW_DIR,'VastDB','event_annotation-Hs2.tsv.gz')
translate_from = 'ENSEMBL'
translate_to = 'GENE'
"""

##### FUNCTIONS #####
def prepare_inputs(genexpr_file, annotation_file, translate_from, translate_to):
    genexpr = pd.read_table(genexpr_file, index_col=0)

    if annotation_file is not None:
        annot = pd.read_table(annotation_file)
        annot = annot[[translate_from, translate_to]].drop_duplicates().dropna()
        genexpr = genexpr.rename(
            index=annot.set_index(translate_from)[translate_to].to_dict()
        )

    return genexpr.loc[~genexpr.index.duplicated()]


def compute_single_adhesion_index(x, adhesion_signature):
    common_genes = set(x.index).intersection(adhesion_signature["gene"].to_list())

    X = x[common_genes]
    W = adhesion_signature.set_index("gene").loc[common_genes, "weight"]

    promotes_adhesion = adhesion_signature.set_index("gene").loc[
        common_genes, "promotes_adhesion"
    ]

    prod = X * W
    index = pd.Series(
        {
            "median_promotes_adhesion": np.median(prod[promotes_adhesion]),
            "median_represses_adhesion": np.median(prod[~promotes_adhesion]),
        }
    )
    index["adhesion_index"] = (
        index["median_promotes_adhesion"] - index["median_represses_adhesion"]
    )

    return index


class TranscriptomeIndices:
    def __init__(self, genexpr):
        self.genexpr_ = genexpr

    def get_mki67(self, genexpr):
        return genexpr.loc["MKI67"]

    def compute_mitotic_index(self, genexpr):
        mitotic_index = genexpr.loc[MITOTIC_INDEX_SYMBOLS].mean()
        mitotic_index.name = "mitotic_index"
        return mitotic_index

    def compute_stemness(self, genexpr):
        """as in TCGAbiolinks::TCGAanalyze_Stemness()"""
        # prepare
        X = genexpr
        w = STEMSIG
        common_genes = set(w.index).intersection(set(X.index))
        X = X.loc[common_genes]
        w = w[common_genes]

        # compute
        s = X.corrwith(w, method="spearman", axis=0)

        # scale scores to be between 0 and 1
        s = s - min(s)
        s = s / max(s)
        s.name = "stemness"
        return s

    def get_tumor_purity(self, samples):
        """
        Requires TCGA barcode up to vial (0 to 16).
        """
        return pd.DataFrame({}, index=samples).join(TUMORPURITY)

    def compute_heterogeneity_score(self, genexpr, thresh_normal=2, only_tumor=False):
        """
        Computed as in DEPTH.
        https://github.com/WangX-Lab/DEPTH/blob/master/DEPTH/R/DEPTH.R
        """
        # pick tumor and normal samples
        tumor_samples = [
            barcode for barcode in genexpr.columns if float(barcode[-3:-1]) < 10
        ]
        normal_samples = [
            barcode for barcode in genexpr.columns if float(barcode[-3:-1]) >= 10
        ]

        # compute DEPTH score
        if (len(normal_samples) > thresh_normal) & (only_tumor == False):
            ## considering normal samples
            score = (
                genexpr.sub(genexpr[normal_samples].mean(axis=1), axis=0) ** 2
            ).std(axis=0)
        else:
            ## only tumor samples
            score = (genexpr.sub(genexpr[tumor_samples].mean(axis=1), axis=0) ** 2).std(
                axis=0
            )
        score.name = "heterogeneity_score"
        return score

    def compute_adhesion_index(self, genexpr):
        result = genexpr.apply(
            lambda x: compute_single_adhesion_index(x, ADHESION_SIGNATURE)
        ).T
        return result

    def compute_transcriptome_indices(self, genexpr):
        # compute transcriptome indices
        mki67 = self.get_mki67(genexpr)
        mitotic_index = self.compute_mitotic_index(genexpr)
        stemness = self.compute_stemness(genexpr)
        tumor_purity = self.get_tumor_purity(list(genexpr.columns))
        adhesion_index = self.compute_adhesion_index(genexpr)

        # prepare output
        transcriptome_indices = pd.concat(
            [mki67, mitotic_index, stemness, tumor_purity, adhesion_index,], axis=1,
        )

        return transcriptome_indices

    def compute(self):
        self.transcriptome_indices_ = self.compute_transcriptome_indices(self.genexpr_)


def parse_args():
    """
    For scripting.
    """
    # required
    parser = argparse.ArgumentParser()
    parser.add_argument("--genexpr_file", type=str)
    parser.add_argument("--annotation_file", type=str, default=None)
    parser.add_argument("--translate_from", type=str, default=None)
    parser.add_argument("--translate_to", type=str, default=None)
    parser.add_argument("--output_file", type=str, default=None)

    args = parser.parse_args()
    return args


def main():
    # parse arguments
    args = parse_args()
    genexpr_file = args.genexpr_file
    annotation_file = args.annotation_file
    translate_from = args.translate_from
    translate_to = args.translate_to
    output_file = args.output_file

    # read and prepare data
    genexpr = prepare_inputs(genexpr_file, annotation_file, translate_from, translate_to)

    # compute transcriptome indices
    ti = TranscriptomeIndices(genexpr)
    ti.compute()

    # save
    ti.transcriptome_indices_.reset_index().to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
