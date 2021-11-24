#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
#
# Script purpose
# --------------
# Preprocess and clean CCLE data.
#
# Outline
# -------
# 1. exchange to official sample identifiers
# 2. Drop rows with all missing values
# 3. Add cancer types to metadata

import argparse
import pandas as pd

# variables
SAVE_PARAMS = {"sep": "\t", "index": False, "compression": "gzip"}

"""
Development
-----------
import os
ROOT = '/home/miquel/projects/publication_splicing_dependency'
RAW_DIR = os.path.join(ROOT,'data','raw')
treatments_file = os.path.join(RAW_DIR,'GDC/TCGA-BRCA/harmonized/Clinical/Clinical_Supplement/1e6b79ff-9787-4cbe-b19d-ebabb6b43589','nationwidechildrens.org_clinical_drug_brca.txt')
drug_correction_file = os.path.join(RAW_DIR,'articles','Spainhour2017','DrugCorrection.csv')
"""

##### FUNCTIONS #####
def load_data(treatments_file, drug_correction_file):
    treatments = pd.read_table(treatments_file, skiprows=range(1, 3))
    drug_correction = pd.read_csv(drug_correction_file)

    return treatments, drug_correction


def preprocess_treatments(treatments, drug_correction):
    df = pd.merge(
        treatments,
        drug_correction.rename(columns={"Correction": "DRUG_NAME"}),
        left_on="pharmaceutical_therapy_drug_name",
        right_on="OldName",
        how="left",
    )

    return df


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--treatments_file", type=str, default="")
    parser.add_argument("--drug_correction_file", type=str, default="")
    parser.add_argument("--output_file", type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_args()
    treatments_file = args.treatments_file
    drug_correction_file = args.drug_correction_file
    output_file = args.output_file

    # load
    treatments, drug_correction = load_data(treatments_file, drug_correction_file)
    output = preprocess_treatments(treatments, drug_correction)

    # save
    output.to_csv(output_file, **SAVE_PARAMS)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
