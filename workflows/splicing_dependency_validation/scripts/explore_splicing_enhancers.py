# Script purpose
# --------------
# Giving the coordinates of the acceptor (5') and donor (3') splice sites of an exon
# and a genome sequence,
#
# Outline
# -------

import argparse
import os
from tqdm import tqdm

import pandas as pd
import numpy as np
import pyfastx
import torch
from pangolin.model import L, W, AR, Pangolin
from pangolin.pangolin import one_hot_encode

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# variables
MODELS_SPLICING_PROB = [0, 2, 4, 6]
MODELS_SPLICING_USAGE = [1, 3, 5, 7]
MODELS_SPLICING = MODELS_SPLICING_PROB
NUCLEOTIDES = ["A", "C", "T", "G"]

"""
Development
-----------
reference_file = "/home/miquel/databases/data/GENCODE/GRCh38.p13.genome.fa.gz"
event_chr = "chr19" 
event_start = 10808569
event_end = 10808580
event_strand = "+"
margin_out = 30
margin_in = 15
output_dir = "."
pangolin_dir = "/home/miquel/repositories/Pangolin"
"""

##### FUNCTIONS #####
def load_models(pangolin_dir):
    models = []
    for i in MODELS_SPLICING:
        for j in range(1, 4):
            model = Pangolin(L, W, AR)
            weights_file = os.path.join(
                pangolin_dir, "pangolin", "models", "final.%s.%s.3.v2"
            ) % (j, i)
            if torch.cuda.is_available():
                weights = torch.load(weights_file)
            else:
                weights = torch.load(weights_file, map_location=torch.device("cpu"))
            model.load_state_dict(weights)
            model.eval()
            models.append(model)
    return models


def prep_sequence(chromosome, position, ref_seq, distance):
    try:
        seq = ref_seq[chromosome][
            position - 5001 - distance : position + 1 + 4999 + distance
        ].seq
    except Exception as e:
        print(e)
        print(
            "[Line %s]" % lnum,
            "WARNING, skipping variant: Could not get sequence, possibly because the variant is too close to chromosome ends. "
            "See error message above.",
        )
        seq = -1

    return seq


def mutagenize(mut_pos, mut_chr, distance, ref_seq):
    """
    Returns alt sequences at a given position.
    """
    # what is the reference nucleotide?
    nucl_ref = ref_seq[mut_chr][mut_pos - 1]  # minus 1 because it is 0 notation
    nucl_alts = list(set(NUCLEOTIDES) - set(nucl_ref))

    # prep reference and mutant sequences for scoring
    input_ref = prep_sequence(mut_chr, mut_pos, ref_seq, distance)
    input_alts = {}
    for idx, alt in enumerate(nucl_alts):
        input_alts["alt%s" % alt] = (
            input_ref[: 5000 + distance] + alt + input_ref[5000 + distance + 1 :]
        )  # we sum by 1 because the reference nucleotide is length one

    return input_ref, input_alts


def compute_score(ref_seq, alt_seq, strand, d, models):
    ref_seq = one_hot_encode(ref_seq, strand).T
    ref_seq = torch.from_numpy(np.expand_dims(ref_seq, axis=0)).float()
    alt_seq = one_hot_encode(alt_seq, strand).T
    alt_seq = torch.from_numpy(np.expand_dims(alt_seq, axis=0)).float()

    if torch.cuda.is_available():
        ref_seq = ref_seq.to(torch.device("cuda"))
        alt_seq = alt_seq.to(torch.device("cuda"))

    pangolin = []
    for j in range(4):
        score = []
        for model in models[3 * j : 3 * j + 3]:
            with torch.no_grad():
                # this outputs only the nucleotides in the distance
                ref = model(ref_seq)[0][[1, 4, 7, 10][j], :].cpu().numpy()
                alt = model(alt_seq)[0][[1, 4, 7, 10][j], :].cpu().numpy()
                if strand == "-":
                    ref = ref[::-1]
                    alt = alt[::-1]
                    
                # the chunk below should not be executed as sequences are the same
                l = 2 * d + 1
                ndiff = np.abs(len(ref) - len(alt))
                if len(ref) > len(alt):
                    alt = np.concatenate(
                        [alt[0 : l // 2 + 1], np.zeros(ndiff), alt[l // 2 + 1 :]]
                    )
                elif len(ref) < len(alt):
                    alt = np.concatenate(
                        [
                            alt[0 : l // 2],
                            np.max(alt[l // 2 : l // 2 + ndiff + 1], keepdims=True),
                            alt[l // 2 + ndiff + 1 :],
                        ]
                    )
                score.append(alt - ref)
        pangolin.append(np.mean(score, axis=0))  # summarize across tissue-based models

    pangolin = np.array(pangolin)
    loss = pangolin[np.argmin(pangolin, axis=0), np.arange(pangolin.shape[1])]
    gain = pangolin[np.argmax(pangolin, axis=0), np.arange(pangolin.shape[1])]
    return loss, gain


def compute_splice_scores(input_ref, input_alts, strand, distance, models):
    # compute splice scores for the corresponding gene strand
    scores_losses = {}
    scores_gains = {}
    for mut_name, input_alt in input_alts.items():
        scores_losses[mut_name], scores_gains[mut_name] = compute_score(
            input_ref, input_alt, strand, distance, models
        )

    return scores_losses, scores_gains


def explore_splicing_enhancers_single(
    position, event_chr, event_start, event_end, event_strand, margin_out, ref_seq, models
):
    # perform point mutagenesis for the positions nearby splice sites
    # choose distance to get information on both splice sites
    if position < event_start:
        # position is upstream of start
        distance = event_end - position + margin_out
    elif position > event_end:
        # position is downstream of end
        distance = position - event_start + margin_out
    else:
        # position is in the exon
        # take distance to the furthest one
        distance = max(position - event_start + margin_out, event_end - position + margin_out)

    # mutagenize reference sequence in that position
    input_ref, input_alts = mutagenize(position, event_chr, distance, ref_seq)

    # score
    scores_losses, scores_gains = compute_splice_scores(
        input_ref, input_alts, event_strand, distance, models
    )

    # prepare output
    ## losses
    loss = pd.DataFrame(scores_losses)
    loss["position"] = list(range(position - distance, position + distance + 1))
    loss = loss.melt(id_vars=["position"], var_name="alt", value_name="score_loss")
    ## gains
    gain = pd.DataFrame(scores_gains)
    gain["position"] = list(range(position - distance, position + distance + 1))
    gain = gain.melt(id_vars=["position"], var_name="alt", value_name="score_gain")

    result = pd.merge(loss, gain, on=["position", "alt"], how="inner")
    result["mut_pos"] = position

    return result


def explore_splicing_enhancers(
    event_chr, event_start, event_end, event_strand, margin_out, ref_seq, models
):
    # positions to mutagenize
    positions = list(range(event_start - margin_out, event_end + margin_out + 1))

    results = []
    for position in tqdm(positions):
        result = explore_splicing_enhancers_single(
            position, event_chr, event_start, event_end, event_strand, margin_out, ref_seq, models
        )
        results.append(result)

    results = pd.concat(results)
    results["event_chr"] = event_chr
    results["event_start"] = event_start
    results["event_end"] = event_end
    results["event_strand"] = event_strand

    return results


def make_plots(result, margin_out=30, output_dir=".", width=5, height=5):
    cm = 1 / 2.54
    
    X = result

    # init
    event = X["event_name"].unique()[0]
    start = X["event_start"].unique()[0]
    end = X["event_end"].unique()[0]
    mut_pos = X[["position", "mut_pos"]].drop_duplicates()

    # plot only scores for region of interest
    idx = (X["position"] >= start - margin_out) & (X["position"] <= end + margin_out)
    X = X.loc[idx]

    ## gain
    fig, ax = plt.subplots(figsize=(width * cm, height * cm))
    ax.scatter(x=X["position"], y=X["score_gain"], alpha=0.5, s=5, zorder=3)
    ax.axvspan(start, end, facecolor='green', alpha=0.25, zorder=1)
    plt.ylim(-1,1)
    plt.title("All Score Gain | %s" % event)

    plt.savefig(os.path.join(output_dir,"all_score_gain.png"), pad_inches=0)
    
    ## loss
    fig, ax = plt.subplots(figsize=(width * cm, height * cm))
    ax.scatter(x=X["position"], y=X["score_loss"], alpha=0.5, s=5, zorder=3)
    ax.axvspan(start, end, facecolor='green', alpha=0.25, zorder=1)
    plt.ylim(-1,1)
    plt.title("All Score Loss | %s" % event)
    
    plt.savefig(os.path.join(output_dir,"all_score_loss.png"), pad_inches=0)
    
    # for each mutant, get the maximum probability for each splice site and compute the 
    # average between the two to measure the change in splicing upon mutagenesis
    ## increase
    X = result
    X = (
        X.loc[X["position"].isin([start, end])]
        .groupby(["mut_pos", "alt", "position"])["score_gain"]
        .max()
        .reset_index()
        .groupby(["mut_pos", "alt"])["score_gain"]
        .mean()
        .reset_index()
    )
    
    fig, ax = plt.subplots(figsize=(width * cm, height * cm))
    ax.scatter(x=X["mut_pos"], y=X["score_gain"], alpha=0.5, s=5, zorder=3)
    ax.axvspan(start, end, facecolor='green', alpha=0.25, zorder=1)
    plt.xlabel("Mutagenesis Position")
    plt.ylabel("Predicted Splice Site Usage Change")
    plt.title("SS Score Gain | %s" % event)
    plt.ylim(-1,1)
    
    plt.savefig(os.path.join(output_dir,"ss_score_gain.png"), pad_inches=0)
        
    ## decrease
    X = result
    X = (
        X.loc[X["position"].isin([start, end])]
        .groupby(["mut_pos", "alt", "position"])["score_loss"]
        .min()
        .reset_index()
        .groupby(["mut_pos", "alt"])["score_loss"]
        .mean()
        .reset_index()
    )
    
    fig, ax = plt.subplots(figsize=(width * cm, height * cm))
    ax.scatter(x=X["mut_pos"], y=X["score_loss"], alpha=0.5, s=5, zorder=3)
    ax.axvspan(start, end, facecolor='green', alpha=0.25, zorder=1)
    plt.xlabel("Mutagenesis Position")
    plt.ylabel("Predicted Splice Site Usage Change")
    plt.title("SS Score Loss | %s" % event)
    plt.ylim(-1,1)
    
    plt.savefig(os.path.join(output_dir,"ss_score_loss.png"), pad_inches=0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--event_name", type=str)
    parser.add_argument("--event_chr", type=str)
    parser.add_argument("--event_start", type=int)
    parser.add_argument("--event_end", type=int)
    parser.add_argument("--event_strand", type=str)
    parser.add_argument("--margin_out", type=int)
    parser.add_argument("--reference_file", type=str)
    parser.add_argument("--pangolin_dir", type=str)
    parser.add_argument("--output_dir", type=str)

    args = parser.parse_args()

    return args


def main():
    # unpack arguments
    args = parse_args()
    event_name = args.event_name
    event_chr = args.event_chr
    event_start = args.event_start
    event_end = args.event_end
    event_strand = args.event_strand
    margin_out = args.margin_out
    reference_file = args.reference_file
    pangolin_dir = args.pangolin_dir
    output_dir = args.output_dir

    # load
    models = load_models(pangolin_dir)
    ref_seq = pyfastx.Fasta(reference_file)

    # analysis
    result = explore_splicing_enhancers(
        event_chr, event_start, event_end, event_strand, margin_out, ref_seq, models
    )
    results["event_name"] = event_name
    
    # save
    os.makedirs(output_dir, exist_ok=True)
    result.to_csv(os.path.join(output_dir,"mutagenesis.tsv.gz"), sep="\t", compression="gzip", index=False)
    
    # viz
    make_plots(result, margin_out, output_dir, 12, 12)
    

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
