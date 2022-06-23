# Script purpose
# --------------
# Giving the coordinates of the acceptor (5') and donor (3') splice sites of an exon
# and a genome sequence,
#
# Outline
# -------

import argparse
import os
import pyfastx

"""
Development
-----------
reference_file = "/home/miquel/databases/data/GENCODE/GRCh38.p13.genome.fa.gz"
event_name = "HsaEX0006970"
event_chr = "chr12" 
event_start = 123751110
event_end = 123751229
event_strand = "+"
margin_out = 100
"""

##### FUNCTIONS #####
def prep_sequence(chromosome, position, ref_seq, distance, size):
    try:
        seq = ref_seq[chromosome][
            position - 5001 - distance : position + size + 4999 + distance
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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--event_name", type=str)
    parser.add_argument("--event_chr", type=str)
    parser.add_argument("--event_start", type=int)
    parser.add_argument("--event_end", type=int)
    parser.add_argument("--event_strand", type=str)
    parser.add_argument("--margin_out", type=int)
    parser.add_argument("--reference_file", type=str)
    parser.add_argument("--output_file", type=str)

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
    output_file = args.output_file

    # load
    seq_ref = pyfastx.Fasta(reference_file)

    # get sequence of interest
    ## subset
    seq_oi = seq_ref[event_chr][
        (event_start - 1 - margin_out) : (event_end + margin_out)
    ]

    full_coords = seq_oi.name
    if event_strand == "+":
        seq = seq_oi.seq
    elif event_strand == "-":
        seq = seq_oi.antisense  # reverse-complement

    # save
    with open(output_file, "w") as f:
        lines = "\n".join(
            [
                ">name={event_name}|event_coords={event_coords}|full_coords={full_coords}".format(
                    event_name=event_name,
                    event_coords=event_chr
                    + ":"
                    + str(event_start)
                    + "-"
                    + str(event_end)
                    + ":"
                    + event_strand,
                    full_coords=full_coords+":"+event_strand,
                ),
                seq[:margin_out],  # margin
                seq[margin_out:][:-margin_out],  # exon
                seq[-margin_out:],  # margin
            ]
        )
        f.writelines(lines)


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
