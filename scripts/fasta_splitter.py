#!/usr/bin/env python3

from Bio import SeqIO
import os
import argparse

def split_fasta(input_file: str, out_dir: str):
    """
    Split a multi-sequence FASTA file into individual files.
    """
    for i, record in enumerate(SeqIO.parse(input_file, 'fasta'), start=1):
        output_filename = f"{record.id}.fasta"
        with open(f"{out_dir}/{output_filename}", 'x') as f:
            SeqIO.write(record, f, 'fasta')
        print(f"Sequence {i} saved to {output_filename}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-f', '--fasta', help='path to multi sequence fasta', required=True)
    parser.add_argument('-o', '--outdir', help='path to output directory', required=True)

    args = parser.parse_args()
    fasta_file = args.fasta
    outdir = args.outdir
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    split_fasta(fasta_file, outdir)
