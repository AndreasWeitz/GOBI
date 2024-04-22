#!/usr/bin/env python
from fileinput import filename
from Bio.Seq import translate
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse as ap
from pathlib import Path
import os


def write_fasta(id_set: set, fasta_file, output_path):
    seq_record_list = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if seq_record.id in id_set:
            seq_record_list.append(seq_record)

    SeqIO.write(seq_record_list, output_path, "fasta")
    return


def get_id_set(id_file):
    id_set = set()

    with open(id_file) as file:
        for line in file:
            id_set.add(line.replace("\n", "").strip())

    return id_set


if __name__=="__main__":

    parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter, description=__doc__)

    
    parser.add_argument('-fasta_in','--input_file_fasta', required=True, help="filename/path for fasta input file")
    parser.add_argument('-fasta_out','--output_file_fasta', help="filename/path for fasta output file")
    parser.add_argument('-id_file',required=True, help="text file containing all ids of fasta file to extract")
    args = parser.parse_args()


    fasta_file_in = args.input_file_fasta
    fasta_file_out = args.output_file_fasta
    id_file = args.id_file
    
    

    if not fasta_file_out:
        filename_fasta =  fasta_file_in.replace(".fasta", "")
        fasta_file_out = f"{filename_fasta}_extract.fasta"
        
    Path(os.path.dirname(fasta_file_out)).mkdir(parents=True, exist_ok=True)

    write_fasta(get_id_set(id_file), fasta_file_in, fasta_file_out)
    

   
    
    
    