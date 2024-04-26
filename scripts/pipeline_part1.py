#!/usr/bin/env python

import subprocess
import argparse as ap
from blast import *
import sys
from pathlib import Path
import os 

# python pipeline_part1.py -dbpath Database/ -qp ../Test_data/Luciferase2.fasta -out test_out_blast -gff data/exon_luc.gff -fasta data/exon_luc.fasta -id data/exons_ids_luc.csv

if __name__=="__main__":
    parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-dbpath', '--db_file_path', help="path to db file", required=True)
    parser.add_argument('-qp', '--query_file_path', help="query, e.g. fasta file", required=True)
    parser.add_argument('-out', '--out_file', help="Output file (and path)", required=True)

    parser.add_argument('-gff','--output_file_gff', help="filename/path for gff output file", default='input_db.gff')
    parser.add_argument('-fasta','--output_file_fasta', help="filename/path for fasta output file", default='input_db.fasta')
    parser.add_argument('-id', '--id_file', help="csv file with 3 columns: id,[start],[end]", default=False)

    parser.add_argument('-e', '--evalue', help="evalue (similarity) for blast query", default=0.001)
    parser.add_argument('-btype', '--blasttype', help="blast type to execute: 'blastn', 'tblastn', 'tblastx'", default='blastn')
    parser.add_argument('-bgff', '--blast_output_gff', help="File name for gff file of blast output", default='output_db_blast.gff')
    
    parser.add_argument('-ftype', '--file_type_blast', action='store_true', help="Output type for blast, if set output can be converted to gff file")

    args = parser.parse_args()
    gff_file = args.output_file_gff
    fasta_file = args.output_file_fasta
    id_file = args.id_file
    db_path = args.db_file_path
    query_path = args.query_file_path
    out_file = args.out_file
    evalue = args.evalue
    btype = args.blasttype
    blast_gff = args.blast_output_gff
    ftype = args.file_type_blast

    
    Path(os.path.dirname(out_file)).mkdir(parents=True, exist_ok=True)
    if id_file:
        Path(os.path.dirname(gff_file)).mkdir(parents=True, exist_ok=True)
        Path(os.path.dirname(fasta_file)).mkdir(parents=True, exist_ok=True)
        Path(os.path.dirname(db_path)).mkdir(parents=True, exist_ok=True)
        subprocess.run(["python", "gff_writer.py", "-gff", gff_file, '-fasta', fasta_file, '-id', id_file])
    else:
        fasta_file = None

        
    if ftype:
        Path(os.path.dirname(blast_gff)).mkdir(parents=True, exist_ok=True)
        if btype == 'blastn':
            blastn(db_path, query_path, out_file, fasta_file, evalue, outtype=6)
        if btype == 'tblastn':
            tblastn(db_path, query_path, out_file, fasta_file, evalue, outtype=6)
        if btype == 'tblastx':
            tblastx(db_path, query_path, out_file, fasta_file, evalue, outtype=6)
        with open(blast_gff, "w") as out_handle:
            subprocess.run(["python", "blast2gff.py", '-b', out_file], stdout=out_handle)
            
    else:
        if btype == 'blastn':
            blastn(db_path, query_path, out_file, fasta_file, evalue)
        if btype == 'tblastn':
            tblastn(db_path, query_path, out_file, fasta_file, evalue)
        if btype == 'tblastx':
            tblastx(db_path, query_path, out_file, fasta_file, evalue)
    
    print("task completed")

