#!/usr/bin/env python3

import os
import subprocess
import argparse
import datetime
from playsound import playsound

def play_audio():
    audio_path = '../data/others/alarm_sound.mp3'
    if os.path.exists(audio_path):
        playsound(audio_path)
    else:
        pass

def multi_blast(fasta_dir: str, blast_db: str, out_dir: str, btype="tblastx", e_value=0.001, num_threads=16, out_fmt=6):
    """
    Query all files in `fasta_dir` against `blast_db` and 
    save non-empty results in `out_dir`.
    """
    start_time = datetime.datetime.now()
    for (root, dirs, files) in os.walk(fasta_dir):
        for file in files:
            query_file_path = os.path.join(root, file)
            out_file_path = os.path.join(out_dir, f"results_{'.'.join(file.split('.')[:-1])}")
            try:
                subprocess.run(f"{btype} -db {blast_db} -query {query_file_path} -out {out_file_path} -evalue {e_value} -num_threads {num_threads} -outfmt {out_fmt}", shell=True, check=True)
            except subprocess.CalledProcessError as err:
                print(f'Error while processing file {query_file_path}')
                print(err)
                continue

            no_hits_found = False
            if out_fmt == 1:
                with open(out_file_path) as file:
                    for line in file:
                        if "***** No hits found *****" in line:
                            no_hits_found = True
                            break
            elif out_fmt == 6:
                no_hits_found = os.stat(out_file_path).st_size == 0
            
            if no_hits_found:
                os.remove(out_file_path)
            
            print(f"Processed {os.path.basename(os.path.join(root, file))}")
    end_time = datetime.datetime.now()
    print(f"Job started at {start_time}\nJob finished at {end_time}\nJob took {end_time - start_time}")
    play_audio()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-q', '--querydir', help='path directory containing queries', required=True, type=str)
    parser.add_argument('-db', '--database', help='name of blast database', required=True, type=str)
    parser.add_argument('-o', '--outdir', help='path to output directory', required=True, type=str)
    parser.add_argument('-b', '--blast', help='blast program, default: tblastx', default='tblastx', type=str)
    parser.add_argument('-e', '--evalue', help='e-value cutoff, default: 0.001', default=0.001, type=float)
    parser.add_argument('-t', '--num_threads', help='number of threads for blast, default: 16', default=16, type=int)
    parser.add_argument('-f', '--outfmt', help='blast output format, default: 6', default=6, type=int)

    args = parser.parse_args()
    querydir = args.querydir
    database = args.database
    out_dir = args.outdir
    btype = args.blast
    e_value = args.evalue
    num_threads = args.num_threads
    out_fmt = args.outfmt

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    multi_blast(
            fasta_dir=querydir, 
            blast_db=database, 
            out_dir=out_dir, 
            btype=btype, 
            e_value=e_value, 
            num_threads=num_threads,
            out_fmt=out_fmt
        )
