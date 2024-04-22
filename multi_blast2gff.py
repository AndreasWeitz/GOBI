#!/usr/bin/env python

import argparse
import os
import subprocess

def multi_blast2gff(in_path, out_path):
    for (root, dirs, files) in os.walk(in_path):
        for file in files:
            in_file_path = os.path.join(root, file)
            out_file_path = os.path.join(out_path, os.path.basename(in_file_path)) + '.gff'
            with open(out_file_path, 'w') as out_file_handle:
                    subprocess.run(f"python blast2gff.py -Q -b {in_file_path}", shell=True, stdout=out_file_handle)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-i', '--inpath', help='Multiblast-output directory', required=True, type=str)
    parser.add_argument('-o', '--outpath', help='Output directory for GFF files', required=True, type=str)

    args = parser.parse_args()
    inpath = args.inpath
    outpath = args.outpath

    if not os.path.exists(outpath):
         os.mkdir(outpath)

    multi_blast2gff(inpath, outpath)
