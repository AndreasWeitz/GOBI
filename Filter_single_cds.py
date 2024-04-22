#!/usr/bin/env python

import argparse as ap
from BCBio import GFF

def filter(gff):
    counter = {}
    with open(gff) as file:
        for line in file:
            desc = line.split("\t")[8]
            id = desc.split(';')[1].split('=')[1]
            if id not in counter.keys():
                counter[id] = 1
            else:
                counter[id] += 1
    
    lines = []
    print(counter)
    
    with open(gff) as file:
        for line in file:
            desc = line.split("\t")[8]
            id = desc.split(';')[1].split('=')[1]
            if counter[id] > 1:
                lines.append(line)
    
    return lines
    

if __name__=="__main__":
    
    parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-gff_in','--gff_file_in', help="matching cds gff file")
    
    args = parser.parse_args()
    
    gff_in = args.gff_file_in    
    seq_records = filter(gff_in)
    
    gff_out = gff_in.replace(".gff", "_filtered.gff")

    with open(gff_out, "w") as file:
        for elem in seq_records:
            file.write(elem)

    print("DONE")
    
    