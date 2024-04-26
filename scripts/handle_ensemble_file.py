#!/usr/bin/env python

import os
from Bio import SeqIO

# id: 
# file|:id

def get_seqRecord_ensembl(id):
    [id_file, id_rec] = id.split("|:")
    
    search_path = "../data/iterations/ensemble_cds_fasta"
    fasta_file = ""
    
    # print(id_file)
    
    for root, dir, files in os.walk(search_path):
        if id_file in files:
            fasta_file = os.path.join(root, id_file)
    
    id_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # print(id_dict[id_rec])
    # print(type(id_dict[id_rec]))
    
    return id_dict[id_rec]

