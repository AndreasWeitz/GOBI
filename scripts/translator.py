#!/usr/bin/env python

from Bio.Seq import translate
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
import argparse as ap
from pathlib import Path
import os

# Translates each record in fasta file
# Writes new protein fasta file



def exon_translate(input_file, out_fasta):

    new_seq_record = {}
    seq_record_protein = []

    for seq_record in SeqIO.parse(input_file, "fasta"):
        id, range = seq_record.id.split("|:")
        description = seq_record.description.split(" ")
        if id in new_seq_record.keys():
            old_seq = new_seq_record.get(id)[0]
            new_seq = old_seq + seq_record.seq
            new_seq_record[id] = [new_seq, description[1:]]
        else:
            new_seq_record[id] = [seq_record.seq, description[1:]]

    for key, value in new_seq_record.items():
        prot_seq = value[0].translate()
        rec = SeqRecord(prot_seq, id=key, description=' '.join(value[1]))
        seq_record_protein.append(rec)

    SeqIO.write(seq_record_protein, out_fasta, "fasta")
    

def cds_translate(input_file, out_fasta):
    seq_record_protein = []
    id_set = set()
    counter = 0
    for seq_record in SeqIO.parse(input_file, "fasta"):
        if "|:" not in seq_record.id:
            prot_seq = seq_record.seq.translate()
            rec = SeqRecord(prot_seq, id=seq_record.id, description=seq_record.description)
            seq_record_protein.append(rec)
        else:
            id, range = seq_record.id.split("|:")
            if id not in id_set:
                print(id)
                id_set.add(id)
                ncbi_id, start, end = id.split("-")
                start = int(start.split(".")[0])
                end = int(end.split(".")[0])
                handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text", seq_start = start, seq_stop = end)
                record = SeqIO.read(handle, "gb")
                handle.close()
                for elem in record.features:
                    if elem.type == "CDS":
                        if "translation" in elem.qualifiers.keys():
                            prot_seq =  elem.qualifiers["translation"][0]
                            rec = SeqRecord(prot_seq, id=id, description= " ".join(seq_record.description.split(" ")[1:]))
                            seq_record_protein.append(rec)
                        else:
                            counter += 1
            else:
                pass
    print(len(id_set))
    print(counter)
    SeqIO.write(seq_record_protein, out_fasta, "fasta")

 
if __name__=="__main__":

    parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter, description=__doc__)

    
    parser.add_argument('-fasta_in','--input_file_fasta',required=True, help="filename/path for fasta input file")
    parser.add_argument('-fasta_out','--output_file_fasta', help="filename/path for fasta output file")
    parser.add_argument('-exon', action='store_true', help="Exon file as input")
    args = parser.parse_args()


    fasta_file_in = args.input_file_fasta
    fasta_file_out = args.output_file_fasta

    Entrez.email = "Ch.Thomas@campus.lmu.de"

    Path(os.path.dirname(fasta_file_out)).mkdir(parents=True, exist_ok=True)
   
    if args.exon:
        exon_translate(fasta_file_in, fasta_file_out)
    else:
        cds_translate(fasta_file_in, fasta_file_out)    
    
    