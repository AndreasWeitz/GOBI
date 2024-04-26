#!/usr/bin/env python

from get_gb import fetch_sequence
from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse as ap
import pandas as pd
import math
from pathlib import Path
import os 


def rewrite_features(sequence_record):
    """
    Generates exons from CDS entry

    Args:
        sequence_record (SeqRecord): A SeqRecord 

    Returns:
        SeqRecord : SeqRecord with exons as additional feature
    """
    
    feature = sequence_record.features
    sub_feature_list = []
    for elem in feature:
        if elem.type == "CDS":
            cds_qualifiers =  elem.qualifiers
            if "db_xref" in cds_qualifiers.keys():
                cds_qualifiers["ID"] = ["".join(["CDS:", cds_qualifiers.get("db_xref")[0].split(":")[1]])]
                cds_qualifiers["locus_tag"] = ["".join(["CDS:", cds_qualifiers.get("db_xref")[0].split(":")[1]])]
            elif "locus_tag" in cds_qualifiers.keys(): 
                cds_qualifiers["ID"] = ["".join(["CDS:", cds_qualifiers.get("locus_tag")[0]])]
            else: 
                cds_qualifiers["ID"]=[f"CDS:{elem.id}"]
                cds_qualifiers["locus_tag"] = ["".join(["CDS:", cds_qualifiers.get("db_xref")[0].split(":")[1]])]           
            
            cds_feature = elem
            locations = cds_feature.location
            location_list = locations.parts
            for loc in location_list:
                sub_feature_list.append(SeqFeature(FeatureLocation(int(loc.start), int(loc.end), strand = int(loc.strand)),
                                                type="exon", qualifiers = {"source": "feature",
                                                                        "gene": cds_qualifiers.get("gene"),
                                                                        "codon_start" : cds_qualifiers.get("codon_start"),
                                                                        "db_xref": cds_qualifiers.get("db_xref")}))
            elem.qualifiers = cds_qualifiers
            new_features = feature + sub_feature_list
            sequence_record.features = new_features
    return sequence_record


def get_seq_record_gff(id_dict):
    """
    Creates sequence record objects, which are written to a gff file

    Args:
        id_dict (dict): _description_

    Returns:
        list: list of sequence records
    """
    seq_record_for_gff = []
    for key, value in id_dict.items(): # Iterate over all entries in dict
        for v in value:
            if not v:  
                # Entry without range
                # Whole sequence is collected
                sequence_record = fetch_sequence(key)
            else:
                (start, end) = v
                # Sequence is collected by range
                sequence_record = fetch_sequence(key, int(start), int(end))
                new_id = f"{sequence_record.id}-{start}-{end}" # Rewrite id
                sequence_record.id = new_id
            # print(rewrite_features(sequence_record))
            seq_record_for_gff.append(rewrite_features(sequence_record))
    return seq_record_for_gff


def get_seq_record_fasta(seq_records):
    """
    

    Args:
        seq_records (_type_): _description_

    Returns:
        _type_: _description_
    """
    id_set = set()
    seq_record_list = []
    for seq_rec in seq_records:
        id = seq_rec.id

        if 'ENS' not in id and 'BRAKER' not in id:
            feature = seq_rec.features
            sequence = seq_rec.seq
            for elem in feature:
                if elem.type == "exon":
                    location = elem.location
                    new_id = f"{id}|:{location.start}-{location.end}"
                    if new_id not in id_set:
                        rec = SeqRecord(Seq(sequence[location.start:location.end]), id=new_id, description=f"{seq_rec.description}")
                        id_set.add(new_id)
                        seq_record_list.append(rec)
        else:
            seq_record_list.append(seq_rec)
    return seq_record_list


def get_id_dict(id_file):
    """
    Read id file to dictionary

    Args:
        id_file (file): _description_

    Returns:
        dict: _description_
    """
    df = pd.read_csv(id_file)
    print(df)
    id_dict = {}
    for i in range(len(df)):
        if not math.isnan(df["start"][i]):
            if str(df["id"][i]) in id_dict.keys():
                id_dict[df["id"][i]].append(tuple([df["start"][i], df["end"][i]]))
            else:
                id_dict[df["id"][i]] = [tuple([df["start"][i], df["end"][i]])]
        else:
            if df["id"][i] in id_dict.keys():
                id_dict[df["id"][i]].append(None)
            else:
                id_dict[df["id"][i]] = [None]
    print(id_dict)
    return id_dict


def get_cds_record(seq_records):
    """
    Rewrites all SeqRecord to get CDS for fasta file

    Args:
        seq_records (list): list of SeqRecord

    Returns:
        list: list of SeqRecord, containing CDS
    """
    seq_record_list = []
    for seq_rec in seq_records:
        id = seq_rec.id

        if 'ENS' not in id and 'BRAKER' not in id:
            feature = seq_rec.features
            sequence = seq_rec.seq
            strand = 0
            cds = ""
            cds_list = []
            cds_counter = 0
        
        
            for entry in feature:
                # start of cds on the + strand
                # end of cds on the - strand
                if entry.type == "CDS" and entry.location.strand == 1:
                    cds_list.append((entry.location.start, entry.location.end))
                    
                elif entry.type == "CDS" and entry.location.strand == -1:
                    cds_list.append((entry.location.end, entry.location.start))

            print(cds_list)
            
            if cds_list != []:
            
                for entry in feature:
                    
                    if entry.type == "exon":
                        
                        location = entry.location
                        strand = location.strand
                        
                        if (location.end == cds_list[cds_counter][0] and strand == -1) or (location.start == cds_list[cds_counter][0] and strand == 1):
                            # current exon is start of new cds
                            if cds != "":
                                # not the first cds
                                strand = location.strand 
                                
                                # if cds is on the "-" strand, the reverse_complement is the cds
                                if strand == 1:
                                    cds = Seq(cds)
                                else:
                                    cds = Seq(cds).reverse_complement() 
                                
                                # Create new SeqRecord for old cds
                                rec = SeqRecord(cds, id=f"{id}|:{cds_list[cds_counter-1][0]}-{cds_list[cds_counter-1][1]}", description=f"{seq_rec.description}")
                                seq_record_list.append(rec)
                            
                            # new cds
                            cds = str(sequence[location.start:location.end]) 
                            cds_counter = min(cds_counter + 1, len(cds_list) - 1)
                        else:
                            # current exon is part of current cds
                            if strand == 1:
                                cds = "".join((cds, str(sequence[location.start:location.end])))
                            else:
                                cds = "".join((str(sequence[location.start:location.end]), cds)) 
                # last cds is added to the list of SeqRecords
                if strand == -1:
                    cds = Seq(cds).reverse_complement()
                else:
                    cds = Seq(cds) 
                rec = SeqRecord(cds, id=f"{id}|:{cds_list[cds_counter][0]}-{cds_list[cds_counter][1]}", description=f"{seq_rec.description}")
                seq_record_list.append(rec)
        else:
            seq_record_list.append(seq_rec)

    return seq_record_list


if __name__ == "__main__":

    parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-gff','--output_file_gff', help="filename/path for gff output file")
    parser.add_argument('-fasta','--output_file_fasta', help="filename/path for fasta output file")
    parser.add_argument('-id', '--id_file', help="csv file with 3 columns: id,[start],[end]")
    args = parser.parse_args()


    gff_file = args.output_file_gff
    fasta_file = args.output_file_fasta
    id_file = args.id_file

    Path(os.path.dirname(gff_file)).mkdir(parents=True, exist_ok=True)
    out_handle_gff = open(gff_file, "w")

    seq_records_final = get_seq_record_gff(get_id_dict(id_file))
    GFF.write(seq_records_final, out_handle_gff)
    #print(seq_records_final)
    SeqIO.write(get_seq_record_fasta(seq_records_final), fasta_file, "fasta")
    
