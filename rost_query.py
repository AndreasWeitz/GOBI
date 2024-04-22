#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
import re

GROUPING = ['[RHK]', '[DE]', '[STNQ]', '[CUGP]', '[AILMFWYV]']

def get_terminus(seq: str, terminus_index: int, group_similar=False) -> str:
    seq = seq.strip('*')
    terminus =  seq[:terminus_index] if terminus_index > 0 else seq[terminus_index:]

    if group_similar:
        terminus_grouped = ''

        for aa in terminus:
            for group in GROUPING:
                if aa in group:
                    terminus_grouped += group

        return terminus_grouped
    else:
        return terminus
    
def contains_subseq(seq: str, sub_seq: str, group_similar=False) -> bool:
    if group_similar:
        grouped_re = r''
        for sub_seq_aa in sub_seq:
            for grouping in GROUPING:
                if sub_seq_aa in grouping:
                    grouped_re += grouping
                    break
        match = re.search(grouped_re, seq)
        return match is not None
    else:
        return sub_seq in seq

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-f', '--fasta', help="fasta file to extract from", type=str, required=True)
    parser.add_argument('-t', '--terminus', help="i.e. 3 for first three letters, -3 for last three", type=int, required=False)
    parser.add_argument('-c', '--contains', help="query whether seq contains substring", type=str, required=False)
    parser.add_argument('-g', '--group_similar', help="group amino acids together based on chemical features", action='store_true', required=False)
    parser.add_argument('-n', '--colname', help="custom colname; however useful colnames are autogenerated", type=str, required=False)
    parser.add_argument('-csv','--csv_file', help="csv file to append col to", type=str, required=True)
    
    args = parser.parse_args()
    fasta = args.fasta
    terminus = args.terminus
    contains = args.contains
    group_similar = args.group_similar
    col_name = args.colname
    csv_file = args.csv_file

    if not col_name and terminus:
        col_name = 'n_term' if terminus > 0 else 'c_term'
    elif not col_name and contains:
        col_name = f"contains_{contains}{'_like' if group_similar else ''}"

    df = pd.read_csv(csv_file)
    
    fasta_records = SeqIO.parse(fasta, 'fasta')
    fasta_record_list = []

    for fasta_record in fasta_records:
        fasta_record.id = fasta_record.id.replace('.', '_')
        fasta_record_list.append((str(fasta_record.id), str(fasta_record.seq)))

    results = []

    for df_row in df.iterrows():
        for (fasta_record_id, fasta_record_seq) in fasta_record_list:
            if fasta_record_id == df_row[1].fasta_header:
                if terminus:
                    result = get_terminus(fasta_record_seq, terminus, group_similar)
                elif contains:
                    result = contains_subseq(fasta_record_seq, contains, group_similar)
                results.append(result)
                fasta_record_list.remove((fasta_record_id, fasta_record_seq))
                break

    df[col_name] = pd.Series(results)
    df.to_csv(csv_file, sep=',', header=True, index=False)
    