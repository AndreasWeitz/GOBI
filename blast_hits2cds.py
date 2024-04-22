#!/usr/bin/env python3

import subprocess
import re
import os
import argparse
import sys

def preprocess_blast_hits_gff(blast_hits_gff: str, scaffold: str) -> str:
    '''
    Helper function for `hits2cds()`. Creates a (temporary) 
    gff file where the ids are exchanged for the id of the 
    scaffold the hit occured on.

    ## Parameters:
    `blast_hits_gff`: Gff file containing blast hits 
    as created by `blast2gff.py`

    ## Returns:
    `tmp_blast_hits_gff`: Path to the newly created gff. 
    It is advised to delete it after usage with:
    ```
    >>> import os
    >>> blast_hits_gff = 'path/to/file.gff'
    >>> tmp_blast_hits_gff = preprocess_blast_hits_gff(blast_hits_gff)
    >>> os.remove(tmp_blast_hits_gff)
    ```
    '''
    temp_file = os.path.join(os.path.dirname(blast_hits_gff), f'tmp_{os.path.basename(blast_hits_gff)}')
    with open(blast_hits_gff) as in_file:
        with open(temp_file, 'w') as out_file:
            for line in in_file:
                if not line.startswith('#'):
                    entries = line.split('\t')
                    entries[0] = scaffold
                    out_file.write("\t".join(entries))

    return temp_file

def filter_features(gff_file: str, features: list) -> str:
    '''
    Helper function for `hits2cds()`. Creates gff file 
    from `gff_file` only containing lines where the feature 
    is in `features`.

    ## Parameters:
    `gff_file`: Path to gff file, most commonly genomic gff
    `features`: List of features to filter for, like `['CDS']`

    ## Returns:
    `feature_gff_file`: Path to newly created gff file
    '''
    feature_gff_file = os.path.join(os.path.dirname(gff_file), f"{'_'.join(features)}_{os.path.basename(gff_file)}")

    with open(gff_file) as in_file:
        with open(feature_gff_file, 'w') as out_file:    
            for line in in_file:
                if not line.startswith('#') and line.split('\t')[2] in features:
                    out_file.write(line)

    return feature_gff_file

def hits2cds(
        blast_hits_gff: str,
        genome_gff: str, 
        scaffold: str,
        out_file=None, 
        no_match_out_file=None,
        overlap=0.5
    ) -> str:
    '''
    Maps blast hits to CDS. Requires a genomic gff file.

    ## Parameters:
    `blast_hits_gff`: Gff file containing blast hits 
    as created by `blast2gff.py`
    `genome_gff`: Gff file of genome
    `out_file` (optional): Path to gff file to write matched CDS to
    `no_match_out_file` (optional): Path to gff file to write 
    blast hits with no CDS matches to. If there are no matchless 
    blast hits, no file is created
    `overlap`: Required percentage overlap of blast hits

    ## Returns:
    `gff_str`: String formatted as gff containing genomic CDS blast hits could be mapped to
    '''
    tmp_blast_hits_gff = preprocess_blast_hits_gff(blast_hits_gff, scaffold)
    cds_genome_gff = filter_features(genome_gff, ['CDS'])
    cmd = f'bedtools intersect -u -F {overlap} -s -a {cds_genome_gff} -b {tmp_blast_hits_gff}'
    
    if no_match_out_file:
        cmd += f' && bedtools intersect -v -a {tmp_blast_hits_gff} -b {cds_genome_gff} > {no_match_out_file}'

    result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    gff_str = result.stdout.decode('utf-8')

    if no_match_out_file and os.stat(no_match_out_file).st_size == 0:
        print('No unannotated potential CDS discovered.', file=sys.stderr)
        os.remove(no_match_out_file)

    os.remove(tmp_blast_hits_gff)
    os.remove(cds_genome_gff)

    if out_file:
        with open(out_file, 'w') as f:
            f.write(gff_str)

    return gff_str

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-bh', '--blast_hits_gff', help='GFF file containing blast hits', required=True, type=str)
    parser.add_argument('-g', '--genome_gff', help='GFF file containing CDS information', required=True, type=str)
    parser.add_argument('-s', '--scaffold', help='Name of scaffold, for example CM070075.1', required=True, type=str)
    parser.add_argument('-o', '--outfile', help='GFF file to write CDS to. If None, print to stdout', required=False, type=str)
    parser.add_argument('-n', '--no_match_outfile', help='GFF file to write hits to that do not map to CDS', required=False, type=str)
    parser.add_argument('-ol', '--overlap', help='', required=False, type=float, default=0.5)

    args = parser.parse_args()
    blast_hits_gff = args.blast_hits_gff
    genome_gff = args.genome_gff
    scaffold = args.scaffold
    out_file = args.outfile
    no_match_out_file = args.no_match_outfile
    overlap = args.overlap

    
    if out_file:
        out_file_dir = os.path.dirname(out_file)
        if not os.path.exists(out_file_dir):
            os.mkdir(out_file_dir)
    
    if no_match_out_file:
        no_match_out_file_dir = os.path.dirname(no_match_out_file)
        if not os.path.exists(no_match_out_file_dir):
            os.mkdir(no_match_out_file_dir)

    gff_str = hits2cds(
        blast_hits_gff=blast_hits_gff,
        genome_gff=genome_gff,
        scaffold=scaffold,
        out_file=out_file,
        no_match_out_file=no_match_out_file,
        overlap=overlap
    )
    
    if not out_file:
        print(gff_str)
        