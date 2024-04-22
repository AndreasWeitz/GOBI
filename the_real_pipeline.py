#!/usr/bin/env python3

import subprocess
import os
import sys
import argparse
from select import select
from playsound import playsound

IS_CHR_SPLIT = {
    1: False,
    2: True,
    3: True,
    4: False,
    5: False,
    6: False,
    7: False,
    8: False,
    9: False,
    10: False,
    11: False,
    12: False,
    13: False,
    14: False,
    15: False,
    16: True,
    17: True,
    18: False,
    19: False,
    20: False,
    21: True,
    22: False
}

GENOME_GFFS = {
    1: '../data/iterations/iteration_1/pyrocoelia_pectoralis_genome/data/GCA_036250285.1/genomic.gff',
    2: '../data/iterations/iteration_2/aquatica_leii_genome/data/GCA_035610365.1/genomic.gff',
    3: '../data/iterations/iteration_3/agriotes_lineatus_genome/Agriotes_lineatus-GCA_940337035.1-2022_08-genes.gff3',
    4: '../data/iterations/iteration_4/rhagonycha_fulva_genome/ncbi_dataset/data/GCA_905340355.1/genomic.gff',
    5: '../data/iterations/iteration_5/agrilus_planipennis_genome/Agrilus_planipennis-GCA_000699045.2-2021_12-genes.gff3',
    6: '../data/iterations/iteration_6/dascillus_cervinus_genome/ncbi_dataset/data/GCA_949768715.1/Dascillus_cervinus-GCA_949768715.1-2023_07-genes.gff3',
    7: '../data/iterations/iteration_7/nicrophorus_investigator_genome/Nicrophorus_investigator-GCA_963457615.1-2023_10-genes.gff3',
    8: '../data/iterations/iteration_8/cetonia_aurata_genome/Cetonia_aurata-GCA_949128085.1-2023_07-genes.gff3',
    9: '../data/iterations/iteration_9/coccinella_septempunctata_genome/Coccinella_septempunctata-GCA_907165205.1-2021_12-genes.gff3',
    10: '../data/iterations/iteration_10/chrysolina_oricalcia_genome/Chrysolina_oricalcia-GCA_944452925.1-2022_08-genes.gff3',
    11: '../data/iterations/iteration_11/anthonomus_grandis_grandis_genome/Anthonomus_grandis_grandis-GCA_022605725.3-2022_09-genes.gff3',
    12: '../data/iterations/iteration_12/brachypterus_glaber_genome/Brachypterus_glaber-GCA_958510825.1-2023_10-genes.gff3',
    13: '../data/iterations/iteration_13/malachius_bipustulatus_genome/Malachius_bipustulatus-GCA_910589415.1-2022_02-genes.gff3',
    14: '../data/iterations/iteration_14/tribolium_castaneum_genome/ncbi_dataset/data/GCF_000002335.3/genomic.gff',
    15: '../data/iterations/iteration_15/carabus_problematicus_genome/Carabus_problematicus-GCA_963422195.1-2023_10-genes.gff3',
    16: '../data/iterations/iteration_16/chrysoperla_carnea_genome/ncbi_dataset/data/GCF_905475395.1/genomic.gff',
    17: '../data/iterations/iteration_17/limnephilus_lunatus_genome/Limnephilus_lunatus-GCA_917563855.2-2022_12-genes.gff3',
    18: '../data/iterations/iteration_18/danaus_plexippus_genome/ncbi_dataset/data/GCF_018135715.1/genomic.gff',
    19: '../data/iterations/iteration_19/drosophila_melanogaster_genome/ncbi_dataset/data/GCF_000001215.4/genomic.gff',
    20: '../data/iterations/iteration_20/apis_mellifera_genome/ncbi_dataset/data/GCF_003254395.2/genomic.gff',
    21: '../data/iterations/iteration_21/schistocerca_americana_genome/ncbi_dataset/data/GCF_021461395.2/genomic.gff',
    22: '../data/iterations/iteration_22/folsomia_candida_genome/ncbi_dataset/data/GCF_002217175.1/genomic.gff'
}

SCAFFOLDS = {
    1: 'CM070075.1',
    2: 'CM069432.1',
    3: '2',
    4: '3',
    5: 'KZ625248.1',
    6: '3',
    7: '2',
    8: '2',
    9: '2',
    10: '10',
    11: '4',
    12: '6',
    13: '3',
    14: 'NC_007422.5',
    15: '6',
    16: 'NC_058339.1',
    17: '11',
    18: 'NC_083541.1',
    19: 'NT_033777.3',
    20: 'NC_037644.1',
    21: 'NC_060119.1',
    22: 'NW_019091213.1'
}

ENSEMBL_NAME = {
    1: '',
    2: '',
    3: 'Agriotes_lineatus-GCA_940337035.1-2022_08-cds.fa',
    4: 'Rhagonycha_fulva-GCA_905340355.1-2021_12-cds.fa',
    5: 'Agrilus_planipennis-GCA_000699045.2-2021_12-cds.fa',
    6: 'Dascillus_cervinus-GCA_949768715.1-2023_07-cds.fa',
    7: 'Nicrophorus_investigator-GCA_963457615.1-2023_10-cds.fa',
    8: 'Cetonia_aurata-GCA_949128085.1-2023_07-cds.fa',
    9: 'Coccinella_septempunctata-GCA_907165205.1-2021_12-cds.fa',
    10: 'Chrysolina_oricalcia-GCA_944452925.1-2022_08-cds.fa',
    11: 'Anthonomus_grandis_grandis-GCA_022605725.3-2022_09-cds.fa',
    12: 'Brachypterus_glaber-GCA_958510825.1-2023_10-cds.fa',
    13: 'Malachius_bipustulatus-GCA_910589415.1-2022_02-cds.fa',
    14: '',
    15: 'Carabus_problematicus-GCA_963422195.1-2023_10-cds.fa',
    16: '',
    17: 'Limnephilus_lunatus-GCA_917563855.2-2022_12-cds.fa',
    18: '',
    19: '',
    20: '',
    21: '',
    22: ''
}

def fastblast(db_path, 
              query_path, 
              blast_out, 
              evalue, 
              blast_type, 
              inputdb_fasta, 
              inputdb_gff, 
              blastout_gff, 
              id_file_db = None, 
              altout_gff = False): 
    
    temp = ""
    if id_file_db != None:
        temp = f"-id {id_file_db}"
        
    if altout_gff:
        temp = f"{temp} -ftype"
        blast_out = f"{blast_out}_alt"
        
    options = f"{temp} -e {evalue} -btype {blast_type} -fasta {inputdb_fasta} -gff {inputdb_gff} -bgff {blastout_gff}"
    
    subprocess.run(f"python3 pipeline_part1.py -dbpath {db_path} -qp {query_path} -out {blast_out} {options}", shell=True, check=True) 

def play_audio():
    audio_path = '../data/others/alarm_sound.mp3'
    if os.path.exists(audio_path):
        playsound(audio_path)
    else:
        pass


def create_blast_db(i: int, base_dir: str):
    id_name = "ppyr_ids"
    id_file = f"../data/{id_name}.csv"
    blast_type = 'tblastx'

    db_path = f"{base_dir}/database/{id_name}/db/{id_name}.db"
    inputdb_fasta = f"{base_dir}/database/{id_name}/db_fasta/{id_name}.fasta"
    inputdb_gff = f"{base_dir}/database/{id_name}/db_gff/{id_name}_db.gff"

    blast_out = f"{base_dir}/blast_output/{id_name}/{blast_type}_{id_name}"
    blastout_gff = f"{base_dir}/blast_output/{id_name}/{blast_type}_{id_name}.gff"

    if i > 1:
        matching_cds_gff = os.path.join(base_dir, f'../iteration_{i-1}/cds_results/matching_cds_filtered.gff')
        cmd = f'python gff_reader.py -gff {matching_cds_gff} -id_file {id_file}'
        ensembl_name = ENSEMBL_NAME[i-1]
        if ensembl_name:
            cmd += f' -ensemble {ensembl_name}'
        subprocess.run(cmd, shell=True, check=True)

    fastblast(
    db_path=db_path, 
    query_path='none', 
    blast_out=blast_out, 
    altout_gff=True, 
    id_file_db=id_file, 
    evalue=0, 
    blast_type=blast_type, 
    inputdb_fasta=inputdb_fasta, 
    inputdb_gff=inputdb_gff, 
    blastout_gff=blastout_gff
    )

def run_blast_search(i: int, base_dir: str):
    query_dir = os.path.join(base_dir, 'query')
    db_path = os.path.join(base_dir, 'database/ppyr_ids/db/ppyr_ids.db')
    out_dir = os.path.join(base_dir, 'multiblast_raw' if IS_CHR_SPLIT[i] else 'blast_raw')
    subprocess.run(
        f'./multi_blast.py -q {query_dir} -db {db_path} -o {out_dir}', 
        shell=True, 
        check=True
    )

def pipeline_after_blast(i: int, base_dir: str):
    genome_gff = GENOME_GFFS[i]
    cds_hits = f"{base_dir}/cds_results/matching_cds.gff"
    hits_without_map_on_cds = f"{base_dir}/cds_results/non_matching_cds.gff"
    cds_hits_filtered = f"{base_dir}/cds_results/matching_cds_filtered.gff"
    scaffold = SCAFFOLDS[i]

    if IS_CHR_SPLIT[i]:
        blast_output = f"{base_dir}/multiblast_raw"
        gffs_folder = f"{base_dir}/multiblast_gff"
        corrected_gffs_folder = f"{base_dir}/correct_ranges"
        concatenated_gff = f"{base_dir}/iteration_{i}_cat.gff"

        scripts = [
        f"multi_blast2gff.py -i {blast_output} -o {gffs_folder}",
        f"correct_ranges_of_blasthits_after_fastasplitter.py -gff {gffs_folder} -out {corrected_gffs_folder}",
        f"concatenate_folder2_1file.py -folder {corrected_gffs_folder} -out {concatenated_gff}",
        f"blast_hits2cds.py -bh {concatenated_gff} -g {genome_gff} -o {cds_hits} -s {scaffold} -n {hits_without_map_on_cds}",
        f"filter_matching_cds.py -inp {cds_hits} -out {cds_hits_filtered}"
        ]
    else:
        blast_output = f"{base_dir}/blast_raw/*"
        gff_folder = f"{base_dir}/blast_gff"
        
        if not os.path.exists(gff_folder):
            os.mkdir(gff_folder)

        hits_gff_file = os.path.join(gff_folder, "blast_hits.gff")

        scripts = [
        f"blast2gff.py -Q -b {blast_output} > {hits_gff_file}",
        f"blast_hits2cds.py -bh {hits_gff_file} -g {genome_gff} -o {cds_hits} -s {scaffold} -n {hits_without_map_on_cds}",
        f"filter_matching_cds.py -inp {cds_hits} -out {cds_hits_filtered}"
        ]

    for script in scripts:
        command = f"python {script}"
        subprocess.run(command, shell=True, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-c', '--checkpoint', help='start of next iteration', required=False, type=int)

    args = parser.parse_args()
    checkpoint = args.checkpoint

    if not checkpoint:
        checkpoint = 1

    for i in range(checkpoint, len(GENOME_GFFS)+1):
        print(f"{'='*25}\nStarting iteration {i}")
        base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"../data/iterations/iteration_{i}")
        
        if i > 1:
            print("Starting blast db build")
            create_blast_db(i, base_dir)
            print("Created blast db")

        print("Starting blast search")
        run_blast_search(i, base_dir)
        print("Completed blast search")

        print("Starting CDS extraction")
        pipeline_after_blast(i, base_dir)
        print(f"Extracted CDS\nIteration {i} done\n{'='*25}")

        play_audio()
        print("Stop process?")
        timeout = 90
        rlist, wlist, xlist = select([sys.stdin], [], [], timeout)

        if rlist:
            print(f"Stopped after iteration {i}")
            exit()
