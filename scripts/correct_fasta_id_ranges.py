#!/usr/bin/env python3

import argparse

GENOME_GFFS = {
    'NW_022170249.1': '../data/photinus_pyralis_genome/data/GCF_008802855.1/genomic.gff',
    'CM070075.1': '../data/iterations/iteration_1/pyrocoelia_pectoralis_genome/data/GCA_036250285.1/genomic.gff',
    'CM069432.1': '../data/iterations/iteration_2/aquatica_leii_genome/data/GCA_035610365.1/genomic.gff',
    'ENSLAO': '../data/iterations/iteration_3/agriotes_lineatus_genome/Agriotes_lineatus-GCA_940337035.1-2022_08-genes.gff3',
    'ENSRFU': '../data/iterations/iteration_4/rhagonycha_fulva_genome/ncbi_dataset/data/GCA_905340355.1/genomic.gff',
    'KZ625248.1': '../data/iterations/iteration_5/agrilus_planipennis_genome/Agrilus_planipennis-GCA_000699045.2-2021_12-genes.gff3',
    'ENSIJZ': '../data/iterations/iteration_6/dascillus_cervinus_genome/ncbi_dataset/data/GCA_949768715.1/Dascillus_cervinus-GCA_949768715.1-2023_07-genes.gff3',
    'ENSMWS': '../data/iterations/iteration_7/nicrophorus_investigator_genome/Nicrophorus_investigator-GCA_963457615.1-2023_10-genes.gff3',
    'ENSUAN': '../data/iterations/iteration_8/cetonia_aurata_genome/Cetonia_aurata-GCA_949128085.1-2023_07-genes.gff3',
    'ENSSIN': '../data/iterations/iteration_9/coccinella_septempunctata_genome/Coccinella_septempunctata-GCA_907165205.1-2021_12-genes.gff3',
    'BRAKEROCH': '../data/iterations/iteration_10/chrysolina_oricalcia_genome/Chrysolina_oricalcia-GCA_944452925.1-2022_08-genes.gff3',
    'ENSAGS': '../data/iterations/iteration_11/anthonomus_grandis_grandis_genome/Anthonomus_grandis_grandis-GCA_022605725.3-2022_09-genes.gff3',
    'BRAKERQTG': '../data/iterations/iteration_12/brachypterus_glaber_genome/Brachypterus_glaber-GCA_958510825.1-2023_10-genes.gff3',
    'BRAKERMBP': '../data/iterations/iteration_13/malachius_bipustulatus_genome/Malachius_bipustulatus-GCA_910589415.1-2022_02-genes.gff3',
    'NC_007422.5': '../data/iterations/iteration_14/tribolium_castaneum_genome/ncbi_dataset/data/GCF_000002335.3/genomic.gff',
    'ENSCZH': '../data/iterations/iteration_15/carabus_problematicus_genome/Carabus_problematicus-GCA_963422195.1-2023_10-genes.gff3',
    'NC_058339.1': '../data/iterations/iteration_16/chrysoperla_carnea_genome/ncbi_dataset/data/GCF_905475395.1/genomic.gff',
    'ENSLLS': '../data/iterations/iteration_17/limnephilus_lunatus_genome/Limnephilus_lunatus-GCA_917563855.2-2022_12-genes.gff3',
    'NC_083541.1': '../data/iterations/iteration_18/danaus_plexippus_genome/ncbi_dataset/data/GCF_018135715.1/genomic.gff',
    'NT_033777.3': '../data/iterations/iteration_19/drosophila_melanogaster_genome/ncbi_dataset/data/GCF_000001215.4/genomic.gff',
    'NC_037644.1': '../data/iterations/iteration_20/apis_mellifera_genome/ncbi_dataset/data/GCF_003254395.2/genomic.gff',
    'NC_060119.1': '../data/iterations/iteration_21/schistocerca_americana_genome/ncbi_dataset/data/GCF_021461395.2/genomic.gff',
    'NW_019091213.1': '../data/iterations/iteration_22/folsomia_candida_genome/ncbi_dataset/data/GCF_002217175.1/genomic.gff'
}

def get_containing_gene_entry(scaffold: str, start: int, end: int) -> str:
    genome_gff = GENOME_GFFS[scaffold]
    with open(genome_gff) as file:
        for line in file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if (fields[0] == scaffold 
                    and fields[2] == 'gene'
                    and start >= int(fields[3])-1
                    and end <= int(fields[4])):
                    return line

def correct_fasta_id_ranges(fasta: str):
    with open(fasta) as in_file:
        for line in in_file:
            if (line.startswith('>') 
                and not line.startswith('>ENS')
                and not line.startswith('>BRAKER')):
                id = line.split(' ')[0]
                rest = ' '.join(line.split(' ')[1:])
                id_parts = id.split('-')
                scaffold = id_parts[0].strip('>')
                start = int(float(id_parts[1]))
                end = int(float(id_parts[2]))

                containing_gene = get_containing_gene_entry(scaffold, start, end)
                fields = containing_gene.strip().split('\t')
                new_start = int(fields[3])
                new_end = int(fields[4])
                new_header = f">{scaffold}-{new_start}-{new_end} {rest}"
                print(new_header.strip())
            else:
                print(line.strip())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-f', '--fasta', help='Fasta containing ids of form scaffold-start-end', required=True, type=str)

    args = parser.parse_args()
    fasta = args.fasta

    correct_fasta_id_ranges(fasta)
