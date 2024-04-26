from Bio import SeqIO
import argparse, os

SCAFFHOLD_ORGANISM = [["NW_022170249.1","Photinus_pyralis", "XXXXXXXX"],
['CM070075.1',"Pyrocoelia_pectoralis", "XXXXXXXX"],
['CM069432.1',"Aquatica_leii", "XXXXXXXX"],
['2',"Agriotes_lineatus", "ENSLAO"],
['3',"Rhagonycha_fulva", "ENSRFU"],
['KZ625248.1',"Agrilus_planipennis", "XXXXXXXX"],
['3',"Dascillus_cervinus", "ENSIJZ"],
['2',"Nicrophorus_investigator", "ENSMWS"],
['2',"Cetonia_aurata", "ENSUAN"],
['2',"Coccinella_septempunctata", "ENSSIN"],
['10',"Chrysolina_oricalcia", "BRAKEROCH"],
['4',"Anthonomus_grandis", "ENSAGS"],
['6',"Brachypterus_glaber", "BRAKERQTG"],
['3',"Malachius_bipustulatus", "BRAKERMBP"],
['NC_"007422.5',"Tribolium_castaneum", "XXXXXXXX"],
['6',"Carabus_problematicus", "ENSCZH"],
['NC_058339.1',"Chrysoperla_carnea", "XXXXXXXX"],
['11',"Limnephilus_lunatus", "ENSLLS"],
['NC_083541.1',"Danaus_plexippus", "XXXXXXXX"],
['NT_033777.3',"Drosophila_melanogaster", "XXXXXXXX"],
['NC_037644.1',"Apis_mellifera", "XXXXXXXX"],
['NC_060119.1',"Schistocerca_americana", "XXXXXXXX"],
['NW_019091213.1',"Folsomia_candida", "XXXXXXXX"]]

def read_bed(file_path):
    data = {}
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if "gene:" in fields[3]:
                data[fields[3].split(':')[1]] = fields[0]
            else:
                data[fields[3]] = fields[0]
    return data


def filter_fasta(fasta_path, bed_path, out):
    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        if "gene:" in record.description:
            sequences[record.description.split('gene:')[1].split(' ')[0]] = str(record.seq)
        else:
            sequences[record.id] = str(record.seq)
    bed = read_bed(bed_path)

    with open(f"{out}/{'_'.join(fasta_path.split('/')[-2:])}", 'w') as writer:
        for sequence in sequences.keys():
            if "ENS" in sequence or "BRAK" in sequence:
                sequence = sequence.split('.')[0]
                if sequence in bed.keys():
                    writer.write(f'>{sequence}\n{sequences[f"{sequence}.1"]}\n')
            else:
                if sequence in bed.keys():
                    writer.write(f'>{sequence}\n{sequences[f"{sequence}"]}\n')


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-f', '--fasta_folder_file', help="file with paths to the fasta files", type=str)
    parser.add_argument('-b', '--bed_folder', help="folder with the bed files", type=str)
    parser.add_argument('-o', '--out', help="folder to safe the filtered fastas", type=str)
    args = parser.parse_args()

    fasta_folder_file = args.fasta_folder_file
    bed_folder = args.bed_folder
    out = args.out

    fasta_files = []
    with open(fasta_folder_file, 'r') as reader:
        for line in reader:
            fasta_files.append(line.rstrip())


    for bed_file_name in os.listdir(bed_folder):
        bed_organism = None
        if ".bed" in bed_file_name:
            bed_file = os.path.join(bed_folder, bed_file_name)
            for i in range(len(SCAFFHOLD_ORGANISM)):
                if SCAFFHOLD_ORGANISM[i][1] in bed_file_name:
                    bed_organism = SCAFFHOLD_ORGANISM[i][1]
        for fasta_file in fasta_files:
            for i in range(len(SCAFFHOLD_ORGANISM)):
                print("here")
                if (SCAFFHOLD_ORGANISM[i][1] in fasta_file or SCAFFHOLD_ORGANISM[i][2] in fasta_file) and SCAFFHOLD_ORGANISM[i][1] == bed_organism:
                    print("here")
                    filter_fasta(fasta_file, bed_file, out)


if __name__ == "__main__":
    main()