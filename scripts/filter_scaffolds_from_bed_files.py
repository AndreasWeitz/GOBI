import argparse, os

SCAFFHOLD_ORGANISM = [["NW_022170249.1","Photinus_pyralis"],
['CM070075.1',"Pyrocoelia_pectoralis"],
['CM069432.1',"Aquatica_leii"],
['2',"Agriotes_lineatus"],
['3',"Rhagonycha_fulva"],
['KZ625248.1',"Agrilus_planipennis"],
['3',"Dascillus_cervinus"],
['2',"Nicrophorus_investigator"],
['2',"Cetonia_aurata"],
['2',"Coccinella_septempunctata"],
['10',"Chrysolina_oricalcia"],
['4',"Anthonomus_grandis"],
['6',"Brachypterus_glaber"],
['3',"Malachius_bipustulatus"],
['NC_007422.5',"Tribolium_castaneum"],
['6',"Carabus_problematicus"],
['NC_058339.1',"Chrysoperla_carnea"],
['11',"Limnephilus_lunatus"],
['NC_083541.1',"Danaus_plexippus"],
['NT_033777.3',"Drosophila_melanogaster"],
['NC_037644.1',"Apis_mellifera"],
['NC_060119.1',"Schistocerca_americana"],
['NW_019091213.1',"Folsomia_candida"]]


def filter_bed(bed_file, scaffhold):
    with open(bed_file, "r") as r:
        lines = r.readlines()

    filtered_lines = [line for line in lines if scaffhold == line.rstrip().split("\t")[0]]

    with open(bed_file, "w") as w:
        w.writelines(filtered_lines)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-b', '--bed_file_folder', help="the folder where the bed files are", type=str)
    args = parser.parse_args()

    bed_file_folder = args.bed_file_folder

    for gff_file_name in os.listdir(bed_file_folder):
        if ".bed" in gff_file_name:
            gff_file = os.path.join(bed_file_folder, gff_file_name)
            for i in range(len(SCAFFHOLD_ORGANISM)):
                 if SCAFFHOLD_ORGANISM[i][1] in gff_file_name:
                      filter_bed(gff_file, SCAFFHOLD_ORGANISM[i][0])
