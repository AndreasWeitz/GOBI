from gff2bed import convert, parse
import pandas as pd
import argparse
import os

def get_cds_gff(gff, out):
    with open(gff, "r") as r:
        lines = r.readlines()

        filtered_lines = []
        for line in lines:
            if "ID=cds-" in line:
                filtered_lines.append(line)

    with open(f"{out}", "w") as w:
        w.writelines(filtered_lines)

def gff_2_bed_function(gff, out_folder, NCBI_info):
    is_NCBI = False
    if "NCBI" in NCBI_info:
        is_NCBI = True
        
    if not is_NCBI:    
        bed_df = pd.DataFrame(convert(parse(gff)))
        bed_df.to_csv(f"{out_folder}/{'_'.join(gff.split('/')[-2:])}.bed", sep='\t', index=False, header=False)
    else:
        get_cds_gff(gff, gff)
        bed = []
        with open(gff, 'r') as reader:
            i = 0
            for line in reader:
                if not line.startswith("#") and ("ID=cds-" in line or "ID=cds:" in line):
                    if "ID=cds-" in line:
                        bed.append([line.split("\t")[0],
                                    line.split("\t")[3],
                                    line.split("\t")[4],
                                    line.split("\t")[6],
                                    line.split("\t")[7],
                                    line.split("\t")[8].split("ID=cds-")[1].split(";")[0]])
                    elif "ID=cds:" in line:
                        bed.append([line.split("\t")[0],
                                    line.split("\t")[3],
                                    line.split("\t")[4],
                                    line.split("\t")[6],
                                    line.split("\t")[7],
                                    line.split("\t")[8].split("ID=cds:")[1].split(";")[0]])
                    i += 1
            bed_new = []
            bed_new.append(bed[0])
            for i in range(len(bed)):
                j = len(bed_new)
                if bed[i][5] != bed_new[j-1][5]:
                    bed_new.append(bed[i])
            bed = bed_new


        with open(f"{out_folder}/{'_'.join(gff.split('/')[-2:])}.bed", 'w') as writer:
            for i in range(len(bed)):
                #if "NW" in bed[i][0] and "grilus" in gff:
                #    bed[i][0] = bed[i][0].replace('NW', 'KZ')
                writer.write(f"{bed[i][0]}\t{bed[i][1]}\t{bed[i][2]}\t{bed[i][5]}\t{bed[i][3]}\t{bed[i][4]}\n")
        


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    #parser.add_argument('-gff_folder', '--gff_file_folder', help="input gff file folder", type=str)
    parser.add_argument('-path_file', '--folder_file', help="file with paths to the gf files", type=str)
    parser.add_argument('-out', '--out_folder', help="path to safe the bed files", type=str)
    args = parser.parse_args()

    #gff_file_folder = args.gff_file_folder
    folder_file = args.folder_file
    out_folder = args.out_folder

    with open(folder_file, 'r') as file:
            for line in file:
                gff_2_bed_function(line.strip().split(" ")[0], out_folder, line.strip().split(" ")[1])

    #for gff_file_name in os.listdir(gff_file_folder):
    #    if ".gff" in gff_file_name:
    #        gff_file = os.path.join(gff_file_folder, gff_file_name)
    #        gff_2_bed_function(gff_file)
