from BCBio import GFF
import pandas as pd
import argparse as ap
import re

# in_file = "blast_output\exon_ids_luc\\tblastx_exon_ids_luc.gff"
# id_file = id_files/exon_ids_luc.csv

if __name__=="__main__":

    parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter, description=__doc__)

    
    parser.add_argument('-gff','--gff_input_file', help="")
    parser.add_argument('-id_file','--csv_id_file', help="")
    parser.add_argument('-ensemble', help="set parameter for ensemble files")
    args = parser.parse_args()

    gff_file = args.gff_input_file
    id_file = args.csv_id_file

    in_handle = open(gff_file)
    loc_list = []
    id_list = set()

    if not args.ensemble:
        for rec in GFF.parse(in_handle):
            for elem in rec.features:
                loc_list.append([rec.id, elem.location.start, elem.location.end])
        in_handle.close()
    else:
        new_id = args.ensemble.replace(".\\", "")
        with open(gff_file) as file:
            for line in file:
                desc = line.split("\t")[8]
                print(desc)
                match = re.search(r'version=(\d+)', desc)
                if match:
                    version = str(match.group(1))
                cds = desc.split(";")[1].split("=")[1].split(":")
                print(cds)
                if cds[0] == "transcript":
                    id = f'{cds[1]}.{version}' if match else cds[1]
                    if id not in id_list:
                        new_id += f"|:{id}"
                        loc_list.append([new_id, None, None])
                        id_list.add(id)
                    new_id = args.ensemble.replace(".\\", "")
    df = pd.DataFrame(loc_list, columns=['ID', 'start', 'end'])
    print(df)
    df.to_csv(id_file, mode='a', index=False, header=False)

    print("Done")