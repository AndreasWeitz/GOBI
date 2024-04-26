import argparse


def filter_gff(gff, out):
    with open(gff, "r") as r:
        lines = r.readlines()

        filtered_lines = []
        for line in lines:
            if "ID=cds-" in line:
                filtered_lines.append(line)

    with open(f"{out}/{'_'.join(gff.split('/')[-2:])}", "w") as w:
        w.writelines(filtered_lines)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-path_file', '--folder_file', help="input gff file folder", type=str)
    parser.add_argument('-out', '--out_folder', help="out folder", type=str)
    args = parser.parse_args()

    folder_file = args.folder_file
    out_folder = args.out_folder

    with open(folder_file, 'r') as file:
            for line in file:
                filter_gff(line.strip(), out_folder)