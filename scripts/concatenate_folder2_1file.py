import os
import argparse as ap


def concatenate_files(folder_path, output):
    output_file_path = output
    
    with open(output_file_path, 'w') as output_file:
        for filename in os.listdir(folder_path):
            if os.path.isdir(os.path.join(folder_path, filename)):
                continue
            with open(os.path.join(folder_path, filename), 'r') as input_file:
                for line in input_file:
                    output_file.write(line)


if __name__ == "__main__":

    parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-folder','--folder', help="folder of gff files")
    parser.add_argument('-out','--out', help="outputfilename/folder")
    args = parser.parse_args()

    folder = args.folder
    output = args.out
    concatenate_files(folder, output)
