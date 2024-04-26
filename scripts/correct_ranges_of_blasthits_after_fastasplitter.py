from BCBio import GFF
import argparse as ap
import re
import os

ONE_MILLION = 1000000

def extract_first_number(string):
    match = re.search(r'_(\d+)B_', string)
    if match:
        return int(match.group(1))
    else:
        return None

def shift_gff_ranges(gff_file, output_file):
    with open(output_file, "w") as out_handle:
        with open(gff_file) as in_handle:
            for rec in GFF.parse(in_handle):
                for feature in rec.features:
                    feature.location = feature.location._shift(extract_first_number(gff_file)*ONE_MILLION)
                GFF.write([rec], out_handle)


if __name__ == "__main__":

    parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter, description=__doc__)

    parser.add_argument('-gff','--gff', help="folder of gff files")
    parser.add_argument('-out', '--out', help="output_folder")
    args = parser.parse_args()

    gff_folder = args.gff
    out_folder = args.out

    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    for filename in os.listdir(gff_folder):
        if filename.endswith(".gff"):
            gff_file = os.path.join(gff_folder, filename)
            output_file = os.path.join(out_folder, "corrected_"+filename)
            shift_gff_ranges(gff_file, output_file)
