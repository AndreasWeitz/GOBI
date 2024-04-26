import argparse
from Bio import SeqIO


def read_blastoutput(filename):
    filtered_entries = []
    seen_prefixes = set()

    with open(filename, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if columns:
                first_column = columns[0]
                prefix = first_column[:7]
                if prefix not in seen_prefixes:
                    seen_prefixes.add(prefix)
                    filtered_entries.append(first_column)

    return filtered_entries


def filter_fasta_with_blastoutput(words_list, input_fasta, output_fasta):
    filtered_entries = []

    # Read FASTA file and filter entries
    with open(input_fasta, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            for word in words_list:
                if word in record.id:
                    filtered_entries.append(record)
                    break  # Move to the next record once a match is found

    # Write filtered entries to new FASTA file
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(filtered_entries, output_handle, "fasta")

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-blast', '--blastoutput', help="", type=str)
    parser.add_argument('-in', '--input_fasta', help="", type=str)
    parser.add_argument('-out', '--output_fasta', help="", type=str)
    args = parser.parse_args()

    blastoutput = args.blastoutput
    input_fasta = args.input_fasta
    output_fasta = args.output_fasta

    filter_fasta_with_blastoutput(read_blastoutput(blastoutput), input_fasta, output_fasta)


if __name__ == '__main__':
    main()
    