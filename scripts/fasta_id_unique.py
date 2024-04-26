import argparse
from Bio import SeqIO

def add_unique_number_to_ids(input_file, output_file):
    with open(input_file, "r") as fasta_file, open(output_file, "w") as output:
        i = 1  # Start index from 1
        for record in SeqIO.parse(fasta_file, "fasta"):
            modified_id = f"{record.id}_{i}"
            modified_record = record.__class__(id=modified_id, description="", seq=record.seq)
            SeqIO.write(modified_record, output, "fasta")
            i += 1

def main():
    parser = argparse.ArgumentParser(description='Add unique numbers to FASTA IDs.')
    parser.add_argument('input_file', type=str, help='Input FASTA file')
    parser.add_argument('output_file', type=str, help='Output file with modified IDs')
    args = parser.parse_args()

    add_unique_number_to_ids(args.input_file, args.output_file)
    print(f"Modified FASTA file with unique IDs written to {args.output_file}")

if __name__ == "__main__":
    main()
