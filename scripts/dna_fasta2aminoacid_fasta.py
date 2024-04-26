import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def translate_dna_to_protein(input_file, output_file):
    with open(input_file, "r") as fasta_file, open(output_file, "w") as protein_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            dna_sequence = record.seq
            protein_sequence = dna_sequence.translate()
            protein_file.write(f">{record.id}\n")
            protein_file.write(f"{protein_sequence}\n")

def main():
    parser = argparse.ArgumentParser(description='Translate DNA sequences to protein sequences.')
    parser.add_argument('input_file', type=str, help='Input FASTA file containing DNA sequences')
    parser.add_argument('output_file', type=str, help='Output file to write translated protein sequences')
    args = parser.parse_args()

    translate_dna_to_protein(args.input_file, args.output_file)
    print(f"Translation complete. Protein sequences written to {args.output_file}")

if __name__ == "__main__":
    main()
