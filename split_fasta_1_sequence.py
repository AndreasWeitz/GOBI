from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math
import argparse as ap

# -file Melanotus_villosus\chromosome2.fasta -n 20000000

def split_fasta(chr2, n, h):
    f = math.ceil(len(chr2.seq)/n)

    for i in range(f):
        if (i+1)*n + h < len(chr2.seq):
            end = (i+1)*n + h
        else:
            end = len(chr2.seq)

         #record = SeqRecord(Seq(chr2.seq[i*n + h:end]), id=f"fragment_{int((i*n + h)/teiler)}{bst}_{int(end/teiler)}{bst}|:{i*n + h}-{end}", description=chr2.description)
        record = SeqRecord(chr2.seq[i*n + h:end], id="", description=chr2.description)
        SeqIO.write(record, f"{file_path.split('.fasta')[0]}_fragment_{int((i*n + h)/teiler)}{bst}_{int(end/teiler)}{bst}.fasta", "fasta")


if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument('-file', required=True, type=str)
    parser.add_argument('-n', required=True, type=int)

    args = parser.parse_args()

    file_path = args.file
    n = args.n

    teiler = 1
    bst = 'H'
    if (n > 999):
        teiler = 1000
        bst = 'K'
        if (n > 999999):
            teiler = 1000000
            bst = 'B'

    chr2 = SeqIO.read(file_path, "fasta")

    split_fasta(chr2, n, 0)
    split_fasta(chr2, n, int(n/2))

