from Bio.Seq import translate
from Bio import SeqIO

'''start counts from 1'''
def ss_translate(input_seq, out_fasta, start, end):
    whole_genome = [(seq_rec.seq, seq_rec.id, seq_rec.description) for seq_rec in SeqIO.parse(open(input_seq), 'fasta')]
    partial_dna = whole_genome[0][0][start-1:end]
    protein = partial_dna.translate()
    
    with open(out_fasta, "w") as out_file:
        out_file.write(f">{whole_genome[0][1]}:{str(start)}-{str(end)} {whole_genome[0][2]}")
        
        i = 0
        for char in protein:
            if i%80 == 0:
                out_file.write("\n")
            i += 1
            out_file.write(char)
            
    