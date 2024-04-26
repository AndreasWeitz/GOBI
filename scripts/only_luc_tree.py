import subprocess
import argparse

def tblastn(db_path, db_fasta, query_fasta, out):
    
    makedb = f"makeblastdb -in {db_fasta} -dbtype prot -out {db_path}"
    subprocess.run(makedb, shell=True)
    
    blast = f"tblastn -db /home/thekiller/gobi/data/only_luc_tree/V2/db/database -query {query_fasta} -out {out} -evalue 0.001"
    subprocess.run(blast, shell=True)

def main():

    #parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    #parser.add_argument('-dp', '--db_path', help="db_translated_file", type=str)
    #parser.add_argument('-df', '--db_fasta', help="MDF for modelfinder", type=str)
    #parser.add_argument('-qf', '--query_fasta', help="db_translated_file", type=str)
    #parser.add_argument('-o', '--out', help="MDF for modelfinder", type=str)
    #args = parser.parse_args()

    #db_path = args.db_path
    #db_fasta = args.db_fasta
    #query_fasta = args.query_fasta
    #out = args.out

    db_path = "/home/thekiller/gobi/data/only_luc_tree/V2/db/database"
    db_fasta = "/home/thekiller/gobi/data/only_luc_tree/V2/db_translated.fasta"
    query_fasta = "/home/thekiller/gobi/data/only_luc_tree/V2/luc_id.fasta"
    out = "/home/thekiller/gobi/data/only_luc_tree/V2/blastoutput"

    tblastn(db_path, db_fasta, query_fasta, out)


if __name__ == '__main__':
    main()
    