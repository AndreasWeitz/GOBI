import subprocess

def blastn(db_path, gene_path, out, db_fasta = None, evalue = 0.001, outtype = None):
    
    '''
    # db_path: path of .db file, or where it should be created
    # gene_path: query, a gene fasta file
    # out: output path & file name
    # db_fasta: to create .db file at db_path from fasta.
    # evalue: define similarity
    '''
    
    if db_fasta != None:
        makedb = f"makeblastdb -in {db_fasta} -dbtype nucl -out {db_path}"
        subprocess.run(makedb, shell=True)
    if outtype:
        blast = f"blastn -db {db_path} -query {gene_path} -out {out} -evalue {str(evalue)} -outfmt {str(outtype)}"
    else:
        blast = f"blastn -db {db_path} -query {gene_path} -out {out} -evalue {str(evalue)}"
    subprocess.run(blast, shell=True)
    
def tblastn(db_path, gene_path, out, db_fasta = None, evalue = 0.001, outtype = None):
    
    '''
    # db_path: path of .db file, or where it should be created
    # gene_path: query, a gene fasta file
    # out: output path & file name
    # db_fasta: to create .db file at db_path from fasta.
    # evalue: define similarity
    '''
    
    if db_fasta != None:
        makedb = f"makeblastdb -in {db_fasta} -dbtype nucl -out {db_path}"
        subprocess.run(makedb, shell=True)

    if outtype:
        blast = f"tblastn -db {db_path} -query {gene_path} -out {out} -evalue {str(evalue)} -outfmt {str(outtype)}"
    else:
        blast = f"tblastn -db {db_path} -query {gene_path} -out {out} -evalue {str(evalue)}"
    subprocess.run(blast, shell=True)
    
def tblastx(db_path, gene_path, out, db_fasta = None, evalue = 0.001, outtype = None):
    
    '''
    # db_path: path of .db file, or where it should be created
    # gene_path: query, a gene fasta file
    # out: output path & file name
    # db_fasta: to create .db file at db_path from fasta.
    # evalue: define similarity
    '''
    
    if db_fasta != None:
        makedb = f"makeblastdb -in {db_fasta} -dbtype nucl -out {db_path}"
        subprocess.run(makedb, shell=True)

    if outtype:
        blast = f"tblastx -db {db_path} -query {gene_path} -out {out} -evalue {str(evalue)} -outfmt {str(outtype)}"
    else:
        blast = f"tblastx -db {db_path} -query {gene_path} -out {out} -evalue {str(evalue)}"
    subprocess.run(blast, shell=True)