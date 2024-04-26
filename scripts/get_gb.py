from Bio import Entrez
from Bio import SeqIO
from handle_ensemble_file import get_seqRecord_ensembl
from urllib.error import HTTPError


def fetch_sequence(id: str, start=-1, end=-1, rettype="gb") -> SeqIO.SeqRecord:
    '''
    Get a nucleotide sequence from NCBI

    ## Parameters:
    `id`: Genbank identifier, e.g. `XM_031473197.1`
    `start`: start index
    `end`: end index
    `rettype`: `"gb"` (default), `"fasta"`

    ## Returns:
    A `Bio.SeqIO.SeqRecord` object of the `id`.
    
    ## Hint:
    To format the record as a string in fasta format, 
    use:

    ```
    record = fetch_sequence("XM_031473197.1", rettype="fasta")
    fasta_str = record.format("fasta")
    ```
    '''
    Entrez.email = "noah.nussbaumer@tum.de"  # Just stick with this address

    # Fetch the sequence from GenBank
    if start < 0 or end < 0:
        if '|:' in id:
            record = get_seqRecord_ensembl(id)
        else:
            handle = Entrez.efetch(db="nucleotide", id=id, rettype=rettype, retmode="text")
            record = SeqIO.read(handle, rettype)
            handle.close()

    else:
        if '|:' in id:
            record = get_seqRecord_ensembl(id)
        else:
            handle = Entrez.efetch(db="nucleotide", id=id, rettype=rettype, retmode="text", seq_start = start, seq_stop = end)
            record = SeqIO.read(handle, rettype)
            handle.close()

    # Parse the record and extract the sequence
    return record

# Example usage:
if __name__ == '__main__':    
    id = "XM_031473197.1"  # Replace with your GenBank ID
    sequence_record = fetch_sequence(id, rettype="fasta")
    print(sequence_record)
