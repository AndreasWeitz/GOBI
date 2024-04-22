import subprocess

results = "results_2/pipeline_output_2"
db_path = "./results_2/db_entries_2/db_entries_2_complete.txt"
out = "results_2/orthologs_new"

for i in range(1,10):
    #subprocess.run(["python", "./find_ortholog_parent.py", "-file_iter", f"./{results}/blast_hits_0{i}.gff", "-file_filt", f"./{results}/matching_cds_filtered_0{i}.gff", "-out", f"./{out}/ortholog_parent_0{i}.tsv"])
    subprocess.run(["python", "./ortholog_parent_2.py", "-file_iter", f"./{results}/blast_hits_0{i}.gff", "-file_filt", f"./{results}/matching_cds_filtered_0{i}.gff", "-db_file", db_path, "-out", f"./{out}/orthologs_0{i}.tsv"])
    
for i in range(10,23):
    #subprocess.run(["python", "./find_ortholog_parent.py", "-file_iter", f"./{results}/blast_hits_{i}.gff", "-file_filt", f"./{results}/matching_cds_filtered_{i}.gff", "-out", f"./{out}/ortholog_parent_{i}.tsv"])
    subprocess.run(["python", "./ortholog_parent_2.py", "-file_iter", f"./{results}/blast_hits_{i}.gff", "-file_filt", f"./{results}/matching_cds_filtered_{i}.gff", "-db_file", db_path, "-out", f"./{out}/orthologs_{i}.tsv"])
