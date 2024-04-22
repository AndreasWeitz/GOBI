import subprocess

# required
blast_output = ""
gffs_folder = ""
corrected_gffs_folder = ""
concatenated_gff = ""
genome_file = ""

# not required
cds_hits = ""
hits_without_map_on_cds = ""

scripts = [
    ("multi_blast2gff.py", ["-i", blast_output, "-o", gffs_folder]),
    ("correct_ranges_of_blasthits_after_fastasplitter.py", ["-gff", gffs_folder, "-out", corrected_gffs_folder]),
    ("concatenate_folder2_1file.py", ["-folder", corrected_gffs_folder, "-out", concatenated_gff]),
    ("blast_hits2cds.py", ["-bh", concatenated_gff, "-g", genome_file, "-o", cds_hits, "-n", hits_without_map_on_cds])
]

def run_scripts(scripts):
    for script, params in scripts:
        command = ["python", script] + params
        subprocess.run(command)

if __name__ == "__main__":
    run_scripts(scripts)
