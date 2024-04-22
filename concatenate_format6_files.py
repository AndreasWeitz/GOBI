result_filename = "../blast-results.alt"
filename = "blastresults_filename_list_file.txt"

with open(result_filename, 'w') as w:
    with open(filename, 'r') as list_r:
        for line in list_r.readlines():
            with open(line.strip(), 'r') as r:
                for line_2 in r.readlines():
                    w.write(line_2)
