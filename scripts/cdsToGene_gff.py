# basename does not include .gff, but existing path!
def cdstogene(basename):
    with open(f"{basename}.gff", "r") as f:
        with open(f"{basename}_togene.gff", "w") as outfile:
            gene_list = []
            gene_id = ""
    
            for line in f:
                columns = line.split("\t")
                attrs = columns[-1].split(";")
                
                for attribute in attrs:
                    cds_attr = attribute.split("=")
                    if cds_attr[0] == "ID":
                        if gene_id != cds_attr[1]:
                            # initialization of the first gene
                            if gene_id != "":
                                outfile.write("\t".join(gene_list))
                            gene_id = cds_attr[1]
                            gene_list = columns[0:-1]
                            gene_list.append(f"{attribute}\n")
                            gene_list[2] = "gene"
                        else:
                            if columns[3] < gene_list[3]:
                                gene_list[3] = columns[3]
                            if columns[4] > gene_list[4]:
                                gene_list[4] = columns[4]
                                
            outfile.write("\t".join(gene_list))
                                