import argparse as ap
import re

if __name__ == "__main__":

    parser = ap.ArgumentParser()
    parser.add_argument('-file_iter', required=True, type=str, help="gff file, iteration output unfiltered")
    parser.add_argument('-file_filt', required=True, type=str, help="gff file, output of iteration filtered")
    parser.add_argument('-out', required=True, type=str, help="output file")
    parser.add_argument('-db_file', required=True, type=str, help="txt file with all db entries in my format")
    parser.add_argument('-cutoff', default=0.9, type=float, help="min overlap percentage of hit sequence and cds")

    args = parser.parse_args()

    file_iter = args.file_iter
    file_filt = args.file_filt
    out_file = args.out
    db_file = args.db_file
    cutoff = args.cutoff

    # read in all database entries and save them in the list of lists entries (list: organism_id/gene_id, gene_range, strand, ggf. ivan_exon_range)
    entries = []
    with open(db_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.split('|') in entries:
                entries.append(line.split('|'))

    # read in filtered file line by line and save infos in the tupel list filt_hits (tupel: organism_id/gene_id, query_start, query_end, query_strand)
    filt_hits = []
    with open(file_filt, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                cont = line.split('\t')
                filt_hits.append((cont[0], cont[3], cont[4], cont[6]))
    filt_hits = list(set(filt_hits))
    
    # read in blast hit file line by line and save lines that are still there after filtering in the tupel list orthologs (tupel: organism_id/gene_id, query_start, query_end, query_strand, ivan_cds_id_ortholog/gene_id_ortholog)
    f_hits = [] # list of found filtered hits
    orthologs = []
    while(cutoff >= 0.9 and len(list(set(f_hits))) != len(filt_hits)):
        f_hits = []
        orthologs = []
        with open(file_iter, 'r') as file:
            for line in file:
                if not line.startswith('#'):
                    cont = line.split("\t")
                    start = int(cont[3])
                    end = int(cont[4])
                    strand = cont[6]

                    for f in filt_hits:
                        s = int(f[1])
                        e = int(f[2])
                        st = f[3]
                        if st == strand:
                            if s <= start and end <= e:
                                f_hits.append(f)
                                orthologs.append((f[0], s, e, st, cont[0]))
                            elif (s <= start and start <= e and (e-start)/(end-start) >= cutoff):
                                f_hits.append(f)
                                orthologs.append((f[0], s, e, st, cont[0]))
                            elif (s <= end and end <= e and (end-s)/(end-start) >= cutoff) :
                                f_hits.append(f)
                                orthologs.append((f[0], s, e, st, cont[0]))
        cutoff -= 0.1
    
    print(cutoff + 0.1)

    # list lines contains all lines to be written into output file
    lines = []
    for (org, start, end, strand, ortho) in list(set(orthologs)):
        query = ""
        queries = []
        ortholog = ""
        # rewrite query
        for [o, g, st, i] in entries:
            if org.split('.')[0] == o.split('.')[0] and strand == st:
                if org.startswith('ENS') or org.startswith('BRAKER') or (int(i.split('-')[0]) == start and int(i.split('-')[1]) == end):
                    if query == "":
                        query = "|".join([o.split('T0')[0], g, st, i])
                        queries.append(query)
                    else:
                        queries.append("|".join([o.split('.')[0], g, st, i]))
                        print("multiple query strings", query, [o, g, st, i])
    
        # rewrite ortholog
        if ortho.startswith('ENS') or ortho.startswith('BRAKER'):
            for [o, g, st, i] in entries:
                if o.split('.')[0] == ortho.split('.')[0]:
                    if ortholog == "":
                        ortholog = "|".join([o.split('T0')[0], g, st, i])
                    else:
                        print("multiple ortholog strings", ortholog, ortho)
        else:
            ort = ortho.split('|:')
            [a, b, c] = ort[0].split('-')
            b = int(float(b))
            c = int(float(c))
            [e, f] = ort[1].split('-')
            e = int(float(e))
            f = int(float(f))
            e += b
            f += b - 1
            for [o, g, st, i] in entries:
                if o == a and g == f"{b}-{c}" and i == f"{e}-{f}":
                    if ortholog == "":
                        ortholog = "|".join([o, g, st, i])
                    else:
                        print("multiple ortholog strings", ortholog, ortho)
        
        if ortholog == "":
            print("no ortholog")
        elif query == "":
            print("no query", (org, start, end, strand, ortho))
        else:
            for q in queries:
                lines.append(f"{ortholog}\t{q}\n")
                #out.write(f"{ortholog}\t{query}\n")
    
    out = open(out_file, 'w')
    for l in list(set(lines)):
        out.write(l)
    out.close()
