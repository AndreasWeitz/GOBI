import argparse as ap

if __name__ == "__main__":
     # -file Melanotus_villosus\chr2_fragment_5.tsv -n 1000
    
     parser = ap.ArgumentParser()
     parser.add_argument('-file_iter', required=True, type=str, help="gff file, iteration output unfiltered")
     parser.add_argument('-file_filt', required=True, type=str, help="gff file, output of iteration filtered")
     parser.add_argument('-out', required=True, type=str, help="output file")
     parser.add_argument('-cutoff', default=0.9, type=float, help="min overlap percentage of hit sequence and cds")

     args = parser.parse_args()

     file_iter = args.file_iter
     file_filt = args.file_filt
     out_file = args.out
     cutoff = args.cutoff

     out = open(out_file, 'w')
     
     # read in not filtered file line by line and save the not comment ones in the dict hits (key: id of new organism in iteration; start match query; end match query, value: list of ortholog ids)
     orth = []
     hits = {}
     with open(file_iter, 'r') as file:
          for line in file:
               if (not (line.startswith('#'))):
                    start = int(line.split('\t')[3])
                    end = int(line.split('\t')[4])
                    feature = line.split('\t')[-1].replace('%7C', '|')

                    #ortholog = feature.split(';')[2].split("=")[1].strip()
                    ortholog = line.split('\t')[0]
                    if (ortholog.startswith("ENS")):
                         orth.append(ortholog)
                    elif (ortholog.startswith("BRAK")):
                         orth.append(ortholog)
                    else:
                         orth.append(ortholog)

                    organism = feature.split(';')[0].split('.m')[0].replace("ID=", "")
                    id = f"{organism};{start};{end}"
                    if (id in hits.keys()):
                         hits[id].append(ortholog)
                    else:
                         hits[id] = [ortholog]
     
     # read in filtered file line by line and save them in the dict last_col (key: id of new organism in iteration; start match query; end match query, value: list of ortholog ids)
     with open(file_filt, 'r') as file:
          for line in file:
               cont = line.split('\t')
               id = f"{cont[0]};{cont[3]};{cont[4]}"

               orthologs = []
               for hit in hits.keys():
                    h = hit.split(';')
                    if (int(cont[3]) <= int(h[1]) and int(h[2]) <= int(cont[4])):
                         for i in hits[hit]:
                              orthologs.append(i)
                    elif((int(cont[3]) <= int(h[1]) and int(h[1]) <= int(cont[4])) or (int(cont[3]) <= int(h[2]) and int(h[2]) <= int(cont[4]))):
                         if (int(cont[3]) <= int(h[1]) and int(h[1]) <= int(cont[4])):
                              if (float(int(cont[4]) - int(h[1])) / float(int(h[2]) - int(h[1])) > cutoff):
                                   for i in hits[hit]:
                                        orthologs.append(i)
                         else:
                              if (float(int(h[2]) - int(cont[3])) / float(int(h[2]) - int(h[1])) > cutoff):
                                   for i in hits[hit]:
                                        orthologs.append(i)
               orthos = list(set(orthologs))
               orthos.sort()

               #print("#", line)
               for o in orthos:
                    out.write(f"{o}\t{cont[0]}:{cont[3]}-{cont[4]}\n")
     
     out.close()
                    
