import argparse as ap
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

if __name__ == "__main__":
     # -file Melanotus_villosus\chr2_fragment_5.tsv -n 1000
    
    parser = ap.ArgumentParser()
    parser.add_argument('-file', required=True, type=str)

    args = parser.parse_args()

    file_path = args.file

    data_dict = {}
    files = []
    with open(file_path, 'r') as file:
        for line in file:
            line_split = line.split('\t')
            id = line_split[8].split(';')[0].split("ID=")[1]
            start = int(line_split[3])
            end = int(line_split[4])
            
            if (id not in data_dict.keys()):
                data_dict[id] = [start+(end-start)/2]
            else:
                data_dict[id].append(start+(end-start)/2)
            files.append(line)

    mu_std = {}
    for key in data_dict.keys():
        mu, std = norm.fit(data_dict[key])
        mu_std[key] = (mu, std)
        
        """
        plt.hist(data_dict[key], bins=25, density=True, color='g')
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)
        plt.plot(x, p, 'k', linewidth=2)
        title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
        plt.title(title)
        plt.axvline(mu)
        plt.axvline(mu-3*std)
        plt.axvline(mu+3*std)
        plt.show()
        """
        
    out_file = open(file_path.split('.gff')[0]+"_std.gff", 'w')
    for line in files:
        line_split = line.split('\t')
        id = line_split[8].split(';')[0].split("ID=")[1]
        start = int(line_split[3])
        end = int(line_split[4])
        middle = start + (end-start)/2

        mu, std = mu_std[id]
        if (middle < mu-3*std or middle > mu+3*std):
            continue
        else:
            out_file.write(line)
    