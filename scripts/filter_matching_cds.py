import argparse


def read_gff(gff_file):

    features = []

    with open(gff_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                feature = line.strip().split('\t')
                seqid = feature[0]
                source = feature[1]
                feature_type = feature[2]
                start = int(feature[3])
                end = int(feature[4])
                score = feature[5]
                strand = feature[6]
                phase = feature[7]
                attributes = feature[8]
                features.append([seqid, source, feature_type, start, end, score, strand, phase, attributes])

    return features


def get_clusters(exons, d):

    clusters = []
    tmp_cluster = []

    for exon in exons:
        if len(tmp_cluster) == 0:
            tmp_cluster.append(exon)
        else:
            if abs(tmp_cluster[-1][4] - exon[3]) < d:
                tmp_cluster.append(exon)
            else:
                clusters.append(tmp_cluster)
                tmp_cluster = [exon]

    if tmp_cluster not in clusters:
        clusters.append(tmp_cluster)
        
    return clusters


def filter_clusters(clusters, number_of_exons, p):

    filtered_clusters = []

    for cluster in clusters:
        if len(cluster) > number_of_exons*p:
            filtered_clusters.append(cluster)

    return filtered_clusters
    

def write_gff(clusters, output_file):

    with open(output_file, 'w') as f:
        for cluster in clusters:
            for feature in cluster:
                f.write('\t'.join([feature[0], feature[1], feature[2], str(feature[3]), str(feature[4]), feature[5], feature[6], feature[7], feature[8]]) + '\n')


def print_cluster_info(clusters, p, number_of_exons):

    clusters.sort(key=lambda x: -len(x))
    
    print("+--------------------------------------------- --- --- --  - -   -")
    print(f"|      >>>   {clusters[0][0][0]}  {clusters[0][0][1]}  <<<   ")
    print("+---------+-------+------------+----------- --- --- --  - -")
    print("| cluster | exons | proceeding |   cluster range ")
    print("+---------+-------+------------+---- --- -- - -  -")

    micro_cluster_counter = 0

    for i, cluster in enumerate(clusters, start=1):
        if len(cluster) > number_of_exons*p:
            cluster_status = "+"
        else:
            cluster_status = "-"
        
        if len(cluster) > 1:
            if len(cluster) > 9 and i < 10:
                print(f"|   {i}     |  {len(cluster)}   |     {cluster_status}      |     {abs(cluster[0][3]-cluster[-1][4])}")
            elif i < 10 and len(cluster) < 10:
                print(f"|   {i}     |  {len(cluster)}    |     {cluster_status}      |     {abs(cluster[0][3]-cluster[-1][4])}")
            elif i > 9 and len(cluster) < 10:
                print(f"|   {i}    |  {len(cluster)}    |     {cluster_status}      |     {abs(cluster[0][3]-cluster[-1][4])}")

        else:
            micro_cluster_counter += 1
    
    if micro_cluster_counter > 0:
        if micro_cluster_counter > 9:
            print(f"|   {micro_cluster_counter}x   |  1    |     -      |")
        else:
            print(f"|   {micro_cluster_counter}x    |  1    |     -      |")
    print("+---------+-------+------------+-- --- -- - -  -")


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-inp', '--input_gff', help="Path to input GFF file", type=str)
    parser.add_argument('-out', '--output_gff', help="Path to output GFF file", type=str)
    parser.add_argument('-d', '--intra_cluster_distance', help="max distance between the 2 nearest elements in a cluster (default=100000)", type=int)
    parser.add_argument('-p', '--percent', help="len(cluster) > number_of_exons*0.3, you can modify the 0.3 with this parameter (to determin which cluster is taken)", type=float)
    args = parser.parse_args()

    input_gff = args.input_gff
    output_gff = args.output_gff
    d = args.intra_cluster_distance
    p = args.percent

    if d is None:
        d = 100000
    if p is None:
        p = 0.1

    print_cluster_info(get_clusters(read_gff(input_gff), d), p, len(read_gff(input_gff)))
    write_gff(filter_clusters(get_clusters(read_gff(input_gff), d), len(read_gff(input_gff)), p), output_gff)


if __name__ == '__main__':
    main()
