def view(iter_file, 
         range_threshhold=2000000, 
         fig_width=150, 
         file_type="png", 
         out_dpi=160, 
         plot_height=3, 
         linewidth=0.8,
         label_size=17, 
         cluster_label=17, 
         show_rglabels=True, 
         tick_size=50, 
         label_rotation=0):
    
    from pygenomeviz import Gff, GenomeViz
    import os
    
    # parse iter-file:
    link_dicts = []
    with open(iter_file, "r") as r:
        for line in r:
            species_tuple = line.rstrip().split("\t")
            species1 = species_tuple[0].split("|")
            species2 = species_tuple[1].split("|")
            link_dicts.append(dict(db_species = species1[0], 
                                   db_gene = species1[1].split("-"), 
                                   db_strand = 1 if species1[2] == "+" else -1, 
                                   db_hit = species1[3],
                                   query_species = species2[0], 
                                   query_gene = species2[1].split("-"), 
                                   query_strand = 1 if species2[2] == "+" else -1, 
                                   query_hit = species2[3]))
            
    # parse "NW_022170249.1" IDs from file
    id_dict = {}
    with open("gene_name_mapping.csv", "r") as r:
        for line in r:
            line = line.rstrip().split(",")
            id_dict[line[0].split("-", 1)[1]] = line[2]
  
    # this is an easter egg! Do not remove, necessary for the code to work ;) 
    def is_nussbaum(noah):
        if noah == "NW_022170249.1":
            return 0
        else:
            return 1

    link_dicts.sort(key = lambda rec : (int(rec["query_gene"][0]), int(rec['query_gene'][1]),
                                        is_nussbaum(rec["db_species"])))
    #for entry in link_dicts:
     #   print((entry["db_species"], entry["query_gene"]))
    
    plot_ranges = []
    start = int(link_dicts[0]["query_gene"][0])
    # end & temp_start must be initialized differently than first record. -> color coding & labeling!
    end = int(link_dicts[0]["query_gene"][1]) + 1
    temp_start = start + 1
    
    overlap_counter = 1
    for record in link_dicts:
        if temp_start == record["query_gene"][0] or end == record["query_gene"][1]:
            record["query_gene"][1] = f"{overlap_counter * "\n"}{record["query_gene"][1]}"
            overlap_counter += 1
        else:
            overlap_counter = 1
          
        # decide plot_ranges
        if temp_start + range_threshhold >= int(record["query_gene"][0]):
            temp_start = int(record["query_gene"][0])
            end = int(record["query_gene"][1])
        else:
            plot_ranges.append([start, end])
            start = int(record["query_gene"][0])
            end = int(record["query_gene"][1])
            temp_start = start
    plot_ranges.append([start, end])
    
    ch = ColorHelper()
    for plot_range in plot_ranges:
        
        # define plot attributes
        gv = GenomeViz(fig_track_height=plot_height, tick_style="axis", fig_width=fig_width, tick_labelsize = tick_size)
        
        # filter for db features in corresponding range of query
        filtered_links = list(filter(lambda r : is_inrange(plt_range=plot_range, record = r), link_dicts))
        plot_dict = {}
        for link_dict in filtered_links:
            # get track ranges of all species in our custom database -> plot_dict[db_species]
            if link_dict["db_species"] in plot_dict.keys():
                if plot_dict[link_dict["db_species"]][0] > int(link_dict["db_gene"][0]):
                    plot_dict[link_dict["db_species"]][0] = int(link_dict["db_gene"][0])
                if plot_dict[link_dict["db_species"]][1] < int(link_dict["db_gene"][1]):
                    plot_dict[link_dict["db_species"]][1] = int(link_dict["db_gene"][1])
            else:
                plot_dict[link_dict["db_species"]] = [int(link_dict["db_gene"][0]), int(link_dict["db_gene"][1])]
                
        # add db tracks (eastern firefly first) & add query track
        track_dict = {}
        if "NW_022170249.1" in plot_dict.keys():
            track_dict["NW_022170249.1"]= gv.add_feature_track(name="NW_022170249.1", 
                                               size=plot_dict["NW_022170249.1"][1] - plot_dict["NW_022170249.1"][0], 
                                               start_pos=plot_dict["NW_022170249.1"][0], 
                                               labelsize=label_size, 
                                               labelmargin=0.005, 
                                               linewidth=30)
        for species in plot_dict.keys():
            if species != "NW_022170249.1":
                track_dict[species]= gv.add_feature_track(name=species, 
                                               size=plot_dict[species][1] - plot_dict[species][0], 
                                               start_pos=plot_dict[species][0], 
                                               labelsize=label_size, 
                                               labelmargin=0.005, 
                                               linewidth=30)
        query_track = gv.add_feature_track(name=link_dicts[0]["query_species"], 
                                           size=plot_range[1] - plot_range[0], 
                                           start_pos=plot_range[0], 
                                           labelsize=label_size, 
                                           labelmargin=0.005, 
                                           linewidth=30)
        
        # add all features to their tracks  
        filtered_links.sort(key = lambda rec : is_nussbaum(rec["db_species"]))
        #for entry in filtered_links:
        #    print((entry["db_species"], entry["query_gene"]))
        for link_dict in filtered_links:
            db_id = "-".join(link_dict["db_gene"])
            query_id = "-".join(link_dict["query_gene"])
            
            # translate photinus pyralis color code to query color code
            if link_dict["db_species"] == "NW_022170249.1":
                ch.add_to_translator(id_dict[db_id], query_id)
               
            track_dict[link_dict["db_species"]].add_feature(start = int(link_dict["db_gene"][0]), 
                                                            end = int(link_dict["db_gene"][1]),
                                                            strand = link_dict["db_strand"],
                                                            plotstyle="arrow", 
                                                            linewidth=linewidth, 
                                                            labelsize=label_size if link_dict["db_species"] != "NW_022170249.1" else cluster_label,
                                                            label = 
                                                            (db_id if show_rglabels else " ") if link_dict["db_species"] != "NW_022170249.1" 
                                                            else id_dict[db_id], 
                                                            facecolor=ch.color_func(query_id, display_type="pergene"), 
                                                            labelrotation = label_rotation)
            
            query_track.add_feature(start = int(link_dict["query_gene"][0]), 
                                    end = int(link_dict["query_gene"][1]), 
                                    strand = link_dict["query_strand"],
                                    plotstyle="arrow", 
                                    linewidth=linewidth, 
                                    labelsize=label_size,
                                    label=query_id if show_rglabels else " ",
                                    facecolor=ch.color_func(query_id, display_type="pergene"), 
                                    labelrotation = label_rotation)
        
        # create plot fo this plot_range
        if not os.path.exists(f"./{iter_file}_th{range_threshhold}/"):
            os.makedirs(f"./{iter_file}_th{range_threshhold}/")
        gv.savefig(f"./{iter_file}_th{range_threshhold}/{iter_file}_{plot_range[0]}-{plot_range[1]}.{file_type}", dpi=out_dpi)
        

def is_inrange(plt_range, record):
    return plt_range[0] <= int(record["query_gene"][0]) and plt_range[1] >= int(record["query_gene"][1])

def idlabel_func(label):
    if isinstance(label, list):
        label = label[0]
    if label[0] == "~":
        return ""
    else:
        parts = label.split(":")
        if len(parts) > 1:
            return parts[1]
        else:
            return parts[0] 
        
        
class ColorHelper:
    
    def __init__(self):
        self.translator = {}
        self.color_helper = [-1, {}]
        self.color_list = ['#FF6347', '#B8860B', '#00FA9A', '#20B2AA', '#48D1CC', '#FFFFFF', '#32CD32', '#800000', '#FF8C00', '#B0E0E6', '#F5DEB3', 
                           '#E6E6FA', '#8B4513', '#FFD700', '#FFC0CB', '#BC8F8F', '#800080', '#008080', '#C0C0C0', '#FF69B4', '#1E90FF', '#00FF00', 
                           '#9932CC', '#000000', '#00BFFF', '#00FFFF', '#7FFFD4', '#D3D3D3', '#B0C4DE', '#6A5ACD', '#FFA07A', '#DAA520', '#9370DB', 
                           '#00FF7F', '#DCDCDC', '#FFFF00', '#F0E68C', '#FF0000', '#DC143C', '#FF1493', '#778899', '#8A2BE2', '#A52A2A', '#5F9EA0', 
                           '#8FBC8F', '#40E0D0', '#FF4500', '#0000FF', '#7B68EE', '#9400D3', '#00CED1', '#008B8B', '#FFA500', '#6495ED', '#808000', 
                           '#FF00FF', '#008000', '#4682B4', '#FFB6C1', '#808080', '#000080']

    def add_to_translator(self, db_id, query_id):
        if db_id in {"4CL1-like_1","LUC-like", "LUC", "4CL1-like_2", "4CL1-like_3"}:
            db_id = "cluster"
        self.translator[query_id] = db_id
        
    def color_func(self, ID, display_type):

        if display_type == "default":
            return "#ffaa0040"
        
        elif display_type == "overlap":
            if ID == "~e":
                return "#23bdff20"
            elif ID == "~s":
                return "#ffaa0020"
            elif ID == "~i":
                return "#67ff7120"
            else:
                return "#fff6d220"
            
        elif display_type == "pergene":
            if ID in self.translator.keys():
                ID = self.translator[ID]
            if ID in self.color_helper[1].keys():
                return self.color_helper[1][ID]
            else:
                self.color_helper[0] += 1
                self.color_helper[1][ID] = f"{self.color_list[self.color_helper[0]]}60"
                return self.color_helper[1][ID]
            