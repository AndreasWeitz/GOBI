library("tidyverse")
library("pheatmap")
library("ape")


setwd("C:/RebEll/Computer/BioinformatikBachelor/Semester5/gobi/data/my_work/folding_proteins/folding_proteins_workflow/results")

tm_out <- read.csv("tm_score.txt", sep = "\t")
colnames(tm_out) <- c("PDBchain1", "PDBchain2", "TM1", "TM2", "RMSD", "ID1", "ID2", "IDali", "L1", "L2", "Lali")

tm_out$PDBchain1 <- gsub("_mod.pdb", "", tm_out$PDBchain1)
tm_out$PDBchain2 <- gsub("_mod.pdb", "", tm_out$PDBchain2)

tm_out$PDBchain1 <- gsub("./", "", tm_out$PDBchain1)
tm_out$PDBchain2 <- gsub("./", "", tm_out$PDBchain2)

tm_out$PDBchain1 <- gsub("pdb_id_1LCI", "Photinus_pyralis_luciferase_(1LCI)", tm_out$PDBchain1)
tm_out$PDBchain2 <- gsub("pdb_id_1LCI", "Photinus_pyralis_luciferase_(1LCI)", tm_out$PDBchain2)

tm_out <- tm_out %>%
  filter(!is.na(TM1))

tm_1 <- tm_out %>%
  select(c("PDBchain1", "PDBchain2", "TM1"))
colnames(tm_1) <- c("PDBchain1", "PDBchain2", "TM")

tm_2 <- tm_out %>%
  select(c("PDBchain2", "PDBchain1", "TM2"))
colnames(tm_2) <- c("PDBchain1", "PDBchain2", "TM")

species <- unique(c(tm_1$PDBchain1,tm_1$PDBchain2))
tm_3 <- data.frame("PDBchain1"=species, "PDBchain2"=species, "TM"=1)

tm_matrix <- rbind(tm_1, tm_2, tm_3)
tm_matrix <- tm_matrix %>%
  #mutate(TM = 1 - TM) %>%
  pivot_wider(names_from = PDBchain1, values_from = TM)

tm_mat <- data.matrix(tm_matrix)
row.names(tm_mat) <- tm_matrix$PDBchain2
tm_mat <- tm_mat[,-c(1)]

sorted <- c("Agrilus_planipennis", "Folsomia_candida", "Danaus_plexippus", "Drosophila_melanogaster", "Chrysoperla_carnea", "Apis_mellifera", "4-coumarate-CoA_ligase_1-like_3", "4-coumarate-CoA_ligase_1-like_2", "Rhagonycha_fulva", "Limnephilus_lunatus", "Cetonia_aurata", "Carabus_problematicus", "4-coumarate-CoA_ligase_1-like_1", "luciferin_4-monooxygenase-like_isoform_X2", "luciferin_4-monooxygenase-like_isoform_X3", "luciferin_4-monooxygenase-like_isoform_X1", "Agriotes_lineatus", "Photinus_pyralis_luciferase_(1LCI)", "Aquatica_leii", "Pyrocoelia_pectoralis", "Photinus_pyralis", "Malachius_bipustulatus", "Dascillus_cervinus", "Brachypterus_glaber", "Nicrophorus_investigator", "Coccinella_septempunctata", "Anthonomus_grandis", "Tribolium_castaneum", "Chrysolina_oricalcia")
tm_mat <- tm_mat[sorted, sorted]

pheatmap(tm_mat, cluster_cols = F, cluster_rows = F)


dist <- dist(tm_mat, diag = TRUE)
hc <- hclust(dist, method="ward.D2")

my_tree <- as.phylo(hc) 
write.tree(phy=my_tree, file="tm_tree_1.newick")


# besserer plot
tm_mat_t <- t(tm_mat)

dist_t <- dist(tm_mat_t, diag = TRUE)
hc_t <- hclust(dist_t, method="ward.D2")

my_tree_t <- as.phylo(hc_t) 
write.tree(phy=my_tree_t, file="tm_tree_2.newick")