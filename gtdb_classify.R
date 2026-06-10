library(ape)
library(ggtree)
library(dplyr)
library(readr)
library(stringr)
library(tidytree)
library(tidyr)
library(ggplot2)
library(magrittr)


#################################
#     specify genome names      #
#################################


rolando_genomes <- read.table("rolando_mag_accessions.txt") %>% pull(V1)
huang_genomes <- read.table("huang_mag_accessions.txt") %>% pull(V1)
osvatic_genomes <- read.table("osvatic_mag_accessions.txt") %>% pull(V1)
giani_genomes <- read.table("giani_mag_accessions.txt") %>% pull(V1)
morel_genomes <- read.table("morel_mag_accessions.txt") %>% pull(V1)
petersen_clam <- c("luna_ont_bin1", "lotti_ont_bin1")
#gtdb_genomes <- read.table("GTDB_subtree_mag_acessions.txt") %>% pull(V1)

ficus_genomes <- c("PIE_T1P1_root_Marsh_MG__3300081559_11",
                   "PIE_T1P1_root_Marsh_MG__3300081559_s3",
                   "PIE_T2P3_root_Marsh_MG__3300077085_6",
                   "PIE_T3P3_root_Marsh_MG__3300077490_s15")







#################################################
#     read in tree and get smaller subtree      #
#################################################

bac120 <- read_tsv("bac120_metadata.tsv")
bac120$accession <- substring(bac120$accession, 4)

# read in tree from gtdbtk
gtdb_tree<-ape::read.tree("gtdbtk.bac120.user_msa.fasta.treefile")


gtdb_tree$tip.label <- ifelse(startsWith(gtdb_tree$tip.label, "GCA"), paste0("GCA_", str_split_i(gtdb_tree$tip.label, "_", 2)), gtdb_tree$tip.label)

# genomes start with a random string so remove these so I can match up with genome accession names
# gtdb_tree$tip.label <- substring(gtdb_tree$tip.label, 4)

# get the ancestor of clam, china and usa tips, then use this node number to tree_subset out a clade to work with a smaller tree, taking a few nodes up for broader context
#getMRCA(gtdb_tree, c("GCA_050305635.1", "GCA_037384565.1"))
#genomes_in_tree <- c(rolando_genomes, huang_genomes, osvatic_genomes, morel_genomes, giani_genomes)[c(rolando_genomes, huang_genomes, osvatic_genomes, morel_genomes, giani_genomes) %in% gtdb_tree$tip.label]
#gtdb_tree_subset <- tree_subset(gtdb_tree, getMRCA(gtdb_tree, genomes_in_tree), 7)


#################################################
#         annotate taxonomy for tree            #
#################################################

# use the total bac120 dataset to get the genus and if available species for each tip
gtdb_tree$genus <- left_join(data.frame(accession=gtdb_tree$tip.label), bac120, by="accession") %>% dplyr::select(accession, ncbi_taxonomy) %>% mutate(ncbi_taxonomy=str_split_i(ncbi_taxonomy, ";", 6)) %>% pull(ncbi_taxonomy) %>% str_replace("g__", "")
gtdb_tree$species <- left_join(data.frame(accession=gtdb_tree$tip.label), bac120, by="accession") %>% dplyr::select(accession, ncbi_taxonomy) %>% mutate(ncbi_taxonomy=str_split_i(ncbi_taxonomy, ";", 7)) %>% pull(ncbi_taxonomy) %>% str_replace("s__", "")
gtdb_tree$species <- ifelse(gtdb_tree$species == "", gtdb_tree$genus, gtdb_tree$species)




# label genomes by which dataset and host they came from
gtdb_tree$genome <- case_when(gtdb_tree$tip.label %in% rolando_genomes ~ "USA Spartina",
                                     gtdb_tree$tip.label %in% huang_genomes ~ "China Spartina",
                                     gtdb_tree$tip.label %in% osvatic_genomes ~ "Clam Osvatic",
                                     gtdb_tree$tip.label %in% giani_genomes ~ "Clam Giani",
                                     gtdb_tree$tip.label %in% morel_genomes ~ "Clam Morel",
                                     gtdb_tree$tip.label %in% petersen_clam ~ "Clam Petersen",
                                     gtdb_tree$tip.label %in% ficus_genomes ~"Ficus Spartina",
                                     #gtdb_tree$tip.label %in% gtdb_genomes ~ "gtdb input",
                                    !(gtdb_tree$tip.label) %in% c(rolando_genomes, huang_genomes, osvatic_genomes) ~ "GTDB")


gtdb_tree$host <- case_when(gtdb_tree$tip.label %in% rolando_genomes ~ "Plant",
                              gtdb_tree$tip.label %in% huang_genomes ~ "Plant",
                              gtdb_tree$tip.label %in% osvatic_genomes ~ "Clam",
                              gtdb_tree$tip.label %in% giani_genomes ~ "Clam",
                              gtdb_tree$tip.label %in% morel_genomes ~ "Clam",
                              gtdb_tree$tip.label %in% petersen_clam ~ "Clam",
                              gtdb_tree$tip.label %in% ficus_genomes ~"Plant")

##########################################################
#       get DTDB genome names to download on lisc        #
##########################################################

#gtdb_tree_subset$tip.label[gtdb_tree_subset$genome == "GTDB"] %>% write.table("GTDB_subtree_mag_acessions.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#################################
#         build tree            #
#################################

dd <- data.frame(taxa=gtdb_tree$tip.label,
                 genome=gtdb_tree$genome,
                 genus=gtdb_tree$genus,
                 species=gtdb_tree$species)

p<-ggtree(gtdb_tree, size=0.3, colour="gray60")


copen <- cophenetic.phylo(gtdb_tree)
hc <- hclust(as.dist(copen), method="complete")
clusters <- cutree(hc, h=0.0005)

keep <- sapply(split(names(clusters), clusters),
               `[`, 1)

gtdb_tree_reduced <- keep.tip(gtdb_tree, keep)


gtdb_tree_reduced
p<-ggtree(gtdb_tree_reduced, size=0.3, colour="gray60")



png("clam_spartina_tree_iqtree_accession_reduced.png", height=2000, width=1000)
p %<+% dd + 
  geom_tippoint(aes(color=genome), size=4) + 
  geom_tiplab(size=5.5) +
  theme(legend.text = element_text(size=20),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=c("seagreen", "dodgerblue", "goldenrod1", "purple", "red", "blue", "gray70", "lightgreen")) +
  #scale_size_manual(values=c(6, 6, 6, 6, 6, 6, 6)) +
  guides(size = guide_legend(override.aes = list(size = 7))) + 
  xlim(NA, 0.15)
dev.off()

png("clam_spartina_tree_iqtree_species_reduced.png", height=2000, width=1000)
p %<+% dd + 
  geom_tippoint(aes(color=genome), size=4) + 
  geom_tiplab(aes(label=species), size=5.5) +
  theme(legend.text = element_text(size=20),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=c("seagreen", "dodgerblue", "goldenrod1", "purple", "red", "blue", "gray70", "lightgreen")) +
  #scale_size_manual(values=c(6, 6, 6, 6, 6, 6, 6)) +
  guides(size = guide_legend(override.aes = list(size = 7))) + 
  xlim(NA, 0.15)
dev.off()


















png("clam_spartina_tree_iqtree_accession.png", height=4300, width=2500)
p %<+% dd + 
  geom_tippoint(aes(color=genome), size=4) + 
  geom_tiplab(size=5.5) +
  theme(legend.text = element_text(size=20),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=c("seagreen", "dodgerblue", "goldenrod1", "purple", "red", "blue", "gray70", "lightgreen")) +
  #scale_size_manual(values=c(6, 6, 6, 6, 6, 6, 6)) +
  guides(size = guide_legend(override.aes = list(size = 7))) + 
  xlim(NA, 0.15)
dev.off()

png("clam_spartina_tree_iqtree_species.png", height=4300, width=2500)
p %<+% dd + 
  geom_tippoint(aes(color=genome), size=4) + 
  geom_tiplab(aes(label=species), size=5.5) +
  theme(legend.text = element_text(size=20),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=c("seagreen", "dodgerblue", "goldenrod1", "purple", "red", "blue", "gray70", "lightgreen")) +
  #scale_size_manual(values=c(6, 6, 6, 6, 6, 6, 6)) +
  guides(size = guide_legend(override.aes = list(size = 7))) + 
  xlim(NA, 0.15)
dev.off()











genome_quality <- read_tsv("quality_report.tsv") %>% 
  dplyr::rename("assembly"="Name") %>%
  dplyr::select(c(assembly, Contamination, Contig_N50, Genome_Size, Completeness)) %>%
  mutate(assembly = ifelse(startsWith(assembly, "GCA"), paste0("GCA_", str_split_i(assembly, "_", 2)), assembly)) %>%
  left_join(data.frame(assembly=gtdb_tree$tip.label, host=gtdb_tree$host, study=gtdb_tree$genome), by="assembly")




genome_quality %>%
  dplyr::select(study, host, Contamination, Completeness, Contig_N50, Genome_Size) %>%
  data.frame() %>% 
  filter(study %in% c("USA Spartina", "China Spartina", "Ficus Spartina")) %>%
  ggplot(aes(x = Completeness, y = Genome_Size)) + 
  geom_point(size=5) +
  ylab("Plant associated Ca. Thio Genome Length (Mb") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25))

genome_all %>%
  dplyr::select(study, host, Contamination, Completeness, n50, length) %>%
  data.frame() %>% 
  filter(study %in% c("rolando", "huang")) %>%
  write_xlsx("thio_genomesize_completeness.xlsx")
  








#png("tree_nodes2.png", height=3000, width=4000)
p %<+% dd + 
  geom_tippoint(aes(color=genome, size=genome)) + 
  scale_colour_manual(values=c("seagreen", "dodgerblue", "goldenrod1", "purple", "red", "gray70", "lightgreen")) + 
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=6)
#dev.off()

p_collapsed <- p %<+% dd + 
  geom_tippoint(aes(color=genome, size=genome)) + 
  scale_colour_manual(values=c("seagreen", "dodgerblue", "goldenrod1", "purple", "red", "gray70", "lightgreen")) + 
  scale_size_manual(values=c(7, 3, 3, 3, 3, 3, 7))


p_collapsed <- scaleClade(p_collapsed, node = 448, scale = 3)
p_collapsed <- collapse(p_collapsed, node = 380, mode = "max", fill="goldenrod1")
p_collapsed <- collapse(p_collapsed, node = 450, mode = "max", fill="goldenrod1")
p_collapsed <- collapse(p_collapsed, node = 583, mode = "max", fill="goldenrod1")

png("test.png", height=1000, width=1000)
p_collapsed
dev.off()

