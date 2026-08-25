library(ape)
library(ggtree)
library(dplyr)
library(readr)
library(stringr)
library(tidytree)
library(tidyr)
library(ggplot2)
library(magrittr)
library(ComplexHeatmap)


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

plant_genomes <- c(rolando_genomes, huang_genomes, ficus_genomes)
clam_genomes <- c(petersen_clam, morel_genomes, giani_genomes, osvatic_genomes)








#################################################
#     read in tree and get smaller subtree      #
#################################################

bac120 <- read_tsv("bac120_metadata.tsv")
bac120$accession <- substring(bac120$accession, 4)


bac120 %>% pull(ncbi_isolation_source)

# read in tree from gtdbtk
gtdb_tree<-ape::read.tree("gtdbtk.bac120.user_msa.fasta.treefile")

gtdb_tree$genome_name <- gtdb_tree$tip.label


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
gtdb_tree$genome <- case_when(gtdb_tree$tip.label %in% rolando_genomes ~ "Spartina Rolando",
                                     gtdb_tree$tip.label %in% huang_genomes ~ "Spartina Huang",
                                     gtdb_tree$tip.label %in% osvatic_genomes ~ "Clam Osvatic",
                                     gtdb_tree$tip.label %in% giani_genomes ~ "Clam Giani",
                                     gtdb_tree$tip.label %in% morel_genomes ~ "Clam Morel",
                                     gtdb_tree$tip.label %in% petersen_clam ~ "Clam Petersen",
                                     gtdb_tree$tip.label %in% ficus_genomes ~"Spartina Ficus",
                                    !(gtdb_tree$tip.label) %in% c(rolando_genomes, huang_genomes, osvatic_genomes) ~ "GTDB")


gtdb_tree$host <- case_when(gtdb_tree$tip.label %in% rolando_genomes ~ "Plant",
                              gtdb_tree$tip.label %in% huang_genomes ~ "Plant",
                              gtdb_tree$tip.label %in% osvatic_genomes ~ "Clam",
                              gtdb_tree$tip.label %in% giani_genomes ~ "Clam",
                              gtdb_tree$tip.label %in% morel_genomes ~ "Clam",
                              gtdb_tree$tip.label %in% petersen_clam ~ "Clam",
                              gtdb_tree$tip.label %in% ficus_genomes ~"Plant")






data.frame(genome_name=ifelse(startsWith(gtdb_tree$genome_name, "GCA"), gtdb_tree$genome_name %>% str_replace("\\.", "_"), gtdb_tree$genome_name),
           host=gtdb_tree$host,
           genome_source=gtdb_tree$genome) %>%
  set_colnames(c("genome", "host", "genome_source")) %>%
  mutate(host=ifelse(is.na(host), "Unknown", host),
         genome=str_replace(genome, "__", "_")) %>%
  write_tsv("anvio_metadata.tsv")



####################################################
#       root tree by most basal sedimenticola      #
####################################################

gtdb_tree <- root(gtdb_tree, "GCA_027068675.1")


#################################
#         build tree            #
#################################

dd <- data.frame(taxa=gtdb_tree$tip.label,
                 genome=gtdb_tree$genome,
                 genus=gtdb_tree$genus,
                 species=gtdb_tree$species) %>%
  mutate(host=case_when(genome == "Clam Osvatic" ~"Clam",
                        genome == "Spartina Rolando" ~"Plant",
                        genome == "Clam Petersen" ~"Clam",
                        genome == "Clam Morel" ~ "Clam",
                        genome == "Clam Giani" ~"Clam",
                        genome == "Spartina Huang" ~"Plant",
                        genome == "Spartina Ficus" ~"Plant",
                        genome == "GTDB" ~ "GTDB"))



p<-ggtree(gtdb_tree, size=0.3, colour="gray60")


p %<+% dd + 
  geom_tippoint(aes(color=genome), size=3) + 
  geom_tiplab(size=3) +
  scale_colour_manual(values=c("#7ec5f4", "#0000C6", "#1D88AB", "#7F00FF", "gray79", "#28652C", "#849E00", "#5CE65C"))



p %<+% dd + 
  geom_tippoint(aes(color=genome), size=3) + 
  geom_tiplab(aes(label=species), size=3) +
  scale_colour_manual(values=c("#7ec5f4", "#0000C6", "#1D88AB", "#7F00FF", "gray79", "#28652C", "#849E00", "#5CE65C"))



png("clam_spartina_tree_iqtree_accession.png", height=2100, width=2500)
p %<+% dd + 
  geom_tippoint(aes(color=genome), size=12) + 
  geom_tiplab(size=15, offset = 0.003) +
  theme(legend.text = element_text(size=45),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=c("#7ec5f4", "#0000C6", "#1D88AB", "#7F00FF", "gray79", "#28652C", "#849E00", "#5CE65C")) +
  guides(size = guide_legend(override.aes = list(size = 7))) + 
  xlim(NA, 0.15)
dev.off()

png("clam_spartina_tree_iqtree_species.png", height=2000, width=2000)
p %<+% dd + 
  geom_tippoint(aes(color=host), size=5.5) + 
  geom_tiplab(aes(label=species), size=7) +
  theme(legend.text = element_text(size=30),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=c("#7ec5f4", "gray79", "#849E00", "#7F00FF", "gray79", "#28652C", "#849E00", "#5CE65C")) +
  #scale_size_manual(values=c(6, 6, 6, 6, 6, 6, 6)) +
  guides(size = guide_legend(override.aes = list(size = 7))) + 
  xlim(NA, 0.15)
dev.off()




# load in AAI values

aai <- read_delim("all_aai.txt") %>% 
  dplyr::select(query, target, AAI_estimate) %>%
  set_colnames(c("genome1", "genome2", "AAI")) %>%
  mutate(genome1 = ifelse(startsWith(genome1, "GCA"), paste0("GCA_", str_split_i(genome1, "_", 2)), genome1),
         genome2 = ifelse(startsWith(genome2, "GCA"), paste0("GCA_", str_split_i(genome2, "_", 2)), genome2)) %>%
  mutate(genome1 = str_replace(genome1, ".fna", ""),
         genome2 = str_replace(genome2, ".fna", "")) %>%
  #filter(genome1 %in% c(plant_genomes, clam_genomes) & genome2 %in% c(plant_genomes, clam_genomes)) %>%
  filter(!(genome1 %in% c("query", "target"))) %>%
  dplyr::select(genome1, genome2, AAI) %>%
  pivot_wider(names_from = genome2,
              values_from = AAI) %>%
  column_to_rownames("genome1")



aai[aai==">90%"]<-"90"
aai <- mutate_all(aai, function(x) as.numeric(as.character(x)))




ha_values <- data.frame(genome=colnames(aai)) %>%
  mutate(Host=case_when(genome %in% clam_genomes ~"Clam",
                        genome %in% plant_genomes ~"Plant",
                        !(genome %in% c(clam_genomes, plant_genomes)) ~"GTDB")) %>%
  pull(Host)

rownames(aai) %<>% str_replace("GCF_963455295.1_xbCteDecu1.Thiodiazotropha_sp._1.1_genomic", "GCF_963455295.1")
colnames(aai) %<>% str_replace("GCF_963455295.1_xbCteDecu1.Thiodiazotropha_sp._1.1_genomic", "GCF_963455295.1")


ha <- HeatmapAnnotation(Host=ha_values, 
                        col = list(Host = c("Clam" = "purple", "Plant" = "seagreen", "GTDB" = "gray78")),
                        annotation_legend_param = list(Host = list(title = "Host"),
                                                       title_gp = gpar(fontsize = 17, fontface = "bold"),
                                                       labels_gp = gpar(fontsize = 23)),
                        gp = gpar(col = "black"),
                        show_annotation_name = FALSE)


colnames(aai) %>% unique() %>% length()

png("AAI_heatmap.png", height=1100, width=1100)
ht <- Heatmap(as.matrix(aai),
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        top_annotation = ha,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 13.5),
        heatmap_legend_param=list(title="AAI",
                                  legend_height = unit(10, "cm"),
                                  grid_height   = unit(1, "cm"),
                                  grid_width    = unit(1, "cm"),
                                  title_gp = gpar(fontsize = 17, fontface = "bold"),
                                  labels_gp = gpar(fontsize = 17)))

draw(ht, padding = unit(c(40, 25, 2, 25), "mm"))
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
  


genome_quality %>%
  filter(Completeness > 85 & Contamination < 5 & host %in% c("Clam", "Plant")) %>%
  dplyr::select(study, host, Contamination, Completeness, Contig_N50, Genome_Size) %>%
  ggplot(aes(x = host, y = Genome_Size)) + 
  geom_boxplot() 






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

