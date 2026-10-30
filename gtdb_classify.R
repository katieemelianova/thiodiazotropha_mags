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
#       read in metadata        #
#################################

anvio_metadata <- read_tsv("anvio_metadata.txt")

#################################################
#     read in tree and get smaller subtree      #
#################################################

gtdb_tree<-ape::read.tree("core-sequences.aln.fasta.treefile")

#################################################
#         annotate taxonomy for tree            #
#################################################

# label genomes by which dataset and host they came from

gtdb_tree$host <- left_join(data.frame(genome=gtdb_tree$tip.label), anvio_metadata) %>% pull(host)
gtdb_tree$host_species <- left_join(data.frame(genome=gtdb_tree$tip.label), anvio_metadata) %>% pull(host_species)
gtdb_tree$genome_name <- left_join(data.frame(genome=gtdb_tree$tip.label), anvio_metadata) %>% pull(genome_name) 

anvio_metadata %>% filter(genome == "GCA_027068675_1_ASM2706867v1_genomic")

####################################################
#       root tree by most basal sedimenticola      #
####################################################

gtdb_tree <- root(gtdb_tree, "GCA_027068675_1_ASM2706867v1_genomic")

#################################
#         build tree            #
#################################

dd <- data.frame(taxa=gtdb_tree$tip.label,
                 genome=gtdb_tree$genome_name,
                 host=gtdb_tree$host,
                 host_species=gtdb_tree$host_species)


p<-ggtree(gtdb_tree, size=1, colour="gray60")


p %<+% dd + 
  geom_tippoint(aes(color=genome), size=2, stroke=2) + 
  geom_tiplab(size=5, offset = 0.003) +
  scale_colour_manual(values=c("deepskyblue1", "steelblue4", "blue2", "cornflowerblue", "purple", "slateblue4", "grey80", "pink",  "chartreuse1", "darkgreen", "olivedrab3", "red", "brown"))+
  xlim(0, max(p$data$x) * 1.5) +
  geom_nodelab(size = 3, hjust = 1.2, vjust = -0.5)

  



colours_plot_accession <- c("darkblue", "#4B3FD1", "#1F4FD8", "#4682B4", "#0091FF", "deepskyblue1", "lightblue2", "grey75",  "palegreen", "#00C853", "#6B8E23", "forestgreen")
colours_plot_host<- c("cornflowerblue", "grey", "palegreen", "olivedrab")


png("clam_spartina_tree_iqtree_accession.png", height=3000, width=2800)
p %<+% dd + 
  geom_tippoint(aes(color=genome), size=10, stroke=4) + 
  geom_tiplab(size=15, offset = 0.003) +
  geom_nodelab(size = 8, hjust = 1.2, vjust = -0.5) +
  theme(legend.text = element_text(size=45),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=colours_plot_accession) + 
  guides(size = guide_legend(override.aes = list(size = 7))) + 
  xlim(NA, 1.1)
dev.off()


png("clam_spartina_tree_iqtree_host.png", height=2200, width=3000)
p %<+% dd + 
  geom_tippoint(aes(color=host_species), size=10, stroke=4) + 
  geom_tiplab(size=15, offset = 0.003) +
  theme(legend.text = element_text(size=45),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=colours_plot_host) + 
  guides(size = guide_legend(override.aes = list(size = 7))) + 
  xlim(NA, 0.16)
dev.off()




png("clam_spartina_tree_iqtree_host.png", height=1300, width=1000)
p %<+% dd + 
  geom_tippoint(aes(color=host), size=9, stroke=0.5) + 
  theme(legend.text = element_text(size=55),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(
    values = c("Clam" = "blue2", "Outgroup" = "gray79", "Plant" = "olivedrab3"),
    breaks = c("Clam", "Plant", "Outgroup")   # legend order
  )


  #scale_shape_manual(values=c(15, 16, 17, 18, 24, 9)) +
  #guides(size = guide_legend(override.aes = list(size = 7))) + 
  #xlim(NA, 0.15)
dev.off()



#################################
#         draw AAI heatmap      #
#################################


# load in AAI values

aai <- read_delim("all_aai.txt", col_names = FALSE) %>% 
  set_colnames(c("query", "target", "avg_jacc_sim", "jacc_SD", "num_shared_SCPs", "poss_shared_SCPs", "AAI_estimate")) %>%
  dplyr::select(query, target, AAI_estimate) %>%
  set_colnames(c("genome1", "genome2", "AAI")) %>%
  mutate(genome1 = ifelse(startsWith(genome1, "GCA"), paste0("GCA_", str_split_i(genome1, "_", 2)), genome1),
         genome2 = ifelse(startsWith(genome2, "GCA"), paste0("GCA_", str_split_i(genome2, "_", 2)), genome2)) %>%
  mutate(genome1 = str_replace(genome1, ".fna", ""),
         genome2 = str_replace(genome2, ".fna", "")) %>%
  mutate(genome1 = str_replace(genome1, ".fa", ""),
         genome2 = str_replace(genome2, ".fa", "")) %>%
  #filter(genome1 %in% c(plant_genomes, clam_genomes) & genome2 %in% c(plant_genomes, clam_genomes)) %>%
  filter(!(genome1 %in% c("query", "target"))) %>%
  dplyr::select(genome1, genome2, AAI) %>%
  pivot_wider(names_from = genome2,
              values_from = AAI) %>%
  tibble::column_to_rownames("genome1")

aai[aai==">90%"]<-"90"
aai <- mutate_all(aai, function(x) as.numeric(as.character(x)))

ha_values <- data.frame(genome=colnames(aai)) %>%
  mutate(Host=case_when(genome %in% clam_genomes ~"Clam",
                        genome %in% plant_genomes ~"Plant",
                        genome %in% c("GCA_027068675.1") ~"Outgroup")) %>%
  pull(Host)

rownames(aai) %<>% str_replace("GCF_963455295.1_xbCteDecu1.Thiodiazotropha_sp._1.1_genomic", "GCF_963455295.1")
colnames(aai) %<>% str_replace("GCF_963455295.1_xbCteDecu1.Thiodiazotropha_sp._1.1_genomic", "GCF_963455295.1")

ha <- HeatmapAnnotation(Host=ha_values, 
                        col = list(Host = c("Clam" = "purple", "Plant" = "seagreen", "Outgroup" = "gray78")),
                        simple_anno_size = unit(1.5, "cm"),
                        annotation_legend_param = list(Host = list(title = "Host"),
                                                       title_gp = gpar(fontsize = 30, fontface = "bold"),
                                                       labels_gp = gpar(fontsize = 30)),
                        gp = gpar(col = "black"),
                        show_annotation_name = FALSE)


colnames(aai) %>% unique() %>% length()

png("AAI_heatmap.png", height=900, width=1100)
ht <- Heatmap(as.matrix(aai),
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        top_annotation = ha,
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_names_gp = gpar(fontsize = 13.5),
        heatmap_legend_param=list(title="AAI",
                                  legend_height = unit(10, "cm"),
                                  grid_height   = unit(1, "cm"),
                                  grid_width    = unit(1, "cm"),
                                  title_gp = gpar(fontsize = 17, fontface = "bold"),
                                  labels_gp = gpar(fontsize = 17)))

draw(ht, padding = unit(c(40, 25, 2, 25), "mm"))
dev.off()


############################################
#         draw genome length boxplots      #
############################################

genome_quality <- read_tsv("quality_report.tsv") %>% 
  dplyr::rename("assembly"="Name") %>%
  dplyr::select(c(assembly, Contamination, Contig_N50, Genome_Size, Completeness)) %>%
  mutate(assembly = ifelse(startsWith(assembly, "GCA"), paste0("GCA_", str_split_i(assembly, "_", 2)), assembly)) %>%
  left_join(data.frame(assembly=gtdb_tree$tip.label, host=gtdb_tree$host, study=gtdb_tree$genome), by="assembly")

png("genome_size.png", height = 500, width=500)
genome_quality %>%
  filter(Completeness > 85 & Contamination < 5 & host %in% c("Clam", "Plant")) %>%
  dplyr::select(study, host, Contamination, Completeness, Contig_N50, Genome_Size) %>%
  ggplot(aes(x = host, y = Genome_Size/1000000, fill=host)) + 
  geom_boxplot() +
  theme(axis.text = element_text(size=30),
        axis.title=element_text(size=30),
        legend.position="none") +
  scale_fill_manual(values=c("dodgerblue", "seagreen")) +
  ylab("Genome Size (Mb)") + 
  xlab("Host")
dev.off()





####################################
####################################
#       END OF FORMAL PLOTS        #
####################################
####################################





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

######################################################################################################
#          alternative to anvio metadata whihc incudes the clade info, not sure if willl use         #
######################################################################################################

data.frame(genome=gtdb_tree$genome_name) %>%
  mutate(genome_join = ifelse(startsWith(genome, "GCA"), paste0(str_split_i(genome, ".1_", 1), ".1"), gtdb_tree$genome_name)) %>%
  dplyr::select(-genome_join) %>%
  mutate(genome=ifelse(startsWith(genome, "GCA"), genome %>% str_replace("\\.", "_"), genome),
         host=gtdb_tree$host,
         genome_source=gtdb_tree$genome) %>%
  mutate(host=ifelse(is.na(host), "Unknown", host),
         genome=str_replace(genome, "__", "_")) %>%
  write_tsv("anvio_metadata.tsv")

