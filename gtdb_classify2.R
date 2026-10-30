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
library(stringr)



rolando_genomes <- read.table("rolando_mag_accessions.txt") %>% pull(V1)
huang_genomes <- read.table("huang_mag_accessions.txt") %>% pull(V1)
seagrass_genomes <- c("S10C43", "S6C31", "S6C37", "metabat2.111.contigs", "metabat2.466.contigs", "metabat2.680.contigs", "metabat2.898.contigs", "metabat2.939.contigs")



osvatic_genomes <- read.table("osvatic_mag_accessions.txt") %>% pull(V1)
giani_genomes <- read.table("giani_mag_accessions.txt") %>% pull(V1)
morel_genomes <- read.table("morel_mag_accessions.txt") %>% pull(V1)
petersen_clam <- c("luna_ont_bin1", "lotti_ont_bin1")
sanger_genomes <- read.table("sanger_mag_accessions.txt") %>% pull(V1)
sudo_genomes <- read.table("sudo_mag_accessions.txt") %>% pull(V1)
petersen_genomes <- read.table("petersen_mag_accessions.txt") %>% pull(V1)
PRJEB85476_genomes <- read.table("PRJEB85476_mag_accessions.txt") %>% pull(V1)
PRJEB87551_genomes <- read.table("PRJEB87551_mag_accessions.txt") %>% pull(V1)
mpi_genomes <- c("GCF_963455295.1_xbCteDecu1.Thiodiazotropha_sp._1.1_genomic")




ficus_genomes <- c("PIE_T1P1_root_Marsh_MG__3300081559_11",
                   "PIE_T1P1_root_Marsh_MG__3300081559_s3",
                   "PIE_T2P3_root_Marsh_MG__3300077085_6",
                   "PIE_T3P3_root_Marsh_MG__3300077490_s15")

sporobolus_genomes <- c(rolando_genomes, huang_genomes, ficus_genomes)

plant_genomes <- c(rolando_genomes, huang_genomes, ficus_genomes, seagrass_genomes)
clam_genomes <- c(PRJEB85476_genomes, PRJEB87551_genomes, petersen_clam, petersen_genomes, morel_genomes, giani_genomes, osvatic_genomes, sanger_genomes, sudo_genomes, mpi_genomes)


#########################
#     read in tree      #
#########################

# read in tree
#gtdb_tree<-ape::read.tree("core_gene_alignment_trimmed.aln.treefile")
gtdb_tree<-ape::read.tree("gtdbtk.bac120.user_msa.fasta.treefile")



# rename tip labels to match genome names in clam and plant genome accessions
gtdb_tree$tip.label <- ifelse(startsWith(gtdb_tree$tip.label, "GCA"), paste0("GCA_", str_split_i(gtdb_tree$tip.label, "_", 2)), gtdb_tree$tip.label)

# set host species
gtdb_tree$host <- case_when(gtdb_tree$tip.label %in% sporobolus_genomes ~ "Sporobolus",
                            gtdb_tree$tip.label %in% seagrass_genomes ~ "Syringodium",
                            gtdb_tree$tip.label %in% clam_genomes ~ "Clam")

####################################################
#       root tree by most basal sedimenticola      #
####################################################

gtdb_tree <- root(gtdb_tree, "GCA_027068675.1")


##########################
#       build tree       #
##########################
dd <- data.frame(taxa=gtdb_tree$tip.label,
                 host=gtdb_tree$host ) %>%
  mutate(host=ifelse(gtdb_tree$tip.label == "GCA_027068675.1", "Outgroup", host))

p<-ggtree(gtdb_tree, size=2, colour="gray4", branch.length = "none")




p2 <- p %<+% dd + 
  geom_tippoint(aes(color=host), size=9, stroke=2) + 
  #geom_tiplab(size=5, offset = 0.003) +
  theme(legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=c("deepskyblue1", "grey80", "darkgreen", "chartreuse1", "purple", "slateblue4", "grey80", "pink",  "chartreuse1", "darkgreen", "olivedrab3")) + 
  guides(size = guide_legend(override.aes = list(size = 7))) #+ 
  #xlim(NA, 0.17)





interest <- c("argP~~~cynR_2", "atoC_4", "atoS_2~~~atoS~~~sasA_2", "ccoP2~~~petJ", 
              "cheB_4~~~rpfG_4", "cobD", "cobD_1~~~hisC_3~~~hisC_2", "cobQ", 
              "ctpM~~~aer~~~aer_3", "cydA_1~~~appC~~~cydA_2~~~cydA", "cydB", "cydC", 
              "cydD", "cynS_2~~~cynS", "dmdC_2~~~dmdC", "fadA", "folD_2", 
              "glgA_3~~~glgA_1~~~glgA_2", "hcaB_1", "ligA_3", "ligJ", "nirQ_1~~~nirQ_3~~~nirQ_4", 
              "nodD2", "norR_2~~~norR", "phsA~~~napA_1", "qoxC~~~ctaE_1", 
              "resA_9~~~resA_8~~~resA_6~~~resA_7~~~resA_1", "rsr", "rtcA", "rtcB_2~~~rtcB_1~~~rtcB1_2", 
              "sasA_27~~~rcsC_10", "sasA_6~~~sasA_23", "slyA_2~~~slyA_1", "ssuB", 
              "torD_2~~~torD_1", "trtA_1~~~biuH", "ttrB_2", "vnfA")


interest %>% data.frame() %>% write.table("genes_of_interest", col.names = FALSE, quote=FALSE, row.names = FALSE)


gpa <- read_tsv("gene_presence_absence.Rtab") 

gpa_interest <- gpa %>% filter(Gene %in% interest) %>% tibble::column_to_rownames("Gene")
colnames(gpa_interest) <- ifelse(startsWith(colnames(gpa_interest), "GCA"), paste0("GCA_", str_split_i(colnames(gpa_interest), "_", 2)), colnames(gpa_interest))


gene_df <- gpa_interest %>% t() %>% data.frame()
gene_df <- gene_df[colSums(gene_df) > 4]
gene_df %<>%  mutate_all(as.character)


colnames(gene_df) %<>% stringr::str_split_i("\\.\\.\\.", 1)

png("tree_and_heatmap.png", width=2700, height=1700)
gheatmap(p2, gene_df,
         offset = 0.05,        # distance from the tips; must be beyond the tip labels
         width = 3.5,         # width of the heatmap relative to the tree
         colnames = TRUE,
         colnames_angle = 90,
         colnames_position = "top",
         colnames_offset_y = 0.5,
         hjust = 0,
         font.size = 15) +
  scale_fill_manual(values = c("red", "black"),
                    labels = c("Absent", "Present")) +
  scale_y_continuous(expand = expansion(add = c(0, 5))) +  # extra room above the heatmap
  coord_cartesian(clip = "off") +                           # don't clip labels at the panel edge
  theme(plot.margin = margin(t = 2, r = 1, b = 1, l = 1, unit = "cm"),
        legend.text = element_text(size=50))
dev.off()







