library(ape)
library(ggtree)
library(dplyr)
library(readr)
library(stringr)
library(tidytree)
library(tidyr)
library(ggplot2)


#################################
#     specify genome names      #
#################################


rolando_genomes <- read.table("rolando_mag_accessions.txt") %>% pull(V1)
huang_genomes <- read.table("huang_mag_accessions.txt") %>% pull(V1)
osvatic_genomes <- read.table("osvatic_mag_accessions.txt") %>% pull(V1)
giani_genomes <- read.table("giani_mag_accessions.txt") %>% pull(V1)
morel_genomes <- read.table("morel_mag_accessions.txt") %>% pull(V1)
petersen_clam <- c("luna_ont_bin1", "lotti_ont_bin1")




#################################################
#     read in tree and get smaller subtree      #
#################################################

bac120 <- read_tsv("bac120_metadata.tsv")
bac120$accession <- substring(bac120$accession, 4)

# read in tree from gtdbtk
gtdb_tree<-ape::read.tree("gtdbtk.bac120.classify.tree.1.tree")
gtdb_tree_subset <- gtdb_tree


#gtdb_tree_subset$tip.label <- ifelse(startsWith(gtdb_tree_subset$tip.label, "GCA"), paste0("GCA_", str_split_i(gtdb_tree_subset$tip.label, "_", 2)), gtdb_tree_subset$tip.label)

# genomes start with a random string so remove these so I can match up with genome accession names
 gtdb_tree$tip.label <- substring(gtdb_tree$tip.label, 4)

# get the ancestor of clam, china and usa tips, then use this node number to tree_subset out a clade to work with a smaller tree, taking a few nodes up for broader context
getMRCA(gtdb_tree, c("GCA_050305635.1", "GCA_037384565.1"))
genomes_in_tree <- c(rolando_genomes, huang_genomes, osvatic_genomes, morel_genomes, giani_genomes)[c(rolando_genomes, huang_genomes, osvatic_genomes, morel_genomes, giani_genomes) %in% gtdb_tree$tip.label]
gtdb_tree_subset <- tree_subset(gtdb_tree, getMRCA(gtdb_tree, genomes_in_tree), 7)



#################################################
#         annotate taxonomy for tree            #
#################################################

# use the total bac120 dataset to get the genus and if available species for each tip
gtdb_tree_subset$genus <- left_join(data.frame(accession=gtdb_tree_subset$tip.label), bac120, by="accession") %>% dplyr::select(accession, ncbi_taxonomy) %>% mutate(ncbi_taxonomy=str_split_i(ncbi_taxonomy, ";", 6)) %>% pull(ncbi_taxonomy) %>% str_replace("g__", "")
gtdb_tree_subset$species <- left_join(data.frame(accession=gtdb_tree_subset$tip.label), bac120, by="accession") %>% dplyr::select(accession, ncbi_taxonomy) %>% mutate(ncbi_taxonomy=str_split_i(ncbi_taxonomy, ";", 7)) %>% pull(ncbi_taxonomy) %>% str_replace("s__", "")
gtdb_tree_subset$species <- ifelse(gtdb_tree_subset$species == "", gtdb_tree_subset$genus, gtdb_tree_subset$species)




# label genomes by which dataset they came from
gtdb_tree_subset$genome <- case_when(gtdb_tree_subset$tip.label %in% rolando_genomes ~ "USA Spartina",
                                     gtdb_tree_subset$tip.label %in% huang_genomes ~ "China Spartina",
                                     gtdb_tree_subset$tip.label %in% osvatic_genomes ~ "Clam Osvatic",
                                     gtdb_tree_subset$tip.label %in% giani_genomes ~ "Clam Giani",
                                     gtdb_tree_subset$tip.label %in% morel_genomes ~ "Clam Morel",
                                     gtdb_tree_subset$tip.label %in% petersen_clam ~ "Clam Petersen",
                                    !(gtdb_tree_subset$tip.label) %in% c(rolando_genomes, huang_genomes, osvatic_genomes) ~ "GTDB")



##########################################################
#       get DTDB genome names to download on lisc        #
##########################################################

gtdb_tree_subset$tip.label[gtdb_tree_subset$genome == "GTDB"] %>% write.table("GTDB_subtree_mag_acessions.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#################################
#         build tree            #
#################################

dd <- data.frame(taxa=gtdb_tree_subset$tip.label,
                 genome=gtdb_tree_subset$genome,
                 genus=gtdb_tree_subset$genus,
                 species=gtdb_tree_subset$species)

p<-ggtree(gtdb_tree_subset, size=0.3, colour="gray60")


png("clam_spartina_tree.png", height=1500, width=1500)
p %<+% dd + 
  geom_tippoint(aes(color=genome, size=genome)) + 
  geom_tiplab(aes(label = species), size=4) +
  theme(legend.text = element_text(size=20),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=c("seagreen", "dodgerblue", "goldenrod1", "purple", "gray70", "lightgreen")) +
  scale_size_manual(values=c(6, 6, 6, 6, 6, 6, 6)) +
  guides(size = guide_legend(override.aes = list(size = 7))) + 
  xlim(NA, 0.15)
dev.off()








###########################################################################################
#     get tip labels from subset tree and from my inout genomes to subset big MSA         #
###########################################################################################


taxonomy <- rbind(read_tsv("giani_taxonomy.txt", skip=1) %>% rename("Assembly"=`# Assembly`) %>% mutate(project="giani"),
      read_tsv("morel_taxonomy.txt", skip=1) %>% rename("Assembly"=`# Assembly`) %>% mutate(project="morel"),
      read_tsv("osvatic_taxonomy.txt", skip=1) %>% rename("Assembly"=`# Assembly`) %>% mutate(project="osvatic"))




taxonomy$Taxonomy




