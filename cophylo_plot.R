

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



################################
#    tip renaming function     #
################################


rename_tips <- function(input_tree){
  input_tree$tip.label <- ifelse(startsWith(input_tree$tip.label, "GCA"), paste0("GCA_", str_split_i(input_tree$tip.label, "_", 2)), input_tree$tip.label)
  return(input_tree)
}


#########################
#     read in tree      #
#########################

# read in trees
gtdb_tree<-ape::read.tree("core_gene_alignment_trimmed.aln.treefile")
#panaroo_tree <- ape::read.tree("core_gene_alignment_trimmed.aln.treefile")
gtdb_tree %<>% rename_tips()

files <- list.files("trees_of_interest", full.names = TRUE)
mytrees <- lapply(files, \(f) ape::read.tree(f) %>% rename_tips())
names(mytrees) <- tools::file_path_sans_ext(basename(files))




# set host species
gtdb_tree$host <- case_when(gtdb_tree$tip.label %in% sporobolus_genomes ~ "Sporobolus",
                            gtdb_tree$tip.label %in% seagrass_genomes ~ "Syringodium",
                            gtdb_tree$tip.label %in% clam_genomes ~ "Clam")

####################################################
#       root tree by most basal sedimenticola      #
####################################################

gtdb_tree <- root(gtdb_tree, "GCA_027068675.1")



####################################################
#       plot cophyloplot      #
####################################################



gene_tree <- mytrees$`cynS_2~~~cynS`
gene_tree <- mytrees$`ttrB_2`

gene_tree <- mytrees$`cheB_4~~~rpfG_4`


host <- setNames(gtdb_tree$host, gtdb_tree$tip.label)

host_cols <- c(Sporobolus = "darkgreen", Syringodium = "steelblue", Clam = "orange")

shared <- intersect(gtdb_tree$tip.label, gene_tree$tip.label)
association <- cbind(shared, shared)
co <- phytools::cophylo(gtdb_tree, gene_tree, assoc = association, rotate = TRUE)


# co$assoc holds the matched pairs in the order links are drawn
link_cols <- unname(host_cols[host[co$assoc[, 1]]])
link_cols[is.na(link_cols)] <- "grey60"   # tips with no host assigned


png("CheB.png", height=1000, width=1500)
plot(co, link.type = "curved", link.lwd = 3, link.col = link_cols,
     fsize = 1, mar = c(2, 4, 2, 4))

phytools::nodelabels.cophylo(text = gene_tree$node.label, which = "right", frame = "none", cex = 1.4)
phytools::nodelabels.cophylo(text = gtdb_tree$node.label, which = "left", frame = "none", cex = 1.4)

legend("bottom",
       legend = c(names(host_cols), "No host assigned"),
       col    = c(host_cols, "grey60"),
       lwd    = 5, lty = 2,          # lty = 2 matches the dashed links in your plot
       bty    = "n", cex = 1.8,
       horiz  = TRUE,                # one row along the bottom
       xpd    = NA, inset = c(0, -0.02))

dev.off()










##############


library(ape)

sp  <- gtdb_tree
gts <- cobq_tree




gt2 <- unroot(keep.tip(gts, intersect(gts$tip.label, sp$tip.label)))
sp2 <- unroot(keep.tip(sp, gt2$tip.label))
test <- comparePhylo(sp2, gt2, plot=TRUE)
splits_in_common <- test$messages[7] %>% str_split_i(" ", 1) %>% as.numeric()
splits_in_common/sp2$Nnode 





###########


library(ape)
library(dplyr)
library(stringr)

sp <- gtdb_tree

results <- lapply(names(mytrees), function(nm) {
  gts <- mytrees[[nm]]
  
  shared <- intersect(gts$tip.label, sp$tip.label)
  
  # Need at least 4 shared tips for an unrooted topology to have informative splits
  if (length(shared) < 4) {
    return(data.frame(tree = nm, n_shared_tips = length(shared),
                      splits_in_common = NA_real_, species_nnode = NA_real_,
                      prop_common = NA_real_))
  }
  
  gt2 <- unroot(keep.tip(gts, shared))
  sp2 <- unroot(keep.tip(sp, gt2$tip.label))
  
  test <- comparePhylo(sp2, gt2, plot = FALSE)
  
  # Pull the "N splits in common" message by pattern rather than by position
  msg <- test$messages[str_detect(test$messages, "splits in common")]
  n_common <- if (length(msg) > 0) {
    as.numeric(str_split_i(msg[1], " ", 1))
  } else {
    NA_real_
  }
  
  data.frame(tree = nm,
             n_shared_tips = length(shared),
             splits_in_common = n_common,
             species_nnode = sp2$Nnode,
             prop_common = n_common / sp2$Nnode)
})

results_df <- bind_rows(results)
results_df %>% filter(prop_common < 0.4)







