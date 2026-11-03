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
gtdb_tree$tip.label.view <- ifelse(startsWith(gtdb_tree$tip.label, "GCA"), paste0("GCA_", str_split_i(gtdb_tree$tip.label, "_", 2)), gtdb_tree$tip.label)


####################################################
#       root tree by most basal sedimenticola      #
####################################################

gtdb_tree <- root(gtdb_tree, "GCA_027068675_1_ASM2706867v1_genomic")


##########################
#       build tree       #
##########################
dd <- data.frame(taxa=gtdb_tree$tip.label,
                 host=gtdb_tree$host ) %>%
  mutate(host=ifelse(gtdb_tree$tip.label == "GCA_027068675_1_ASM2706867v1_genomic", "Outgroup", host))



p<-ggtree(gtdb_tree, size=2, colour="gray4", branch.length = "none")



p2 <- p %<+% dd + 
  geom_tippoint(aes(color=host), size=9, stroke=2) + 
  #geom_tiplab(size=5, offset = 0.003) +
  theme(legend.text = element_text(size=10),
        legend.title = element_blank(),
        plot.margin = margin(1,5,1,5, "cm")) +
  scale_colour_manual(values=c("#0A6CFF", "grey80", "#00D1B2", "chartreuse1", "purple", "slateblue4", "grey80", "pink",  "chartreuse1", "darkgreen", "olivedrab3")) + 
  guides(size = guide_legend(override.aes = list(size = 7))) #+ 
  #xlim(NA, 0.17)


gpa <- read_tsv("gene_clusters_long.txt")
result_annotated <- read_tsv("result_annotated.tsv")
interest <- result_annotated %>% filter(`fq*ep` < 0.001) %>% pull(Gene)


gpa_interest <- gpa %>% filter(Gene %in% interest) %>% tibble::column_to_rownames("Gene")

gene_df <- gpa_interest %>% t() %>% data.frame()
gene_df <- gene_df[colSums(gene_df) > 4]
gene_df %<>%  mutate_all(as.character)


# instead og gene cluster IDs have gene annotations
colnames(gene_df) <- data.frame(Gene=colnames(gene_df)) %>% left_join(result_annotated, by="Gene") %>% pull(COG24_FUNCTION)
gene_df <- gene_df[, !is.na(colnames(gene_df)), drop = FALSE]

# hierarchically cluster the matrix columsn so it plots more similar clusters together
colnames(gene_df) <- make.unique(colnames(gene_df), sep = "_")
col_dist <- dist(t(gene_df), method = "binary")
col_hc <- hclust(col_dist, method = "average")   # or "ward.D2", "complete"
gene_df <- gene_df[, col_hc$order, drop = FALSE]



png("tree_and_heatmap_anvio.png", width=3500, height=1700)
gheatmap(p2, gene_df,
         offset = 0.7,         # distance from the tips; must be beyond the tip labels
         width = 3.5,         # width of the heatmap relative to the tree
         colnames = TRUE,
         colnames_angle = 40,
         colnames_position = "top",
         colnames_offset_y = 0.5,
         hjust = 0,
         font.size = 6) +
  scale_fill_manual(values = c("#FBEFB8", "#A63D2F"),
                    labels = c("Absent", "Present")) +
  scale_y_continuous(expand = expansion(add = c(0, 5))) +  # extra room above the heatmap
  coord_cartesian(clip = "off") +                           # don't clip labels at the panel edge
  scale_y_continuous(expand = expansion(add = c(0, 25))) +   # more room at the top
  theme(plot.margin = margin(t = 2, r = 8, b = 1, l = 1, unit = "cm"),
        legend.text = element_text(size=50))
dev.off()







