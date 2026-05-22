library(ape)
library(ggtree)
library(tidytree)


bac120 <- read_tsv("bac120_metadata.tsv")
bac120$accession <- substring(bac120$accession, 4)

gtdb_tree<-ape::read.tree("gtdbtk.bac120.classify.tree.1.tree")

# genomes start with a random string so remove these so I can match up with genome accession names
gtdb_tree$tip.label <- substring(gtdb_tree$tip.label, 4)

# get the ancestor of a china and usa tip, then use this node number to tree_subset out a clade to work with a smaller tree
getMRCA(gtdb_tree, c("GCA_050305635.1", "GCA_037384565.1"))
gtdb_tree_subset <- tree_subset(gtdb_tree, 32064, 14)


# use the total bac120 dataset to get the genus and if available species for each tip
gtdb_tree_subset$genus <- left_join(data.frame(accession=gtdb_tree_subset$tip.label), bac120, by="accession") %>% dplyr::select(accession, ncbi_taxonomy) %>% mutate(ncbi_taxonomy=str_split_i(ncbi_taxonomy, ";", 6)) %>% pull(ncbi_taxonomy) %>% str_replace("g__", "")
gtdb_tree_subset$species <- left_join(data.frame(accession=gtdb_tree_subset$tip.label), bac120, by="accession") %>% dplyr::select(accession, ncbi_taxonomy) %>% mutate(ncbi_taxonomy=str_split_i(ncbi_taxonomy, ";", 7)) %>% pull(ncbi_taxonomy) %>% str_replace("s__", "")



rolando_genomes <- c("GCA_037384565.1", "GCA_037384545.1", "GCA_037384465.1",
                     "GCA_037384445.1", "GCA_037384425.1", "GCA_037384385.1",
                     "GCA_037381645.1", "GCA_037382205.1", "GCA_037382185.1", "GCA_037382065.1")


huang_genomes <- c("GCA_050302455.1", "GCA_050303695.1", "GCA_050305775.1",
                   "GCA_050305635.1", "GCA_050304775.1", "GCA_050291245.1",
                   "GCA_050291065.1", "GCA_050292945.1", "GCA_050292405.1", "GCA_050292275.1")

osvatic_genomes <- c("GCA_026042565.1", "GCA_026042595.1", "GCA_026042555.1", "GCA_026042335.1", "GCA_026042345.1",
                     "GCA_026042155.1", "GCA_026042195.1",  "GCA_026041995.1", "GCA_026041825.1",  "GCA_026041815.1",
                     "GCA_026041775.1", "GCA_026041615.1", "GCA_026041535.1", "GCA_026041445.1", "GCA_016842515.1",
                     "GCA_016842625.1",  "GCA_016842965.1",  "GCA_016843295.1", "GCA_026041415.1", "GCA_026041375.1",
                     "GCA_026041365.1",  "GCA_026041355.1",  "GCA_026041325.1", "GCA_026041295.1",  "GCA_026041315.1",
                     "GCA_026041435.1", "GCA_016848925.1", "GCA_016842425.1", "GCA_016842285.1", "GCA_016842265.1",
                     "GCA_016842245.1",  "GCA_016842225.1", "GCA_016842205.1", "GCA_016842165.1",  "GCA_016842175.1",
                     "GCA_016842145.1", "GCA_016842125.1",  "GCA_016842325.1", "GCA_016842065.1", "GCA_026041475.1",
                     "GCA_026041495.1", "GCA_026041675.1",  "GCA_026041735.1", "GCA_026041835.1",  "GCA_026042075.1",
                     "GCA_026042035.1", "GCA_026041915.1", "GCA_026041985.1", "GCA_026042315.1", "GCA_026042205.1",
                     "GCA_026042475.1", "GCA_026042415.1", "GCA_026042425.1", "GCA_016842295.1",  "GCA_016842345.1",
                     "GCA_016842075.1",  "GCA_016842395.1", "GCA_016842505.1", "GCA_016842485.1", "GCA_016842465.1",
                     "GCA_016843285.1", "GCA_016843265.1", "GCA_016843245.1", "GCA_016843225.1", "GCA_016843205.1",
                     "GCA_016843185.1", "GCA_016843165.1", "GCA_016843125.1", "GCA_016843135.1", "GCA_016843105.1",
                     "GCA_016843075.1", "GCA_016843065.1",  "GCA_016843045.1", "GCA_016843025.1", "GCA_016842995.1",
                     "GCA_016842985.1", "GCA_016842945.1", "GCA_016842925.1", "GCA_016842895.1", "GCA_016842885.1",
                     "GCA_016842865.1",  "GCA_016842845.1",  "GCA_016842825.1", "GCA_016842805.1", "GCA_016842785.1",
                     "GCA_016842765.1",  "GCA_016842745.1", "GCA_016842725.1", "GCA_016842705.1", "GCA_016842685.1", 
                     "GCA_016842655.1", "GCA_016842645.1", "GCA_016842605.1", "GCA_016842585.1", "GCA_016842565.1", 
                     "GCA_016842545.1",  "GCA_026041225.1", "GCA_026041245.1", "GCA_016842445.1",  "GCA_016842385.1", "GCA_016842365.1")


# label genomes by which dataset they came from
gtdb_tree_subset$genome <- case_when(gtdb_tree_subset$tip.label %in% rolando_genomes ~ "USA",
                                     gtdb_tree_subset$tip.label %in% huang_genomes ~ "China",
                                    !(gtdb_tree_subset$tip.label) %in% c(rolando_genomes, huang_genomes) ~ "GTDB")





dd <- data.frame(taxa=gtdb_tree_subset$tip.label,
                 genome=gtdb_tree_subset$genome,
                 genus=gtdb_tree_subset$genus,
                 species=gtdb_tree_subset$species)

p<-ggtree(gtdb_tree_subset, size=0.3, colour="gray60")
p %<+% dd + 
  geom_tippoint(aes(color=genome, size=genome)) + 
  geom_tiplab(aes(label = genus)) +
  theme(legend.title = element_text(size=25),
        legend.text = element_text(size=20)) +
  scale_colour_manual(values=c("deeppink2", "gray80", "green")) +
  scale_size_manual(values=c(4, 2, 4)) +
  guides(size = guide_legend(override.aes = list(size = 7)))














# check how many of each genome got into the final tree
table(gtdb_tree$genome)



ggtree(gtdb_tree)  + geom_tiplab()




tree_all_sedi <- ggtree::read.tree("sedimenticola_asv_seqs_all tree.newick")
tree_all_sedi$range <- case_when(startsWith(tree_all_sedi$tip.label, "usa") ~ "USA",
                                 startsWith(tree_all_sedi$tip.label, "france") ~ "France")

ggtree(tree_all_sedi)  + geom_tiplab()

dd <- data.frame(taxa=tree_all_sedi$tip.label,
                 Sedimenticola=tree_all_sedi$range)
p<-ggtree(tree_all_sedi, size=0.3, colour="gray60")
p %<+% dd + geom_tippoint(aes(color=Sedimenticola, size=Sedimenticola)) + 
  theme(legend.title = element_text(size=25),
        legend.text = element_text(size=20)) +
  scale_colour_manual(values=c("deeppink2", "darkolivegreen4")) +
  scale_size_manual(values=c(4, 1.2)) +
  guides(size = guide_legend(override.aes = list(size = 7)))


