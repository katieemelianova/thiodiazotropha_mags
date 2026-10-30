library(microViz)
library(DESeq2)
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyverse)
library(DESeq2)
library(phyloseq)
library(patchwork)

#####################################################
#           Read in RDS and sample data             # 
#####################################################

phylo_rennes <- readRDS("/Users/katieemelianova/Desktop/Spartina/JMF_results/JMF-2508-06_16S_raw_phyloseq.rds")
sample_info <- read_tsv("/Users/katieemelianova/Desktop/Spartina/Spartina2025//Rennes_Sampling.tsv") %>% 
  dplyr::select(Locality, Species, `Sample number`)


##############################################
#    filter on prevalence, set metadata      # 
##############################################

# abundance and prevalence filter
phylo_rennes <- tax_filter(
  phylo_rennes,
  min_prevalence = 10,
  prev_detection_threshold = 5,
  min_total_abundance = 5,
  min_sample_abundance = 5)

# sample names currently in format "1058 r", split this and set sample ID as the number 
# and compartment as the letter translated to the compartment name
phylo_rennes@sam_data
phylo_rennes@sam_data$User_sample_ID_number <- phylo_rennes@sam_data %>% data.frame() %>% pull(User_sample_ID) %>% str_split_i(" ", 1)
phylo_rennes@sam_data$compartment <- phylo_rennes@sam_data %>% data.frame() %>% pull(User_sample_ID) %>% str_split_i(" ", 2) %>% 
  str_replace("root", "r") %>%
  str_replace("rhizome", "z") %>%
  str_replace("soil", "s") %>%
  str_replace("r", "Root") %>%
  str_replace("z", "Rhizome") %>%
  str_replace("s", "Rhizosphere") %>% 
  replace_na('Unknown')

# set locality and species as metadata slots
phylo_rennes@sam_data$Locality <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Locality) %>% replace_na('Unknown')
phylo_rennes@sam_data$Species <- phylo_rennes@sam_data %>% data.frame() %>% left_join(sample_info, by = c("User_sample_ID_number" = "Sample number")) %>% pull(Species) %>% replace_na('Unknown')
phylo_rennes@sam_data$Species %<>% 
  str_replace("Spartina alternifllora", "Sporobolus alterniflorus") %>%
  str_replace("Spartina maritima", "Sporobolus maritimus") %>%
  str_replace("Spartina anglica", "Sporobolus anglicus")



#############################################################################################
#   filter plastids & low coverage, rarefy and transform abundances to be proportional      # 
#############################################################################################

# remove chloroplast and mitochonria and unknowns
phylo_rennes <- subset_samples(phylo_rennes, compartment != "Unknown")
phylo_rennes <- subset_samples(phylo_rennes, Species != "Unknown")
phylo_rennes <- subset_taxa(phylo_rennes, !(Family %in% c("Mitochondria", "Chloroplast"))) %>% subset_taxa(!(Order %in% c("Mitochondria", "Chloroplast")))

# rmeove samples with less than 2000 total reads
phylo_rennes <- prune_samples(sample_sums(phylo_rennes) >= 2000, phylo_rennes)

# rarefy to lowest depth sample
set.seed = 1
phylo_rennes_rarefied <- rarefy_even_depth(phylo_rennes, sample.size = min(sample_sums(phylo_rennes)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# transform rarefied counts to proportional
phylo_rennes_rarefied_prop <- transform_sample_counts(phylo_rennes_rarefied, function(otu) otu/sum(otu))


phylo_rennes_rarefied %>% 
  subset_taxa(Genus == "Candidatus Thiodiazotropha") %>%
  plot_richness(x = "JMF_sample_ID", measures = c("Observed")) + 
  geom_boxplot(aes(fill = Species)) +
  facet_wrap(~Species, scales="free") +
  ylim(0,10)

phylo_rennes_rarefied_prop %>%
  subset_taxa(Genus == "Candidatus Thiodiazotropha") %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(~sample_Species+compartment, ncol=3, scales="free") +
  ylim(0,.25) +
  ylab("Relative Abundance")


