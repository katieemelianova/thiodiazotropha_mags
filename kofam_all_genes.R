library(dplyr)
library(readr)
library(stringr)
library(tidytree)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vroom)



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


f = list.files(path = "/Users/katieemelianova/Desktop/Spartina/thiodiazotropha_mags/kofam", full.names = TRUE)

kofam <- vroom(f, delim = "\t", id = "assembly") %>%
  set_colnames(c("assembly", "significant", "gene_name", "KO", "threshold", "score", "evalue", "KO_definition")) %>%
  mutate(assembly = str_split_i(assembly, "/", 8) %>% str_replace("_kofam", ""))

kofam$assembly <- ifelse(startsWith(kofam$assembly, "GCA"), paste0("GCA_", str_split_i(kofam$assembly, "_", 2)), kofam$assembly)

kofam$study <- case_when(kofam$assembly %in% rolando_genomes ~ "Spartina Rolando",
                         kofam$assembly %in% huang_genomes ~ "Spartina Huang",
                         kofam$assembly %in% osvatic_genomes ~ "Clam Osvatic",
                         kofam$assembly %in% giani_genomes ~ "Clam Giani",
                         kofam$assembly %in% morel_genomes ~ "Clam Morel",
                         kofam$assembly %in% petersen_clam ~ "Clam Petersen",
                         kofam$assembly %in% ficus_genomes ~"Spartina Ficus",
                         !(kofam$assembly) %in% c(rolando_genomes, huang_genomes, osvatic_genomes) ~ "GTDB")

kofam$host <- case_when(kofam$assembly %in% rolando_genomes ~ "Plant",
                         kofam$assembly %in% huang_genomes ~ "Plant",
                         kofam$assembly %in% osvatic_genomes ~ "Clam",
                         kofam$assembly %in% giani_genomes ~ "Clam",
                         kofam$assembly %in% morel_genomes ~ "Clam",
                         kofam$assembly %in% petersen_clam ~ "Clam",
                         kofam$assembly %in% ficus_genomes ~"Plant",
                         !(kofam$assembly) %in% c(rolando_genomes, huang_genomes, osvatic_genomes) ~ "GTDB")

kofam <- kofam %>% filter(significant == "*" & study != "GTDB")
kofam %>% group_by(study) %>% summarise(length(unique(assembly)))
kofam %>% group_by(host) %>% summarise(length(unique(assembly)))


########################################
#       CLAM more than PLANT           #
########################################


kofam_clam_mt_plant <- distinct_at(kofam, vars(assembly, KO_definition), .keep_all=TRUE) %>%
  group_by(KO_definition, host) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = host, values_from = test, values_fill = 0) %>%
  mutate(clam_percent = (Clam/11) * 100,
         plant_percent = (Plant/9) * 100) %>%
  filter(clam_percent > 80 & plant_percent < 20) 




########################################
#       PLANT more than CLAM           #
########################################

kofam_plant_mt_clam <- distinct_at(kofam, vars(assembly, KO_definition), .keep_all=TRUE) %>%
  group_by(KO_definition, host) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = host, values_from = test, values_fill = 0) %>%
  mutate(clam_percent = (Clam/11) * 100,
         plant_percent = (Plant/9) * 100) %>%
  filter(clam_percent < 20 & plant_percent > 80)

kofam_clam_mt_plant %>% writexl::write_xlsx("clam_morethan_plant_kofam.xlsx")
kofam_plant_mt_clam %>% writexl::write_xlsx("plant_morethan_clam_kofam.xlsx")

