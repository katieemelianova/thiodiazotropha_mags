library(dplyr)
library(readr)
library(stringr)
library(tidytree)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vroom)



f = list.files(path = "/Users/katieemelianova/Desktop/Spartina/thiodiazotropha_mags/kofam", full.names = TRUE)

kofam <- vroom(f, delim = "\t", id = "assembly") %>%
  set_colnames(c("assembly", "significant", "gene_name", "KO", "threshold", "score", "evalue", "KO_definition")) %>%
  mutate(assembly = str_split_i(assembly, "/", 8) %>% str_replace("_kofam", ""))

kofam$assembly <- ifelse(startsWith(kofam$assembly, "GCA"), paste0("GCA_", str_split_i(kofam$assembly, "_", 2)) %>% paste0(".1"), kofam$assembly)




kofam_genome <- left_join(kofam, genome_quality, by="assembly") %>% filter(significant == "*" & study != "GTDB" & Completeness > 85)
kofam_genome <- kofam_genome %>% mutate(study = case_when(study == "China Spartina" ~ "huang",
                                          study == "Clam Giani" ~ "giani",
                                          study == "Clam Morel" ~ "morel",
                                          study == "Clam Osvatic" ~ "osvatic",
                                          study == "Ficus Spartina" ~ "ficus",
                                          study == "USA Spartina" ~ "rolando"))

kofam_genome %>% group_by(study) %>% summarise(length(unique(assembly)))



########################################
#       CLAM more than PLANT           #
########################################

kofam_clam_mt_plant <- kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  mutate(huang_percent = (huang/4) * 100,
         rolando_percent = (rolando/7) * 100,
         osvatic_percent = (osvatic/81) * 100,
         morel_percent = (morel/156) * 100,
         giani_percent = (giani/15) * 100) %>%
  filter(huang_percent < 10 & rolando_percent < 10 & osvatic_percent > 80 & giani_percent > 80 & morel > 80)

test <- kofam_clam_mt_plant %>% pull(KO_definition)
clam_genes <- kofam_genome %>% filter(KO_definition %in% test) %>% pull(gene_name)

kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  mutate(huang_percent = (huang/4) * 100,
         rolando_percent = (rolando/7) * 100,
         osvatic_percent = (osvatic/81) * 100,
         morel_percent = (morel/156) * 100,
         giani_percent = (giani/15) * 100) %>%
  filter(huang_percent < 10 & rolando_percent < 10 & osvatic_percent > 80 & giani_percent > 80 & morel > 80) %>%
  data.frame() %>%
  writexl::write_xlsx("clam_morethan_plant.xlsx")


########################################
#       PLANT more than CLAM           #
########################################

kofam_plant_mt_clam <- kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  mutate(huang_percent = (huang/4) * 100,
         rolando_percent = (rolando/7) * 100,
         osvatic_percent = (osvatic/81) * 100,
         morel_percent = (morel/156) * 100,
         giani_percent = (giani/15) * 100) %>%
  filter(huang_percent > 80 & rolando_percent > 80 & osvatic_percent < 5 & giani_percent < 5 & morel_percent < 5)


test <- kofam_plant_mt_clam %>% pull(KO_definition)
plant_genes <- kofam_genome %>% filter(KO_definition %in% test) %>% pull(gene_name)

kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  mutate(huang_percent = (huang/4) * 100,
         rolando_percent = (rolando/7) * 100,
         osvatic_percent = (osvatic/81) * 100,
         morel_percent = (morel/156) * 100,
         giani_percent = (giani/15) * 100) %>%
  filter(huang_percent > 80 & rolando_percent > 80 & osvatic_percent < 5 & giani_percent < 5 & morel_percent < 5) %>%
  data.frame() %>%
  writexl::write_xlsx("plant_morethan_clam.xlsx")





kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  filter(KO_definition == "MFS transporter, DHA2 family, tetracycline/oxytetracycline resistance protein")


# KO definition most similar to Otr which was found to be 
#"MFS transporter, DHA2 family, tetracycline/oxytetracycline resistance protein"

