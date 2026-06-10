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

kofam %<>% mutate(assembly=ifelse(startsWith(assembly, "GCA"), paste0(assembly, ".1"), assembly))



kofam_genome <- left_join(kofam, genome_quality, by="assembly") %>% filter(significant == "*" & study != "GTDB" & Completeness > 85)

# check how many of each there are to calculate percentage in next step
kofam_genome %>% group_by(study) %>% summarise(length(unique(assembly)))

########################################
#       CLAM more than PLANT           #
########################################

kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  mutate(huang_percent = (`China Spartina`/4) * 100,
         rolando_percent = (`USA Spartina`/7) * 100,
         osvatic_percent = (`Clam Osvatic`/66) * 100,
         morel_percent = (`Clam Morel`/114) * 100,
         giani_percent = (`Clam Giani`/14) * 100,
         ficus_percent = (`Ficus Spartina`/4) * 100) %>%
  filter(huang_percent < 10 & rolando_percent < 10 & ficus_percent < 10 & osvatic_percent > 80 & giani_percent > 80 & morel_percent > 80) %>%
  data.frame() %>%
  writexl::write_xlsx("clam_morethan_plant.xlsx")


########################################
#       PLANT more than CLAM           #
########################################

kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  mutate(huang_percent = (`China Spartina`/4) * 100,
         rolando_percent = (`USA Spartina`/7) * 100,
         osvatic_percent = (`Clam Osvatic`/66) * 100,
         morel_percent = (`Clam Morel`/114) * 100,
         giani_percent = (`Clam Giani`/14) * 100,
         ficus_percent = (`Ficus Spartina`/4) * 100) %>%
  filter(huang_percent > 80 & rolando_percent > 80 & ficus_percent > 80 & osvatic_percent < 5 & giani_percent < 5 & morel_percent < 5) %>%
  data.frame() %>%
  writexl::write_xlsx("plant_morethan_clam.xlsx")





kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  filter(KO_definition == "MFS transporter, DHA2 family, tetracycline/oxytetracycline resistance protein")


# KO definition most similar to Otr which was found to be 
#"MFS transporter, DHA2 family, tetracycline/oxytetracycline resistance protein"

