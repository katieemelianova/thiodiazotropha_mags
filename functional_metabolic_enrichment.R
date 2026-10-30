

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(vroom)

metadata <- read_delim("/Users/katieemelianova/Desktop/Spartina/thiodiazotropha_mags/anvio_metadata.tsv") %>% dplyr::rename("genome_name" = "genome" )
files <- list.files("/Users/katieemelianova/Desktop/Spartina/thiodiazotropha_mags/anvio_estimate_metabolism/", pattern = "*_modules.txt", full.names = TRUE)
metabolism <- vroom(files)
metabolism_sub <- metabolism %>% dplyr::select(module, module_name, genome_name, module_subcategory, stepwise_module_completeness) %>% left_join(metadata, by="genome_name")


metabolism_sub %>% 
  group_by(module_name, host) %>%
  summarise(avg=mean(stepwise_module_completeness)) %>%
  ungroup %>%
  pivot_wider(names_from = "host", values_from="avg") %>%
  mutate(delta = abs(Clam-Plant)) %>%
  filter(delta > 0.3)


read_tsv("/Users/katieemelianova/Desktop/Spartina/thiodiazotropha_mags/test.txt") %>%
  filter(unadjusted_p_value < 0.01 & p_Clam > 0.9 & p_Plant < 0.1) %>%
  pull(COG24_FUNCTION)


read_tsv("/Users/katieemelianova/Desktop/Spartina/thiodiazotropha_mags/test.txt") %>%
  filter(unadjusted_p_value < 0.0001 & p_Plant > 0.9 & p_Clam < 0.1) %>%
  pull(COG24_FUNCTION)



