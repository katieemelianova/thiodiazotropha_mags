

f = list.files(path = "/Users/katieemelianova/Desktop/Spartina/thiodiazotropha_mags/kofam", full.names = TRUE)

kofam <- vroom(f, delim = "\t", id = "assembly") %>%
  set_colnames(c("assembly", "significant", "gene_name", "KO", "threshold", "score", "evalue", "KO_definition")) %>%
  mutate(assembly = str_split_i(assembly, "/", 8) %>% str_replace("_kofam", "") %>% paste0(".1"))


kofam_genome <- left_join(kofam, genome_all, by="assembly") %>% filter(significant == "*" & study != "gtdb")


kofam_genome %>% group_by(study) %>% summarise(length(unique(assembly)))



########################################
#       CLAM more than PLANT           #
########################################

kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  mutate(huang_percent = (huang/10) * 100,
         rolando_percent = (rolando/10) * 100,
         osvatic_percent = (osvatic/29) * 100,
         giani_percent = (giani/3) * 100) %>%
  filter(huang_percent < 10 & rolando_percent < 10 & osvatic_percent > 70 & giani_percent > 70) %>%
  data.frame()


########################################
#       PLANT more than CLAM           #
########################################

kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  mutate(huang_percent = (huang/10) * 100,
         rolando_percent = (rolando/10) * 100,
         osvatic_percent = (osvatic/29) * 100,
         giani_percent = (giani/3) * 100) %>%
  filter(huang_percent > 80 & rolando_percent > 80 & osvatic_percent < 5 & giani_percent < 5) %>%
  data.frame()





kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  filter(KO_definition == "MFS transporter, DHA2 family, tetracycline/oxytetracycline resistance protein")


# KO definition most similar to Otr which was found to be 
#"MFS transporter, DHA2 family, tetracycline/oxytetracycline resistance protein"

