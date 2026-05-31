

f = list.files(path = "/Users/katieemelianova/Desktop/Spartina/thiodiazotropha_mags/kofam", full.names = TRUE)

kofam <- vroom(f, delim = "\t", id = "assembly") %>%
  set_colnames(c("assembly", "significant", "gene_name", "KO", "threshold", "score", "evalue", "KO_definition")) %>%
  mutate(assembly = str_split_i(assembly, "/", 8) %>% str_replace("_kofam", "") %>% paste0(".1"))


kofam_genome <- left_join(kofam, genome_all, by="assembly") 
#%>% filter(significant == "*")


kofam_genome %>% group_by(study) %>% summarise(length(unique(assembly)))


kofam_genome$KO_definition


kofam_genome %>%
  group_by(KO_definition, study) %>%
  summarise(test=n()) %>%
  pivot_wider(names_from = study, values_from = test, values_fill = 0) %>%
  mutate(gtdb_percent = (gtdb/13) * 100,
         huang_percent = (huang/10) * 100,
         rolando_percent = (rolando/10) * 100,
         osvatic_percent = (osvatic/3) * 100) %>%
  filter(huang_percent < 10 & rolando_percent < 10 & osvatic_percent > 70) %>%
  data.frame()



