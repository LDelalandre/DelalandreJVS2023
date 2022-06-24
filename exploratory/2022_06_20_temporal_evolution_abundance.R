
# Question : les variations interannuelles d'abondance sont-elles plus fortes chez les annuelles ?
# Cela corrèle-t-il avec les variations environnementales ?
# Fertilized treatment ####
ABUNDANCE <- read.csv2("data/abundance/pooled_abundance_data.csv") %>%  
  mutate(treatment = case_when(grepl("C",paddock) ~"Fer",
                               grepl("P",paddock) ~"Nat",
                               grepl("N",paddock) ~"Nat",
                               grepl("T",paddock) ~"Tem"))


ABUNDANCE
