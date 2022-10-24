library(tidyverse)
try_fage <- read.csv2("data/traits/TRY.fage.csv")
sp <- read.csv2("data/annual_species_experiment.csv")
sp_names <- read.csv2("data/species_names_lifehistory.csv") %>% 
  rename(code_sp = Code_Sp)
species <- merge(sp,sp_names)

try_grime <- try_fage %>% 
  filter(SpeciesName %in% species$Species | SpeciesName == "Draba Verna") %>% 
  filter(TraitName == "Species strategy type according to Grime") %>% 
  select(SpeciesName,OrigValueStr,Reference)

write.csv2(try_grime,"outputs/data/Species_strategy_type_according_to_Grime.csv",row.names=F)
