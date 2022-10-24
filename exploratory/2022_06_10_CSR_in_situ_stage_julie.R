library(tidyverse)
source("scripts/Data_traits.R")
ann_sp <- read.csv2("data/annual_species_experiment.csv") 

LM_csr_julie <- LeafMorpho %>% 
  filter(LifeForm1 == "The") %>% 
  # filter(Code_Sp %in% ann_sp$code_sp) %>%
  # filter(measurementDeterminedBy=="Léo Delalandre") %>%
  filter(Treatment%in% c("Nat_Sab","Fer_Clc")) %>% 
  mutate(origin = if_else(Treatment == "Nat_Sab","Nat","Fer")) %>% 
  select(Species,Code_Sp,LifeForm1,origin,Plot,measurementDeterminedBy,Rep,L_Area,LDMC,SLA) %>%
  rename(species = Species, code_sp = Code_Sp,rep = Rep) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA))

write.csv2(LM_csr_julie,"outputs/data/data_for_csr_julie_annuals.csv",row.names = F)


sp_in_situ <- LM_csr_julie %>% 
  pull(code_sp) %>% 
  unique()

setdiff(ann_sp%>% 
            pull(code_sp),sp_in_situ)

setdiff(sp_in_situ,ann_sp%>% 
          pull(code_sp))

