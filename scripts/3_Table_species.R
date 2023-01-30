library(tidyverse)

# trait measurement
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv") %>%
  filter(!(species== "Geranium dissectum - pétiole")) %>% 
  mutate(log_LA = log(L_Area)) %>% 
  unique()

table_MEAN <- MEAN %>% 
  select(species,LifeForm1,LifeHistory,treatment) %>% 
  mutate(measured = 1) %>% 
  spread(key = treatment,value = measured) %>% 
  replace(is.na(.),0) %>% 
  rename(Intensive = Fer, Extensive = Nat) %>% 
  arrange(LifeHistory,LifeForm1,species)

write.csv2(table_MEAN,"outputs/data/Table_species.csv")




# trait measurement

# LeafMorpho
# LeafCN
# Leaf13C
# Biovolume
# Pheno
# Seed


leo <- Seed %>% 
  filter(measurementDeterminedBy == "Léo Delalandre") %>% 
  select(Species,Code_Sp,LifeForm1,Treatment) %>% 
  unique()
dim(leo)
leo %>% select(-Treatment) %>% unique() %>% dim()

fer <- leo %>% 
  filter(Treatment=="Fer_Clc") %>% 
  pull(Code_Sp)
nat <- leo %>% 
  filter(Treatment=="Nat_Sab") %>% 
  pull(Code_Sp)
intersect(fer,nat)
