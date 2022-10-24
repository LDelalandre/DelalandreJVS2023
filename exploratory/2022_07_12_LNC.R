library(tidyverse)
source("scripts/Data_traits.R")

LeafCN_pour_Lila <- LeafCN %>% 
  filter(LifeForm1=="The") %>% 
  select(Species,Plot,Treatment,Code_Sp,Treatment,LNC,measurementDeterminedBy)

write.csv2(LeafCN_pour_Lila,"outputs/data/LeafCN_in_situ.csv")
