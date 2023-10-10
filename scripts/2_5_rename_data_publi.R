library(tidyverse)

# Rename trait datasets ####

# Traits averaged at the management regime level
MEAN_multivar <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv") %>% 
  select(-Form) %>% 
  select(-c(Nb_Lf,Nb_Lflet,Sple_Area,Sple_FM,Sple_DM,L_Length, L_Width,LTmes,Hveg,Dmax,Dmin,Flo,Mat_Per,Mat,LPC )) %>% 
  select(species,treatment,Hrepro,Ldelta13C) %>% 
  rename(Hrepro_multivar = Hrepro, Ldelta13C_multivar = Ldelta13C)

MEAN_univar <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM.csv") %>% 
  select(-Form) %>% 
  select(-c(Nb_Lf,Nb_Lflet,Sple_Area,Sple_FM,Sple_DM,L_Length, L_Width,LTmes,Hveg,Dmax,Dmin,Flo,Mat_Per,Mat,LPC ))

MEAN_pop <- merge(MEAN_univar,MEAN_multivar)

write.csv(MEAN_pop,"outputs/data/traits_management_level.csv",row.names=F)

# Traits averaged at the site level
MEAN_site <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_site_level.csv")%>% 
  select(-Form) %>% 
  select(-c(Nb_Lf,Nb_Lflet,Sple_Area,Sple_FM,Sple_DM,L_Length, L_Width,LTmes,Hveg,Dmax,Dmin,Flo,Mat_Per,Mat,LPC ))
write.csv(MEAN_site,"outputs/data/traits_site_level.csv",row.names=F)

