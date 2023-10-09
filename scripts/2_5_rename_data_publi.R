library(tidyverse)

# Rename trait datasets ####
MEAN_multivar <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv") %>% 
  select(-Form) %>% 
  select(-c(Nb_Lf,Nb_Lflet,Sple_Area,Sple_FM,Sple_DM,L_Length, L_Width,LTmes,Hveg,Dmax,Dmin,Flo,Mat_Per,Mat,LPC ))
write.csv(MEAN_multivar,"outputs/data/traits_multivariate.csv",row.names=F)

MEAN_univar <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM.csv") %>% 
  select(-Form) %>% 
  select(-c(Nb_Lf,Nb_Lflet,Sple_Area,Sple_FM,Sple_DM,L_Length, L_Width,LTmes,Hveg,Dmax,Dmin,Flo,Mat_Per,Mat,LPC ))
write.csv(MEAN_univar,"outputs/data/traits_univariate.csv",row.names=F)
           
MEAN_site <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_site_level.csv")%>% 
  select(-Form) %>% 
  select(-c(Nb_Lf,Nb_Lflet,Sple_Area,Sple_FM,Sple_DM,L_Length, L_Width,LTmes,Hveg,Dmax,Dmin,Flo,Mat_Per,Mat,LPC ))
write.csv(MEAN_site,"outputs/data/traits_site_level.csv",row.names=F)

