library(tidyverse)

MEAN_multivar <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv")
write.csv2(MEAN_multivar,"outputs/data/traits_multivariate.csv",row.names=F)

MEAN_univar <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM.csv")
write.csv2(MEAN_univar,"outputs/data/traits_univariate.csv",row.names=F)
           
MEAN_site <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_site_level.csv")
write.csv2(MEAN_site,"outputs/data/traits_site_level.csv",row.names=F)
