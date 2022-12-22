library(tidyverse)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab.csv")%>%
  filter(!is.na(SLA)) %>% 
  filter(!is.na(L_Area)) 

# units for Segrestin's the shiny app
MEAN_shiny <- MEAN %>% 
  mutate(H = Hrepro/100 ) %>%  # from cm to m
  mutate(LA=L_Area*10)  %>% # to change unit from cm² to mm²
  mutate(SLA_m2_g = 0.001 * SLA) %>%  # change SLA from mm²/mg to m²/g
  mutate(LMA = 1/SLA_m2_g) %>%  # LMA en g/m²
  rename(Nmass = LNC) %>% # LNC in mg/g 
  rename(SM = SeedMass) %>% # SeedMass in mg

  select(code_sp,treatment,LifeHistory,
         H,LA,LMA,Nmass,SM)

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

data_fer <- MEAN_shiny %>%
  filter(code_sp %in% ab_fer$code_sp & treatment == "Fer")
data_nat <- MEAN_shiny %>%
  filter(code_sp %in% ab_nat$code_sp & treatment == "Nat")
MEAN_shiny2 <- rbind(data_fer,data_nat)

write.csv2(MEAN_shiny,"outputs/data/global spectrum/data_shiny.csv",row.names=F)
write.csv2(MEAN_shiny2,"outputs/data/global spectrum/data_shiny_in_ab.csv",row.names=F)

# Predict SSD from LDMC ####

# We thus used three different equations to predict
# SSD for 1963 herbaceous species for which LDMC values were available in TRY:
#   one for monocotyledons (SSD = 0.888 × LDMC + 2.69), one for Leguminosae
# (SSD = 0.692 × LDMC + 47.65), and a third one for other non-monocotyledons
# (SSD = 0.524 × LDMC + 95.87).


phenospace <- read.csv2("outputs/data/Global Spectrum/PhenoSpace_3traits.csv",dec=".")
phenospace <- read.csv2("outputs/data/Global Spectrum/PhenoSpace_in_ab.csv",dec=".")

phenospace %>% 
  ggplot(aes(x=PC1,y=PC2,color = LifeHistory)) +
  geom_point() +
  facet_wrap(~treatment)
