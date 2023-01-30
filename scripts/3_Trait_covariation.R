library(tidyverse)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_completed_seed_mass_flore.csv") %>%
  filter(!is.na(SLA)) %>% 
  filter(!(species== "Geranium dissectum - pétiole"))

# I case I want to perform the analyses with species both in the trait and abundance data:

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(line = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") %>% 
  filter(depth == "S")

data_fer <- MEAN %>%
  filter(code_sp %in% ab_fer$code_sp & treatment == "Fer")
data_nat <- MEAN %>%
  filter(code_sp %in% ab_nat$code_sp & treatment == "Nat")
MEAN_intersect <- rbind(data_fer,data_nat)


traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "H_FLORE",#"Hrepro"   , "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "FLO_FLORE", #Disp",#"Mat_Per", #"Mat","Flo",
            "SeedMass"
)


MEAN %>% 
  ggplot(aes(x=SLA,y=LDMC,color=LifeHistory)) +
  geom_point() 

# 1) SLA-LNC ####
MEAN %>% 
  ggplot(aes(x=SLA,y=LNC,color=LifeHistory)) +
  geom_point() +
  facet_wrap(~treatment) +
  geom_smooth(method="lm")

# 2) Height-SeedMass ####
MEAN %>% 
  ggplot(aes(x=H_FLORE,y=SeedMass,color=LifeHistory)) +
  geom_point() +
  facet_wrap(~treatment) +
  geom_smooth(method="lm")

# 3) LDMC-Ldelta13C ####
MEAN %>% 
  ggplot(aes(x=LDMC,y=Ldelta13C,color=LifeHistory)) +
  geom_point() +
  facet_wrap(~treatment) +
  geom_smooth(method="lm")

# 4) Height-flowering date ####
MEAN %>% 
  ggplot(aes(x=H_FLORE,y=FLO_FLORE,color=LifeHistory)) +
  geom_point() +
  facet_wrap(~treatment) +
  geom_smooth(method="lm")








