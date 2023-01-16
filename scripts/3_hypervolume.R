source("scripts/Packages.R")
library(FactoMineR)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library("missMDA")
# require("hypervolume")


MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv")

code_sp_lifeform <- read.csv2("data/species_names_lifehistory.csv")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , #"Dmax"  , #    "Dmin" ,"Hveg"  , "H_FLORE",#
            "Disp", #"Mat_Per", #"Mat","Flo","FLO_FLORE", #
            "SeedMass"
)

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(line = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") 

#_______________________________________________________________________________

MEAN_traits <- MEAN %>% 
  select(code_sp,LifeHistory,treatment,all_of(traits)) %>% 
  select(-Disp) %>% 
  na.omit()
MEAN_traits %>% 
  group_by(LifeHistory,treatment) %>% 
  summarize(n = n())

MEAN %>% 
  filter(treatment=="Nat") %>%
  filter(LifeHistory=="annual") %>% 
  column_to_rownames("code_sp") %>% 
  select(all_of(traits)) %>% View
  # select(-c(SeedMass,Ldelta13C)) %>% 
  # filter(!is.na(LDMC))
  # filter(!is.na(Ldelta13C))
  # # select(LDMC,SLA,L_Area,LCC,LNC) %>% 
  # # filter(!is.na(SLA)) %>% 
  na.omit()

MEAN  %>% 
  left_join(code_sp_lifeform)

MEAN %>% 
  select(code_sp,LifeHistory,treatment,all_of(traits)) %>%  
  na.omit()

in_traits <- MEAN_traits %>% 
  rownames()

in_ab <- ab_nat %>% 
  # filter(!(LifeForm1=="The")) %>% 
  pull(code_sp) %>% 
  unique()

intersect(in_traits,in_ab)
setdiff(in_traits,in_ab)
setdiff(in_ab,in_traits)


# coverage in terms of abundance ####
relat_ab_nat <- ab_nat %>%
  filter(depth == "S") %>% 
  group_by(species,code_sp) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance)) %>% 
  left_join(sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial"))

relat_ab_fer <- ab_fer %>% 
  group_by(species,code_sp) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance))%>% 
  left_join(sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial")) %>% 
  #add missing life history (NB: not missing lifeform!!)
  mutate(LifeHistory = case_when(
    is.na(LifeHistory) & !(species %in% c("Vicia sativa ssp. sativa","Trifolium stellatum"))~"perennial",
    is.na(LifeHistory) & species %in% c("Vicia sativa ssp. sativa","Trifolium stellatum") ~ "annual",
    TRUE ~ LifeHistory))


choice_life_history <- "annual"

trait_available <- MEAN_traits %>% 
  rownames_to_column("code_sp") %>% 
  left_join(code_sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial")) %>% 
  mutate(sp_in_trait = 1) %>% 
  select(species,code_sp,LifeHistory,sp_in_trait) %>% 
  # full_join(relat_ab_nat,by=c("species","code_sp","LifeHistory")) %>%
  full_join(relat_ab_fer,by=c("species","code_sp","LifeHistory")) %>%
  filter(LifeHistory==choice_life_history) %>%
  
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance,na.rm=T)) %>% 
  select(species,sp_in_trait,code_sp,sp_relat_abundance,LifeHistory) %>% 
  replace(is.na(.),0)

# cumulated abundance for species for which we have the trait
info_coverage <- trait_available %>%
  group_by(sp_in_trait) %>% 
  summarise(abundance_covered = sum(sp_relat_abundance)) %>% 
  spread(key = sp_in_trait, value = abundance_covered)
info_coverage


# missMDA ####
# Check the strength of imputations
nbdim <- estim_ncpPCA(MEAN_traits,ncp.min = 1)
res.comp <- MIPCA(MEAN_traits, ncp = nbdim$ncp, nboot = 1000)
plot(res.comp)
