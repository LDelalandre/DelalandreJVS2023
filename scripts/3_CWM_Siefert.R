library(tidyverse)


# NB A REFAIRE (COMME TOUT LE RESTE) UNE FOIS QUE LES COQUILLES DE FORME DE VIE SERONT CORRIGEES AVEC LE DATA PAPER

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , #"Dmax"  , #    "Dmin" ,"Hveg"  , "H_FLORE",#
            "Disp", #"Mat_Per", #"Mat","Flo","FLO_FLORE", #
            "SeedMass"
)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv")
MEAN_site <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_site_level.csv")



ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(transect = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") %>% 
  filter(depth == "S") %>% 
  select(-depth) %>% 
  rename(transect = Ligne)

# Join abundance datasets
ab_fer2 <- ab_fer %>% 
  select(any_of(colnames(ab_nat))) %>% 
  mutate(treatment = "Fer")
ab_nat2 <- ab_nat %>% 
  mutate(treatment = "Nat") %>% 
  select(all_of(colnames(ab_fer2)))
ab <- rbind(ab_fer2,ab_nat2) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial")) %>% 
  select(species,code_sp,LifeForm1,LifeHistory,treatment,paddock,transect,abundance)



# faire une fonction pour calculer les cwm des traits:
# avec des traits moyennés à un niveau spatial


# pour une forme de vie
LH <- "annual"


# regarder les deux traitements


trait_averaging <- "treatment" # "site" or "treatment"

if (trait_averaging == "treatment"){
  fMEAN <- MEAN
} else {
  fMEAN <- MEAN_site
}

fMEAN <- fMEAN %>% 
  select("species","code_sp","LifeForm1","LifeHistory","treatment",
         all_of(traits))

ab_traits <- ab %>% 
  left_join(MEAN,by = c("species","code_sp","LifeForm1","LifeHistory","treatment")) %>% 
  group_by(transect,LifeHistory,treatment) %>% 
  mutate(relat_ab = abundance/sum(abundance))

data_cwm <- ab_traits %>% 
  mutate_at(vars(all_of(traits)),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique() %>% 
  select(LifeHistory,treatment,paddock, transect,starts_with("CWM")) %>% 
  unique() %>% 
  rename_at( vars( contains( "CWM_") ), list( ~ gsub("CWM_", "", .) ) )

