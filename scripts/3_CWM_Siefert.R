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




compute_cwm <- function(ab,MEAN,MEAN_site,level_of_trait_averaging){
  
  if (level_of_trait_averaging == "treatment"){
    fMEAN <- MEAN %>% 
      select("species","code_sp","LifeForm1","LifeHistory","treatment",
             all_of(traits))
  } else if (level_of_trait_averaging == "site") {
    fMEAN <- MEAN_site %>% 
      select("species","code_sp","LifeForm1","LifeHistory",
             all_of(traits))
  }
  
  ab_traits <- ab %>% 
    left_join(fMEAN) %>% 
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
  
  data_cwm
}

trait_averaging <- "site" # "site" or "treatment"
compute_cwm(ab,MEAN,MEAN_site,level_of_trait_averaging = trait_averaging)
# pb when working on site level

# Compute the different CWM ####

# CWM with both species replacement and intrasp variability  
CWM <- compute_cwm(ab,MEAN,MEAN_site,level_of_trait_averaging = "treatment")

# CWM with species replacement but no intraspecific variability
CWMfixed <- compute_cwm(ab,MEAN,MEAN_site,level_of_trait_averaging = "site")

# CWM intra  = CWM - CWMfixed pour chaque trait


# Explain variation at different levels
ftrait <- "SLA"
# make a data frame with variation at these levels
fCWM <- CWM %>% 
  select(LifeHistory,treatment,paddock,transect, ftrait) %>% 
  mutate(CWM = get(ftrait)) %>% 
  select(-ftrait)
fCWMfixed <- CWMfixed %>% 
  select(LifeHistory,treatment,paddock,transect, ftrait) %>% 
  mutate(CWMfixed = get(ftrait)) %>% 
  select(-ftrait)

# CWM data for the focal trait, at the three levels
fCWMall <- full_join(fCWM,fCWMfixed) %>% 
  mutate(CWMintra = CWM - CWMfixed)

fLH <- "perennial"

fdata <- fCWMall %>% 
  filter(LifeHistory == fLH)

mod <- lm(CWM ~ 1, data = fdata)
ano <- anova(mod)
SS <- ano$`Sum Sq`

modfixed <- lm(CWMfixed ~ 1, data = fdata)
anofixed <- anova(modfixed)
SSfixed <- anofixed$`Sum Sq`

modintra <- lm(CWMintra ~ 1, data = fdata)
anointra <- anova(modintra)
SSintra <- anointra$`Sum Sq`


# aITV (Siefert), à faire pour annuelles et pérennes. Et voir l'effet du traitement (comment le coder ?)
aITV <- log(SSintra/SSfixed)
