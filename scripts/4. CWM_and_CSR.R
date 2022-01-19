source("scripts/1. Packages.R")

# Load data ####
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")

# PCA ####
# Compute PCA on traits from all species measured in each treatment
MEAN %>% 
  filter(treatment == "Fer")

MEAN %>% 
  filter(treatment == "Nat")

# Compute CSR scores ####
# Species level
MEAN_goodunit <- MEAN %>% 
  # select(species,code_sp,LifeHistory,LifeForm1,treatment,L_Area,LDMC,SLA) %>% 
  relocate(species,treatment,code_sp,Form,LifeForm1,LifeHistory,L_Area,LDMC,SLA) %>%
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% # change from mg/g to %
  filter(!is.na(L_Area)) %>% 
  filter(!is.na(SLA))
write.csv2(MEAN_goodunit,"outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt.csv" ,row.names=F)
  
MEAN_CSR <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_completed.csv",dec=",") %>% 
  # merge(name_LH,by="Code_Sp") %>% 
  # relocate(C,S,R) %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric())

# Compute CWM of traits and CSR scores ####
# I can compute CWM of CSR in the two ways (CWM of traits and CSR, and CSR and CWM of these index).

# i) Fer ####
ab_traits_fer <- ab_fer %>% 
  left_join(MEAN_CSR %>% filter(treatment == "Fer"),
            by = c("species","code_sp","LifeForm1","treatment")) 
  
# I keep it with species level I case I want to compute distance to CWM for some species
CWM_sp_fer <- ab_traits_fer %>% 
  ungroup() %>% 
  group_by(id_transect_quadrat) %>% 
  mutate_at(vars(L_Area:R),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique()

CWM_fer <- CWM_sp_fer %>% 
  select(paddock, id_transect_quadrat,starts_with("CWM")) %>% 
  unique()

# Compute CSR scores on CWM
write.csv2(CWM_fer,"outputs/data/Pierce CSR/Traits_CWM_fer.csv" ,row.names=F)
# NB : CWM_S = CWM of S scores
# and S_CWM = S computed on CWM of L_Area, SLA, and LDMC scores.
CWM2_fer <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_fer_completed.csv")

# ii) Nat ####
ab_traits_nat <- ab_nat %>% 
  left_join(MEAN_CSR %>% filter(treatment == "Nat"),
            by = c("species","code_sp")) 
# /!\ pb: Anthyllis vulneraria has two LifeForm1 : it adds lines as new species. To be corrected.

# I keep it with species level I case I want to compute distance to CWM for some species
CWM_sp_nat <- ab_traits_nat %>% 
  ungroup() %>% 
  group_by(depth,paddock) %>% 
  mutate_at(vars(L_Area:R),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique()

CWM_nat <- CWM_sp_nat %>% 
  select(PC1score,depth,paddock,line,starts_with("CWM")) %>% 
  unique()
# Compute CSR scores on CWM
write.csv2(CWM_nat,"outputs/data/Pierce CSR/Traits_CWM_nat.csv" ,row.names=F)
# NB : CWM_S = CWM of S scores
# and S_CWM = S computed on CWM of L_Area, SLA, and LDMC scores.
CWM2_nat <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_nat_completed.csv")



