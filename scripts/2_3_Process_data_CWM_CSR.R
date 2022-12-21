# This script computes CWM of CSR scores (using Pierce's spreadsheet outside of R).

source("scripts/Packages.R")

subset_gt_nat <- T # traits measured in the GUs and GUi, not all GU

# Load data ####
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv")
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")

# Compute CSR scores ####
# Species level

# -- EXCEL SPREADSHEET -- (manually import and transform data)

# Compute CWM of traits and CSR scores ####
# I can compute CWM of CSR in the two ways (CWM of traits and CSR, and CSR and CWM of these index).
if(subset_gt_nat == F){
  MEAN_CSR <- read.csv2("outputs/data/mean_attribute_per_treatment_completed.csv")
}else{
  MEAN_CSR <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_completed.csv")%>%
    filter(!is.na(SLA))
}


MEAN_CSR <- MEAN_CSR %>% 
  mutate(C=str_replace(C,",",".") %>% as.numeric())%>% 
  mutate(S=str_replace(S,",",".") %>% as.numeric())%>% 
  mutate(R=str_replace(R,",",".") %>% as.numeric())

# Add Ellenberg info ####
ellenberg_fage <- read.csv2("data/traits/ellenberg_la_fage.csv") %>% 
  rename(species = AccSpeciesName)
MEAN_CSR_ellenberg <- MEAN_CSR %>% 
  left_join(ellenberg_fage,by="species")

# 1) CWM in all the communities ####
# i) Fer ####
ab_traits_fer <- ab_fer %>% 
  left_join(MEAN_CSR_ellenberg %>% filter(treatment == "Fer"),
            by = c("species","code_sp","LifeForm1","treatment")) 
  
# I keep it with species level I case I want to compute distance to CWM for some species
CWM_sp_fer <- ab_traits_fer %>% 
  ungroup() %>% 
  group_by(id_transect_quadrat) %>% 
  mutate_at(vars(L_Area:temperature),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique()

CWM_fer <- CWM_sp_fer %>% 
  select(paddock, id_transect_quadrat,starts_with("CWM")) %>% 
  unique()

write.csv2(CWM_fer,"outputs/data/Pierce CSR/Traits_CWM_fer.csv" ,row.names=F)

# Optional: Compute CSR scores on CWM.
# -- EXCEL SPREADSHEET -- (manually import and transform data)
# NB : CWM_S = CWM of S scores
# and S_CWM = S computed on CWM of L_Area, SLA, and LDMC scores.
CWM2_fer <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_fer_completed.csv")



# ii) Nat ####
ab_traits_nat <- ab_nat %>% 
  left_join(MEAN_CSR_ellenberg %>% filter(treatment == "Nat"),
            by = c("species","code_sp","LifeHistory")) 
# /!\ pb: Anthyllis vulneraria has two LifeForm1 : it adds lines as new species. To be corrected.

# I keep it with species level I case I want to compute distance to CWM for some species
CWM_sp_nat <- ab_traits_nat %>% 
  ungroup() %>% 
  group_by(depth,paddock) %>% 
  mutate_at(vars(L_Area:temperature),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique()

CWM_nat <- CWM_sp_nat %>% 
  select(PC1score,depth,paddock,line,starts_with("CWM")) %>% 
  unique()

if(subset_gt_nat == F){
  write.csv2(CWM_nat,"outputs/data/Pierce CSR/Traits_CWM_nat.csv" ,row.names=F)
}else{
  write.csv2(CWM_nat,"outputs/data/Pierce CSR/Traits_CWM_subset_nat_sab.csv" ,row.names=F)
}


# Optional: compute CSR scores on CWM
# -- EXCEL SPREADSHEET -- (manually import and transform data)
# NB : CWM_S = CWM of S scores
# and S_CWM = S computed on CWM of L_Area, SLA, and LDMC scores.
CWM2_nat <- read.csv2("outputs/data/Pierce CSR/Traits_CWM_nat_completed.csv")


# 2) CWM in the guild of annuals ####
# i) Fer ####
CWM_annuals_fer <- ab_traits_fer %>% 
  filter(LifeHistory == "annual") %>% 
  group_by(id_transect_quadrat) %>% 
  mutate(relat_ab = abundance/sum(abundance)) %>% 
  mutate_at(vars(L_Area:R),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique() %>% 
  select(paddock, id_transect_quadrat,starts_with("CWM")) %>% 
  unique() %>% 
  rename_at( vars( contains( "CWM_") ), list( ~ gsub("CWM_", "", .) ) )

write.csv2(CWM_annuals_fer,"outputs/data/CWM_annuals_fer.csv" ,row.names=F)

# ii) Nat ####
CWM_annuals_nat <- ab_traits_nat %>% 
  filter(LifeHistory == "annual") %>% 
  group_by(paddock,line,depth) %>% 
  mutate(relat_ab = abundance/sum(abundance)) %>% 
  mutate_at(vars(L_Area:R),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique() %>% 
  select(depth,paddock,line, starts_with("CWM")) %>% 
  unique() %>% 
  rename_at( vars( contains( "CWM_") ), list( ~ gsub("CWM_", "", .) ) )

write.csv2(CWM_annuals_nat,"outputs/data/CWM_annuals_nat.csv" ,row.names=F)


# 3) CWM in the guild of perennials ####
# i) Fer ####
CWM_perennials_fer <- ab_traits_fer %>% 
  filter(LifeHistory == "perennial") %>% 
  group_by(id_transect_quadrat) %>% 
  mutate(relat_ab = abundance/sum(abundance)) %>% 
  mutate_at(vars(L_Area:R),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique() %>% 
  select(paddock, id_transect_quadrat,starts_with("CWM")) %>% 
  unique() %>% 
  rename_at( vars( contains( "CWM_") ), list( ~ gsub("CWM_", "", .) ) )

write.csv2(CWM_perennials_fer,"outputs/data/CWM_perennials_fer.csv" ,row.names=F)

# ii) Nat ####
CWM_perennials_nat <- ab_traits_nat %>% 
  filter(LifeHistory == "perennial") %>% 
  group_by(paddock,line,depth) %>% 
  mutate(relat_ab = abundance/sum(abundance)) %>% 
  mutate_at(vars(L_Area:R),
            .funs = list(CWM = ~ weighted.mean(.,relat_ab,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique() %>% 
  select(depth,paddock,line, starts_with("CWM")) %>% 
  unique() %>% 
  rename_at( vars( contains( "CWM_") ), list( ~ gsub("CWM_", "", .) ) )

write.csv2(CWM_perennials_nat,"outputs/data/CWM_perennials_nat.csv" ,row.names=F)

