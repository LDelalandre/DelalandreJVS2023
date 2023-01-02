library(tidyverse)
library(dplyr)
# detach("package:MASS")    

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_completed_seed_mass_flore.csv") %>%
  filter(!is.na(SLA)) %>% 
  filter(!(species== "Geranium dissectum - pétiole"))

sp_lifeform <- read.csv2("data/species_names_lifehistory.csv") %>% 
  rename(code_sp = Code_Sp, species= Species)

trait_unit <- read.csv2("data/trait_unit.csv",encoding = "latin1")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "H_FLORE",#"Hrepro"   , "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "FLO_FLORE", #Disp",#"Mat_Per", #"Mat","Flo",
            "SeedMass"
)

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(line = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") 

ab_nat %>% group_by(paddock,depth,line) %>% 
  summarize(n=n())

ab_nat %>% filter(depth == "S") %>% pull(Ligne) %>% unique() %>% length

ab_fer_ann <- ab_fer %>% 
  filter(LifeForm1 == "The")
ab_nat_ann <- ab_nat %>% 
  filter(LifeHistory == "annual")


#___________________________________________________________
# Species in the abundance, but not in the trait, data? ####
# lifehistory == "annual"

FTRAIT <- c("LDMC","LCC","Ldelta13C","H_FLORE","SeedMass")
FDF <- NULL

for (ftrait in FTRAIT){
  
  MEAN_ftrait <- MEAN %>% 
    filter(!is.na(get(ftrait)))
  
  per_ab_fer <- ab_fer %>% 
    filter(!(LifeForm1=="The")) %>%
    arrange(code_sp) %>% 
    pull(code_sp) %>% 
    unique() 
  ann_ab_fer <- ab_fer %>% 
    filter((LifeForm1=="The")) %>% 
    pull(code_sp) %>% 
    unique()
  
  per_trait_fer <- MEAN_ftrait %>% 
    filter(!(LifeForm1=="The")) %>% 
    filter(treatment == "Fer") %>%  
    pull(code_sp) %>% unique()
  ann_trait_fer <- MEAN_ftrait %>% 
    filter((LifeForm1=="The")) %>% 
    filter(treatment == "Fer") %>%  
    pull(code_sp) %>% unique()
  
  A <- per_ab_fer
  B <- per_trait_fer
  
  A2 <- ann_ab_fer
  B2 <- ann_trait_fer
  
  fer_ab_onlyP <- setdiff(A,B) %>% length() # dans le premier mais pas dans le deuxième
  fer_trait_onlyP <- setdiff(B,A) %>% length()
  fer_instersectP <- intersect(A,B) %>% length()
  
  fer_ab_onlyA <- setdiff(A2,B2) %>% length() # dans le premier mais pas dans le deuxième
  fer_trait_onlyA <- setdiff(B2,A2) %>% length()
  fer_instersectA <- intersect(A2,B2) %>% length()
  
  per_ab_nat <- ab_nat %>%   
    filter(!(LifeForm1=="The")) %>% 
    filter(depth=="S") %>% 
    arrange(code_sp) %>%
    pull(code_sp) %>% unique()
  ann_ab_nat <- ab_nat %>%   
    filter((LifeForm1=="The")) %>% 
    filter(depth=="S") %>% 
    pull(code_sp) %>% unique()
  
  per_trait_nat <- MEAN_ftrait %>%   
    filter(!(LifeForm1=="The")) %>% 
    filter(treatment == "Nat") %>%  
    pull(code_sp) %>% unique()
  ann_trait_nat <- MEAN_ftrait %>%   
    filter((LifeForm1=="The")) %>% 
    filter(treatment == "Nat") %>%  
    pull(code_sp) %>% unique()
  
  C <- per_ab_nat
  D<- per_trait_nat
  
  C2 <- ann_ab_nat
  D2 <- ann_trait_nat
  

  
  nat_ab_onlyP <- setdiff(C,D) %>% length() # dans le premier mais pas dans le deuxième
  nat_trait_onlyP <- setdiff(D,C) %>% length()
  nat_instersectP <- intersect(C,D) %>% length()
  
  nat_ab_onlyA <- setdiff(C2,D2) %>% length() # dans le premier mais pas dans le deuxième
  nat_trait_onlyA <- setdiff(D2,C2) %>% length()
  nat_instersectA <- intersect(C2,D2) %>% length()
  
  
  fdf <- data.frame(trait = ftrait,
    treatment = c(rep("Fer",2),rep("Nat",2)),
                    LifeForm1 = c(rep("Per",1),rep("Ann",1),rep("Per",1),rep("Ann",1)),
                    in_ab_only = c(fer_ab_onlyP,fer_ab_onlyA,nat_ab_onlyP,nat_ab_onlyA) ,
                    in_trait_only = c(fer_trait_onlyP,fer_trait_onlyA,nat_trait_onlyP,nat_trait_onlyA),
                    intersect = c(fer_instersectP,fer_instersectA,nat_instersectP,nat_instersectA) ) %>% 
    mutate(in_traits = in_trait_only + intersect)%>% 
    mutate(trait_coverage = intersect/(intersect+in_ab_only))
  
  FDF <- rbind(FDF,fdf)

}

FDF 

# table à faire
table_trait_abundance_overlap <- FDF %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Trait", "Abbr.", "Unit",
                                   "Abundance in G+F",
                                   "Abundance in GUs",
                                   "Abundance in G+F",
                                   "Abundance in GUs")) %>%
  kableExtra::kable_styling("hover", full_width = F) %>% 
  kableExtra::add_header_above(c(" " = 3,"Annuals" = 2,"Perennials" = 2))


#________________________________________________
# Abundance and trait coverage info together ####

# The regional relative abundance of species (RA) 
# is represented in percent of total pin-point contacts observed across all 
# 60 lines of vegetation survey. 
# Here I do it within annual species and all species

# Fer
tot_contact_fer <- ab_fer %>% 
  pull(abundance) %>% 
  sum()
tot_contact_fer_ann <- ab_fer_ann %>% 
  pull(abundance) %>% 
  sum()
tot_contact_fer_ann/tot_contact_fer # percent of contacts represented by annuals

# Regional relative abundance of each species 
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

# Regional relative abundance of each annual species within the annual guild
relat_ab_fer_ann <- ab_fer_ann %>% 
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

# Nat
# Regional relative abundance in GUS
tot_contact_nat <- ab_nat %>%
  filter(depth == "S") %>% 
  pull(abundance) %>% 
  sum()
tot_contact_nat_ann <- ab_nat_ann %>% 
  filter(depth == "S") %>% 
  pull(abundance) %>% 
  sum()
tot_contact_nat_ann/tot_contact_nat

# Regional relative abundance of each species 
relat_ab_nat <- ab_nat %>%
  filter(depth == "S") %>% 
  group_by(species,code_sp) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance)) %>% 
  left_join(sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial"))

# Regional relative abundance of each annual species within the annual guild
relat_ab_nat_ann <- ab_nat_ann %>%
  filter(depth == "S") %>% 
  group_by(species,code_sp) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance))%>% 
  left_join(sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial"))

choice_life_history <- "annual"

TRAITS_ABUNDANCE_COVERAGE <- list()
i <- 0
for (choice_life_history in c("annual","perennial")){
  i <- i+1
  traits_abundance_coverage <- NULL
  for (ftrait in traits){
    # Natif
    trait_available_nat <- MEAN %>% # or MEAN_annual
      filter(treatment=="Nat") %>% 
      mutate(trait = get(ftrait)) %>%  # choose the trait here
      # LES SP AVEC NA POUR LA LIFEFORM ONT LIFEHISTORY = NA !!!!!
      mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial")) %>% 
      select(species,code_sp,trait,LifeHistory) %>% 
      full_join(relat_ab_nat,by=c("species","code_sp","LifeHistory")) %>%  # or relat_ab_nat_ann
      # IL FAUT RAJOUTER LIFEHISTORY ICI POUR LES ESPECES AJOUTEES DEPUIS relat_ab_nat
      # NB: LE FAIRE DIRECTEMENT DANS LE FICHIER D'ABONDANCE ?
      # Et NA --> 0 en abundance
      # filter(!(LifeHistory==choice_life_history)) %>% # /!\ only for relative abundance of annuals !!
      filter(LifeHistory==choice_life_history) %>%
      
      mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance,na.rm=T)) %>% 
      
      mutate(trait_available = if_else(is.na(trait),"no","yes")) %>%
      # mutate(trait_available = if_else(trait==0,"no","yes")) %>% 
      select(species,code_sp,trait_available,sp_relat_abundance,LifeHistory) %>% 
      replace(is.na(.),0)
    # 
    # trait_available_nat %>%
    #   ggplot(aes(x=trait_available,y=sp_relat_ab_nat)) +
    #   geom_histogram(position='stack',stat="identity") +
    #   ggtitle(paste("Natif",ftrait))
    
    
    # cumulated abundance for species for which we have the trait
    Tnat <- trait_available_nat %>%
      group_by(trait_available) %>% 
      summarise(abundance_covered = sum(sp_relat_abundance)) %>% 
      spread(key = trait_available, value = abundance_covered) %>% 
      mutate(trait = ftrait) %>% 
      mutate(treatment = "Nat")
    
    if (is.null(Tnat$yes) ){ 
      Tnat$yes <- 0
    }
    if (is.null(Tnat$no) ){ 
      Tnat$no <- 0
    }
    
    Tnat <- Tnat %>% select(all_of(c("trait","treatment","yes","no")))
    
    
    # Fertile
    trait_available_fer <- MEAN %>% # or MEAN_annual
      filter(treatment=="Fer") %>% 
      mutate(trait = get(ftrait)) %>%  # choose the trait here
      # LES SP AVEC NA POUR LA LIFEFORM ONT LIFEHISTORY = NA !!!!!
      mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial")) %>% 
      select(species,LifeHistory,code_sp,trait) %>% 
      full_join(relat_ab_fer,by=c("species","code_sp","LifeHistory")) %>% # or relat_ab_fer_ann
      
      filter(!(LifeHistory==choice_life_history)) %>% # /!\ only for relative abundance of annuals !!
      
      mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance,na.rm=T)) %>%
      
      mutate(trait_available = if_else(is.na(trait),"no","yes")) %>% 
      select(species,code_sp,trait_available,sp_relat_abundance)%>% 
      replace(is.na(.),0)
    # 
    # trait_available_fer %>% 
    #   ggplot(aes(x=trait_available,y=sp_relat_ab_fer)) +
    #   geom_histogram(position='stack',stat="identity") +
    #   ggtitle("Fertile")
    
    # cumulated abundance for species for which we have the trait
    Tfer <- trait_available_fer %>% 
      group_by(trait_available) %>% 
      summarise(abundance_covered = sum(sp_relat_abundance)) %>% 
      spread(key = trait_available, value = abundance_covered) %>% 
      mutate(trait = ftrait) %>% 
      mutate(treatment = "Fer")
    
    if (is.null(Tfer$no) ){ 
      Tfer$no <- 0
    }
    
    Tfer <- Tfer %>% select(all_of(c("trait","treatment","yes","no")))
    
    traits_abundance_coverage <- rbind(traits_abundance_coverage, Tnat )
    traits_abundance_coverage <- rbind(traits_abundance_coverage, Tfer )
  }
  TRAITS_ABUNDANCE_COVERAGE[[i]] <- traits_abundance_coverage
}

cover_ann <- TRAITS_ABUNDANCE_COVERAGE[[1]] %>% 
  mutate(yes = round(yes,3)) %>% 
  mutate(no = round(no,3)) %>% 
  select(-no) %>% 
  spread(key=treatment,value=yes) %>% 
  arrange(factor(trait,levels = traits))

cover_per <- TRAITS_ABUNDANCE_COVERAGE[[2]]%>% 
  mutate(yes = round(yes,3)) %>% 
  mutate(no = round(no,3)) %>% 
  select(-no) %>% 
  spread(key=treatment,value=yes) %>% 
  arrange(factor(trait,levels = traits))

cover <- merge(cover_ann,cover_per,by="trait")%>% 
  merge(trait_unit %>% rename(trait = Trait)) %>% 
  arrange(factor(trait,levels = traits)) %>% 
  select(Full_name,trait,Unit,everything())

# traits_abundance_coverage %>% 
#   ggplot(aes(x=trait,y=yes))+
#   geom_point() +
#   facet_wrap(~treatment)


table_trait_coverage <- cover %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Trait", "Abbr.", "Unit",
                                   "Abundance in G+F",
                                   "Abundance in GUs",
                                   "Abundance in G+F",
                                   "Abundance in GUs")) %>%
  kableExtra::kable_styling("hover", full_width = F) %>% 
  kableExtra::add_header_above(c(" " = 3,"Annuals" = 2,"Perennials" = 2))

table_trait_coverage 


cat(table_trait_coverage, file = "draft/cover_traits.doc")
# cat(table_trait_coverage, file = "draft/cover_traits_annuals.doc")
# cat(table_trait_coverage, file = "draft/cover_traits_perennials_approx_SM.doc")

