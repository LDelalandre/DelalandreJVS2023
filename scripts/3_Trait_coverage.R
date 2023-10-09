# trait coverage

library(tidyverse)

MEAN <- read.csv("outputs/data/traits_univariate.csv")
MEAN_multivar <- read.csv("outputs/data/traits_multivariate.csv") %>% 
  select(-Disp)

code_sp_lifeform <- read.csv2("data/species_names_lifehistory.csv")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , #"Dmax"  , #    "Dmin" ,"Hveg"  , "H_FLORE",#
            "Disp", #"Mat_Per", #"Mat","Flo","FLO_FLORE", #
            "SeedMass"
)

data_abundance <- read.csv("outputs/data/data_abundance.csv")
ab_fer <- data_abundance %>% filter(treatment == "Int")
ab_nat <- data_abundance %>% filter(treatment == "Ext") 


# coverage (abundance) ####
relat_ab_nat <- ab_nat %>%
  filter(depth == "S") %>% 
  group_by(species,code_sp,LifeHistory) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance)) %>% 
  left_join(code_sp_lifeform) 

relat_ab_fer <- ab_fer %>% 
  group_by(species,code_sp,LifeHistory) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance))%>% 
  left_join(code_sp_lifeform) 



INFO_COVERAGE <- NULL
NB_SP <- NULL

for(trtmt in c("Fer","Nat")){
  if (trtmt == "Nat"){
    file_abundance <- relat_ab_nat
  }else{
    file_abundance <- relat_ab_fer
  }
  
  for (ftrait2 in c(traits,"Multivariate") ){
    
    if (ftrait2 == "Multivariate"){
      ftrait <- c("LDMC"    ,  "SLA"    ,   "L_Area"   , "LCC"  ,     "LNC"    ,   "Ldelta13C" ,"Hrepro"   ,     "SeedMass")
    } else{
      ftrait <- ftrait2
    }
    
    trait_available <- MEAN %>% 
      filter(treatment==trtmt) %>%
      group_by(LifeHistory) %>%
      select(species,code_sp,LifeHistory,all_of(ftrait)) %>% 
      na.omit() %>% 
      mutate(sp_in_trait = 1) %>% 
      select(species,code_sp,LifeHistory,sp_in_trait) %>% 
      full_join(file_abundance,by=c("species","code_sp","LifeHistory")) %>%
      # filter(LifeHistory==choice_life_history) %>%
      select(-LifeForm1) %>% 
      select(-code_sp) %>% 
      replace(is.na(.),0) %>%
      group_by(LifeHistory) %>%
      mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance,na.rm=T)) %>% 
      select(species,sp_in_trait,sp_relat_abundance,LifeHistory)
    
    # if (length(ftrait)>1  ){
    #   ftrait2 <- "Multivariate"
    # }else{
    #   ftrait2 <- ftrait
    # }
    
    # cumulated abundance for species for which we have the trait
    info_coverage <- trait_available %>%
      group_by(sp_in_trait,LifeHistory) %>% 
      summarise(abundance_covered = sum(sp_relat_abundance)) %>% 
      spread(key = sp_in_trait, value = abundance_covered) %>% 
      mutate(treatment = trtmt) %>% 
      rename(no = '0', yes = '1') %>% 
      mutate(trait = ftrait2)
    INFO_COVERAGE <- rbind(INFO_COVERAGE,info_coverage)
    
    # proportion of species
    nb_sp <- trait_available %>%
      mutate(sp_in_ab = case_when(sp_relat_abundance >0 ~ 1,
                                  TRUE ~0)) %>%
      group_by(sp_in_ab,sp_in_trait,LifeHistory) %>% 
      summarise(n = n()) %>% 
      mutate(ab_trait = paste(sp_in_ab,sp_in_trait)) %>%
      ungroup() %>% 
      select(LifeHistory,ab_trait,n) %>% 
      spread(key = ab_trait, value = n) %>% 
      rename(in_trait = '0 1',in_ab = '1 0',in_both = '1 1') %>% 
      mutate(treatment = trtmt) %>% 
      mutate(trait = ftrait2)
    NB_SP <- rbind(NB_SP,nb_sp)
    
  }
}


# abundance ####
cover <- INFO_COVERAGE %>% 
  select(-no) %>% 
  unique() %>% 
  mutate(yes = round(yes,2)) %>% 
  mutate(categ = paste(LifeHistory,treatment,sep = "_")) %>%
  select(-c(LifeHistory,treatment)) %>% 
  spread(key = categ,value = yes) %>% 
  arrange(factor(trait,levels = c(traits,"Multivariate"))) 

table_trait_coverage <- cover %>% 
  filter(!(trait == "Multivariate")) %>%
  kableExtra::kable( escape = F,
                     col.names = c("Trait",
                                   "Intensive",
                                   "Extensive",
                                   "Intensive",
                                   "Extensive")) %>%
  kableExtra::kable_styling("hover", full_width = F)  %>% 
  kableExtra::add_header_above(c(" " = 1,"Annuals" = 2,"Perennials" = 2))

table_trait_coverage 

cat(table_trait_coverage, file = "draft/table_trait_coverage_abundance.doc")
               


# proportion of species ####
FTRAIT <- traits
FDF <- NULL

for (ftrait in c(FTRAIT,"multivariate")){
  
  if(ftrait == "multivariate" ){
    MEAN_ftrait <- MEAN_multivar %>% 
      select(species,treatment,code_sp,LifeHistory,
             any_of(traits)) %>% 
      na.omit()
  }else{
    MEAN_ftrait <- MEAN %>% 
      filter(!is.na(get(ftrait)))
  }
  

  
  per_ab_fer <- ab_fer %>% 
    filter(!(LifeHistory=="annual")) %>%
    arrange(code_sp) %>% 
    pull(code_sp) %>% 
    unique() 
  ann_ab_fer <- ab_fer %>% 
    filter((LifeHistory=="annual")) %>% 
    pull(code_sp) %>% 
    unique()
  
  per_trait_fer <- MEAN_ftrait %>% 
    filter(!(LifeHistory=="annual")) %>% 
    filter(treatment == "Fer") %>%  
    pull(code_sp) %>% unique()
  ann_trait_fer <- MEAN_ftrait %>% 
    filter((LifeHistory=="annual")) %>% 
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
    filter(!(LifeHistory=="annual")) %>% 
    filter(depth=="S") %>% 
    arrange(code_sp) %>%
    pull(code_sp) %>% unique()
  ann_ab_nat <- ab_nat %>%   
    filter((LifeHistory=="annual")) %>% 
    filter(depth=="S") %>% 
    pull(code_sp) %>% unique()
  
  per_trait_nat <- MEAN_ftrait %>%   
    filter(!(LifeHistory=="annual")) %>% 
    filter(treatment == "Nat") %>%  
    pull(code_sp) %>% unique()
  ann_trait_nat <- MEAN_ftrait %>%   
    filter((LifeHistory=="annual")) %>% 
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
                    LifeHistory = c(rep("Per",1),rep("Ann",1),rep("Per",1),rep("Ann",1)),
                    in_ab_only = c(fer_ab_onlyP,fer_ab_onlyA,nat_ab_onlyP,nat_ab_onlyA) ,
                    in_trait_only = c(fer_trait_onlyP,fer_trait_onlyA,nat_trait_onlyP,nat_trait_onlyA),
                    intersect = c(fer_instersectP,fer_instersectA,nat_instersectP,nat_instersectA) ) %>% 
    mutate(in_traits = in_trait_only + intersect)%>% 
    mutate(trait_coverage = intersect/(intersect+in_ab_only))
  
  FDF <- rbind(FDF,fdf)
  
}

FDF 

# table à faire

table_trait_coverage_nb_sp <- FDF %>% 
  mutate(LH_tr = paste(LifeHistory,treatment)) %>% 
  mutate(trait_coverage = round(trait_coverage,2)) %>% 
  select(trait,LH_tr,trait_coverage) %>% 
  spread(key = LH_tr,value = trait_coverage) %>% 
  arrange(trait = factor (trait, levels = traits)) %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Trait",
                                   "Intensive",
                                   "Extensive",
                                   "Intensive",
                                   "Extensive")) %>%
  kableExtra::kable_styling("hover", full_width = F) %>% 
  # kableExtra::add_header_above(c(" " =1,"Management" = 4)) %>% 
  kableExtra::add_header_above(c(" " = 1,"Annuals" = 2,"Perennials" = 2))

table_trait_coverage_nb_sp
cat(table_trait_coverage_nb_sp, file = "draft/table_trait_coverage_richness.doc")



# nb of species ####
# NB: too complex for the paper, but useful working tool
details <- F
if (details == T){
  
  table_nb_sp <- NB_SP %>% 
    select(LifeHistory,treatment,trait,in_ab,in_trait,in_both) %>% 
    arrange(LifeHistory,treatment,factor(trait,levels = c(traits, "Multivariate"))) %>% 
    mutate(trait = case_when(trait %in% c("SLA","LDMC","L_Area") ~ "Leaf morphological traits",
                             trait %in% c("LCC","LNC")~"Leaf chemical traits",
                             TRUE ~trait)) %>%
    unique() %>% 
    mutate(treatment = if_else(treatment == "Fer","Int.","Ext.")) %>% 
    mutate(percent_sp = round(in_both/(in_both+in_ab),digits =2))
  
  table_trait_coverage_nb_sp <- table_nb_sp %>% 
    kableExtra::kable( escape = F,
                       col.names = c("Life History ",
                                     "Management regime",
                                     "Trait",
                                     "in transects only",
                                     "in traits only",
                                     "in both",
                                     "Percentage of species covered")) %>%
    kableExtra::kable_styling("hover", full_width = F)  %>% 
    kableExtra::add_header_above(c(" " = 3,"Number of species" = 3," "=1))
  
  table_trait_coverage_nb_sp 
  
  cat(table_trait_coverage_nb_sp, file = "draft/table_trait_coverage_richness.doc")
}







  

  
#   
#   
#   
# # test sur ldmc et nb sp ok
# sp_nat_ab <- ab_nat %>% 
#   filter(depth == "S") %>% 
#   group_by(species,code_sp,LifeHistory) %>% 
#   summarize(sp_abundance = sum(abundance)) %>% 
#   filter(LifeHistory=="annual") %>% 
#   pull(code_sp) %>% 
#   unique()
# sp_nat_traits <- MEAN %>% 
#   filter(treatment == "Nat" & LifeHistory == "annual") %>% 
#   filter(!is.na(LDMC)) %>% 
#   pull(code_sp) %>% 
#   unique()
# 
# int <- intersect(sp_nat_traits,sp_nat_ab)
# diff <- setdiff(sp_nat_ab,sp_nat_traits)
# length(int) /(length(int)+length(diff))
