library(tidyverse)
library(dplyr)
detach("package:MASS")    

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_completed_seed_mass.csv") %>%
  filter(!is.na(SLA)) %>% 
  filter(!(species== "Geranium dissectum - pétiole"))

trait_unit <- read.csv2("data/trait_unit.csv",encoding = "latin1")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "Disp",#"Mat_Per", #"Mat","Flo",
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
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance))

# Regional relative abundance of each annual species within the annual guild
relat_ab_fer_ann <- ab_fer_ann %>% 
  group_by(species,code_sp) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance))

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
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance))

# Regional relative abundance of each annual species within the annual guild
relat_ab_nat_ann <- ab_nat_ann %>%
  filter(depth == "S") %>% 
  group_by(species,code_sp) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance))

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
      select(species,code_sp,trait,LifeHistory) %>% 
      full_join(relat_ab_nat,by=c("species","code_sp")) %>%  # or relat_ab_nat_ann
      
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
      select(species,LifeHistory,code_sp,trait) %>% 
      full_join(relat_ab_fer,by=c("species","code_sp")) %>% # or relat_ab_fer_ann
      
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


# Species in the abundance, but not in the trait, data? ####

sp_ab_fer <- ab_fer %>% pull(code_sp) %>% unique()
sp_trait_fer <- MEAN %>% filter(treatment == "Fer") %>%  pull(code_sp) %>% unique()

setdiff(sp_ab_fer,sp_trait_fer) # dans le premier mais pas dans le deuxième
setdiff(sp_trait_fer,sp_ab_fer)


#_______________________________________________________________________________
# Mass grave ?

# Coverage for each trait ####
# For each trait: proportion of the annual species occurring in the trait data for which we have measurement
trait_unit <- read.csv2("data/trait_unit.csv")

MEAN_annual <- MEAN %>% 
  filter(LifeHistory=="annual")

MEAN_annuals_coverage_nat <- MEAN_annual %>% 
  filter(treatment=="Nat") %>% 
  summarise_each(funs(100*mean(is.na(.)))) %>% # counts the proportion of lines with NAs for each column (= each trait)
  select(-c(species,treatment,code_sp,Form,LifeHistory,LifeForm1)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("trait") %>% 
  mutate(coverage_nat=100 - V1) %>% # 100 - proportion of NAs
  select(-V1) 

MEAN_annuals_coverage_fer <- MEAN_annual %>% 
  filter(treatment=="Fer") %>% 
  summarise_each(funs(100*mean(is.na(.)))) %>% # counts the proportion of lines with NAs for each column (= each trait)
  select(-c(species,treatment,code_sp,Form,LifeHistory,LifeForm1)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("trait") %>% 
  mutate(coverage_fer=100 - V1) %>% 
  select(-V1)

coverage <- merge(MEAN_annuals_coverage_fer,MEAN_annuals_coverage_nat,by="trait") %>% 
  arrange(desc(coverage_nat)) %>% 
  filter(trait %in% traits) 
write.csv2(coverage,"outputs/data/trait_coverage_annuals.csv",row.names=F)


# Species*trtmt not measured ####
ftrait <- "Hveg"

MEAN_annual %>% 
  filter(treatment=="Nat") %>% 
  filter(is.na(get(ftrait))) %>% 
  arrange(code_sp) %>% 
  select(code_sp) %>% 
  unique()
# Je peux récupérer les données de masse des graines pour 11 espèces du natif. Top !

# Pour les traits phéno et de hauteur, regarder dans quelle mesure on a les annuelles les plus abondantes
MEAN_annual %>% 
  filter(treatment=="Fer") %>% 
  filter(!is.na(Disp)) %>% 
  arrange(code_sp) %>% 
  select(code_sp) %>% 
  unique()

# Relative abundance of annuals is computed for all the transects by dividing the sum of 
# annuals' abundance by the sum of all abundances.
# Then, we divide a species' cumulated relative abundance by the cumulated relative abundance of annuals
relat_ab_fer <- ab_fer_ann %>%
  mutate(tot_relat_ab = sum(relat_ab)) %>%  # total relative abundance of annuals
  group_by(species,code_sp) %>% 
  summarize(sp_relat_ab_fer = sum(relat_ab)/tot_relat_ab) %>% 
  unique()

relat_ab_nat <- ab_nat_ann %>%
  mutate(tot_relat_ab = sum(relat_ab)) %>%  # total relative abundance of annuals
  group_by(species,code_sp) %>% 
  summarize(sp_relat_ab_nat = sum(relat_ab)/tot_relat_ab) %>% 
  unique()

# relat_abundances <- full_join(relat_ab_fer,relat_ab_nat,by=c("species","code_sp"))


