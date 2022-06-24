library(tidyverse)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab.csv")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , "Dmax"  , #    "Dmin" ,"Hveg"  , 
            "Disp",#"Mat_Per", #"Mat","Flo",
            "SeedMass"
)

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(line = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") 

ab_fer_ann <- ab_fer %>% 
  filter(LifeForm1 == "The")
ab_nat_ann <- ab_nat %>% 
  filter(LifeHistory == "annual")

# Coverage for each trait ####
# For each trait: proportion of the annual species occurring in the trait data for which we have measurement
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


# Abundance and trait coverage info together ####
# Abundance is computed as 

traits_abundance_coverage <- NULL
for (ftrait in traits){
  # Natif
  trait_available_nat <- MEAN_annual %>% 
    filter(treatment=="Nat") %>% 
    mutate(trait = get(ftrait)) %>%  # choose the trait here
    select(species,code_sp,trait) %>% 
    full_join(relat_ab_nat,by=c("species","code_sp")) %>% 
    mutate(trait_available = if_else(is.na(trait),"no","yes")) %>%
    # mutate(trait_available = if_else(trait==0,"no","yes")) %>% 
    select(species,code_sp,trait_available,sp_relat_ab_nat) %>% 
    replace(is.na(.),0)
  # 
  # trait_available_nat %>%
  #   ggplot(aes(x=trait_available,y=sp_relat_ab_nat)) +
  #   geom_histogram(position='stack',stat="identity") +
  #   ggtitle(paste("Natif",ftrait))
  
  
  # cumulated abundance for species for which we have the trait
  Tnat <- trait_available_nat %>% 
    group_by(trait_available) %>% 
    summarise(abundance_covered = sum(sp_relat_ab_nat)) %>% 
    spread(key = trait_available, value = abundance_covered) %>% 
    mutate(trait = ftrait) %>% 
    mutate(treatment = "Nat")
  
  if (is.null(Tnat$yes) ){ 
    Tnat$yes <- 0
  }
  
  Tnat <- Tnat %>% select(all_of(c("trait","treatment","yes","no")))
  
  
  # Fertile
  trait_available_fer <- MEAN_annual %>% 
    filter(treatment=="Fer") %>% 
    mutate(trait = get(ftrait)) %>%  # choose the trait here
    select(species,code_sp,trait) %>% 
    full_join(relat_ab_fer,by=c("species","code_sp")) %>% 
    mutate(trait_available = if_else(is.na(trait),"no","yes")) %>% 
    select(species,code_sp,trait_available,sp_relat_ab_fer)%>% 
    replace(is.na(.),0)
  # 
  # trait_available_fer %>% 
  #   ggplot(aes(x=trait_available,y=sp_relat_ab_fer)) +
  #   geom_histogram(position='stack',stat="identity") +
  #   ggtitle("Fertile")
  
  # cumulated abundance for species for which we have the trait
  Tfer <- trait_available_fer %>% 
    group_by(trait_available) %>% 
    summarise(abundance_covered = sum(sp_relat_ab_fer)) %>% 
    spread(key = trait_available, value = abundance_covered) %>% 
    mutate(trait = ftrait) %>% 
    mutate(treatment = "Fer")
  
  Tfer <- Tfer %>% select(all_of(c("trait","treatment","yes","no")))
  
  traits_abundance_coverage <- rbind(traits_abundance_coverage, Tnat )
  traits_abundance_coverage <- rbind(traits_abundance_coverage, Tfer )
}


traits_abundance_coverage %>% 
  ggplot(aes(x=trait,y=yes))+
  geom_point() +
  facet_wrap(~treatment)


traits_abundance_coverage %>% 
  select(-no) %>% 
  spread(key=treatment,value=yes) %>% 
  arrange(factor(trait,levels = traits))


table_mod <- df_mod %>%
  merge(trait_unit,by="Trait") %>% 
  select(Trait,Unit,Intercept,Estimate,p.value) %>% 
  arrange(factor(Trait,levels = vec_traits)) %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Trait",
                                   "Unit",
                                   "Value in G+F",
                                   "Value in GUs",
                                   "p.value")) %>%
  kableExtra::kable_styling("hover", full_width = F)

table_mod

cat(table_mod, file = "draft/mixed_model_intra_annual.doc")
