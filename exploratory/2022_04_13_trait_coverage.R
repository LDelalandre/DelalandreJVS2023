library(tidyverse)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C","LPC",
            "Hveg"  ,    "Hrepro"   , "Dmax"  , #    "Dmin" ,
            "Flo","Disp","Mat_Per", #"Mat",
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

MEAN_annual <- MEAN %>% 
  filter(LifeHistory=="annual")

MEAN_annuals_coverage_nat <- MEAN_annual %>% 
  filter(treatment=="Nat") %>% 
  summarise_each(funs(100*mean(is.na(.)))) %>%
  select(-c(species,treatment,code_sp,Form,LifeHistory,LifeForm1)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("trait") %>% 
  mutate(coverage_nat=100 - V1) %>% 
  select(-V1) 
  
MEAN_annuals_coverage_fer <- MEAN_annual %>% 
  filter(treatment=="Fer") %>% 
  summarise_each(funs(100*mean(is.na(.)))) %>%
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
ftrait <- "SeedMass"

MEAN_annual %>% 
  filter(treatment=="Nat") %>% 
  filter(is.na(SeedMass)) %>% 
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

# Natif
trait_available_nat <- MEAN_annual %>% 
  filter(treatment=="Nat") %>% 
  mutate(trait = Hrepro) %>%  # choose the trait here
  select(species,code_sp,trait) %>% 
  full_join(relat_ab_nat,by=c("species","code_sp")) %>% 
  mutate(trait_available = if_else(is.na(trait),"no","yes")) %>% 
  select(species,code_sp,trait_available,sp_relat_ab_nat)

trait_available_nat %>% 
  ggplot(aes(x=trait_available,y=sp_relat_ab_nat)) +
  geom_histogram(position='stack',stat="identity") +
  ggtitle("Natif")



# Fertile
trait_available_fer <- MEAN_annual %>% 
  filter(treatment=="Fer") %>% 
  mutate(trait = Hrepro) %>%  # choose the trait here
  select(species,code_sp,trait) %>% 
  full_join(relat_ab_fer,by=c("species","code_sp")) %>% 
  mutate(trait_available = if_else(is.na(trait),"no","yes")) %>% 
  select(species,code_sp,trait_available,sp_relat_ab_fer)

trait_available_fer %>% 
  ggplot(aes(x=trait_available,y=sp_relat_ab_fer)) +
  geom_histogram(position='stack',stat="identity") +
  ggtitle("Fertile")




