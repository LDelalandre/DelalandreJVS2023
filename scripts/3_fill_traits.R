# Take the trait data averaged by species*origin
# Complete missing trait values

library(tidyverse)

# Trait values computed in the G+F and GU(S+I) conditions, i.e. removing GUd for the latter
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_completed.csv")%>%
  filter(!is.na(SLA)) %>% 
  filter(!(species== "Geranium dissectum - pétiole")) %>% 
  filter(!species == "Geranium dissectum - pÃ©tiole")

# Trait values computed in the G+F and GU conditions
MEAN_no_subset <- read.csv2("outputs/data/mean_attribute_per_treatment.csv",encoding = "latin1") %>%
  filter(!is.na(SLA)) %>%
  filter(!(LifeForm1 %in% c("DPh","EPh")))%>% 
  filter(!(species== "Geranium dissectum - pétiole"))%>% 
  filter(!species == "Geranium dissectum - pÃ©tiole") 



#___________________________________________________
# Seed Mass ####

# NB: GERADISS et MYOSRAMO: on a la seed mass
# mais myosramo-ramo: pb de code ?
# GERADISS: refaire refaire tourner le script pour moyenner traits (et refaire CSR Pierce)

# j'ajoute la masse des graines des mêmes espèces mesurées dans tout le natif
# aux traits mesurés dans le GUs (avec l'approx. que la masse individuelle des graines est peu plastique... à tester!!)
MEAN_no_subset_seed_mass <-  MEAN_no_subset %>% 
  select(code_sp,treatment,LifeHistory,SeedMass) %>% 
  group_by(code_sp)

# add seed mass measured in the rest of GU to missing seed mass in GU
MEAN_fill <- MEAN %>%
  filter(is.na(SeedMass)) %>% 
  select(-SeedMass) %>% 
  merge(MEAN_no_subset_seed_mass, by = c("code_sp","treatment","LifeHistory"))

MEAN_completed <- MEAN %>%
  filter(!is.na(SeedMass)) %>%
  rbind(MEAN_fill)

# complete seed mass from fer
# mean seed mass in fer
MEAN_fer <- MEAN_completed %>%
  filter( treatment=="Fer") %>% 
  select(code_sp,SeedMass)

# add seed mass in fer to missing seed mass in nat
sp_missing_seed_mass_nat <- MEAN_completed %>% 
  filter(treatment == "Nat") %>% 
  filter(is.na(SeedMass)) %>% 
  pull(code_sp)

MEAN_nat <- MEAN_completed %>%
  filter( treatment=="Nat") %>% 
  filter(code_sp %in% sp_missing_seed_mass_nat) %>% 
  select(-SeedMass) %>% 
  left_join(MEAN_fer, by = c("code_sp"))

MEAN_completed2_nat <- MEAN_completed %>% 
  filter(treatment=="Nat") %>% 
  filter(!(code_sp %in% MEAN_nat$code_sp)) %>% 
  rbind(MEAN_nat)

MEAN_completed2 <- MEAN_completed %>% 
  filter(treatment %in% "Fer") %>% 
  rbind(MEAN_completed2_nat)

#___________________________________________________
# Phenology and height ####
traits_flore <- read.table("data/traits/flores/TRAIT_ESP_FLORE.txt",header=T) %>% 
  select(CODE_ESP,FLO_FLORE,FRU_FLORE,H_FLORE) %>% 
  rename(code_sp = CODE_ESP)

MEAN_completed2_flore <- left_join(MEAN_completed2,traits_flore)

# check accuracy of estimation by data from flora
MEAN_completed2_flore %>% 
  ggplot(aes(x=Hrepro,y=H_FLORE)) + 
  geom_point()

MEAN_completed2_flore %>% 
  ggplot(aes(x=disp,y=FRU_FLORE)) + 
  geom_point()


#___________________________________________________
# test seed mass low plasticity
ftrait <- "SeedMass"
seed_fer_nat <- MEAN_no_subset %>% 
  select(code_sp,treatment,all_of(ftrait)) %>% 
  spread(key = treatment,value = ftrait) 

seed_fer_nat %>% 
  ggplot(aes(x=log(Fer),y=log(Nat))) +
  # ggplot(aes(x=Fer,y=Nat))+
  geom_point() +
  geom_abline(slope = 1,intercept=0) +
  geom_smooth(method="lm")

mod <- lm(log(Nat) ~ log(Fer), data = seed_fer_nat)
# plot(mod)
anova(mod)
sum <- summary(mod)

mod$coefficients
sum$adj.r.squared
# " OU BIEN: regarder la masse des graines comme une fonction des espèces et du traitement
# pour montrer que la donnée espèce explique une part de variance bien plus grande


write.csv2(MEAN_completed2_flore,
           "outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_completed_seed_mass_flore.csv",
           row.names=F)
