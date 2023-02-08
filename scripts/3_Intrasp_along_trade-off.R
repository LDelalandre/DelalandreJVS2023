library(tidyverse)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv")
MEAN_site <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_site_level.csv")
traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , #"Dmax"  , #    "Dmin" ,"Hveg"  , "H_FLORE",#
            "Disp", #"Mat_Per", #"Mat","Flo","FLO_FLORE", #
            "SeedMass"
)


# reasons for outliers ?
source("scripts/Data_traits.R") # Load traits per group of traits (e.g. LeafMorpho)
fdata <- LeafCN %>%
  filter(Treatment %in% c("Fer_Clc","Fer_Dlm","Nat_Sab","Nat_Int")) %>% 
  mutate(treatment = str_sub(Treatment,1,3) ) 

ftrait <- "LCC"
fdata %>% 
  group_by(treatment,Code_Sp) %>% 
  summarise(n = n(),
            mean = mean(get(ftrait),na.rm=T),
            sd = sd(get(ftrait),na.rm=T),
            cv = sd/mean) %>% 
  arrange(Code_Sp) %>% 
  View()

fdata %>% 
  filter(Code_Sp=="SESEMONT") %>%
  group_by(Code_Sp,treatment) %>% 
  ggplot(aes_string(x=ftrait)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~treatment)

# ATTENTION, IL Y A DES OUTLIERS SUR CETTE VAR INTRA!
ftrait <- "LCC"
# for (ftrait in traits){
  # compute trait difference and ratio across the two treatments
  intrasp_var <- MEAN %>% 
    # filter(LifeHistory=="annual") %>% 
    select(species,code_sp,LifeHistory,treatment,all_of(ftrait)) %>% 
    spread(key=treatment,value=ftrait) %>% 
    mutate(trait = ftrait) %>% 
    na.omit() %>%  # keep only sp measured in the 2 treatments
    mutate(ratio = Nat / Fer) %>% 
    mutate(diff = Nat-Fer) %>% 
    mutate(RDPI = diff/Fer)
  
  intrasp_toplot <- MEAN_site %>% 
    select(code_sp,SLA) %>% 
    merge(intrasp_var)
  intrasp_toplot %>% 
    ggplot(aes(x=SLA,y=diff)) +
    geom_point()

# }