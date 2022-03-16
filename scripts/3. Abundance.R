source("scripts/1. Packages.R")

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")
names_LH <- MEAN %>% # correspondence name, short name, lifehistory
  select(code_sp,species,LifeHistory) %>% 
  unique()

# Fertilized treatment ####
ABUNDANCE <- read.csv2("data/abundance/pooled_abundance_data.csv") %>%  
  mutate(treatment = case_when(grepl("C",paddock) ~"Fer",
                               grepl("P",paddock) ~"Nat",
                               grepl("N",paddock) ~"Nat",
                               grepl("T",paddock) ~"Tem"))
# I use abondance data from the 2005 survey
ab_diachro_2005 <- ABUNDANCE %>% 
  filter(dataset == "Diachro") %>% 
  filter(year == 2005) %>%  # /!\ No measure available after 2005! Should I keep 2006, or several years?
  filter(treatment == "Fer") %>%
  # add relative abundance
  group_by(id_transect_quadrat) %>% 
  mutate(relat_ab = abundance/sum(abundance))

write.csv2(ab_diachro_2005,"outputs/data/abundance_fertile.csv",row.names=F)


# Regarder dans diachro nombre points contacts f(traitement)
ABUNDANCE %>% 
  filter(dataset == "Diachro") %>% 
  # filter(year == 2004) %>%
  group_by(id_transect_quadrat,year) %>% 
  mutate(sumab = sum(abundance)) %>% 
  ggplot(aes(x=treatment,y=sumab))+
  geom_boxplot() +
  facet_wrap(~year) 
  # ggsave("outputs/plots/comparison_intercept_nat_fer.png",width = 10,height=10)

# Regarder la météo
meteo <- read.csv2("data/environment/Meteo_LaFage_1973-2006.csv")
meteo %>% 
  group_by(AN) %>% 
  summarize(sum_RR = sum(RR), mean_TX = mean(TX), mean_RGC = mean(RGC)) %>% 
  ggplot(aes(x=AN,y=mean_TX))+
  geom_line()

# Unfertilized treatment #### 
# (Maud's relevés)
soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))
# par ligne
ab_maud <- read.xlsx("data/abundance/maud_Relevés d'abondance La Fage Juin 2009.xlsx", 
                     sheet = "abondance par ligne", 
                     startRow = 1, colNames = TRUE, rowNames = F) %>% 
  remove_rownames() %>%  
  gather(species,abundance,-Ligne) %>%
  mutate(species=str_replace(species,"_"," ")) %>% 
  mutate(depth = str_sub(Ligne,start = -2L,end=-2L)) %>% 
  mutate(line = str_sub(Ligne,start = -1L,end=-1L)) %>% 
  mutate(paddock = str_sub(Ligne,start = 1L,end=-3L)) %>%
  # select(-Ligne) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  filter(abundance >0 ) %>% 
  merge(names_LH, by="species") %>% 
  # add relative abundance
  group_by(line,depth,paddock) %>% 
  filter(!code_sp =="STRIFSCAB")%>% # /!\ à corriger à la source !!
  mutate(relat_ab = abundance/sum(abundance))

write.csv2(ab_maud,"outputs/data/abundance_natif.csv",row.names=F)
