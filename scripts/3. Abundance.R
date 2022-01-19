source("scripts/1. Packages.R")

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")
names_LH <- MEAN %>% 
  select(code_sp,species,LifeHistory) %>% 
  unique()

# Fertilized treatment ####
ABUNDANCE <- read.csv2("data/abundance/pooled_abundance_data.csv") %>%  
  mutate(treatment = case_when(grepl("C",paddock) ~"Fer",
                               grepl("P",paddock) ~"Nat",
                               grepl("N",paddock) ~"Nat",
                               grepl("T",paddock) ~"Tem"))

ab_diachro_2005 <- ABUNDANCE %>% 
  filter(dataset == "Diachro") %>% 
  filter(year == 2005) %>%  # /!\ No measure available after 2005! Should I keep 2006, or several years?
  filter(treatment == "Fer") %>% 
  # add relative abundance
  group_by(id_transect_quadrat) %>% 
  mutate(relat_ab = abundance/sum(abundance))

write.csv2(ab_diachro_2005,"outputs/data/abundance_fertile.csv",row.names=F)

# Unfertilized treatment #### 
# (Maud's relevés)
# NB : I should take the sheet "abondance par ligne" to be homogeneous with data from diachro relevés.
# Sinon il faut agréger au niveau du diachro, mais je n'aurai plus assez d'observations... à moins de prendre plusieurs années.
soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))


# # Par parcelle
# ab_maud <- read.xlsx("data/abundance/maud_Relevés d'abondance La Fage Juin 2009.xlsx",
#                      sheet = "abondances par parcelle",
#                      startRow = 1, colNames = TRUE, rowNames = F) %>%
#   remove_rownames() %>%
#   gather(Species,abundance,-plot) %>%
#   mutate(Species=str_replace(Species,"_"," ")) %>%
#   mutate(depth = str_sub(plot,start = 1L,end=1L)) %>%
#   mutate(paddock = str_sub(plot,start = 3L,end=-1L)) %>%
#   select(-plot) %>%
#   full_join(soil_Maud,.,by=c("paddock","depth")) %>%
#   filter(abundance >0 ) %>%
#   merge(names_LH, by="Species") %>%
#   # add relative abundance
#   group_by(depth,paddock) %>%
#   mutate(relat_ab = abundance/sum(abundance)) %>%
#   dplyr::rename(species = Species, code_sp = Code_Sp)

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
  select(-Ligne) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  filter(abundance >0 ) %>% 
  merge(names_LH, by="species") %>% 
  # add relative abundance
  group_by(depth,paddock) %>% 
  mutate(relat_ab = abundance/sum(abundance))

write.csv2(ab_maud,"outputs/data/abundance_natif.csv",row.names=F)






# Merge to traits (do it in another file) ####
# To be cleaned from now on.
# Duplicated rows containing the same species
# Better to add info about lifeform after... and to merge traits and 
ab_Maud_traits <- ab_maud %>% 
  full_join(MEAN,by="Species") %>% 
  filter(Trtmt == "Nat") %>% 
  filter(abundance >0 )


  dplyr::rename(species=Species,code_sp=Code_Sp) %>% 
  filter(!is.na(abundance)) %>% 
  mutate(depth = str_sub(plot,start = 1L,end=1L)) %>% 
  mutate(paddock = str_sub(plot,start = 3L,end=-1L)) %>%
  select(-plot) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  filter(Trtmt == "Nat") %>% 
  filter(abundance >0 )



ab_Maud_traits2 <- read.xlsx("data/abundance/maud_Relevés d'abondance La Fage Juin 2009.xlsx", sheet = "abondances par parcelle", 
                            startRow = 1, colNames = TRUE, rowNames = F) %>% 
  remove_rownames() %>%  
  gather(Species,abundance,-plot) %>%
  mutate(Species=str_replace(Species,"_"," ")) %>% 
  full_join(MEAN,by="Species") %>% 
  dplyr::rename(species=Species,code_sp=Code_Sp) %>% 
  filter(!is.na(abundance)) %>% 
  mutate(depth = str_sub(plot,start = 1L,end=1L)) %>% 
  mutate(paddock = str_sub(plot,start = 3L,end=-1L)) %>%
  select(-plot) %>% 
  full_join(soil_Maud,.,by=c("paddock","depth")) %>% 
  filter(Trtmt == "Nat") %>% 
  filter(abundance >0 )
write.csv2(ab_Maud_traits,"outputs/data/ab_Maud_traits.csv",row.names = F)

ab_Maud_traits %>% 
  select(code_sp,SLA,LDMC,L_Area) %>% 
  unique() %>% 
  filter(!(is.na(SLA))) %>%
  relocate(L_Area,LDMC,SLA) %>% 
  mutate(L_Area=L_Area*100) %>% # to change unit from cm² to mm²
  mutate(LDMC=LDMC/1000*100) %>% 
  write.csv2("outputs/data/Pierce CSR/Traits_Maud.csv",row.names=F)