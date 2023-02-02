# This script imports trait values measured before my PhD 
# and includes data that I measured in spring 2021.

library(openxlsx)
library(tidyverse)
# Characteristics of the species  ####
names_LH <- read.csv2("data/species_names_lifehistory.csv") %>% 
  rename(Code_Sp = code_sp)
names_family <- read.csv2("data/traits/names and family.csv")


data_file <- "LaFage_PlantTraitsDP_vp.xlsx"

# Leaf Morphological traits ####
LeafMorpho1 <-  read.xlsx(paste0("data/traits/",data_file), sheet = "LeafMorpho_traits", startRow = 1, colNames = TRUE)  
LeafMorpho_leo <- read.csv2("data/traits/leafMorpho_3.csv",encoding = "latin1") %>% 
  filter(!(Code_Sp %in% c("EROPVERN","STELMEDI"))) # senescent leaves
LeafMorpho <- rbind(LeafMorpho1,LeafMorpho_leo)

species_measured <- LeafMorpho %>% 
  filter(Treatment=="Nat_Sab") %>% 
  filter(LifeForm1=="The") %>% 
  arrange(Code_Sp) %>% 
  select(Code_Sp) %>% 
  unique()


# Leaf chemical traits ####
# LeafDimensions <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafDimensions (àsupprimer)", startRow = 1, colNames = TRUE)
LeafCN1 <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafC&N", startRow = 1, colNames = TRUE) 
LeafP <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafP", startRow = 1, colNames = TRUE) 
Leaf13C1 <- read.xlsx(paste0("data/traits/",data_file), sheet = "Leaf13C", startRow = 1, colNames = TRUE)
# measurements C, N, 13C spring 2021
Leafchim_leo <- read.csv2(paste0("data/traits/leafchim_leo.csv"),encoding ="latin1") %>% 
  filter(!(Code_Sp %in% c("EROPVERN","STELMEDI"))) # senescent leaves
# Species name and LifeForm
L2 <- Leafchim_leo %>% 
  merge(names_LH %>% rename(Species=species),by="Code_Sp") %>% 
  mutate(LifeForm2=NA)
# Day of measurement
L3 <- LeafMorpho_leo %>% 
  select(Code_Sp,Day,Plot,Treatment) %>% 
  unique() %>% 
  merge(L2,by=c("Code_Sp","Plot"))
# Family
fam <- read.csv2("data/traits/names and family.csv")
L4 <- merge(L3,fam)
# Rep
L5 <- L4 %>% 
  group_by(Species,Plot) %>% 
  mutate(Rep = paste("RCN",Rep,sep=""))

LeafCN_leo <- L5[,colnames(LeafCN1)] %>% 
  mutate(LCC = 10*LCC) %>% 
  mutate(LNC = 10*LNC) # to change the unit from % to mg/g.
Leaf13C_leo <- L5[,colnames(Leaf13C1)]

LeafCN <- rbind(LeafCN1,LeafCN_leo) %>% 
  mutate(LifeForm1 = if_else(Code_Sp == "SANGMINO","Hem",LifeForm1)) %>% 
  mutate(LifeForm1 = if_else(Code_Sp == "ANTHVULN","Hem",LifeForm1)) %>%
  mutate(LifeForm1 = if_else(Code_Sp == "LOTUCORN","Hem",LifeForm1))
Leaf13C <- rbind(Leaf13C1,Leaf13C_leo)

# Biovolume ####
Biovolume1 <- read.xlsx(paste0("data/traits/",data_file), sheet = "Biovolume", startRow = 1, colNames = TRUE) %>% 
  mutate(Hrepro = as.numeric(Hrepro))

Biovolume_leo <- read.csv2("data/traits/Biovolume_leo.csv",encoding = "latin1" )
# Sort min and max diameter
nrows <- dim(Biovolume_leo)[1]
for (i in c(1:nrows)){
  d1 <- Biovolume_leo[i,]$Dmax
  d2 <- Biovolume_leo[i,]$Dmin
  
  Biovolume_leo[i,]$Dmax <- max(d1,d2)
  Biovolume_leo[i,]$Dmin <- min(d1,d2)
}

Biovolume <- rbind(Biovolume1 , Biovolume_leo)

# Phenology ####
Pheno1 <- read.xlsx(paste0("data/traits/",data_file), sheet = "Pheno", startRow = 1, colNames = TRUE) %>% 
  mutate(Rep = "None") %>% 
  mutate(Day = "None")

day_month <- read.table("data/phenology/day_month.txt",header=T)
get_day_of_year <- function(Month,Day2){
  day_month %>% 
    filter(month == Month) %>% 
    filter(day_of_month == Day2) %>% 
    pull(day_of_year)
}
colnames <- colnames(Pheno1)

pheno_leo <- read.xlsx("data/phenology/Pheno_leo.xlsx", sheet = "rawdata_seeds", startRow = 1, colNames = TRUE) %>% 
  mutate(Day2 = as.Date(date - 25569, origin = "1970-01-01")) %>% 
  mutate(Day = str_replace_all(Day2,"-","")) %>% 
  filter(!(plot=="C1")) # wrong pheno estimation in the fertile treatment

pheno_leo_DP <- pheno_leo %>% 
  mutate(Code_Sp = toupper(code_sp)) %>%
  filter(!(Code_Sp == "rhinpumi")) %>% # mesuré trop tard
  mutate(Site= "La Fage", Block = "None",Plot = plot, Treatment = "Nat_Sab",Year = 2021,
  ) %>% 
  left_join(names_LH,by="Code_Sp") %>% 
  left_join(names_family %>% rename(species = Species),by="species") %>% 
  separate(col = Day2,into = c("Year","Month","Day")) %>% 
  mutate(Year = as.numeric(Year), Month = as.numeric(Month), Day = as.numeric(Day)) %>% 
  mutate(Disp = map2_dbl(Month,Day,get_day_of_year)) %>% 
  select(-c(code_sp,date,plot)) %>% 
  mutate(LifeForm2 = NA) %>%  # A compléter pour le data paper
  mutate(Flo = NA, Mat_Per = NA, nameOfProject = "Annuals",measurementDeterminedBy = "Léo Delalandre") %>% 
  mutate(Rep = "None") %>% 
  group_by(species) %>% 
  filter(Disp == min (Disp)) %>%  # NB: une valeur de phéno par espèce*traitment dans la BDD ; je prends la min !
  rename(Species = species) %>% 
  select(all_of(colnames)) %>% 
  mutate(Species = if_else(Code_Sp=="SAXITRIDA","Saxifraga tridactylites",Species))

Pheno <- rbind(Pheno1,pheno_leo_DP)
# Virer filago pyramidata, pour lequel j'ai peut-être estimé trop tard la dispersion (et elle avait aussi lieu
# dans le fertile!!)

# Seed mass ####
Seed1 <- read.xlsx(paste0("data/traits/",data_file), sheet = "Seed", startRow = 1, colNames = TRUE) 

Seed_leo <- read.csv2("data/traits/Seed_leo_lila.csv",encoding="latin1") %>%
  mutate(Rep = paste0("Rsd",Rep)) %>% 
  mutate(Site = "La Fage",Block = "None",Year = 2022, 
         LifeForm2 = "None",
         SeedMass = TotSeedMass.mg./SeedNb, Mat = NA, 
         nameOfProject = "Annuals",measurementDeterminedBy = "Léo Delalandre")  %>% 
  merge(names_LH , by = "Code_Sp") %>%
  rename(Species = species) %>% 
  merge(names_family,by = "Species") %>% 
  mutate(Treatment = if_else(Plot %in% c("P8","P1","P10"), "Nat_Sab","Fer_Clc")) %>% 
  select(all_of(colnames(Seed1)))

Seed <- rbind(Seed1,Seed_leo) 
  # filter(!(Code_Sp == "GERADISS")) # truc chelou avec les masses de graines, à voir



# Data paper ####
# checker données pour être sûr
# exporter en csv
# les copier_coller dans la dernière version des fichiers de La Fage après avoir demandé à Eric

LeafMorpho_leo %>%
  mutate(Species = case_when(
    Species == "Myosostis ramosissima subsp. ramosissima" ~ "Myosotis ramosissima subsp. ramosissima", # s en trop
    Species == "Vicia sativa sativa" ~ "Vicia sativa subsp. sativa",
    TRUE ~ Species)) %>% 
  filter(!(Species == "Medicago rigidula")) # je ne suis pas sûr de l'identification
# corriger idem ensuite
LeafCN_leo
Leaf13C_leo
Biovolume_leo
pheno_leo_DP
Seed_leo

# species Léo
data.frame( Species = c(
  LeafMorpho_leo$Species,
  LeafCN_leo $Species,
  Leaf13C_leo$Species,
  Biovolume_leo$Species,
  pheno_leo_DP$Species,
  Seed_leo$Species) %>% unique() ) %>% 
  arrange(Species) %>% 
  write.csv2("outputs/data/species_leo.csv",row.names=F)
