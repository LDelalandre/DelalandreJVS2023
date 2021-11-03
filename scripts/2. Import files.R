source("scripts/1. Packages.R")

# Plant traits ####
data_file <- "LaFage_PlantTraitsDP_vp.xlsx"

LeafMorpho1 <-  read.xlsx(paste0("data/traits/",data_file), sheet = "LeafMorpho_traits", startRow = 1, colNames = TRUE)  
LeafMorpho_leo <- read.csv2("data/traits/leafMorpho_3.csv")
LeafMorpho <- rbind(LeafMorpho,LeafMorpho_leo)
# LeafDimensions <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafDimensions (àsupprimer)", startRow = 1, colNames = TRUE)
LeafCN <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafC&N", startRow = 1, colNames = TRUE) 
LeafP <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafP", startRow = 1, colNames = TRUE) 
Leaf13C <- read.xlsx(paste0("data/traits/",data_file), sheet = "Leaf13C", startRow = 1, colNames = TRUE) 
Biovolume <- read.xlsx(paste0("data/traits/",data_file), sheet = "Biovolume", startRow = 1, colNames = TRUE) %>% 
  mutate(Hrepro = as.numeric(Hrepro))
Pheno <- read.xlsx(paste0("data/traits/",data_file), sheet = "Pheno", startRow = 1, colNames = TRUE) %>% 
  mutate(Rep = "None") %>% 
  mutate(Day = "None")
Seed <- read.xlsx(paste0("data/traits/",data_file), sheet = "Seed", startRow = 1, colNames = TRUE) 


# Abundance ####
ABUNDANCE <- read.csv2("data/abundance/pooled_abundance_data.csv") %>%  # NB checker comment j'ai construit ce jeu de données.
# Notamment, est-ce que j'ai bien les mêmes données d'abondance que Maud dans son papier de 2012 ? (ça pourrait expliquer mes résultats différents).
  mutate(treatment = case_when(grepl("C",paddock) ~"Fer",
                               grepl("P",paddock) ~"Nat",
                               grepl("N",paddock) ~"Nat",
                               grepl("T",paddock) ~"Tem"))
