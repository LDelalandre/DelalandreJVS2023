source("scripts/1. Packages.R")

# Plant traits ####
data_file <- "LaFage_PlantTraitsDP_vp+KB.xlsx"

LeafMorpho <-  read.xlsx(paste0("data/",data_file), sheet = "LeafMorpho_traits", startRow = 1, colNames = TRUE)  
LeafDimensions <- read.xlsx(paste0("data/",data_file), sheet = "LeafDimensions (àsupprimer)", startRow = 1, colNames = TRUE)
LeafCN <- read.xlsx(paste0("data/",data_file), sheet = "LeafC&N", startRow = 1, colNames = TRUE) 
LeafP <- read.xlsx(paste0("data/",data_file), sheet = "LeafP", startRow = 1, colNames = TRUE) 
Biovolume <- read.xlsx(paste0("data/",data_file), sheet = "Biovolume", startRow = 1, colNames = TRUE)
Seed <- read.xlsx(paste0("data/",data_file), sheet = "Seed", startRow = 1, colNames = TRUE) 
