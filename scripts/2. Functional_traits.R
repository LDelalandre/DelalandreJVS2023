source("scripts/1. Packages.R")


# 1) Import plant trait values ####
data_file <- "LaFage_PlantTraitsDP_vp.xlsx"

LeafMorpho1 <-  read.xlsx(paste0("data/traits/",data_file), sheet = "LeafMorpho_traits", startRow = 1, colNames = TRUE)  
LeafMorpho_leo <- read.csv2("data/traits/leafMorpho_3.csv")
LeafMorpho <- rbind(LeafMorpho1,LeafMorpho_leo)

# LeafDimensions <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafDimensions (àsupprimer)", startRow = 1, colNames = TRUE)
LeafCN1 <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafC&N", startRow = 1, colNames = TRUE) 
LeafP <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafP", startRow = 1, colNames = TRUE) 
Leaf13C1 <- read.xlsx(paste0("data/traits/",data_file), sheet = "Leaf13C", startRow = 1, colNames = TRUE) 

# measurements C, N, 13C spring
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(species == "Geranium dissectum - limbe")) %>% 
  filter(!(species == "Geranium dissectum - pétiole")) %>% 
  filter(!(species == "Carex humilis?"))

names_LH <- MEAN %>% 
  select(code_sp,species,LifeForm1) %>% 
  unique() %>% 
  rename(Code_Sp = code_sp,Species=species)


Leafchim_leo <- read.csv2(paste0("data/traits/leafchim_leo.csv"))
# Species name and LifeForm
L2 <- Leafchim_leo %>% 
  merge(names_LH,by="Code_Sp") %>% 
  mutate(LifeForm2=NA)
# Day of measurement
L3 <- LeafMorpho_leo %>% 
  select(Code_Sp,Day,Plot,Treatment) %>% 
  unique() %>% 
  merge(L2,by=c("Code_Sp","Plot"))
# Family
fam <- read.csv2("data/traits/names and family.csv")
L4 <- merge(L3,fam,by="Species")
# Rep
L5 <- L4 %>% 
  group_by(Species,Plot) %>% 
  mutate(Rep = paste("RCN",Rep,sep=""))

LeafCN_leo <- L5[,colnames(LeafCN)]
Leaf13C_leo <- L5[,colnames(Leaf13C)]

LeafCN <- rbind(LeafCN1,LeafCN_leo)
Leaf13C <- rbind(Leaf13C1,Leaf13C_leo)

Biovolume <- read.xlsx(paste0("data/traits/",data_file), sheet = "Biovolume", startRow = 1, colNames = TRUE) %>% 
  mutate(Hrepro = as.numeric(Hrepro))
Pheno <- read.xlsx(paste0("data/traits/",data_file), sheet = "Pheno", startRow = 1, colNames = TRUE) %>% 
  mutate(Rep = "None") %>% 
  mutate(Day = "None")
Seed <- read.xlsx(paste0("data/traits/",data_file), sheet = "Seed", startRow = 1, colNames = TRUE) 

# 2) Mean trait value per species * treatment ####
mean_attribute_per_species <- function(dataset){
  dataset2 <- dataset %>% 
    mutate(Trtmt = str_sub(Treatment,1,3) ) %>% 
    filter(Trtmt %in% c("Fer",'Nat','Tem','Che')) %>% 
    mutate(LifeHistory = if_else(LifeForm1=="The","annual","perennial")) %>% 
    mutate(Form = case_when(Family == "Poaceae" ~ "Grass",
                            Family =="Juncaceae" ~ "Rush",
                            Family == "Cyperaceae" ~ "Sedge",
                            TRUE ~ "Other")) %>% 
    # NB : /!\ I sould keep family and lifeform info, but only when these are variables are clean in the dataset.
    # (sinon, ça découple une espèce en deux artificiellement dans le calcul de la moyenne).
    select(-c(nameOfProject,measurementDeterminedBy,Rep))
  
  # Variables on which we want to summarize:
  vars <- dataset2 %>% 
    select(!c(Site,Block,Plot,Treatment,Year,Day,Species,Code_Sp,Family,LifeForm1,LifeForm2,
              Trtmt,LifeHistory,Form)) %>% 
    colnames()
  
  dataset2 %>%     
    group_by(Species,Trtmt,Code_Sp ,    Form, LifeHistory, LifeForm1) %>% 
    summarise_at( .vars = vars , mean,na.rm=T)
}


# Compile mean of all the traits in one dataset
MEAN_list <- list( mean_attribute_per_species(LeafMorpho),
                   mean_attribute_per_species(LeafCN),
                   mean_attribute_per_species(LeafP),
                   mean_attribute_per_species(Leaf13C),
                   mean_attribute_per_species(Biovolume),
                   mean_attribute_per_species(Pheno),
                   mean_attribute_per_species(Seed) )
MEAN <- MEAN_list[[1]]
MEAN$Species <- recode(MEAN$Species,"Cirsium acaulon" = "Cirsium acaule")
for (i in 2:length(MEAN_list)){
  MEAN <- full_join(MEAN,MEAN_list[[i]], by = c('Species','Trtmt','Code_Sp','LifeHistory','LifeForm1','Form'))
  MEAN$Species <- recode(MEAN$Species,"Cirsium acaulon" = "Cirsium acaule")
}

# Clean (uniformize) data
MEAN$Species <- recode(MEAN$Species,"Festuca christiani-bernardii" = "Festuca christianii-bernardii")
MEAN2 <- MEAN %>% 
  dplyr::rename(species = Species, code_sp = Code_Sp, treatment = Trtmt) %>% 
  filter(!(species %in% c("Carex humilis?","Carex sp.","Geranium dissectum - petiole","Geranium dissectum - limbe"))) %>% 
  filter(!(treatment %in% c("Tem","Che")))


write.csv2(MEAN2,"outputs/data/mean_attribute_per_treatment.csv",row.names=F)



