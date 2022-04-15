source("scripts/1. Packages.R")

# List of species
names_LH <- read.csv2("data/species_names_lifehistory.csv")
names_family <- read.csv2("data/traits/names and family.csv")

# 1) Import plant trait values ####
data_file <- "LaFage_PlantTraitsDP_vp.xlsx"

LeafMorpho1 <-  read.xlsx(paste0("data/traits/",data_file), sheet = "LeafMorpho_traits", startRow = 1, colNames = TRUE)  
LeafMorpho_leo <- read.csv2("data/traits/leafMorpho_3.csv") %>% 
  filter(!(Code_Sp %in% c("EROPVERN","STELMEDI"))) # senescent leaves
LeafMorpho <- rbind(LeafMorpho1,LeafMorpho_leo)

# LeafDimensions <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafDimensions (àsupprimer)", startRow = 1, colNames = TRUE)
LeafCN1 <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafC&N", startRow = 1, colNames = TRUE) 
LeafP <- read.xlsx(paste0("data/traits/",data_file), sheet = "LeafP", startRow = 1, colNames = TRUE) 
Leaf13C1 <- read.xlsx(paste0("data/traits/",data_file), sheet = "Leaf13C", startRow = 1, colNames = TRUE) 

# measurements C, N, 13C spring 2021
Leafchim_leo <- read.csv2(paste0("data/traits/leafchim_leo.csv")) %>% 
  filter(!(Code_Sp %in% c("EROPVERN","STELMEDI"))) # senescent leaves
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

LeafCN_leo <- L5[,colnames(LeafCN1)] %>% 
  mutate(LCC = 10*LCC) %>% 
  mutate(LNC = 10*LNC) # to change the unit from % to mg/g.
Leaf13C_leo <- L5[,colnames(Leaf13C1)]

LeafCN <- rbind(LeafCN1,LeafCN_leo)
Leaf13C <- rbind(Leaf13C1,Leaf13C_leo)

Biovolume <- read.xlsx(paste0("data/traits/",data_file), sheet = "Biovolume", startRow = 1, colNames = TRUE) %>% 
  mutate(Hrepro = as.numeric(Hrepro))

# Pheno ####
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
  left_join(names_family,by="Species") %>% 
  separate(col = Day2,into = c("Year","Month","Day")) %>% 
  mutate(Year = as.numeric(Year), Month = as.numeric(Month), Day = as.numeric(Day)) %>% 
  mutate(Disp = map2_dbl(Month,Day,get_day_of_year)) %>% 
  select(-c(code_sp,date,plot)) %>% 
  mutate(LifeForm2 = NA) %>%  # A compléter pour le data paper
  mutate(Flo = NA, Mat_Per = NA, nameOfProject = "Annuals",measurementDeterminedBy = "Léo Delalandre") %>% 
  mutate(Rep = "None") %>% 
  group_by(Species) %>% 
  filter(Disp == min (Disp)) %>%  # NB: une valeur de phéno par espèce*traitment dans la BDD ; je prends la min !
  select(all_of(colnames))

Pheno <- rbind(Pheno1,pheno_leo_DP)

# Seed
Seed <- read.xlsx(paste0("data/traits/",data_file), sheet = "Seed", startRow = 1, colNames = TRUE) 

# 2) Mean trait value per species * treatment ####
mean_attribute_per_species <- function(dataset,subset_gt_nat = F){
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
  if (subset_gt_nat == T){
    dataset2 <- dataset2 %>% 
      filter(Treatment %in% c("Fer_Clc","Fer_Dlm","Nat_Sab"))
  }
  
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
take_nat_sab_only <- T

if(take_nat_sab_only == F){
  MEAN_list <- list( mean_attribute_per_species(LeafMorpho),
                     mean_attribute_per_species(LeafCN),
                     mean_attribute_per_species(LeafP),
                     mean_attribute_per_species(Leaf13C),
                     mean_attribute_per_species(Biovolume),
                     mean_attribute_per_species(Pheno),
                     mean_attribute_per_species(Seed) )
} else{
  MEAN_list <- list( mean_attribute_per_species(LeafMorpho, subset_gt_nat = T),
                     mean_attribute_per_species(LeafCN, subset_gt_nat = T),
                     mean_attribute_per_species(LeafP, subset_gt_nat = T),
                     mean_attribute_per_species(Leaf13C, subset_gt_nat = T),
                     mean_attribute_per_species(Biovolume, subset_gt_nat = T),
                     mean_attribute_per_species(Pheno, subset_gt_nat = T),
                     mean_attribute_per_species(Seed, subset_gt_nat = T) )
}

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



if(take_nat_sab_only == F){
  write.csv2(MEAN2,"outputs/data/mean_attribute_per_treatment.csv",row.names=F)
} else{
  write.csv2(MEAN2,"outputs/data/mean_attribute_per_treatment_subset_nat_sab.csv",row.names=F)
}

