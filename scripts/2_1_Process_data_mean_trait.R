source("scripts/Packages.R")
source("scripts/Data_traits.R") # Load traits per group of traits (e.g. LeafMorpho)

# 2) Mean trait value per species * treatment ####

# Compute mean trait value per species*treatment.
# treatment is either Fer against Nat, 
# or G+F (= Fer) against GU-S (= Nat_Sab).

mean_attribute_per_species <- function(dataset,subset_gt_nat = F){
  # Computes mean trait values for a given dataset corresponding
  # to a category of traits (e.g. LeafMorpho)
  # Decide whether to measure mean trait values in the Natif, or Nat_Sab, treatment
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
  
  e <- dataset2 %>%     
    group_by(Species,Trtmt,Code_Sp ,    Form, LifeHistory, LifeForm1) %>% 
    summarise_at( .vars = vars , mean,na.rm=T)
}


# The code below uses the function 'mean_attribute_per_species' to compute mean attributes
# across all the traits and merges it into one data frame

take_nat_sab_only <-  T

# Generate a list of mean attributes per sp*treatment
MEAN_list <- list( mean_attribute_per_species(LeafMorpho, subset_gt_nat = take_nat_sab_only),
                   mean_attribute_per_species(LeafCN, subset_gt_nat = take_nat_sab_only),
                   mean_attribute_per_species(LeafP, subset_gt_nat = take_nat_sab_only),
                   mean_attribute_per_species(Leaf13C, subset_gt_nat = take_nat_sab_only),
                   mean_attribute_per_species(Biovolume , subset_gt_nat = take_nat_sab_only), 
                   mean_attribute_per_species(Pheno, subset_gt_nat = take_nat_sab_only),
                   mean_attribute_per_species(Seed, subset_gt_nat = take_nat_sab_only) )



# Convert the list into a data frame
MEAN <- MEAN_list[[1]]
MEAN$Species <- recode(MEAN$Species,"Cirsium acaulon" = "Cirsium acaule")
for (i in 2:length(MEAN_list)){
  MEAN <- full_join(MEAN,MEAN_list[[i]], by = c('Species','Trtmt','Code_Sp','LifeHistory','LifeForm1','Form'))
  MEAN$Species <- recode(MEAN$Species,"Cirsium acaulon" = "Cirsium acaule")
}

# Clean (uniformize) the data
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
