source("scripts/1. Packages.R")
source("scripts/2. Import files.R")


#_______________________________________________________________________________
# 1) Mean trait value per species * treatment ####

mean_attribute_per_species <- function(dataset){
  dataset2 <- dataset %>% 
    mutate(Trtmt = str_sub(Treatment,1,3) ) %>% 
    filter(Trtmt %in% c("Fer",'Nat')) %>% 
    group_by(Species,Trtmt,Code_Sp ,    Family, LifeForm1, LifeForm2) %>% 
    select(-c(nameOfProject,measurementDeterminedBy,Rep))
  col <- dim(dataset2)[2] # number of columns
  # last column is Trtmt. Take all columns from the 12th (first trait measured) to the col-1 th (last trait measured).
  dataset2 %>% 
    summarise_at( .vars = colnames(.)[12:(col-1)] , mean,na.rm=T)
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
for (i in 2:length(MEAN_list)){
  MEAN <- full_join(MEAN,MEAN_list[[i]], by = c('Species','Trtmt','Code_Sp','Family','LifeForm1','LifeForm2'))
}

# Clean (uniformize) data
MEAN$Species <- recode(MEAN$Species,"Festuca christiani-bernardii" = "Festuca christianii-bernardii")

write.csv2(MEAN,"outputs/data/mean_attribute_per_treatment.csv",row.names=F)
#_______________________________________________________________________________
# 2) Add mean species attribute per treatment to abundance data per transect/quadrat ####

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")

# Trait values sampled in the same treatment
get_trait_value <- function(sp,trtmt,trait){
  if (trtmt != "Tem"){
    trait_value <- MEAN %>% 
      filter(Species==sp & Trtmt==trtmt) %>% 
      pull(trait) %>% 
      first()
    
    if (length(trait_value==1)){
      trait_value
    } else{
      NA
    }
  } else {
    NA
  }
}

# Traits available
traits <- c("LDMC","SLA","LCC","LNC","Ldelta13C","LPC","SeedMass","Mat","Flo","Disp","Mat_Per","L_Area" ,  
            "Hveg"  ,    "Hrepro"   , "Dmax"  ,    "Dmin" )


# Get the trait values of each species in each treatment, for each trait
get_trait <- function(AB,trait){
  # AB is the ABUNDANCE dataset
  varname <- paste0(trait)
  AB <- mutate(AB, !!varname:= map2_dbl(species,treatment,get_trait_value,trait) )
  AB
}

ABUNDANCE_traits <- ABUNDANCE
for(trait in traits){
  ABUNDANCE_traits <- get_trait(ABUNDANCE_traits,trait)
}
ABUNDANCE_traits <- ABUNDANCE_traits %>% 
  mutate(annual = case_when(LifeForm1 == "The" ~ T,
                            TRUE ~ F)) # add info: annual or not

write.csv2(ABUNDANCE_traits,"outputs/data/pooled_abundance_and_traits.csv",row.names=F)




#_______________________________________________________________________________
# CWM ####
ABUNDANCE_traits <- read.csv2("outputs/data/pooled_abundance_and_traits.csv")


ABUNDANCE_CWM_treatment <- ABUNDANCE_traits %>% 
  group_by(treatment) %>% # NB choose the level at which to compute moments
  mutate_at(vars(LDMC:Mat_Per),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) )
  
# Keep one row per species*treatment, with its trait value and the CWM
attribute_CWM_treatment <- ABUNDANCE_CWM_treatment %>% 
  distinct(species, treatment,.keep_all=T) %>% 
  select(-c("dataset" ,   "method" ,             "year"    ,            "paddock"  ,           "id_transect_quadrat",
            "abundance"    ,       "line_length"    , "soil"       ,         "depth"      ,         "point")) # remove variables that


compute_moment <- function(){
  
}

# Maud
deep_maud <- ABUNDANCE_traits %>%
  filter(dataset=="Maud") %>% 
  filter(depth=="D") %>% 
  mutate_at(vars(LDMC:Mat_Per),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) )
deep_maud %>% 
  select(CWM_LDMC,CWM_SLA)%>% 
  unique()

intermediary_maud <- ABUNDANCE_traits %>%
  filter(dataset=="Maud") %>% 
  filter(depth=="I") %>% 
  mutate_at(vars(LDMC:Mat_Per),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) )
intermediary_maud %>% 
  select(CWM_LDMC,CWM_SLA)%>% 
  unique()

superficial_maud <- ABUNDANCE_traits %>%
  filter(dataset=="Maud") %>% 
  filter(depth=="S") %>% 
  mutate_at(vars(LDMC:Mat_Per),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) )
superficial_maud %>% 
  select(CWM_LDMC,CWM_SLA)%>% 
  unique()


ggplot(deep_maud %>% mutate(LifeHistory = if_else(LifeForm1=="The","Annual","Perennial")),
       aes(x=LifeHistory, y=LDMC))+
  geom_boxplot()
ggplot(intermediary_maud %>% mutate(LifeHistory = if_else(LifeForm1=="The","Annual","Perennial")),
       aes(x=LifeHistory, y=LDMC))+
  geom_boxplot()
ggplot(superficial_maud %>% mutate(LifeHistory = if_else(LifeForm1=="The","Annual","Perennial")),
       aes(x=LifeHistory, y=LDMC))+
  geom_boxplot()




#_______________________________________________________________________________
# Abundance ####
# Spread abundance data
ab_maud <- ABUNDANCE_traits %>% 
  filter(dataset=="Maud") %>% 
  select(paddock,id_transect_quadrat,depth,code_sp,abundance) %>% 
  spread(key = code_sp,value = abundance, fill = 0)
# NB : il me faut surtout ses données par paddock * depth, par par ligne !!
# Je peux le faire dans mes gros data frame pooled.

weighted.mean(ABUNDANCE_traits %>% filter() ,
               ab_maud[1,c(4:dim(ab_maud)[2])])
