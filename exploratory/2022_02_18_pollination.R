source("scripts/1. Packages.R")

# Importer aussi indices d'Ellenberg

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv")
baseflor <- read.csv2("data/traits/baseflor.csv")


remove_subsp_name <- function(text){
  # text is the scientific name, e.g., Cymodocea nodosa (Ucria) Asch.
  #  the output is Cymodocea nodosa
  spacestwice <- str_locate_all(text," ")
  spaces <- spacestwice[[1]][,1]
  if (length(spaces)>1){
    text <- str_sub(text,start=0,end=spaces[2]-1)
  }else{
  }
  str_replace(text," "," ")
  
}

sp_list <- MEAN %>% 
  select(species,code_sp,LifeForm1,LifeHistory)

info_baseflor <- baseflor %>% 
  mutate(species_1 = map_chr(NOM_SCIENTIFIQUE,remove_subsp_name)) %>% 
  mutate(species_2 = map_chr(nomH,remove_subsp_name)) %>% 
  select(species_1,species_2,everything()) %>% 
  filter( (species_1 %in% sp_list$species) | (species_2 %in% sp_list$species)) %>% 
  select(inflorescence,sexualité,ordre_maturation,pollinisation,fruit,dissémination,floraison,
         CARACTERISATION_ECOLOGIQUE_(HABITAT_OPTIMAL),
         )
