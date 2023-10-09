library(tidyverse)

## Nitrogen and Phosporus nutrition indices ####
INraw <- read.xlsx("data/environment/IN La Fage.xlsx", sheet = "resultat", startRow = 1, colNames = TRUE) %>% 
  mutate()
nom <- read.xlsx("data/environment/IN La Fage.xlsx", sheet = "nom", startRow = 1, colNames = TRUE)  %>% 
  filter(!is.na(plot))
IN <- merge(nom,INraw,by="nature.echantillon") %>% 
  filter(!(treatment =="temoin")) %>% 
  mutate(soil = if_else(treatment=="fertilise","Fer",plot)) %>% 
  select(soil,INN,INP) %>% 
  filter(soil %in% c("Fer","Sable"))

## Biomass produced ####
biomass <- read.xlsx("data/environment/Biomasses et indices La Fage.xlsx", 
                     sheet = "2009", 
                     startRow = 1, colNames = TRUE, rowNames = F) %>% 
  mutate(Dates = as.Date(Dates- 25569, origin = "1970-01-01"))

biomass_may <- biomass %>% 
  filter(Parcs %in% c("C1","C2","1","6","8","10")) %>% 
  filter(Dates == "2009-05-01") %>% 
  mutate(rdt.T.ha = as.numeric(rdt.T.ha)) %>% 
  mutate(Position = str_sub(Position,1,1)) %>% 
  mutate(Position = case_when(is.na(Position) ~ "Fer",
                              TRUE ~ Position)) %>% 
  select(-INN)
biomass_may$Position <- factor(biomass_may$Position , levels = c("Fer","D","I","S"))
biomass_may_reduced <- biomass_may %>% 
  filter(Position %in% c("Fer","S"))
biomass_may_reduced$Position <- factor(biomass_may_reduced$Position , levels = c("Fer","S"))

## Vegetation consumption ####
disturbance <- read.table("data/environment/Disturbance_DivHerbe.txt",header=T,sep="\t",dec=",") %>% 
  select(-c(CODE_SITE,TRAITEMENT,ReTimDist,DistInt))

# Group the environmental data in one
# data utilisable pr figures
IN_gathered <- IN %>% 
  rename(management = soil) %>% 
  mutate(management = if_else(management=="Fer","Int","Ext")) %>% 
  na.omit() %>% 
  gather(key = variable, value = value, -management) %>% 
  na.omit() %>% 
  mutate(unit = "%")

biomass_gathered <- biomass_may_reduced %>% 
  rename(management = Position) %>% 
  mutate(management = if_else(management=="Fer","Int","Ext")) %>% 
  select(-c(Dates,Parcs,Repet)) %>% 
  mutate(variable = "Biomass", unit = "t/ha") %>% 
  rename(value = rdt.T.ha)


disturbance_gathered <- disturbance %>% 
  rename(management = Trtmt) %>% 
  filter(!(management == "Tem")) %>% 
  mutate(management = if_else(management=="Fer","Int","Ext")) %>% 
  select(-c(PLOT,REPETITION)) %>% 
  mutate(variable = "Disturbance", unit = "%") %>% 
  rename(value = Tx_CalcPic)

env_data <- rbind(IN_gathered,biomass_gathered,disturbance_gathered) %>% 
  mutate(variable = case_when(variable == "INN" ~ "NNI",
                              variable == "INP" ~ "PNI",
                              TRUE ~variable))
env_data$management = factor(env_data$management, levels = c("Int", "Ext"))

write.csv(env_data,"outputs/data/env_data.csv",row.names=F)
