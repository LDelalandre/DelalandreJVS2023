library(tidyverse)

soil <- read.csv2("data/environment/Soil_maud_adeline.csv")

soil %>% 
  filter(treatment == "Nat") %>% 
  ggplot(aes(x=measured_by,y=Humidité))+
  geom_point()

soil %>% 
  filter(treatment == "Nat") %>% 
  ggplot(aes(x=measured_by,y=C.N))+
  geom_point()

# Pas de concordance franche entre les mesures d'adeline et de Maud dans le 
# natif --> je ne peux pas combiner ces données.
