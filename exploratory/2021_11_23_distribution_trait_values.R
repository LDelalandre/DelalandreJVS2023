source("scripts/1. Packages.R")
source("scripts/2. Import files.R")

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  filter(!(Species == "Carex humilis?")) %>% 
  filter(!(Species == "Cirsium acaule")) # il faudra le réintégrer

MEAN_treatment <- MEAN %>% filter(Trtmt == "Nat")

trait <- "LPC"
ggplot(MEAN_treatment,aes_string(x=trait))+
  geom_density() +
  geom_point(data = MEAN %>% filter(LifeHistory=="annual"), aes_string(x=trait,y=0))
