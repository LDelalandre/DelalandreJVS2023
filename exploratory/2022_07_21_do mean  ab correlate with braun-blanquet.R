library(tidyverse)


# 4) Abundance ####
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(line = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") 

soil_Maud <- data.frame(PC1score = c(-3.08,-2.85,-2.52,-1.78,-1.60,-1.56,-0.03,0.16,1.97,2.66,4.05,4.58),
                        depth = c("S","S","S","S","I","I","I","I","D","D","D","D" ),
                        paddock = c("P8","P10","P6","P1","P6","P8","P10","P1","P10","P1","P6","P8"))

FER <- ab_fer %>%
  filter(LifeForm1=="The") %>% 
  group_by(code_sp) %>% 
  summarize(mean_ab_fer = mean(relat_ab))

NAT <- ab_nat %>% 
  filter(LifeForm1=="The") %>% 
  group_by(code_sp) %>% 
  summarize(mean_ab_nat = mean(relat_ab))

comp_ab <- full_join(FER,NAT,by="code_sp")
