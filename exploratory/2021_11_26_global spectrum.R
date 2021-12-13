library(tidyverse)
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  filter(!(Species == "Carex humilis?")) %>% 
  filter(!(Species == "Cirsium acaule")) # il faudra le réintégrer

data_shiny <- MEAN %>% 
  select(Code_Sp,Trtmt,LifeHistory,Hrepro,L_Area,SLA,LNC,SeedMass) %>% 
  filter(!(is.na(SLA | Hrepro | SeedMass))) %>% 
  transmute(Code_Sp = Code_Sp,
            Trtmt = Trtmt,
            LifeHistory = LifeHistory,
            H=Hrepro * 0.01, # from cm to m
            LA = 100*L_Area, # from cm² to mm²
            LMA = 1/SLA*1000, # from m²/kg to g/m²
            Nmass = LNC,
            SM = SeedMass)

write.csv2(data_shiny %>% 
             filter(LifeHistory == "annual") ,"outputs/data/Global spectrum/data_shiny_annuals.csv",row.names = F) 

# shiny app always disconnect from server...

write.csv2(data_shiny %>% 
             filter(Trtmt=="Fer") %>% 
             select(-Trtmt),"outputs/data/Global spectrum/data_shiny_sp_fer.csv",row.names = F)

write.csv2(data_shiny %>% 
             filter(Trtmt=="Nat") %>% 
             select(-Trtmt),"outputs/data/Global spectrum/data_shiny_sp_nat.csv",row.names = F)
