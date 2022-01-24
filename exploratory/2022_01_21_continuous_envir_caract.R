source("scripts/1. Packages.R")

biomass <- read.xlsx("data/environment/Biomasses et indices La Fage.xlsx", 
          sheet = "2009", 
          startRow = 1, colNames = TRUE, rowNames = F)%>% 
  mutate(Dates = as.Date(Dates- 25569, origin = "1970-01-01"))

biomass_may <- biomass %>% 
  filter(Parcs %in% c("C1","C2","1","6","8","10")) %>% 
  filter(Dates == "2009-05-01") %>% 
  mutate(rdt.T.ha = as.numeric(rdt.T.ha))

ggplot(biomass_may,aes_string(x="Position",y="rdt.T.ha"))+
  geom_boxplot()

ggplot(biomass_may,aes_string(x="Position",y="INN"))+
  geom_boxplot()
