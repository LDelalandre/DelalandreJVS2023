library("tidyverse")
source("scripts/Data_traits.R")

Seed %>% 
  filter(LifeForm1=="The") %>% 
  filter(SeedMass < 1) %>% 
  ggplot(aes(x=measurementDeterminedBy,y=SeedMass)) +
  geom_point() +
  facet_wrap(~Code_Sp) +
  geom_jitter()

Seed %>% 
  filter(measurementDeterminedBy=="Leo Delalandre") %>% 
  ggplot(aes(x=Treatment, y=SeedMass)) +
  geom_boxplot()

fdata %>% 
  filter(measurementDeterminedBy=="Leo Delalandre") %>%
  ggplot(aes(x=Treatment, y=SeedMass,label = Code_Sp)) +
  geom_jitter()+
  facet_wrap(~Code_Sp)

fdata %>% 
  filter(measurementDeterminedBy=="Leo Delalandre") %>%
  filter(Code_Sp == "TRIFCAMP") %>% 
  ggplot(aes(x=Treatment, y=SeedMass)) +
  # geom_boxplot() +
  geom_point()
  
fdata %>% 
  filter(measurementDeterminedBy=="Leo Delalandre") %>% 
  pull(Code_Sp) %>% 
  unique()
