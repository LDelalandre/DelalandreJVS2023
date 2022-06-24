library(tidyverse)
seed <- read.csv2("data/traits/Seed_leo.csv")

seed %>% 
  mutate(indiv_seed_mass = TotSeedMass.mg./SeedNb) %>% 
  filter(Code_Sp=="TRIFCAMP") %>%
  mutate(sp_plot = paste(Code_Sp,Plot,sep="_")) %>% 
  ggplot(aes(x=sp_plot,y=indiv_seed_mass)) +
  geom_boxplot() +
  geom_point() 

# masse des graines de la manip à corréler avec les taux de germination

seed_focus <- seed %>% 
  mutate(indiv_seed_mass = TotSeedMass.mg./SeedNb) %>% 
  filter(Code_Sp=="MEDIMINI")

t.test(indiv_seed_mass ~ Plot, data = seed_erop)


seed %>% 
  mutate(indiv_seed_mass = TotSeedMass.mg./SeedNb) %>% 
  filter(Code_Sp=="CERAGLOM") %>%
  mutate(sp_plot = paste(Code_Sp,Plot,sep="_")) %>% 
  ggplot(aes(x=sp_plot,y=indiv_seed_mass)) +
  geom_point()

seed_aren <- seed %>% 
  mutate(indiv_seed_mass = TotSeedMass.mg./SeedNb) %>% 
  filter(Code_Sp=="ARENSERP")

t.test(indiv_seed_mass ~ Plot, data = seed_aren)

 