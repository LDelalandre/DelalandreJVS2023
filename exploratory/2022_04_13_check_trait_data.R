library(tidyverse)

# /!\ Importer données de scripts/2. Functional_traits.R /!\

who_measures <- LeafCN %>% 
  filter(grepl("Nat",Treatment)) %>% 
  arrange(Code_Sp) %>% 
  filter(LifeForm1 == "The")
 

who_measures_sp <- who_measures %>% 
  pull(Code_Sp) %>% 
  unique()

who_measures %>%
  # filter(Code_Sp %in% who_measures_sp[1:20]) %>%
  # filter(Code_Sp %in% who_measures_sp[-(1:20)]) %>%
  # filter(Code_Sp %in% who_measures_sp[21:40]) %>%
  # filter(Code_Sp %in% who_measures_sp[41:60]) %>%
  # filter(Code_Sp %in% who_measures_sp[61:80]) %>%
  # filter(Code_Sp %in% who_measures_sp[-(1:80)]) %>%
  ggplot(aes(x=Code_Sp,y= LNC ,color = measurementDeterminedBy)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))



# Focus sur certaines espèces et certains traits ####
# alyssum alyssoides SLA LDMC
LeafMorpho %>% 
  filter(Code_Sp == "ALYSALYS") %>% 
  ggplot(aes(x=Treatment,y=SLA)) +
  geom_boxplot()+
  facet_wrap(~measurementDeterminedBy)

LeafMorpho %>% 
  filter(Code_Sp == "TRIFINCA") %>% 
  ggplot(aes(x=Treatment,y=LDMC)) +
  geom_boxplot()+
  facet_wrap(~measurementDeterminedBy)

# Biovolume ####
Biovolume %>% 
  filter(LifeForm1=="The") %>% 
  filter(Treatment == "Fer_Clc") %>%  # traitement commun avec Adeline
  ggplot(aes(x=Code_Sp,y= Hveg ,color = measurementDeterminedBy)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))
# Je mesure Hveg plus petit chez GERADISS et CERAGLOM

# Hrepro
Biovolume %>% 
  filter(LifeForm1=="The") %>%
  filter(Treatment == "Fer_Clc") %>%  # traitement commun avec Adeline
  ggplot(aes(x=Code_Sp,y= Hrepro ,color = measurementDeterminedBy)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

Biovolume %>% 
  filter(LifeForm1=="The") %>%
  filter(grepl("Nat",Treatment)) %>%  # Je compare dans le natif, pas nécessairement Nat_Sab (où les autres n'ont pas mesuré)
  ggplot(aes(x=Code_Sp,y= Hrepro ,color = measurementDeterminedBy)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

Biovolume_leo %>% 
  ggplot(aes(x = Dmax))+
  geom_histogram(stst="identity")
