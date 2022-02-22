library(tidyverse)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(species == "Geranium dissectum - limbe")) %>% 
  filter(!(species == "Geranium dissectum - pétiole")) %>% 
  filter(!(species == "Carex humilis?"))

names_LH <- MEAN %>% 
  select(code_sp,species,LifeHistory) %>% 
  unique()
annuals <- names_LH %>% filter(LifeHistory=="annual") %>% 
  pull(species)


# Eric ####
tempo_var_eric <- read.table("data/abundance/Decreasers increasers Garnier et al 2018/Supplementary Garnier 2018.txt",sep="\t",header=T)
andecinc <- tempo_var_eric %>% 
  filter(Life.cycle == "Annual" )


# Adeline ####
tempo_var_adeline <- read.csv2("data/abundance/Decreasers increasers Adeline/temporal_variation_adeline.CSV")



tempo_var_adeline_annuals <- tempo_var_adeline %>% 
  filter(species %in% annuals)

# Diachro ####
diachro <- read.csv2("data/abundance/diachro_releves_tidy2.csv")
# 33 transects, 16 dans le fertile, 17 dans le natif.

diachro2 <- diachro %>% 
  filter(fertilized==T) %>% 
  group_by(line, year) %>% 
  mutate(rel_ab = abundance/sum(abundance)) %>% 
  group_by(species,year) %>% 
  mutate(delta_year = year - 1978, mean_ab = mean(abundance)) %>% 
  filter(LifeForm1=="The")

# Plus ou moins d'annuelles dans le natif au cours du temps?
diachro %>% 
  filter(fertilized == FALSE) %>% 
  select(year,line) %>% 
  unique() %>% 
  arrange(year,line)

diachro_nat<- diachro %>% 
  filter(fertilized==F) %>% 
  group_by(line, year) %>% 
  mutate(rel_ab = abundance/sum(abundance)) %>% 
  group_by(species,year) %>% 
  mutate(delta_year = year - 1978, mean_ab = mean(abundance)) %>% 
  filter(LifeForm1=="The") %>% 
  ungroup()

diachro_nat_richness <- diachro_nat %>% 
  group_by(year) %>% 
  filter(abundance >0) %>% 
  select(species,year,fertilized) %>% 
  unique() %>% 
  summarize(richness = n())

ggplot(diachro_nat_richness,aes(x=year,y=richness))+
  geom_point()

ggplot(diachro_nat %>% filter(abundance>0),aes(x=year,y=abundance))+
  geom_point() +
  facet_wrap(~species)


ggplot(diachro2 %>% filter(species == "Arenaria serpyllifolia"),aes(x=year,y=mean_ab))+
  geom_point()

diachro3 <- diachro2 %>% 
  ungroup() %>% 
  group_by(species) %>%
  filter(LifeForm1=="The") %>% 
  summarize(Rs = cor(delta_year,mean_ab,method = "spearman"),
            pval = cor.test(delta_year, mean_ab, method = "spearman")$p.value)


ggplot(diachro2 %>% filter(species == "Bupleurum baldense"),aes(x=year,y=abundance))+
  geom_point() 

ggplot(diachro2 %>% filter(species == "Arenaria serpyllifolia"),aes(x=year,y=abundance))+
  geom_point() +
  # facet_wrap(~line) + 
  geom_smooth(method = "lm")


sp_manip <- c("BUPLBALD",
          "ALYSALYS","CAPSBURS","EROPVERN","HORNPETR",
          "ARENSERP","CERAGLOM","CERAPUMI","MINUHYBR",
          "MEDIMINI","TRIFCAMP",
          "GERADISS",
          "BROMHORD","VULPMYUR",
          "SHERARVE",
          "SAXITRIDA",
          "VEROARVE",
          "MYOSRAMOS",
          "FILAPYRAM")

diachro4 <- diachro3 %>% 
  merge(names_LH,by="species") %>% 
  filter(code_sp %in% sp_manip)


# Lien à abondance ####
ab_nat <- read.csv2("outputs/data/abundance_natif.csv")
ab_nat_ann <- ab_nat %>% filter(LifeHistory=="annual")
ab_nat_ann %>% 
  group_by(species,code_sp) %>% 
  summarize(mean_ab = mean(abundance)) %>% 
  merge(diachro3,by="species") %>% 
  ggplot(aes(x=Rs,y=mean_ab,label=code_sp))+
  geom_point() +
  geom_text()
  # geom_smooth(method="lm")

  
  

# Natif ####
diachro_tot <- openxlsx::read.xlsx("data/abundance/DATA_DIACHRO_AB.xlsx", sheet = "AB_ALL", startRow = 1, colNames = TRUE)  

annuals <- names_LH %>% 
  filter(LifeHistory=="annual") %>% 
  pull(code_sp) %>% 
  unique()

diachro_nat <- diachro_tot %>% 
  filter(GESTION == "NATIF50") %>% 
  group_by(ANNEE,LIGNE_STANDARD) %>% 
  mutate(AB_relat = AB/sum(AB)) %>% 
  filter(!(ANNEE %in% c(1972,1976)))


# evolution de la richesse des annuelles
diachro_nat %>% 
  filter(CODE_ESP %in% annuals) %>% 
  select(CODE_ESP,ANNEE) %>% 
  group_by(ANNEE) %>% 
  unique() %>% 
  summarize(richness = n())
  
# Nombre de transects au cours des années
nb_transects <- diachro_nat %>% 
  select(ANNEE,LIGNE_STANDARD) %>% 
  unique() %>% 
  group_by(ANNEE) %>% 
  summarize(NB_TRANSECTS=n())
# prendre les 63 transects de 1995 sur chaque année. 
# Regarder l'évolution par rapport à la position spatiale dans ces transects ?
  
# evolution de la richesse des annuelles corrigée pour nombre de relevés
# NB c'est pas propre, et le mieux c'est de regarder en proportion du nombre de points contacts, comme a fait Adeline
diachro_nat %>% 
  merge(nb_transects, by = "ANNEE") %>% 
  select(CODE_ESP,ANNEE,NB_TRANSECTS) %>% 
  group_by(ANNEE) %>% 
  unique() %>% 
  geom_point()

# Evolution de l'abondance de chaque espèce
diachro_nat %>% 
  filter(CODE_ESP %in% annuals) %>% 
  filter(!(ANNEE %in% c(1972,1976))) %>% 
  merge(nb_transects, by = "ANNEE") %>%
  mutate(AB_corr = AB/NB_TRANSECTS) %>% 
  ggplot(aes(x=ANNEE,y=AB))+
  geom_point() +
  facet_wrap(~CODE_ESP) +
  geom_smooth(method="lm")
  