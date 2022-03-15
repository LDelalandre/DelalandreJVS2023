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
  merge(names_LH,by="species") 

diachro5 <- diachro4 %>% 
  filter(code_sp %in% sp_manip)

write.csv2(diachro4,"outputs/data/temporal_evolution_fer.csv",row.names = F)


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

lines_kept <- diachro_tot %>% 
  filter(GESTION == "NATIF50") %>% 
  filter(ANNEE == 1995) %>% 
  pull(LIGNE_STANDARD) %>% 
  unique()

diachro_nat <- diachro_tot %>% 
  filter(GESTION == "NATIF50") %>% 
  group_by(ANNEE,LIGNE_STANDARD) %>% 
  filter(LIGNE_STANDARD %in% lines_kept) %>% 
  mutate(AB_relat = AB/sum(AB)) %>% 
  filter(!(ANNEE %in% c(1972,1976))) %>% 
  ungroup()

# ajouter les abondances nulles
diachro_nat_zeros <- diachro_nat %>% 
  expand(ANNEE,LIGNE_STANDARD,CODE_ESP) %>% 
  full_join(diachro_nat,by=c("ANNEE","LIGNE_STANDARD","CODE_ESP")) %>% 
  select(-GESTION) %>% 
  replace(is.na(.), 0)


# evolution de la richesse des annuelles
diachro_nat %>% 
  filter(CODE_ESP %in% annuals) %>% 
  select(CODE_ESP,ANNEE) %>% 
  group_by(ANNEE) %>% 
  unique() %>% 
  summarize(richness = n())

# Evolution de l'abondance de chaque espèce
diachro_nat %>% 
  filter(CODE_ESP %in% annuals) %>% 
  ggplot(aes(x=ANNEE,y=AB_relat))+
  geom_point()  +
  facet_wrap(~CODE_ESP) +
  geom_smooth(method="lm")

diachro_nat_zeros %>% 
  filter(CODE_ESP %in% annuals) %>% 
  group_by(ANNEE,CODE_ESP) %>% 
  summarize(mean_AB = mean(AB)) %>% 
  ggplot(aes(x=ANNEE,y=mean_AB))+
  geom_point()  +
  facet_wrap(~CODE_ESP) +
  geom_smooth(method="lm")


# correlation abondance temps
tete <- diachro_nat_zeros %>% 
  mutate(delta_year = ANNEE - 1978) %>% 
  group_by(CODE_ESP,ANNEE,delta_year) %>% 
  summarize(mean_ab = mean(AB)) %>% 
  summarize(Rs = cor(delta_year,mean_ab,method = "spearman"),
            pval = cor.test(delta_year, mean_ab, method = "spearman")$p.value)

# 
# A faire
# 1) Mesurer increasers et dec du natif
# 2) Reproduire figure Adeline p. 104
  