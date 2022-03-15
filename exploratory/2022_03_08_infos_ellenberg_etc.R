library(tidyverse)

ellenberg <- read.csv2("outputs/data/ellenberg_annuals.csv") %>% 
  rename(species = AccSpeciesName)
evolution <- read.csv2("outputs/data/temporal_evolution_fer.csv") 
baseflor <- read.csv2("outputs/data/baseflor_annuals.csv") %>% 
  rename(species = species_1) %>% 
  rename(caract_ecol = "CARACTERISATION_ECOLOGIQUE_.HABITAT_OPTIMAL.") %>% 
  select(-species_2) %>% 
  unique()

ellenberg_evolution <- full_join(ellenberg,evolution,by="species")

jointot <- full_join(ellenberg_evolution,baseflor,by="species")


MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(species == "Geranium dissectum - limbe")) %>% 
  filter(!(species == "Geranium dissectum - pétiole")) %>% 
  filter(!(species == "Carex humilis?"))

names_LH <- MEAN %>% 
  select(code_sp,species,LifeHistory) %>% 
  unique()

sp_manip <- c("BUPLBALD",
              "FILAPYRA",
              "ALYSALYS","CAPSBURS","EROPVERN","HORNPETR",
              "ARENSERP","CERAGLOM","CERAPUMI","MINUHYBR",
              "MEDIMINI","TRIFCAMP",
              "GERADISS",
              "BROMHORD","VULPMYUR",
              "SHERARVE",
              "SAXITRID",
              "VEROARVE",
              "MYOSRAMOS")
names_LH_manip <- names_LH %>% 
  filter(code_sp %in% sp_manip)

# 1) Classification des espèces selon indices d'Ellenberg ####
ellenberg %>% 
  full_join(names_LH_manip,by="species") %>% 
  filter(code_sp %in% sp_manip) %>%  
  arrange(nitrogen) %>% 
  select(species,nitrogen) %>% 
  write.csv2("outputs/data/manip_ellenberg_nitrogen.csv",row.names=F)
  
  
  # ggplot(aes(y=nitrogen)) +
  # geom_histogram(stat="identity")


# baseflor_evoltempo ####
baseflor_evolution <- merge(baseflor,evolution %>% filter(species %in% names_LH_manip$species),
                            by="species") %>% 
  filter(!is.na(dissémination))
baseflor_evolution %>%
  ggplot(aes(x=as.factor(dissémination),y=Rs,label = species))+
  geom_point() +
  ggrepel::geom_text_repel()
  
M_dissem <- lm(data = baseflor_evolution, Rs~dissémination)
anova(M_dissem)
summary(M_dissem)

# J'ai les indices d'azote pour 16 sur les 19 espèces de ma manip. Pas mal !
ellenberg %>% 
  filter(species %in% names_LH_manip$species) %>% 
  pull(species) %>% 
  unique()

# Pollinisation et dissémination ####
baseflor_manip <- names_LH %>% 
  full_join(baseflor,by = "species") %>% 
  filter(code_sp %in% sp_manip) %>% 
  select(species,code_sp,pollinisation,dissémination,caract_ecol)

baseflor_manip %>% 
  select(species,dissémination) %>% 
  unique() %>% 
  arrange(dissémination)
# manque saxitrida. Pas grave, je vais possiblement pas la garder... ?

env <- baseflor_manip %>% 
  select(species,caract_ecol) %>% 
  unique() %>% 
  arrange(caract_ecol)

# ellenberg et évol tempo ####

# modèle
mod <- lm(Rs ~ light, data = ellenberg_evolution)
mod <- lm(Rs ~ light, data = jointot %>% filter(!is.na(Rs)))
anova(mod)
summary(mod)

# azote très structurant
ellenberg_evolution %>% 
  filter(!is.na(nitrogen)) %>% 
  ggplot(aes(y=Rs,x=as.factor(nitrogen),label = code_sp))+
  geom_point() +
  # geom_smooth(method="lm")
  ggrepel::geom_text_repel()
  
NRs <- lm(data = ellenberg_evolution, Rs ~ nitrogen )
anova(NRs)
summary(NRs)

ellenberg_evolution %>% 
  filter(!is.na(nitrogen)) %>% 
  ggplot(aes(x=as.factor(nitrogen),y=Rs,label = code_sp))+
  geom_boxplot()

# sécheresse moins structurant. Logique.
ellenberg_evolution %>% 
  ggplot(aes(x=Rs,y=moisture,label = code_sp))+
  geom_point() +
  ggrepel::geom_text_repel()
# En revanche, la sécheresse peut être structurante au niveau du gradient de Maud, à l'échelle communauté.
# Bonne idée pour identifier les facteurs structurants, en plus du CSR.

# Lumière, donc compétition, structurante
ellenberg_evolution %>% 
  ggplot(aes(y=Rs,x=light,label = code_sp))+
  geom_point() 
  ggrepel::geom_text_repel()
# Stylé ! ont besoin de peu de compétition pour s'en sortir dans le natif.
ellenberg_evolution %>% 
  filter(!is.na(light)) %>% 
  ggplot(aes(x=as.factor(light),y=Rs,label = code_sp))+
  geom_boxplot()

# Température: pas assez de données
ellenberg_evolution %>% 
  ggplot(aes(x=Rs,y=temperature,label = code_sp))+
  geom_point() +
  ggrepel::geom_text_repel()

# pH
ellenberg_evolution %>% 
  ggplot(aes(x=Rs,y=pH,label = code_sp))+
  geom_point() +
  ggrepel::geom_text_repel()

# abondances maud ####
# NB faire pareil sur les abondances du fertile... mais c'est tautologique, donc non.
ab_maud_ann <- read.csv2("outputs/data/abundance_natif.csv") %>% 
  filter(depth=="S") %>%
  filter(LifeHistory=="annual") %>% 
  group_by(species) 

ab_ellenberg_evolution <- full_join(ellenberg_evolution,ab_maud_ann,by="species")
ab_ellenberg_evolution %>% 
  # filter(species %in% names_LH_manip$species) %>% 
  ggplot(aes(x=Rs, y = abundance))+
  geom_point()

# azote et abondance: intéressant
ab_ellenberg_evolution %>%
  filter(!is.na(abundance)) %>% 
  filter(!is.na(nitrogen)) %>% 
  ggplot(aes(x= as.factor(nitrogen), y = abundance))+
  geom_boxplot()
# les plus abondantes sont soit les moins nitrophiles, soit les plus nitrophiles.
# Lié à la perturbation par les brebis (crottes très localement...) et on n'a pas l'info.

# lumière et abondance: ne dit rien... Mais il faudrait comparer avec le fertile, ça doit bien marcher!
# = les plus abondantes du fertile vont avoir peu de besoins en lumière.
# Je peux intégrer ellenberg sur les ACP en comparaison entre annuelles (ou le faire en boxplot, parce 
# que ça me donne une valeur par espèce, pas par pop, contrairement aux ACP).
# Les ACP en comparaison entre annuelles auraient un sens avec des valeurs par sp, pour écarter les effets de plasticité!
ab_ellenberg_evolution %>%
  filter(!is.na(abundance)) %>% 
  filter(!is.na(light)) %>% 
  ggplot(aes(x= as.factor(light), y = abundance))+
  geom_boxplot()


# présence et évoluton temporelle ####
sp_maud <- ab_maud_ann %>% 
  filter(abundance > 0) %>% 
  pull(code_sp) %>% 
  unique()

evol_maud <- evolution %>% 
  mutate(in_maud = if_else(code_sp %in% sp_maud,T,F))

evol_maud %>% 
  ggplot(aes(x=in_maud,y=Rs,label=code_sp))+
  geom_point() +
  ggrepel::geom_text_repel()


# abondantes increaser et dec dans mêmes plots ? ####
ab_maud_ann %>%
  mutate(plot = paste(depth,paddock,line,sep="_")) %>% 
  ggplot(aes(x=as.factor(plot), y=abundance, label = code_sp))+
  geom_point() 
  ggrepel::geom_text_repel()

ab_maud_ann %>% filter(species)
