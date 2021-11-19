source("scripts/1. Packages.R")
source("scripts/2. Import files.R")
library("vegan")

name_lifeform <- ABUNDANCE %>% 
  select(species,code_sp,LifeForm1) %>% 
  unique() %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The", "annual", "perennial"))

IN_fer <- read.csv2("data/abundance/IN_fertile.csv") %>% 
  dplyr::rename(cage = plot) %>% 
  select(parc,cage,INN,INP,group) %>% 
  # filter(Traitement == "P+F+") %>% 
  mutate(plot = paste(parc,cage,sep=""))

ab_iris <- read.xlsx("data/abundance/iris_Relevés bota.xlsx",sheet="Fertile") %>% 
  mutate(cage = str_sub(Cage,1,1)) %>% 
  dplyr::rename(abundance = 'Présence.(tot)') %>% 
  mutate(plot = paste(Parc,cage,sep=""))

ab_iris_IN <- full_join(ab_iris,IN_fer,by="plot") %>% 
  mutate(species=Espece) %>% 
  select(-Espece) %>% 
  full_join(name_lifeform,by="species") %>% 
  filter(!is.na(abundance))




ggplot(ab_iris_IN,aes(x=group,y=abundance))+
  geom_boxplot() +
  facet_wrap(~LifeHistory)

# CWM
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  rename(species = Species, code_sp = Code_Sp, treatment = Trtmt)

CWM_iris <- MEAN %>% 
  filter(treatment=="Fer") %>% 
  select(-LifeHistory,code_sp) %>% 
  full_join(ab_iris_IN,by="species") %>% 
  group_by(group) %>% # NB choose the level at which to compute moments
  mutate_at(vars(Nb_Lf:Mat),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) )


CWM <- CWM_iris %>% 
  filter(!is.na(group)) %>% 
  select(starts_with("CWM")) %>% 
  unique() %>% 
  arrange(group)
# Le CWM de LNC ne suit pas les INN. C'est donc que ce n'est pas du remplacement 
# d'espèces qu'on a ici.

comm <- ab_iris_IN %>% 
  filter(!(Cage=="A2")) %>% 
  select(species,abundance,plot) %>% 
  unique() %>% 
  spread(plot,abundance) %>% 
  column_to_rownames("species") %>% 
  replace(is.na(.),0) %>% 
  relocate(C2C)

distance <- vegdist(t(comm),na.rm=T,method = "bray")
# REF : C2C: le plus faible INN
# compute sorensen similarity index

reshape::melt(as.matrix(distance)) %>% 
  dplyr::rename(distance=value)

distance_IN <- distance %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  select(C2C) %>% 
  rownames_to_column("plot") %>% 
  merge(IN_fer,by="plot") %>% 
  arrange(INN)

order_plot <- distance_IN %>% 
  pull(plot)

distance_IN <- distance_IN %>% 
  arrange(factor(plot,levels=order_plot) ) %>% 
  mutate(plot=factor(plot,levels=plot))

ggplot(distance_IN,aes(x=INN,y=C2C))+
  geom_point() +
  xlab("Indice de nutrition azotée (INN)")+
  ylab("Distance de Jaccard au plot de plus faible INN")
