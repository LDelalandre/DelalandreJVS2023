library(tidyverse)
library(ggpubr)
library(rstatix)


#______________________________________________________________________________
# Comparaison annuelles spring - norme de réaction ####

# /!\ Importer donnees de scripts/2. Functional_traits.R 

# Paired tests to compare trait values of annuals 
# (trait by trait comparison on the traits I measured in spring).

# Annuals in both treatments
ann_fer <- LeafMorpho %>% 
  filter(measurementDeterminedBy == "Léo Delalandre") %>% 
  filter(Treatment == "Fer_Clc") %>% 
  pull(Code_Sp) %>% 
  unique()

ann_nat <- LeafMorpho %>% 
  filter(measurementDeterminedBy == "Léo Delalandre") %>% 
  filter(Treatment == "Nat_Sab") %>% 
  pull(Code_Sp) %>% 
  unique()

ann_both <- intersect(ann_fer,ann_nat)

LeafMorpho_ann_both <- LeafMorpho %>% 
  filter(measurementDeterminedBy == "Léo Delalandre") %>%
  filter(Code_Sp %in% ann_both)
ggplot(LeafMorpho_ann_both,aes(x=Treatment,y= SLA ,label = Code_Sp)) +
  geom_point() +
  # geom_boxplot() +
  facet_wrap(~ Code_Sp) +
  theme(axis.text.x = element_text(angle = 45)) +
  ggsignif::geom_signif( comparisons = list(c("Fer_Clc", "Nat_Sab")) , map_signif_level = TRUE)

tempo_evol <- read.csv2("outputs/data/temporal_evolution_fer.csv") %>% 
  rename(Code_Sp = code_sp)

logratio <- LeafMorpho_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  summarize(mean = mean(SLA)) %>%
  spread(key = Treatment,value = mean) %>% 
  mutate(LR = log(Fer_Clc/Nat_Sab))

LR_tempo <- left_join(logratio,tempo_evol,by="Code_Sp")

ggplot(LR_tempo,aes(x=Rs,y=LR,label=Code_Sp))+
  geom_point() + 
  geom_label()
  # geom_smooth(method="lm")

# SLA in fer and nat
ggplot(logratio,aes(x=Fer_Clc,y=Nat_Sab,label= Code_Sp))+
  geom_point() +
  xlab("SLA in Fer")+
  ylab("SLA in Nat") +
  # ggrepel::geom_text_repel() +
  geom_abline(slope=1,intercept=0) +
  geom_smooth(method = "lm")

# SLA and RS
# aa <- LR_tempo %>% 
#   gather(key=treatment,value=SLA,-Code_Sp,-Rs) %>% 
#   mutate(SLA=as.double(SLA)) 
#   # ggplot(aes(x=Rs,y=SLA)) +
#   # geom_point()

ggplot(LR_tempo,aes(x=Rs,y=Fer_Clc))+
  geom_point()

ggplot(LR_tempo,aes(x=Rs,y=Nat_Sab))+
  geom_point() 
  geom_smooth(method="lm")
  



# regarder ça sur toute la base de données (pérennes comprises), pas que les annuelles dans les deux

# phéno
list_ann_nat <- Pheno %>% 
  filter(LifeForm1=="The") %>% 
  filter(Treatment == "Nat_Sab") %>% 
  pull(Code_Sp)

list_ann_fer <- Pheno %>% 
  filter(LifeForm1=="The") %>% 
  filter(Treatment %in%  c("Fer_Clc","Fer_Dlm")) %>% 
  pull(Code_Sp)

common <- intersect(list_ann_fer,list_ann_nat)
Pheno %>% 
  unique() %>% 
  filter(Code_Sp %in% common) %>%
  filter(Treatment%in% c("Nat_Sab","Fer_Clc","Fer_Dlm")) %>% 
  spread(key = Treatment,value = Disp)
  

  
#______________________________________________________________________________
# Moyennes par espèce
# C'est moins intéressant de travailler sur des moyennes par espèce :
# on perd en puissance statistique.