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

LeafCN_ann_both <- LeafCN %>% 
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

logratio_LeafMorpho <- LeafMorpho_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  mutate(LMA = 1/SLA) %>% 
  summarize(mean = mean(LMA)) %>%
  spread(key = Treatment,value = mean) %>% 
  mutate(LR = log(Fer_Clc/Nat_Sab)) %>% 
  mutate(difference = Fer_Clc - Nat_Sab)

logratio_LeafCN <- LeafCN_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  summarize(mean = mean(LNC)) %>%
  spread(key = Treatment,value = mean) %>% 
  mutate(LR = log(Fer_Clc/Nat_Sab)) %>% 
  mutate(difference = Fer_Clc - Nat_Sab)

LR_tempo_LeafMorpho <- left_join(logratio_LeafMorpho,tempo_evol,by="Code_Sp")
LR_tempo_LeafCN <- left_join(logratio_LeafCN,tempo_evol,by="Code_Sp")

# Plots Eric ####
# i) 3 espèces ####

# LMA
fig1_LMA <- LeafMorpho_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  mutate(LMA = 1/SLA) %>% 
  summarize(mean = mean(LMA)) %>% 
  filter(Code_Sp %in% c("ALYSALYS","SHERARVE","GERAMOLL")) %>% 
  ggplot(aes(x=Treatment,y=mean,color=Code_Sp,group=Code_Sp)) +
  geom_point() +
  geom_line() +
  ylab("LMA") +
  theme(axis.title.x = element_blank())+
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  scale_colour_discrete(labels = c('Alyssum alyssoides', 'Geranium molle',"Sherardia arvensis")) +
  scale_x_discrete( labels=c("G+F", "GU")) 

# LNC
fig1_LNC <- LeafCN_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  summarize(mean = mean(LNC)) %>% 
  filter(Code_Sp %in% c("ALYSALYS","SHERARVE","GERAMOLL")) %>% 
  ggplot(aes(x=Treatment,y=mean,color=Code_Sp,group=Code_Sp)) +
  geom_point() +
  geom_line() +
  ylab("LNC") +
  theme(axis.title.x = element_blank())+
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  scale_colour_discrete(labels = c('Alyssum alyssoides', 'Geranium molle',"Sherardia arvensis")) +
  scale_x_discrete( labels=c("G+F", "GU"))

# ii) Différence ####
# LMA
fig2_LMA <- LR_tempo_LeafMorpho %>% 
  ggplot(aes(x=reorder(Code_Sp,-difference),y=difference,label=Code_Sp))+
  geom_point() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_abline(intercept = 0,slope = 0)

fig2_LMA_Rs <- ggplot(LR_tempo_LeafMorpho,aes(x=Rs,y=difference,label=Code_Sp))+
  geom_point() 

# LNC
fig2_LNC <- LR_tempo_LeafCN %>% 
  ggplot(aes(x=reorder(Code_Sp,-difference),y=difference,label=Code_Sp))+
  geom_point() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_abline(intercept = 0,slope = 0)

fig2_LNC_Rs <- ggplot(LR_tempo_LeafCN,aes(x=Rs,y=difference,label=Code_Sp))+
  geom_point() 


# iii) Lien SLA LNC ####
U <- LeafMorpho_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  mutate(LMA = 1/SLA) %>% 
  summarize(LMA = mean(LMA))

V <- LeafCN_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  summarize(LNC = mean(LNC))

fig_3_LMA_LNC <- merge(U,V,by=c("Code_Sp","Treatment")) %>% 
  ggplot(aes(x=LMA,y=LNC)) +
  facet_wrap(~Treatment) +
  geom_point()


ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig1_LMA.png",fig1_LMA)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig1_LNC.png",fig1_LNC)

ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig2_LMA.png",fig2_LMA)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig2_LNC.png",fig2_LNC)

ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig2_LMA_Rs.png",fig2_LMA_Rs)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig2_LNC_Rs.png",fig2_LNC_Rs)

ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig3_LMA_LNC.png",fig_3_LMA_LNC)

# Essais ####
ggplot(LR_tempo,aes(x=Rs,y=LR,label=Code_Sp))+
  geom_point() + 
  # geom_label() +
  geom_smooth(method="lm") +
  ylab("LR = (mean attribute in Fer)/(mean attribute in Nat)") +
  xlab("Rs: Increasers towards positive values")

A <- logratio_LeafMorpho %>% 
  rename(SLA_Fer = Fer_Clc, SLA_Nat = Nat_Sab,LR_SLA = LR)

B <- logratio_LeafCN %>% 
  rename(LNC_Fer = Fer_Clc, LNC_Nat = Nat_Sab,LR_LNC = LR)

sla_lnc <- merge(A,B,by=c("Code_Sp")) %>% 
  filter(!is.na(LNC_Nat))

ggplot(sla_lnc,aes(x=SLA_Nat,y=LNC_Nat,label = Code_Sp))+
  geom_point() +
  ggrepel::geom_text_repel()

ggplot(sla_lnc,aes(x=SLA_Nat,y=LNC_Nat))+
  geom_point()


# SLA in fer and nat
ggplot(LR_tempo %>% filter(!is.na(Rs)),aes(x=Fer_Clc,y=Nat_Sab,label= Code_Sp,color = Rs))+
  geom_point() +
  xlab("LDMC in Fer")+
  ylab("LDMC in Nat") +
  # ggrepel::geom_text_repel() +
  geom_abline(slope=1,intercept=0) +
  scale_color_gradient(low = "blue",high = "orange")
  # geom_smooth(method = "lm")

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
  # geom_smooth(method="lm")
  
  
  



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


