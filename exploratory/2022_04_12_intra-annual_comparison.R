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
  mutate(LMA = 1/SLA*100) %>% 
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

# i) Rs
histogram_Rs <- tempo_evol %>% 
  ggplot(aes(x=Rs)) +
  geom_histogram(binwidth = 0.1, col = "black",fill = "grey") +
  theme_classic() +
  ylab("Species number") +
  theme(text=element_text(size=15))

histogram_Rs_ann_both <- tempo_evol %>% 
  filter(Code_Sp %in% ann_both) %>% 
  ggplot(aes(x=Rs)) +
  geom_histogram(binwidth = 0.1, col = "black",fill = "grey") +
  theme_classic() +
  ylab("Species number") +
  theme(text=element_text(size=15))

# ii) 3 espèces ####

# LMA
fspecies <- c("ALYSALYS","BROMHORD","SHERARVE","FILAPYRA")

fig2_LMA <- LeafMorpho_ann_both %>% 
  group_by(Species,Code_Sp,Treatment) %>% 
  mutate(LMA = 1/SLA*100) %>% 
  summarize(mean = mean(LMA)) %>% 
  filter(Code_Sp %in% fspecies) %>%
  arrange(factor(Treatment,levels = c("Nat_Sab","Fer_Clc"))) %>% 
  mutate(Treatment = factor(Treatment,levels = c("Nat_Sab","Fer_Clc"))) %>%
  mutate(Evolution = if_else(Code_Sp %in% c("ALYSALYS","FILAPYRA"), "decreaser","increaser")) %>% 
  
  ggplot(aes(x=Treatment,y=mean,color=Evolution,group=Code_Sp, shape = Species)) +
  theme_classic() +
  geom_point(size= 4) +
  geom_line() +
  ylab("LMA") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "right") +
  scale_x_discrete( labels=c( "GU","G+F")) +
  
  theme(text=element_text(size=15)) #change font size of all text
          # axis.text=element_text(size=10), #change font size of axis text
          # # axis.title=element_text(size=20), #change font size of axis titles
          # # plot.title=element_text(size=20), #change font size of plot title
          # legend.text=element_text(size=10), #change font size of legend text
          # legend.title=element_text(size=10)) #change font size of legend title 
fig2_LMA


# LNC
fig2_LNC <- LeafCN_ann_both %>% 
  group_by(Species,Code_Sp,Treatment) %>% 
  summarize(mean = mean(LNC)) %>% 
  filter(Code_Sp %in% fspecies) %>%
  arrange(factor(Treatment,levels = c("Nat_Sab","Fer_Clc"))) %>% 
  mutate(Treatment = factor(Treatment,levels = c("Nat_Sab","Fer_Clc"))) %>%
  mutate(Evolution = if_else(Code_Sp %in% c("ALYSALYS","FILAPYRA"), "decreaser","increaser")) %>% 
  
  ggplot(aes(x=Treatment,y=mean,color=Evolution,group=Code_Sp, shape = Species)) +
  theme_classic() +
  geom_point(size= 4) +
  geom_line() +
  ylab("LNC") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "right") +
  scale_x_discrete( labels=c( "GU","G+F")) +
  
  theme(text=element_text(size=15))
fig2_LNC 

# iii) Différence ####
# LMA

lm_to_plot <- LR_tempo_LeafMorpho %>% 
  mutate(Evolution = if_else(Rs< 0, "decreaser","increaser")) %>% 
  filter(!is.na(Evolution))

fig3_LMA_Rs <- 
  lm_to_plot  %>% 
  ggplot( ) +
  ggrepel::geom_label_repel(data = lm_to_plot %>% 
                              filter(Code_Sp %in% fspecies),
                            # mutate(sp = c("A. alyssoides","B. hordeaceus","F. pyramidata","S. arvensis") ),
                            aes(x=Rs,y=difference,label=species,color = Evolution)) +
  geom_point(aes(x=Rs,y=difference,label=Code_Sp, color = Evolution)) +
  
  theme_classic() +
  ylab("LMA(G+F) - LMA(GU)")+
  # ggtitle("LMA") +
  theme(legend.position="none")+
  theme(text=element_text(size= 15 ))
fig3_LMA_Rs
# ggplot shape for subset of points

# LNC
lnc_to_plot <- LR_tempo_LeafCN %>% 
  mutate(Evolution = if_else(Rs< 0, "decreaser","increaser")) %>% 
  filter(!is.na(Evolution))

fig3_LNC_Rs <- 
  lnc_to_plot  %>% 
  ggplot( ) +
  ggrepel::geom_label_repel(data = lnc_to_plot %>% 
                              filter(Code_Sp %in% fspecies),
                            # mutate(sp = c("A. alyssoides","B. hordeaceus","F. pyramidata","S. arvensis") ),
                            aes(x=Rs,y=difference,label=species,color = Evolution)) +
  geom_point(aes(x=Rs,y=difference,label=Code_Sp, color = Evolution)) +
  
  theme_classic() +
  ylab("LNC(G+F) - LNC(GU+F)")+
  # ggtitle("LMA") +
  theme(legend.position="none")+
  theme(text=element_text(size= 15 ))
fig3_LNC_Rs


# iii) Lien SLA LNC et photosynthèse ####
# log Amass = 0.74 log Nmass –  0.57 log LMA + 2.96
U <- LeafMorpho_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  mutate(LMA = 1/SLA) %>% 
  summarize(LMA = mean(LMA))

V <- LeafCN_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  summarize(LNC = mean(LNC))

PS_evol <- merge(U,V,by=c("Code_Sp","Treatment")) %>% 
  mutate(logAmass = 0.74 * log(LNC) - 0.57*log(LMA)) %>% 
  merge(tempo_evol, by = "Code_Sp") %>%
  select(Code_Sp,logAmass,Rs,Treatment) %>% 
  spread(key = Treatment,value = logAmass) %>% 
  mutate(diff_logAmass = Fer_Clc - Nat_Sab)

ggplot(PS_evol, aes(x=Rs,y=diff_logAmass,label=Code_Sp))+
  geom_point() +
  ggrepel::geom_label_repel() +
  theme_classic() +
  ggtitle("photosynthesis difference")

ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig1_hist_Rs.png",histogram_Rs)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig1_hist_Rs_annuelles_2_trtmts.png",histogram_Rs_ann_both)

ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig2_LMA.png",fig2_LMA)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig2_LNC.png",fig2_LNC)

ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig3_LMA_Rs.png",fig3_LMA_Rs)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig3_LNC_Rs.png",fig3_LNC_Rs)


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


