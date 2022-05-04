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
  summarize(mean = mean(SLA)) %>%
  spread(key = Treatment,value = mean) %>% 
  mutate(ratio = Fer_Clc/Nat_Sab) %>% 
  mutate(LR = log(Fer_Clc/Nat_Sab)) %>% 
  mutate(difference = Fer_Clc - Nat_Sab)

logratio_LeafCN <- LeafCN_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  summarize(mean = mean(LNC)) %>%
  spread(key = Treatment,value = mean) %>% 
  mutate(ratio = Fer_Clc/Nat_Sab) %>% 
  mutate(LR = log(Fer_Clc/Nat_Sab)) %>% 
  mutate(difference = Fer_Clc - Nat_Sab)

LR_tempo_LeafMorpho <- left_join(logratio_LeafMorpho,tempo_evol,by="Code_Sp")
LR_tempo_LeafCN <- left_join(logratio_LeafCN,tempo_evol,by="Code_Sp")

# Plots Eric ####

# i) Rs ####
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
fspecies <- c("ALYSALYS","CERAPUMI","SHERARVE","ERODCICU")

fig2_LMA <- LeafMorpho_ann_both %>% 
  group_by(Species,Code_Sp,Treatment) %>% 
  # mutate(LMA = 1/SLA*1000) %>% 
  summarize(mean = mean(SLA)) %>% 
  filter(Code_Sp %in% fspecies) %>%
  arrange(factor(Treatment,levels = c("Nat_Sab","Fer_Clc"))) %>% 
  mutate(Treatment = factor(Treatment,levels = c("Nat_Sab","Fer_Clc"))) %>%
  # mutate(Evolution = if_else(Code_Sp %in% c("ALYSALYS","FILAPYRA"), "decreaser","increaser")) %>% 
  
  ggplot(aes(x=Treatment,y=mean,group=Code_Sp, shape = Species, color = Species)) +
  theme_classic() +
  geom_point(size= 3) +
  geom_line() +
  ylab("SLA (mm²/mg)") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "right") +
  scale_x_discrete( labels=c( "GU","G+F")) +
  
  theme(text=element_text(size=15)) + #change font size of all text
          # axis.text=element_text(size=10), #change font size of axis text
          # # axis.title=element_text(size=20), #change font size of axis titles
          # # plot.title=element_text(size=20), #change font size of plot title
          # legend.text=element_text(size=10), #change font size of legend text
          # legend.title=element_text(size=10)) #change font size of legend title 
  scale_color_manual(values = c("#F8766D","#39B600","#00BFC4","#00BFC4")) 
fig2_LMA


# LNC
fig2_LNC <- LeafCN_ann_both %>% 
  group_by(Species,Code_Sp,Treatment) %>% 
  summarize(mean = mean(LNC)) %>% 
  filter(Code_Sp %in% fspecies) %>%
  arrange(factor(Treatment,levels = c("Nat_Sab","Fer_Clc"))) %>% 
  mutate(Treatment = factor(Treatment,levels = c("Nat_Sab","Fer_Clc"))) %>%
  # mutate(Evolution = if_else(Code_Sp %in% c("ALYSALYS","FILAPYRA"), "decreaser","increaser")) %>% 
  
  ggplot(aes(x=Treatment,y=mean,group=Code_Sp, shape = Species, color = Species)) +
  theme_classic() +
  geom_point(size= 3) +
  geom_line() +
  ylab("LNC (mg/g)") +
  theme(axis.title.x = element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "right") +
  scale_x_discrete( labels=c( "GU","G+F")) +
  
  theme(text=element_text(size=15)) +
  scale_color_manual(values = c("#F8766D","#39B600","#00BFC4","#00BFC4")) 
fig2_LNC 

# iii) Différence ####
# LMA

lm_to_plot <- LR_tempo_LeafMorpho %>% 
  # filter(!is.na(Rs)) %>% 
  mutate(Evolution = case_when(
                    Rs <= -0.1 ~  "Decreaser",
                     Rs >= 0.3  ~ "Increaser",
                     TRUE ~ "Stable")) %>% 
  filter(!is.na(Evolution)) %>% 
  mutate(Evolution = factor(Evolution, levels = c("Decreaser","Stable","Increaser")))

lm_to_plot  %>%
  # ggplot( aes(x= Rs,y= difference,label = Code_Sp))+ geom_point()
  ggplot( aes(x= Fer_Clc,y= Nat_Sab,label = Code_Sp))+ 
  ggrepel::geom_label_repel()

# pearson
lm_cor <- lm_to_plot %>% 
  na.omit()

lm_test <- cor.test(lm_cor$Rs , lm_cor$ratio , method = "pearson")
lm_test$estimate
lm_test$p.value

fig3_LMA_Rs <- 
  lm_to_plot  %>% 
  ggplot( aes(x=Rs,y=ratio,label=Code_Sp, color = Evolution)) +
  geom_point() +
  scale_color_manual(values = c("#F8766D","#39B600","#00BFC4"),
                      breaks = c("Decreaser","Stable","Increaser")) +
  ggrepel::geom_label_repel(data = lm_to_plot %>%
                              filter(Code_Sp %in% fspecies),
                            aes(x=Rs,y=ratio,label=species,color = Evolution,label.size = 0.1)) +
  theme_classic() +
  ylab("SLA(G+F) / SLA(GU) ")+ #(g/m²)
  # ggtitle("LMA") +
  theme(legend.position="none")+
  theme(text=element_text(size= 13 )) +
  geom_abline(slope = 0 , intercept = 1)

fig3_LMA_Rs
# ggplot shape for subset of points

# LNC
lnc_to_plot <- LR_tempo_LeafCN %>% 
  mutate(Evolution = case_when(
    Rs <= -0.1 ~  "Decreaser",
    Rs >= 0.3  ~ "Increaser",
    TRUE ~ "Stable")) %>% 
  mutate(Evolution = factor(Evolution, levels = c("Decreaser","Stable","Increaser"))) %>% 
  arrange(factor(Evolution, levels = c("Decreaser","Stable","Increaser"))) %>% 
  filter(!is.na(Evolution))

# pearson
lnc_cor <- lnc_to_plot %>% 
  na.omit()

lnc_test <- cor.test(lnc_cor$Rs , lnc_cor$ratio , method = "pearson")
lnc_test$estimate
lnc_test$p.value

fig3_LNC_Rs <- 
  lnc_to_plot  %>% 
  ggplot( ) +
  geom_point(aes(x=Rs,y=ratio,label=Code_Sp, color = Evolution)) +

  ggrepel::geom_label_repel(data = lnc_to_plot %>% 
                              filter(Code_Sp %in% fspecies),
                            # mutate(sp = c("A. alyssoides","B. hordeaceus","F. pyramidata","S. arvensis") ),
                            aes(x=Rs,y=ratio,label=species,color = Evolution,label.size = 0.1)) +
  scale_color_manual(values = c("#F8766D","#39B600","#00BFC4"),
                     breaks = c("Decreaser","Stable","Increaser")) +
  theme_classic() +
  ylab("LNC(G+F) / LNC(GU)")+ # (mg/g)
  # ggtitle("LMA") +
  theme(legend.position="none") +
  theme(text=element_text(size= 13 ))
fig3_LNC_Rs


# iv) Lien SLA LNC et photosynthèse ####
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
  mutate(diff_logAmass = Fer_Clc - Nat_Sab) %>% 
  mutate(ratio_logAmass = Fer_Clc / Nat_Sab) %>% 
  mutate(Evolution = case_when(
    Rs <= -0.1 ~  "Decreaser",
    Rs >= 0.3  ~ "Increaser",
    TRUE ~ "Stable"))

fig_logAmass_Rs <- ggplot(PS_evol, aes(x=Rs,y=ratio_logAmass,label=Code_Sp,color = Evolution))+
  geom_point() +
  ggrepel::geom_label_repel() +
  theme_classic() +
  ggtitle("photosynthesis difference")+
  theme(legend.position="none")+scale_color_manual(values = c("#F8766D","#39B600","#00BFC4"),
                                                   breaks = c("Decreaser","Stable","Increaser"))

logAmass_cor <- PS_evol %>% 
  na.omit()
cor(logAmass_cor$Rs , logAmass_cor$ratio_logAmass , method = "pearson")
lnc_test <- cor.test(logAmass_cor$Rs , logAmass_cor$ratio_logAmass , method = "pearson")
lnc_test$estimate
lnc_test$p.value

ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig1_hist_Rs.png",histogram_Rs)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig1_hist_Rs_annuelles_2_trtmts.png",histogram_Rs_ann_both)

ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig2_LMA.png",fig2_LMA)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig2_LNC.png",fig2_LNC)

ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig3_SLA_Rs.png",fig3_LMA_Rs)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig3_LNC_Rs.png",fig3_LNC_Rs)
ggsave("notebook/2022_04_26_Eric_rep_tcheque/fig3_logAmass_Rs.png",fig_logAmass_Rs)


# Pérennes ####
intra_var_per <- read.csv2("data/traits/instrasp_var.csv")
toplot <- intra_var_per %>% 
  select(verbatimScientificName , trtmt, Rs, LMA) %>% 
  filter(!(trtmt == "Tem")) %>% 
  spread(key = trtmt,value = LMA) %>% 
  na.omit() %>% 
  mutate(difference = Fer - Nat) %>% 
  rename(species = verbatimScientificName) %>% 
  mutate(LifeHistory = "perennial")

toplot %>% 
  ggplot(aes(x=Rs, y = difference))+
  geom_point()


toplot_both <- lm_to_plot %>% 
  ungroup() %>% 
  select(species,Rs,Fer_Clc,Nat_Sab,difference) %>% 
  rename(Fer = Fer_Clc,Nat = Nat_Sab) %>% 
  mutate(LifeHistory = "annual") %>% 
  rbind(toplot)

toplot_both %>% 
  ggplot(aes(x=Rs, y = difference, color = LifeHistory))+
  geom_point() +
  facet_wrap(~LifeHistory) +
  ylab("LMA(Fer) - LMA(Ntat)")
#_______________________________
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


