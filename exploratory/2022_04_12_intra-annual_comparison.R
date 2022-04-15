library(tidyverse)
library(ggpubr)
library(rstatix)


#______________________________________________________________________________
# Comparaison annuelles spring - norme de réaction ####

# /!\ Importer donnees de scripts/2. Functional_traits.R 

# Paired tests to compare trait values of annuals 
# (trait by trait comparison on the traits I measured in spring).

# Annuals in both treatments
ann_fer <- Leaf13C %>% 
  filter(measurementDeterminedBy == "Léo Delalandre") %>% 
  filter(Treatment == "Fer_Clc") %>% 
  pull(Code_Sp) %>% 
  unique()

ann_nat <- Leaf13C %>% 
  filter(measurementDeterminedBy == "Léo Delalandre") %>% 
  filter(Treatment == "Nat_Sab") %>% 
  pull(Code_Sp) %>% 
  unique()

ann_both <- intersect(ann_fer,ann_nat)

LeafMorpho_ann_both <- Pheno %>% 
  filter(measurementDeterminedBy == "Léo Delalandre") %>%
  filter(Code_Sp %in% ann_both)
ggplot(LeafMorpho_ann_both,aes(x=Treatment,y= Disp ,label = Code_Sp)) +
  geom_point() +
  # geom_boxplot() +
  facet_wrap(~ Code_Sp) +
  theme(axis.text.x = element_text(angle = 45)) +
  ggsignif::geom_signif( comparisons = list(c("Fer_Clc", "Nat_Sab")) , map_signif_level = TRUE)

tempo_evol <- read.csv2("outputs/data/temporal_evolution_fer.csv") %>% 
  rename(Code_Sp = code_sp)

logratio <- LeafMorpho_ann_both %>% 
  group_by(Code_Sp,Treatment) %>% 
  summarize(mean = mean(LDMC)) %>%
  spread(key = Treatment,value = mean) %>% 
  mutate(LR = log(Fer_Clc/Nat_Sab))

LR_tempo <- left_join(logratio,tempo_evol,by="Code_Sp")

ggplot(LR_tempo,aes(x=Rs,y=LR,label=Code_Sp))+
  geom_point() +
  geom_label()
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
  
# Models, intra-annual comparison ####
  
# Pour le LCC: la comparaison n'est jamais significative

# Pour le carbone 13, le SLA :
# The comparison is significant when comparing all the annuals that I measured,
# both not when subsetting on annuals measured in both treatments
# --> Il y a bien des communautés d'annuelles différentes entre les deux traitements,
# et ça semble se jouer au niveau spécifique
# ... Mais pour le SLA, ça joue beaucoup en sens contraire aux attendus !

# Pour le LDMC, ça semble différent dans les deux cas de figure : une composante de 
# variabilité intraspécifique malgré tout ! (adaptation locale ? Plasticité ? Regarder pour quelles 
# espèces c'est valable : Alyssum, Cerapumi, Filago, Geramoll...et Geradiss en sens contraire aux attendus)
# Et idem  pour le leaf_area et le LNC (pas étonnant pour ces derniers)

  
  
# Tests appariés 
# Faire un modèle mixte ?

# LDMC
data_LDMC <- LeafMorpho %>%
  filter(measurementDeterminedBy == "Léo Delalandre") %>%
  filter(Code_Sp %in% ann_both)
data_LDMC %>% 
  group_by(Treatment) %>% 
  get_summary_stats(LDMC, type = "mean_sd")


data_LDMC_mean <- data_LDMC %>% 
  group_by(Code_Sp,Treatment) %>% 
  summarize(LDMC = mean(LDMC)) %>% 
  ungroup() %>% 
  group_by(Code_Sp)
bxp <- ggpaired(data_LDMC_mean, x = "Treatment", y = "LDMC", 
                order = c("Fer_Clc", "Nat_Sab"),
                ylab = "LDMC", xlab = "Treatment")
bxp


library("lmerTest")
mod0_ldmc <- lme4::lmer(L_Area~ (1|Code_Sp),data = data_LDMC)
mod_ldmc <- lme4::lmer(L_Area ~ Treatment + (1|Code_Sp),data = data_LDMC)
anova(mod0_ldmc,mod_ldmc)
summary(mod_ldmc)
# Dans le fertile, augmentation du LDMC, diminution de la surface 
# foliaire, et SLA inchangé. Regarder en scores CSR.

# --> Calculer les scores CSR par individu !

m0_ldmc <- lm(LDMC ~ Code_Sp * Treatment,data=data_LDMC)
shapiro.test(residuals(m0_ldmc)) # normality of residuals
lmtest::bptest(m0_ldmc) # no homoscedasticity
lmtest::dwtest(m0_ldmc) # autocorrelation
summary(m0_ldmc)
anova(m0_ldmc)
  
#______________________________________________________________________________
# Moyennes par espèce
# C'est moins intéressant de travailler sur des moyennes par espèce :
# on perd en puissance statistique.