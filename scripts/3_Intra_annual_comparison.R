library(tidyverse)
source("scripts/Data_traits.R")

# importer traits
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
# Faire un modèle mixte : voir l'effet

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
bxp <- ggplot(data_LDMC_mean, aes(x = Treatment, y = LDMC)) + 
            geom_boxplot()

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



# Abundances ####
ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  filter(LifeForm1=="The") %>% 
  rename(transect = id_transect_quadrat) %>% 
  mutate(code_sp = case_when(species == "Crepis vesicaria ssp. haenseleri" ~ "CREPVESI",
                             species == "Vicia sativa ssp. sativa" ~ "VICISATI",
                             TRUE ~ code_sp))

ab_nat <- read.csv2("outputs/data/abundance_natif.csv")%>% 
  filter(LifeHistory == "annual") %>% 
  mutate(transect=paste(paddock,depth,line,sep="_"))

Nat <- ab_nat %>% 
  group_by(code_sp) %>% 
  summarize(abundance_nat=mean(relat_ab))
Fer <- ab_fer %>% 
  group_by(code_sp) %>% 
  summarize(abundance_fer=mean(relat_ab))
NatFer <- full_join(Nat,Fer,by="code_sp") %>% 
  mutate(abundance_fer = abundance_fer %>%   replace(is.na(.), 0)) %>% 
  mutate(abundance_nat = abundance_nat %>%   replace(is.na(.), 0))

abundance_annuals_both_treatments <- ggplot(NatFer,aes(x=abundance_nat,abundance_fer,label=code_sp)) +
  geom_point() +
  # ggrepel::geom_text_repel() +
  geom_abline(slope = 1,intercept=0)
abundance_annuals_both_treatments

cor.test(NatFer$abundance_nat,NatFer$abundance_fer)

ggsave("outputs/figures/Appendix/3_Correlation_abundance_both_treatments.jpg")
