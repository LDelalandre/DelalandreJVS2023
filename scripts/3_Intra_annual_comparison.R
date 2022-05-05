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

ab_nat %>% 
  pull(paddock) %>% 
  unique()

library("vegan")
x_abundance <- ab_fer %>%
  select(code_sp,transect,relat_ab) %>% 
  rbind(ab_nat %>%  
          filter(depth == "S") %>% # INTRA NAT SUP ONLY
          select(code_sp,transect,relat_ab) ) %>% 
  spread(key = code_sp,value = relat_ab) %>% 
  column_to_rownames("transect") %>% 
  replace(is.na(.),0) %>% 
  as.matrix()

x_abundance_paddock <- ab_fer %>%
  ungroup() %>% 
  group_by(paddock) %>% 
  mutate(sum = sum(abundance)) %>% 
  group_by(code_sp,paddock) %>% 
  summarize(relat_ab = sum(abundance)/sum) %>%
  unique() %>% 
  select(code_sp,paddock,relat_ab) %>% 
  
  rbind(ab_nat %>%   
          group_by(paddock) %>% 
          mutate(sum = sum(abundance)) %>% 
          group_by(code_sp,paddock) %>% 
          summarize(relat_ab = sum(abundance)/sum) %>%
          unique() %>% 
          select(code_sp,paddock,relat_ab) ) %>% 
  spread(key = code_sp,value = relat_ab) %>% 
  column_to_rownames("paddock") %>% 
  replace(is.na(.),0) %>% 
  as.matrix()

x_presence <- ab_fer %>%
  mutate(presence = if_else(relat_ab == 0,0,1)) %>%
  select(code_sp,transect,presence) %>% 
  rbind(ab_nat %>%   
          mutate(presence = if_else(relat_ab == 0,0,1)) %>%
          select(code_sp,transect,presence) ) %>% 
  spread(key = code_sp,value = presence) %>% 
  column_to_rownames("transect") %>% 
  replace(is.na(.),0) %>% 
  as.matrix()


# comparer présence dans chaque traitement
list_sp_each <- ab_fer %>%
  mutate(presence = if_else(relat_ab == 0,0,1)) %>%
  select(code_sp,transect,presence) %>% 
  rbind(ab_nat %>%   
          filter(depth == "S") %>% # INTRA NAT SUP ONLY
          mutate(presence = if_else(relat_ab == 0,0,1)) %>%
          select(code_sp,transect,presence) ) %>% 
  mutate(treatment = if_else(str_detect(transect ,"F"),"Fer","Nat_Sab")) %>% 
  select(-transect) %>%
  unique() %>% 
  spread(key = code_sp,value = presence) %>%
  column_to_rownames("treatment") %>% 
  replace(is.na(.),0) %>% 
  as.matrix()

distance <- vegdist(x = list_sp_each,method="jaccard")
distance # 50% of species replacement 

# Distance for all transects
distance <- vegdist(x = x_abundance,method="bray")

transects <- x_abundance %>% rownames()
df_transects <- data.frame(transect = transects, number = c(1:length(transects)))
df_transects_row <- df_transects %>% 
  mutate(row = number) %>% 
  mutate(transect1 = transect) %>% 
  select(row,transect1)
df_transects_col <- df_transects %>% 
  mutate(col = number) %>% 
  mutate(transect2 = transect) %>% 
  select(col,transect2)


# Transform it into a data frame
dist.to.df <- function(d){
  size <- attr(d, "Size")
  return(
    data.frame(
      subset(expand.grid(row=2:size, col=1:(size-1)), row > col),
      distance=as.numeric(d),
      row.names = NULL
    )
  )
}

df_dist <- dist.to.df(distance) 

df_dist2 <- df_dist %>% 
  merge(df_transects_row,by="row")%>% 
  merge(df_transects_col,by="col")


df_dist3 <- df_dist2 %>%
  mutate(comparison = case_when(
    str_detect(transect1 ,"F") & str_detect(transect2 ,"F") ~ "intra_fer",
    str_detect(transect1 ,"F") & str_detect(transect2 ,"P") ~ "inter",
    str_detect(transect1 ,"P") & str_detect(transect2 ,"F") ~ "inter",
    str_detect(transect1 ,"P") & str_detect(transect2 ,"P") ~ "intra_nat"
    ))

# Si on bosse sur les paddocks
# df_dist3 <- df_dist2 %>% 
#   mutate(comparison = case_when( 
#     str_detect(transect1 ,"C") & str_detect(transect2 ,"C") ~ "intra_fer",
#     str_detect(transect1 ,"C") & str_detect(transect2 ,"P") ~ "inter",
#     str_detect(transect1 ,"P") & str_detect(transect2 ,"C") ~ "inter",
#     str_detect(transect1 ,"P") & str_detect(transect2 ,"P") ~ "intra_nat"
#   ))

df_dist3 %>% 
  ggplot(aes(x=distance))+
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~comparison) +
  theme_classic()
