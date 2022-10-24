library(tidyverse)
source("scripts/Data_traits.R")

MEAN <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_subset_nat_sab_completed.csv") %>%
  filter(!is.na(SLA)) %>% 
  filter(!species == "Geranium dissectum - pétiole")


# MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab.csv") %>%
#   filter(!is.na(SLA))
MEAN_annuals <- MEAN %>% 
  filter(LifeHistory == "annual")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",# "LPC",
            # "Hveg"  ,    "Hrepro"   , "Dmax"  , #    "Dmin" ,
            "Disp",#,"Mat_Per", #"Mat","Flo",
            # "SeedMass"
            "C","S","R"
)

# boxplots traits annuelles ####
# Comme d'hab, le faire sur les points bruts, et pas moyennés par espèce (et avoir un effet espèce dans le modèle) ?

# 1) Quelles annuelles ont été mesurées dans chaque traitement ? 
# Combien des abondances des annuelles cela représente-t-il ?
sp_measured <- MEAN_annuals %>% 
  select(code_sp,treatment) %>% 
  unique() %>%
  mutate(present = 1) %>% 
  spread(key=treatment,value = present) %>% 
  replace(is.na(.),0) %>% 
  arrange(Fer)

sp_measured %>% 
  summarize(Fer=sum(Fer),Nat=sum(Nat))

sp_measured_both <- sp_measured %>% 
  filter(Fer==1 & Nat==1) %>% 
  pull(code_sp)

sp_measured_both %>% length()

# 2) Comparer trait par trait


PLOT <- list()
i <- 0
for (trait in traits){
  i <- i+1
  PLOT[[i]] <- MEAN_annuals %>% 
    # filter(!(code_sp %in% sp_measured_both)) %>%
    filter(code_sp %in% sp_measured_both) %>%
    ggplot(aes_string(x="treatment",y=trait,label="code_sp")) +
    geom_boxplot() +
    geom_point() +
    geom_line(aes(group = code_sp)) +
    # ggrepel::geom_text_repel() +
    theme_classic() 
}

gridExtra::grid.arrange( PLOT[[1]],PLOT[[2]], PLOT[[3]],
                         PLOT[[4]], PLOT[[5]],PLOT[[6]],
                         PLOT[[7]],
                         ncol = 4, nrow = 2)

gridExtra::grid.arrange( PLOT[[8]],PLOT[[9]], PLOT[[10]],
                         ncol = 3, nrow = 2)

# espèces mesurées dans les deux
MEAN_annuals %>% 
  select(code_sp,treatment,Ldelta13C) %>% 
  spread(key = treatment,value = Ldelta13C) %>% 
  ggplot(aes(x=Fer,y=Nat)) +
  geom_point() +
  geom_abline(slope=1,intercept=0)

# ggpubr::ggarrange(PLOT[[3]],PLOT[[2]],
#                   labels = c("A","B"),
#                   ncol = 2,nrow=1
# )






# 3) stats 
mmod <- lme4::lmer(SLA ~ treatment + (1| code_sp), data = MEAN_annuals)
mmod0 <- lme4::lmer(SLA ~ 1 + (1| code_sp), data = MEAN_annuals)
anova(mmod0,mmod)
summary(mmod)

qqnorm(residuals(mmod),main="")
plot(fitted(mmod),residuals(mmod),xlab="Fitted",ylab="Residuals") ; abline(h=0)

# importer traits ####
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

# Hauteur à mesurer à La Fage 
LeafMorpho
Biovolume %>% 
  filter(LifeForm1=="The") %>% 
  filter(grepl("Nat",Treatment)) %>% 
  filter(Treatment)
  filter(!is.na(Hrepro)) %>% 
  pull(Code_Sp) %>% 
  unique()

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


ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  filter(LifeForm1=="The") %>% 
  rename(transect = id_transect_quadrat) %>% 
  mutate(code_sp = case_when(species == "Crepis vesicaria ssp. haenseleri" ~ "CREPVESI",
                             species == "Vicia sativa ssp. sativa" ~ "VICISATI",
                             TRUE ~ code_sp))

ab_nat <- read.csv2("outputs/data/abundance_natif.csv")%>% 
  filter(LifeHistory == "annual") %>% 
  mutate(transect=paste(paddock,depth,line,sep="_"))

library("vegan")
# Sp abundance per transect
x_abundance <- ab_fer %>%
  select(code_sp,transect,relat_ab) %>% 
  rbind(ab_nat %>%  
          filter(depth == "S") %>% # INTRA NAT SUP ONLY
          select(code_sp,transect,relat_ab) ) %>% 
  spread(key = code_sp,value = relat_ab) %>% 
  column_to_rownames("transect") %>% 
  replace(is.na(.),0) %>% 
  as.matrix()

# Sp abundance per paddock
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

# comparer abundance dans chaque traitement
# Faire la somme des abondances au niveau du traitement a du sens ou pas ?

# Sp presence per transect
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

# abondance moyenne et bootstrap ####
matrix_ab_each <- NatFer %>%
  column_to_rownames("code_sp") %>% 
  t() %>% 
  as.matrix()

distance <- vegdist(x = matrix_ab_each,method="bray") %>% 
  as.numeric
distance

set.seed(10052022)
DIST <- NULL
n_bootstrap <- 10000
for (i in c(1:n_bootstrap)){
  # randomized abundances
  s1 <- matrix_ab_each[1,] %>% 
    sample() %>% 
    as.vector()
  s2 <- matrix_ab_each[2,] %>% 
    sample() %>% 
    as.vector()
  
  matrix_random <- rbind(s1,s2) %>% 
    as.data.frame()
  dist <- vegdist(x = matrix_random,method="bray") %>% 
    as.numeric
  DIST <- c(DIST,dist)
}

hist(DIST)
abline(v=distance,col = "red")
# random <- data.frame(d_random = DIST) %>% 
#   mutate(diff = d_random - distance)

# ça converge assez bien. A 10000 et à 100000, mêmes résultats.
inf <- length(which(DIST < distance))
sup <- length(which(DIST >= distance))
inf/(sup+inf) # 94 % des rendom sont en-dessous
p.val = sup/(sup+inf) # 6 % des random sont au-dessus

data.frame(distance = distance,
           mean_random_distance = mean(DIST),
           sd_random_distance = sd(DIST),
           pval = p.val)

# comparer présence dans chaque traitement ####
matrix_sp_each <- ab_fer %>%
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

distance <- vegdist(x = matrix_ab_each,method="bray")
distance # 65% of species replacement 


DIST <- NULL
n_bootstrap <- 10000
for (i in c(1:n_bootstrap)){
  sp1 <- sum(matrix_sp_each[1,]) # Nb sp in Fer
  sp2 <- sum(matrix_sp_each[2,]) # Nb sp in Nat_sab
  sptot <- dim(matrix_sp_each)[2]
  vec1 <- c(rep(1,sp1),rep(0,sptot-sp1))
  vec2 <- c(rep(1,sp2),rep(0,sptot-sp2))
  
  s1 <- sample(vec1)
  s2 <- sample(vec2)
  matrix_random <- rbind(s1,s2) %>% 
    as.data.frame()
  dist <- vegdist(x = matrix_random,method="bray") %>% 
    as.numeric
  DIST <- c(DIST,dist)
}

hist(DIST)
random <- data.frame(distance = DIST)

# Distance for all transects





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

# NB: jaccard: presence-absence data
# https://rdrr.io/cran/abdiv/man/jaccard.html

dist_to_df_global <- function(x_abundance){
  # abundance dataset
  # compute distance between pairs of transects
  d <- vegdist(x = x_abundance,method="bray")
  
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
  
  df_dist <- dist.to.df(d) 
  
  df_dist2 <- df_dist %>% 
    merge(df_transects_row,by="row")%>% 
    merge(df_transects_col,by="col")
  
  # Si on bosse sur les transects
  df_dist3 <- df_dist2 %>%
    mutate(comparison = case_when(
      str_detect(transect1 ,"F") & str_detect(transect2 ,"F") ~ "intra_fer",
      str_detect(transect1 ,"F") & str_detect(transect2 ,"P") ~ "inter",
      str_detect(transect1 ,"P") & str_detect(transect2 ,"F") ~ "inter",
      str_detect(transect1 ,"P") & str_detect(transect2 ,"P") ~ "intra_nat"
    ))
  
  # # Si on bosse sur les paddocks
  # df_dist3 <- df_dist2 %>%
  #   mutate(comparison = case_when(
  #     str_detect(transect1 ,"C") & str_detect(transect2 ,"C") ~ "intra_fer",
  #     str_detect(transect1 ,"C") & str_detect(transect2 ,"P") ~ "inter",
  #     str_detect(transect1 ,"P") & str_detect(transect2 ,"C") ~ "inter",
  #     str_detect(transect1 ,"P") & str_detect(transect2 ,"P") ~ "intra_nat"
  #   ))
}

# Distance niveau transect ####
df_dist3 <- dist_to_df_global(x_abundance)

df_dist3 %>% 
  ggplot(aes(x= distance  ))+
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~comparison) +
  theme_classic()

df_dist3 %>% 
  ggplot(aes(x=comparison, y = distance))+
  geom_boxplot()

mod <- lm( distance^2 ~ comparison , data=df_dist3)
summary(mod)

par(mfrow=c(2,2)) ; plot(mod) # diagnostic_graphs
par(mfrow= c(1,1)) ; plot(density(residuals(mod))) # normality_graph

shapiro.test(residuals(mod))
lmtest::bptest(mod) # homoscedasticity
lmtest::dwtest(mod) # independence



# distance inter
x_abundance
distance_ab <- vegdist(x = x_abundance,method="bray")
df_dist3 <- dist_to_df_global(distance_ab) 

df_mean <- df_dist3 %>% 
  group_by(comparison) %>% 
  summarize(mean=mean(distance)) %>% 
  mutate(simul = "original")
# average interpoint dissimilarities (Anderson et al., 2011, Ecology Letters)

df_mean %>% 
  ggplot(aes(x=comparison,y=mean)) +
  geom_point()


# randomize abundances 
nb_random <- 100
DF_MEAN <- df_mean

for (j in c(1:nb_random)){
  x_abundance_random <- x_abundance
  for (i in c(1: dim(x_abundance)[1])){
    x_abundance_random[i,] <- x_abundance[i,] %>% sample()
  }
  
  distance_ab_random <- vegdist(x = x_abundance_random,method="bray")
  df_dist_random <- dist_to_df_global(distance_ab_random) 
  
  df_mean_random <- df_dist_random %>% 
    group_by(comparison) %>% 
    summarize(mean=mean(distance)) %>% 
    mutate(simul = paste0("random_",j))
  DF_MEAN <- rbind(DF_MEAN,df_mean_random)
}

DF_MEAN %>% 
  filter(!(simul == "original")) %>% 
  ggplot(aes(x=comparison,y=mean)) +
  geom_boxplot()


DF_MEAN %>% 
  filter(!(simul == "original")) 
# spread(key = "comparison", value = "mean" )

fcomp <- "inter"

boot <- DF_MEAN %>% 
  filter(!(simul == "original")) %>% 
  filter(comparison == fcomp) %>% 
  pull(mean)
real <- df_mean %>% 
  filter(comparison == fcomp) %>% 
  pull(mean)



inf <- length(which( boot < real))
sup <- length(which( boot >= real))

p.val = inf/(sup+inf) # proportion des average interponit dissimilarities
# qui sont inférieures en vrai qu'en bootstrapant
# C'est donc la proba d'observer nos données de distance entre transects si les 
# abondances sont réparties aléatoirement au sein d'un transect
# C'est pas ce que je veux... Je veux voir si le regroupement est bon.
p.val








# Abundance and trait coverage ####
TRAITS <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",# "LPC",
            "Hveg"  ,    "Hrepro"   , "Dmax"  , "Dmin" ,
            "Disp",#,"Mat_Per", #"Mat","Flo",
            "SeedMass"
)
COVER_FER <- c()
COVER_NAT <- c()

# ftrait <- "SLA"
for (ftrait in TRAITS){
  sp_measured <- MEAN_annuals %>% 
    filter(!is.na(UQ(sym(ftrait)))) %>% # species for which the focal trait was measured
    select(code_sp,treatment) %>%
    unique() %>%
    mutate(present = 1) %>% 
    spread(key=treatment,value = present) %>% 
    replace(is.na(.),0) %>% 
    arrange(Fer)
  
  ab_traits <- NatFer %>% 
    full_join(sp_measured,by="code_sp") %>% 
    replace(is.na(.),0)
  
  
  # Coverage in Nat = plants whose traits were measured (Nat ==1) represent 0.826% of the abundance of the com
  # NB: le fait trait par trait !!
  # natif
  if (dim(ab_traits)[2]==5){ # = if a column "Nat" exists
    cover_nat <- ab_traits %>%
      mutate(relat_ab_nat = abundance_nat/sum(abundance_nat),
             relat_ab_fer = abundance_fer/sum(abundance_fer)) %>% 
      filter(Nat==1) %>% 
      summarize(coverage = sum(relat_ab_nat)) %>% 
      pull(coverage)
  } else {
    cover_nat <- 0
  }
  
  
  # fertile
  cover_fer <- ab_traits %>%
    mutate(relat_ab_nat = abundance_nat/sum(abundance_nat),
           relat_ab_fer = abundance_fer/sum(abundance_fer)) %>% 
    filter(Fer==1) %>% 
    summarize(coverage = sum(relat_ab_fer)) %>% 
    pull(coverage)
  # In addition, we measured traits of some species absent from abunndance data.
  # See the sensibility of the results to their removal!!
  
  COVER_FER <- c(COVER_FER,cover_fer)
  COVER_NAT <- c(COVER_NAT,cover_nat)
}

cover <- data.frame(trait = TRAITS, cover_fer=COVER_FER,cover_nat=COVER_NAT)

