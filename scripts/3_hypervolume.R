source("scripts/Packages.R")
library("hypervolume")
library("FactoMineR")


MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv")

code_sp_lifeform <- read.csv2("data/species_names_lifehistory.csv")

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , #"Dmax"  , #    "Dmin" ,"Hveg"  , "H_FLORE",#
            "Disp", #"Mat_Per", #"Mat","Flo","FLO_FLORE", #
            "SeedMass"
)

ab_fer <- read.csv2("outputs/data/abundance_fertile.csv") %>% 
  rename(line = id_transect_quadrat)
ab_nat <- read.csv2("outputs/data/abundance_natif.csv") 


# functions 
# compute euclidean distances between centroids
# (I can also compute distance in each dimension and see to which trait they are linked)
euclidean <- function(a, b) sqrt(sum((a - b)^2)) # compute euclidean distance between two points

#_______________________________________________________________________________


# Coverage (nb of species covered) ####
MEAN_traits <- MEAN %>% 
  select(code_sp,LifeHistory,treatment,all_of(traits)) %>% 
  select(-Disp) %>%
  # select(-Hrepro) %>% 
  na.omit()
MEAN_traits %>% 
  group_by(LifeHistory,treatment) %>% 
  summarize(n = n())

MEAN_to_use <- MEAN %>% 
  select(code_sp,LifeHistory,treatment,all_of(traits)) %>% 
  select(-Disp) %>% 
  na.omit()

in_traits <- MEAN_to_use %>% 
  filter(treatment=="Nat") %>% 
  pull(code_sp)

in_ab <- ab_nat %>% 
  # filter(!(LifeForm1=="The")) %>% 
  pull(code_sp) %>% 
  unique()

intersect(in_traits,in_ab)
setdiff(in_traits,in_ab)
setdiff(in_ab,in_traits)


# coverage (abundance) ####
relat_ab_nat <- ab_nat %>%
  filter(depth == "S") %>% 
  group_by(species,code_sp) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance)) %>% 
  left_join(code_sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial"))

relat_ab_fer <- ab_fer %>% 
  group_by(species,code_sp) %>% 
  summarize(sp_abundance = sum(abundance)) %>% 
  ungroup() %>% 
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance))%>% 
  left_join(code_sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial")) %>% 
  #add missing life history (NB: not missing lifeform!!)
  mutate(LifeHistory = case_when(
    is.na(LifeHistory) & !(species %in% c("Vicia sativa ssp. sativa","Trifolium stellatum"))~"perennial",
    is.na(LifeHistory) & species %in% c("Vicia sativa ssp. sativa","Trifolium stellatum") ~ "annual",
    TRUE ~ LifeHistory))


trtmt = "Nat"
if (trtmt == "Nat"){
  file_abundance <- relat_ab_nat
}else{
  file_abundance <- relat_ab_fer
}


trait_available <- MEAN %>% 
  filter(treatment==trtmt) %>%
  # filter(LifeHistory==choice_life_history) %>% 
  group_by(LifeHistory) %>%
  
  select(species,code_sp,LifeHistory,all_of(traits)) %>% 
  select(-Disp) %>% 
  na.omit() %>% 
  mutate(sp_in_trait = 1) %>% 
  select(species,code_sp,LifeHistory,sp_in_trait) %>% 
  full_join(file_abundance,by=c("species","code_sp","LifeHistory")) %>%
  # filter(LifeHistory==choice_life_history) %>%
  select(-LifeForm1) %>% 
  select(-code_sp) %>% 
  replace(is.na(.),0) %>%
  group_by(LifeHistory) %>%
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance,na.rm=T)) %>% 
  select(species,sp_in_trait,sp_relat_abundance,LifeHistory)

# cumulated abundance for species for which we have the trait
info_coverage <- trait_available %>%
  group_by(sp_in_trait,LifeHistory) %>% 
  summarise(abundance_covered = sum(sp_relat_abundance)) %>% 
  spread(key = sp_in_trait, value = abundance_covered) %>% 
  mutate(treatment = trtmt) %>% 
  rename(no = '0', yes = '1')
info_coverage

trait_available %>%
  group_by(sp_in_trait,LifeHistory) %>% 
  summarise(n = n()) %>% 
  spread(key = sp_in_trait, value = n) %>% 
  rename(absent = '0',present = '1') %>% 
  mutate(proportion_sp = present/(present+absent))

#_______________________________________________________________________________
# Hypervolume ####
# Steps follow Blonder et al., 2017, Methods in Ecology and Evolution, TABLE 2

fMEAN <- MEAN %>% 
  select(code_sp,LifeHistory,treatment,all_of(traits)) %>% 
  select(-Disp) %>% 
  na.omit() %>% 
  mutate(sp_trt = paste(code_sp,treatment,sep="_"))
rownames(fMEAN) <- NULL

data_hypervolume <- fMEAN %>% 
  column_to_rownames("sp_trt") %>% 
  select(-c(code_sp,LifeHistory,treatment)) 

# intrasp:
sp_fer <- fMEAN %>% 
  filter(treatment == "Fer") %>% 
  pull(code_sp)
sp_nat <- fMEAN %>% 
  filter(treatment == "Nat") %>% 
  pull(code_sp)
  
# ann_common <- intersect( rownames(coord_AF), rownames(coord_AN))
# per_common <- intersect( rownames(coord_PF), rownames(coord_PN))
sp_common <- intersect(sp_fer,sp_nat)


data_hypervolume_intrasp <- fMEAN %>% 
  filter(code_sp %in% sp_common) %>% 
  column_to_rownames("sp_trt") %>% 
  select(-c(code_sp,LifeHistory,treatment))


#_______________________________________________________________________________
## Perform a PCA: ####
# 1) Rescale data
# 2) reduce dimensionality:
#   i) traits are highly correlated
#   ii) there are few observations


PCA_hypervolume <- PCA(data_hypervolume,scale.unit=TRUE)
coord_ind <- PCA_hypervolume$ind$coord %>% 
  as.data.frame() %>% 
  rownames_to_column("sp_trt") %>% 
  separate("sp_trt",into = c("code_sp","treatment"),sep="_") %>% 
  left_join(code_sp_lifeform) %>% 
  mutate(LifeHistory = if_else(LifeForm1 == "The","annual","perennial"))


# Plot PCA
coord_ind %>% 
  ggplot(aes(x=Dim.1,y=Dim.2,color = treatment,shape=LifeHistory))+
  geom_point() +
  facet_wrap(~LifeHistory)

# Plot intrasp variation
coord_ind %>% 
  filter(code_sp %in% sp_common) %>% 
  group_by(code_sp) %>% 
  ggplot(aes(x=Dim.1,y=Dim.2,color = treatment,shape=LifeHistory))+
  geom_point() +
  facet_wrap(~LifeHistory) +
  geom_line(aes(group = code_sp))


# I will work on 40 data points for each "bootstrap" 
# (because I only have all traits for 10 annual species in the Nat treatment)
# I want log(nb_obs) > n_dim (natural logarithm)
# with log(40) = 3.69, I can use maximum 3 dimensions
# If I want 4 hypervolumes with these 40 points, I need n_dim < log(10) = 2.3
# So n_dim = 2
# and for intrasp var, log(8) = 2.07 > 2
percent_var <- factoextra::fviz_eig(PCA_hypervolume, addlabels = TRUE)
# no real inflexion point

# I keep 3 dimensions
dim_to_keep <-c("Dim.1","Dim.2")#,"Dim.3")


# coordinates of each of 4 groups
coord_AF <- coord_ind %>% # Annuals in Fertile
  filter(LifeHistory == "annual" & treatment == "Fer") %>% 
  column_to_rownames("code_sp") %>% 
  select(all_of(dim_to_keep))

coord_AN <- coord_ind %>% # Annuals in Natif
  filter(LifeHistory == "annual" & treatment == "Nat") %>% 
  column_to_rownames("code_sp") %>% 
  select(all_of(dim_to_keep))

coord_PF <- coord_ind %>% # Perennials in Fertile
  filter(LifeHistory == "perennial" & treatment == "Fer") %>% 
  column_to_rownames("code_sp") %>% 
  select(all_of(dim_to_keep))

coord_PN <- coord_ind %>% # Perennials in Natif
  filter(LifeHistory == "perennial" & treatment == "Nat") %>% 
  column_to_rownames("code_sp") %>% 
  select(all_of(dim_to_keep))

#_______________________________________________________________________________


## NB: add bandwidth and threshold ####
estimate_bandwidth(coord_AF, method = "silverman")
estimate_bandwidth(coord_AN, method = "silverman")
estimate_bandwidth(coord_PF, method = "silverman")
estimate_bandwidth(coord_PN, method = "silverman")

# compute hypervolumes and extract centroids
nb_replicates <- 1000

intrasp <- F
col <- c("AF","AN","PF","PN")
CENTROID <- NULL
CENTROID_intrasp <- NULL
for(i in c(1:nb_replicates)) {
  
  if(intrasp==F){ 
    # work on all the species
    # I take subsets of 10 species for each category
    scoord_AF <- sample_n(coord_AF,10,replace = F)
    scoord_AN <- sample_n(coord_AN,10,replace = F)
    scoord_PF <- sample_n(coord_PF,10,replace = F)
    scoord_PN <- sample_n(coord_PN,10,replace = F)
    
    bandwidth = c(0.6,0.6) # VOIR LA SENSIBILITE A CE TRUC!
    
  }else if (intrasp == T){ 
    # only species with intraspecific variability
    # real bootstraps: I resample with replacement, keeping the same size of dataset
    scoord_AF <- sample_n(coord_AF %>% filter(rownames(.) %in% sp_common),
                          8,replace = T)
    scoord_AN <- sample_n(coord_AN %>% filter(rownames(.) %in% sp_common),
                          8,replace = T)
    scoord_PF <- sample_n(coord_PF %>% filter(rownames(.) %in% sp_common),
                                              8,replace = T)
    scoord_PN <- sample_n(coord_PN %>% filter(rownames(.) %in% sp_common),
                                              8,replace = T)
    
    bandwidth = c(0.3,0.3)
  }

  H_AF <- hypervolume_gaussian(scoord_AF, 
                               kde.bandwidth = estimate_bandwidth(scoord_AF), # default
                               quantile.requested = 0.95, # default
                               quantile.requested.type = "probability") # default
  H_AN <- hypervolume_gaussian(scoord_AN)
  H_PF <- hypervolume_gaussian(scoord_PF)
  H_PN <- hypervolume_gaussian(scoord_PN) #,kde.bandwidth = bandwidth
  
  C_AF <- get_centroid(H_AF) %>% t() %>% t() %>% as.data.frame() %>% mutate(AF = V1)
  C_AN <- get_centroid(H_AN) %>% t() %>% t() %>% as.data.frame() %>% mutate(AN = V1)
  C_PF <- get_centroid(H_PF) %>% t() %>% t() %>% as.data.frame() %>% mutate(PF = V1)
  C_PN <- get_centroid(H_PN) %>% t() %>% t() %>% as.data.frame() %>% mutate(PN = V1)
  
  centroid <- cbind(C_AF,C_AN,C_PF,C_PN) %>% 
    rownames_to_column("dim") %>% 
    mutate(num = i)  %>% 
    select(num,dim,all_of(col))
  
  if (intrasp == F){
    CENTROID <- rbind(CENTROID,centroid)
  }else if (intrasp == T){
    CENTROID_intrasp <- rbind(CENTROID_intrasp,centroid)
  }
}



write.csv2(CENTROID,"outputs/data/centroid_interspecific.csv",row.names=F)

# write.csv2(CENTROID_intrasp,"outputs/data/centroid_intraspecific_small_PCA.csv",row.names=F)


# Fonctions utiles
# get_centroid()
# get_centroid_weighted()
?hypervolume_set_n_intersection()


#______________________________________________________________________________
# Plot Intersp ####
## PCA axis ####
PCA_hypervolume <- PCA(data_hypervolume,scale.unit=TRUE)

coord_axes <- PCA_hypervolume$var$coord %>%
  as.data.frame() %>% 
  rownames_to_column("trait") %>% 
  mutate(trait = case_when(trait == "L_Area" ~ "LA",
                           trait == "Hrepro" ~ "H",
                           TRUE ~ trait))

plot_axis_pca <- ggplot(coord_axes) +
  geom_segment( aes(x=0, y=0, xend=Dim.1, yend=Dim.2), 
                arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") +
  ggrepel::geom_label_repel( aes(x=Dim.1, Dim.2, label=trait), size = 4, vjust=1, color="black")  +
  theme_classic() +
  xlab("Dim.1") + ylab("Dim.2")
plot_axis_pca

## hypervolume ####

# One example of hypervolume:
# 10 1000 10000 
# Pas trop mal (per nat en deux ensembles) : 1000, 
# pas mal 11, 12
# keep seed = 80
# set.seed(50)
set.seed(80)

scoord_AF <- sample_n(coord_AF,10,replace = F)
scoord_AN <- sample_n(coord_AN,10,replace = F)
scoord_PF <- sample_n(coord_PF,10,replace = F)
scoord_PN <- sample_n(coord_PN,10,replace = F)

H_AF <- hypervolume_gaussian(scoord_AF, 
                             kde.bandwidth = estimate_bandwidth(scoord_AF), # default
                             quantile.requested = 0.95, # default
                             quantile.requested.type = "probability") # default
H_AN <- hypervolume_gaussian(scoord_AN)
H_PF <- hypervolume_gaussian(scoord_PF)
H_PN <- hypervolume_gaussian(scoord_PN) #,kde.bandwidth = bandwidth


df_AF <- hypervolume_to_data_frame(H_AF) %>% 
  mutate(LifeHistory = "annual",
         treatment = "Fer")
df_AN <- hypervolume_to_data_frame(H_AN) %>% 
  mutate(LifeHistory = "annual",
         treatment = "Nat")
df_PF <- hypervolume_to_data_frame(H_PF) %>% 
  mutate(LifeHistory = "perennial",
         treatment = "Fer")
df_PN <- hypervolume_to_data_frame(H_PN) %>% 
  mutate(LifeHistory = "perennial",
         treatment = "Nat")

df_hypervolumes <- rbind(df_AF,df_AN,df_PF,df_PN)

example_hypervolume <- ggplot(df_hypervolumes,aes(x=Dim.1,y=Dim.2,color = treatment)) +
  geom_point(alpha = 0.1)+
  xlim(c(-5,9.5))+
  ylim(c(-5,8))+
  facet_wrap(~LifeHistory) +
  theme_classic()
# ggsave("draft/plot_hypervolume.png",example_hypervolume)

## centroid shift ####
CENTROID <- read.csv2("outputs/data/centroid_interspecific.csv")  

# compute euclidean distances in the two dimensions
euclid_dist <- CENTROID %>% 
  group_by(num) %>% 
  summarize(within_A = euclidean(AF,AN),
            within_P = euclidean(PF,PN),
            within_F = euclidean(AF,PF),
            within_N = euclidean(AN,PN))


## stats
euclid_dist_separate_dim <- CENTROID %>% 
  ungroup() %>% 
  mutate(within_A = map2_dbl(AF,AN,euclidean),
         within_P = map2_dbl(PF,PN,euclidean),
         within_F = map2_dbl(AF,PF,euclidean),
         within_N = map2_dbl(AN,PN,euclidean))

euclid_dist2 <- euclid_dist_separate_dim %>% 
  select(-c(AF,AN,PF,PN)) %>% 
  gather(key = comparison,value = distance,-c(num,dim)) %>% 
  filter(comparison %in% c("within_A","within_P"))


mod <- lm(distance~ comparison * dim, data = euclid_dist2)
table_anova <- anova(mod)
summary(mod)

rownames(table_anova) <- c("Life History","Dimension","Interaction","Residuals")
table_anova_centroid <- table_anova %>% 
  kableExtra::kable( escape = F,
                     col.names = c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")) %>%
  kableExtra::kable_styling("hover", full_width = F) 

cat(table_anova_centroid, file = "draft/table_anova_centroid_intersp.doc")


posthoc <- multcomp::cld(emmeans::emmeans(mod, specs = c("comparison","dim"),  type = "response",
                                          adjust = "tuckey"),
                         Letters = "abcdefghi", details = T)
comp <- as.data.frame(posthoc$emmeans)  %>% 
  arrange(dim,comparison) %>% 
  mutate(comparison = case_when(comparison == "within_A"~"Annuals",
                                TRUE ~ "Perennials"))
  # mutate(treatment = factor(treatment,levels = c("Fer","Nat"))) %>% 
  # mutate(LifeHistory = factor(LifeHistory, levels = c("annual","perennial"))) %>% 
  # arrange(LifeHistory,treatment)

## plot
# euclid_dist_2dim <- euclid_dist %>% 
#   gather(key = comparison,value = distance,-num) %>% 
#   filter(comparison %in% c("within_A","within_P")) %>% 
#   ggplot(aes(x=factor(comparison, levels = c("within_A","within_P","within_F","within_N")),
#              y=distance)) +
#   geom_boxplot() +
#   xlab("comparison") +
#   theme_classic() 

# also compute distances in each dimension! (dim2 = LDMC - SLA)

  
maxy <- euclid_dist2$distance %>% max()
plot_euclid_dist2 <- 
  euclid_dist2 %>% 
  mutate(comparison = case_when(comparison == "within_A"~"Annuals",
                                comparison == "within_P"~"Perennials",
                                TRUE ~ comparison)) %>%
  ggplot(aes(x=factor(comparison, levels = c("Annuals","Perennials","within_F","within_N")),
             y=distance)) +
  geom_boxplot() +
  # geom_point()+
  xlab("comparison") +
  facet_wrap(~dim) +
  ylim(c(0,maxy+0.1*maxy))+
  theme_classic() +
  geom_text(
    data    = comp,
    aes(y=maxy + 0.1*maxy,label = .group),
    position = position_dodge(width = .75)) + # ne pas oublier ça pour les lettres du test post hoc: POSITION DODGE
  theme(axis.title.x=element_blank()) +
  ylab("Euclidean distance")
plot_euclid_dist2 

  # ggsignif::geom_signif(comparisons = list(
  #   c("within_A","within_P")),
  #   map_signif_level = T)

plot_euclid_dist2


## complete plot ####
first_row = cowplot::plot_grid(example_hypervolume, labels = c('A'))
second_row = cowplot::plot_grid(plot_axis_pca, plot_euclid_dist2, labels = c('B', 'C'), nrow = 1)
plot_hypervolume = cowplot::plot_grid(first_row, second_row, labels=c('', ''), ncol=1)

ggsave("draft/plot_hypervolume.png",plot_hypervolume,
       width = 30,height = 30,units='cm')





# Intraspecific ####
CENTROID_intrasp <- read.csv2("outputs/data/centroid_intraspecific_small_PCA.csv")

# compute euclidean distances in the two dimensions
euclid_dist <- CENTROID_intrasp %>% 
  group_by(num) %>% 
  summarize(within_A = euclidean(AF,AN),
            within_P = euclidean(PF,PN))

euclid_dist %>% 
  gather(key = comparison,value = distance,-num) %>% 
  ggplot(aes(x=factor(comparison, levels = c("within_A","within_P")),
             y=distance)) +
  geom_boxplot() +
  xlab("comparison") +
  geom_point()

# also compute distances in each dimension! (dim2 = LDMC - SLA)
euclid_dist_separate_dim <- CENTROID_intrasp %>% 
  ungroup() %>% 
  mutate(within_A = map2_dbl(AF,AN,euclidean),
         within_P = map2_dbl(PF,PN,euclidean)) %>% 
  select(-c(AF,AN,PF,PN)) %>% 
  gather(key = comparison,value = distance,-c(num,dim)) %>% 
  filter(comparison %in% c("within_A","within_P"))

euclid_dist_separate_dim  %>% 
  ggplot(aes(x=factor(comparison, levels = c("within_A","within_P")),
             y=distance)) +
  geom_boxplot() +
  xlab("comparison") +
  facet_wrap(~dim) 


mod_intra <- lm(distance ~ comparison*dim , data = euclid_dist_separate_dim)
anova(mod_intra)
# plot(mod)
summary(mod)
posthoc_intra <- multcomp::cld(emmeans::emmeans(mod_intra, specs = c("comparison","dim"),  type = "response",
                                          adjust = "tuckey"),
                         Letters = "abcdefghi", details = T)
comp_intra <- as.data.frame(posthoc_intra$emmeans)  %>% 
  arrange(dim,comparison) %>% 
  mutate(comparison = case_when(comparison == "within_A"~"Annuals",
                                TRUE ~ "Perennials"))
