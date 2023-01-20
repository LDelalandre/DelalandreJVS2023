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


choice_life_history <- "perennial"

trait_available <- MEAN %>% 
  filter(treatment=="Fer") %>%
  filter(LifeHistory==choice_life_history) %>% 
  select(species,code_sp,LifeHistory,all_of(traits)) %>% 
  select(-Disp) %>% 
  na.omit() %>% 
  mutate(sp_in_trait = 1) %>% 
  select(species,code_sp,LifeHistory,sp_in_trait) %>% 
  # full_join(relat_ab_nat,by=c("species","code_sp","LifeHistory")) %>%
  full_join(relat_ab_fer,by=c("species","code_sp","LifeHistory")) %>%
  filter(LifeHistory==choice_life_history) %>%
  
  mutate(sp_relat_abundance = sp_abundance/sum(sp_abundance,na.rm=T)) %>% 
  select(species,sp_in_trait,code_sp,sp_relat_abundance,LifeHistory) %>% 
  replace(is.na(.),0)

# cumulated abundance for species for which we have the trait
info_coverage <- trait_available %>%
  group_by(sp_in_trait) %>% 
  summarise(abundance_covered = sum(sp_relat_abundance)) %>% 
  spread(key = sp_in_trait, value = abundance_covered)
info_coverage



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

axes_toplot <- PCA_hypervolume$var$coord %>%
  as.data.frame() %>% 
  rownames_to_column("trait") 
ggplot(axes_toplot) +
  ggrepel::geom_text_repel( aes(x=Dim.1, Dim.2, label=trait), size = 6, vjust=1, color="black")  +
  geom_segment( aes(x=0, y=0, xend=Dim.1, yend=Dim.2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black") 

  


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
nb_replicates <- 100

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

ggplot(df_hypervolumes,aes(x=Dim.1,y=Dim.2,color = treatment)) +
  geom_point(alpha = 0.1)+
  xlim(c(-5,6))+
  ylim(c(-4,5))+
  facet_wrap(~LifeHistory) +
  theme_classic()

# write.csv2(CENTROID,"outputs/data/centroid_interspecific.csv",row.names=F)
CENTROID <- read.csv2("outputs/data/centroid_interspecific.csv")  

# compute euclidean distances in the two dimensions
euclid_dist <- CENTROID %>% 
  group_by(num) %>% 
  summarize(within_A = euclidean(AF,AN),
         within_P = euclidean(PF,PN),
         within_F = euclidean(AF,PF),
         within_N = euclidean(AN,PN))
  
euclid_dist %>% 
  gather(key = comparison,value = distance,-num) %>% 
  filter(comparison %in% c("within_A","within_P")) %>% 
  ggplot(aes(x=factor(comparison, levels = c("within_A","within_P","within_F","within_N")),
             y=distance)) +
  geom_boxplot() +
  xlab("comparison") +
  theme_classic()

# also compute distances in each dimension! (dim2 = LDMC - SLA)
euclid_dist_separate_dim <- CENTROID %>% 
  ungroup() %>% 
  mutate(within_A = map2_dbl(AF,AN,euclidean),
            within_P = map2_dbl(PF,PN,euclidean),
            within_F = map2_dbl(AF,PF,euclidean),
            within_N = map2_dbl(AN,PN,euclidean))

euclid_dist_separate_dim %>% 
  select(-c(AF,AN,PF,PN)) %>% 
  gather(key = comparison,value = distance,-c(num,dim)) %>% 
  filter(comparison %in% c("within_A","within_P")) %>% 
  ggplot(aes(x=factor(comparison, levels = c("within_A","within_P","within_F","within_N")),
             y=distance)) +
  geom_boxplot() +
  xlab("comparison") +
  facet_wrap(~dim) +
  geom_point() +
  theme_classic()


## intraspecific ####
# write.csv2(CENTROID_intrasp,"outputs/data/centroid_intraspecific_small_PCA.csv",row.names=F)
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
  facet_wrap(~dim) +
  geom_point()


mod <- lm(distance ~ comparison*dim , data = euclid_dist_separate_dim)
anova(mod)
# plot(mod)
summary(mod)


# Fonctions utiles
# get_centroid()
# get_centroid_weighted()
?hypervolume_set_n_intersection()

