source("scripts/1. Packages.R")
source("scripts/2. Import files.R")
library("mFD")

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  filter(!(Species == "Carex humilis?")) %>% 
  filter(!(Species == "Cirsium acaule")) # il faudra le réintégrer


abundance <- read.table("data/abundance/Adeline.txt",header=T) %>% 
  filter(METHOD == "BIOMASS")
# " faire correspondre nom d'espèce d'adeline et non"


Adeline_ab_tr <- read.csv2("outputs/data/pooled_abundance_and_traits.csv") %>% 
  filter(dataset == "Adeline") %>% 
  group_by(paddock,id_transect_quadrat) # NB choose the level at which to compute moments. group, or  plot...

CWM_adeline0 <- Adeline_ab_tr %>% 
  mutate_at(vars(Nb_Lf:Mat),
            .funs = list(CWM = ~ weighted.mean(.,abundance,na.rm=T) )) %>% 
  rename_at( vars( contains( "_CWM") ), list( ~paste("CWM", gsub("_CWM", "", .), sep = "_") ) ) %>% 
  unique()%>% 
  filter(Hrepro<100)

# CWM_adeline <- CWM_adeline0 %>% 
#   select(contains("CWM")) %>% 
#   unique() %>% 
#   rename_at( vars( contains( "CWM_") ), list( ~paste(gsub("CWM_", "", .), sep = "_") ) ) %>% 
#   mutate(Trtmt = case_when(
#     str_detect(paddock, "C") ~ "Fer",
#     str_detect(paddock, "N") ~ "Nat",
#     str_detect(paddock, "T") ~ "Tem")) 


trait <- "LDMC"

distance_CWM <- CWM_adeline0 %>% 
  group_by(code_sp,species,treatment,LifeHistory,id_transect_quadrat) %>% 
  select(!!trait,!!paste0("CWM_",trait)) %>% 
  mutate(diff_trait = UQ(sym(paste0("CWM_",trait))) - UQ(sym(trait))    ) %>% 
  filter(!is.na(LifeHistory)) %>% 
  filter(! (treatment == "Tem")) 

ackerly <- ggplot(distance_CWM,aes_string(x=paste0("CWM_",trait),y=trait,color = "LifeHistory",shape = "treatment"))+
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  ggtitle(trait)
ggsave(paste0("outputs/figures/figure_Ackerly_Cornwell_",trait,".png"),ackerly,height = 9,width = 9)




# Comparison distance to CWM for annuals v.s. perennials ####
ggplot(distance_CWM,aes(x=LifeHistory,y=diff_trait,color = LifeHistory))+
  geom_boxplot() +
  facet_wrap(~treatment)+
  ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                        map_signif_level = TRUE,vjust = 0.5,col="black") +
  ggtitle(trait)

ggplot(distance_CWM,aes(x=LifeHistory,y=abs(diff_trait),color = LifeHistory))+
  geom_boxplot() + 
  facet_wrap(~treatment) +
  ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                                                        map_signif_level = TRUE,vjust = 0.5,col="black") +
  ggtitle(trait)

ggplot(distance_CWM,aes(x=abs(diff_trait),color = LifeHistory))+
  geom_density() +
  facet_wrap(~treatment)



# Proportion of annuals



# Indices diversité fonctionnelle ####
# Indices de diversité fonctionnelle ####
fer_status <- "Nat" # choix du traitement dans lequel on veut travailler
traits_selected <- c("Hrepro","SLA","SeedMass","LDMC","Mat")

# Traits
traits_studied <- MEAN %>% 
  filter(Trtmt == fer_status) %>% # je sélectionne les mesures faites dans le traitment voulu
  select(Species,all_of(traits_selected)) %>% # je garde les colonnes correspondant aux noms d'espèces, etaux traits sélectionnés
  drop_na() # j'enlève les espèces dont les valeurs de traits n'ont pas été mesurées

# abundance
abundance <- Adeline_ab_tr %>% 
  select(species,code_sp,paddock,id_transect_quadrat,abundance,treatment,LifeHistory)
  

abundance_studied <- abundance %>% 
  filter(treatment == fer_status) %>% 
  mutate(trait_measured = if_else(species %in% traits_studied$Species,"yes","no")) # j'indique si les valeurs de traits ont été mesurées pour cette espèce


# trait_coverage <- abundance_studied %>%
#   group_by(line,trait_measured) %>%
#   summarize(cumul_ab=sum(abundance)) %>%
#   mutate(trait_coverage = cumul_ab/sum(cumul_ab)) %>% 
#   select(-cumul_ab) %>% 
#   spread(trait_measured,-line)


ab_sp <- abundance_studied %>%   
  ungroup() %>% 
  select(id_transect_quadrat,abundance,species) %>% 
  spread(key = species,abundance) %>% 
  # column_to_rownames("line") %>%
  replace(is.na(.),0) 

abundance_studied[c(56,511),]

sp_common <- intersect( colnames(ab_sp) , traits_studied$species )
# Regarder pour quelle proportion des espèces on a les traits !!

data_traits <- traits_studied %>% 
  filter(species %in% sp_common) %>% 
  column_to_rownames("species") 

data_abundance <- ab_sp[,sp_common] %>% 
  as.matrix()


# Analyse
fspace <- tr.cont.fspace(
  sp_tr        = data_traits, 
  pca          = TRUE, 
  nb_dim       = 2, 
  scaling      = "scale_center",
  compute_corr = "pearson")


# fspace$"quality_metrics" 
# fspace$"eigenvalues_percentage_var"
# dist_mat <- as.matrix(fspace$sp_dist_multidim$"2D")
# head(dist_mat)
# fspace$"tr_correl"


funct.space.plot(fspace$sp_faxes_coord)





alpha_fd <- alpha.fd.multidim(fspace$sp_faxes_coord, asb_sp_w = data_abundance)
# alpha_fd$functional_diversity_indices



indices <- alpha_fd$functional_diversity_indices

ggplot(indices,aes(x=fric))+
  geom_density()

indices$sp_richn # NB : on ne peut calculer les indices de diversité fonctionnelle... qu'à richesse constante


alpha.multidim.plot(alpha_fd,"F12",
                    ind_nm = c("fric"))



