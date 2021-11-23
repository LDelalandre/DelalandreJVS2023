source("scripts/1. Packages.R")
source("scripts/2. Import files.R")

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment.csv") %>% 
  filter(!(Species == "Geranium dissectum - limbe")) %>% 
  filter(!(Species == "Geranium dissectum - pétiole")) %>% 
  filter(!(Species == "Carex humilis?")) %>% 
  filter(!(Species == "Cirsium acaule")) # il faudra le réintégrer

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


trait <- "SeedMass"

distance_CWM <- CWM_adeline0 %>% 
  group_by(code_sp,species,treatment,LifeHistory,id_transect_quadrat) %>% 
  select(!!trait,!!paste0("CWM_",trait)) %>% 
  mutate(diff_trait = UQ(sym(paste0("CWM_",trait))) - UQ(sym(trait))    ) %>% 
  filter(!is.na(LifeHistory)) %>% 
  filter(! (treatment == "Tem")) 

ackerly <- ggplot(distance_CWM,aes_string(x=paste0("CWM_",trait),y=trait,color = "LifeHistory",shape = "treatment"))+
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
ggsave(paste0("outputs/figures/figure_Ackerly_Cornwell_",trait,".png"),ackerly,height = 9,width = 9)


# distance
# ggplot(distance_CWM,aes(x=LifeHistory,y=diff_trait,color = LifeHistory))+
#   geom_point() + 
#   facet_wrap(~treatment)

ggplot(distance_CWM,aes(x=LifeHistory,y=abs(diff_trait),color = LifeHistory))+
  geom_boxplot() + 
  facet_wrap(~treatment) +
  ggsignif::geom_signif(comparisons = list(c("annual", "perennial")),
                                                        map_signif_level = TRUE,vjust = 0.5,col="black") +
  ggtitle(trait)

ggplot(distance_CWM,aes(x=diff_trait,color = LifeHistory))+
  geom_density() +
  facet_wrap(~treatment)



# Proportion of annuals

