library(tidyverse)

MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int.csv")
MEAN_site <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_site_level.csv")
traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",#"LPC",
            "Hrepro"   , #"Dmax"  , #    "Dmin" ,"Hveg"  , "H_FLORE",#
            "Disp", #"Mat_Per", #"Mat","Flo","FLO_FLORE", #
            "SeedMass"
)


# reasons for outliers ?
source("scripts/Data_traits.R") # Load traits per group of traits (e.g. LeafMorpho)
fdata <- LeafMorpho %>%
  filter(Treatment %in% c("Fer_Clc","Fer_Dlm","Nat_Sab","Nat_Int")) %>% 
  mutate(treatment = str_sub(Treatment,1,3) ) 

# ftrait <- "L_Area"
# fdata %>% 
#   group_by(treatment,Code_Sp) %>% 
#   summarise(n = n(),
#             mean = mean(get(ftrait),na.rm=T),
#             sd = sd(get(ftrait),na.rm=T),
#             cv = sd/mean) %>% 
#   arrange(Code_Sp) %>% 
#   View()

# fdata %>% 
#   filter(Code_Sp=="ERYNCAMP") %>%
#   group_by(Code_Sp,treatment) %>% 
#   ggplot(aes_string(x=ftrait)) +
#   geom_histogram(binwidth = 1) +
#   facet_wrap(~treatment)


# ATTENTION, IL Y A DES OUTLIERS SUR CETTE VAR INTRA!
ftrait <- "L_Area"
PLOT <- NULL
i <- 0
for (ftrait in traits){
  i <- i+1
  # compute trait difference and ratio across the two treatments
  intrasp_var <- MEAN %>% 
    # filter(LifeHistory=="annual") %>% 
    select(species,code_sp,LifeHistory,treatment,all_of(ftrait)) %>% 
    spread(key=treatment,value=ftrait) %>% 
    mutate(trait = ftrait) %>% 
    na.omit() %>%  # keep only sp measured in the 2 treatments
    mutate(ratio = Nat / Fer) %>% 
    mutate(diff = Nat-Fer) %>% 
    mutate(RDPI = diff/Fer)
  
  intrasp_toplot <- MEAN_site %>% 
    select(code_sp,SLA) %>%
    merge(intrasp_var)
  
  plot <- intrasp_toplot %>% 
    ggplot(aes(x=SLA,y=RDPI,label = code_sp)) +
    geom_point(aes(color = LifeHistory)) +
    geom_hline(yintercept = 0) +
    ggtitle(ftrait) +
    geom_smooth(method = "lm")
  
  mod <- lm(RDPI~SLA,data = intrasp_toplot)
  anova(mod)
  # ggrepel::geom_label_repel()
  
  # intrasp_toplot %>% 
  #   ggplot(aes(x=LifeHistory,y=RDPI,label = code_sp,color = LifeHistory)) +
  #   geom_boxplot() +
  #   geom_point() +
  #   geom_hline(yintercept = 0) +
  #   ggtitle(ftrait) 
    
  PLOT[[i]] <- plot
}

# SeedMass: pas assez de couverture (que des annuelles)
boxplot <- intrasp_var %>% 
  ggplot(aes(x=LifeHistory,y=SLA,color = LifeHistory)) +
  geom_boxplot() +
  geom_point() +
  ggsignif::geom_signif(comparisons = list(c("annual", "perennial")), 
                        map_signif_level=TRUE)
  
traits
PLOT[[9]]
PLOT

# Var in intra trop forte:
  # sur heitht et LA pour GALICORR
  # sur LCC pour SESEMONT
  # Sur LDMC pour HIERPILO