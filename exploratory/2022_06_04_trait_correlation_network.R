library(tidyverse)
library(igraph) 
source("scripts/Data_traits.R")

# igraphs https://kateto.net/netscix2016.html
# messier used kamada-kawai method

traits <- c("LDMC","SLA","L_Area",
            "LCC","LNC","Ldelta13C",# "LPC",
            "Hrepro"   , "Dmax"  , #    "Dmin" ,"Hveg"  ,
            "Disp"#,"Mat_Per", #"Mat","Flo",
            # "SeedMass"
)

MEAN <- read.csv2("outputs/data/Pierce CSR/Traits_mean_sp_per_trtmt_subset_nat_sab_completed.csv") %>%
  filter(!is.na(SLA)) %>% 
  filter(!species == "Geranium dissectum - pétiole")

MEAN %>% 
  select(all_of(traits)) %>% 
  cor(
    use = "na.or.complete", # équivalent à na.rm=T
    method = "spearman")



# igraph

net.bg <- sample_pa(80) 


