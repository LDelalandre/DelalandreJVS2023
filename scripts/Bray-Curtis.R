# BOOTSTRAP final, à garder ####
# Données dans la section # Abundances
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

# 16 transects dans le fertile et 18 dans le natif
# randomize position of the transect (fer ou nat)

df_dist3 <- dist_to_df_global(x_abundance)

df_mean <- df_dist3 %>% 
  group_by(comparison) %>% 
  summarize(mean=mean(distance)) %>% 
  mutate(simul = "original")
# average interpoint dissimilarities (Anderson et al., 2011, Ecology Letters)



nb_random <- 100 #00
DF_MEAN <- df_mean

for (j in c(1:nb_random)){
  x_abundance_random <- x_abundance
  rownames(x_abundance_random) <-  sample(rownames(x_abundance))
  
  # distance_ab_random <- vegdist(x = x_abundance_random,method="bray")
  df_dist_random <- dist_to_df_global(x_abundance_random) 
  
  df_mean_random <- df_dist_random %>% 
    group_by(comparison) %>% 
    summarize(mean=mean(distance)) %>% 
    mutate(simul = paste0("random_",j))
  DF_MEAN <- rbind(DF_MEAN,df_mean_random)
}

true_distance <- df_mean %>% filter(comparison == "inter") %>% pull(mean)

distance_bray <- DF_MEAN %>% 
  filter(!(simul=="original")) %>% 
  filter(comparison=="inter") %>% 
  mutate(diff = true_distance - mean) %>% 
  ggplot(aes(x= diff  ))+
  geom_histogram(binwidth = 0.001) +
  theme_classic() +
  xlim(c(0,0.12)) +
  xlab("Difference in Bray-Curtis distance (observed - random)")

ggsave("draft/bray_curtis_bootstrap.png",distance_bray)

