library(tidyverse)


# NB: I use trait value scomputed at the SITE scale (otherwise, two tdifferent trait values for species from both treatments)

# NB A REFAIRE (COMME TOUT LE RESTE) UNE FOIS QUE LES COQUILLES DE FORME DE VIE SERONT CORRIGEES AVEC LE DATA PAPER


# Load data ####

# trait data
MEAN <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_SM_H_13C.csv")
MEAN_site <- read.csv2("outputs/data/mean_attribute_per_treatment_subset_nat_sab_int_site_level.csv")

ftrait <- "SLA"



traits_names <- c("Leaf Dry Matter Content (mg/g)", "Specific Leaf Area (mm²/mg)"," Leaf Area (cm²) (log scale)",
                  "Leaf Carbon Content (mg/g)","Leaf Nitrogen Content (mg/g)", "Leaf delta 13C (part per thousand)",
                  "Reproductive Height (cm)", 
                  "Date of first dispersal (Julian day)",
                  "Seed Mass (mg)")

i <- 0
PLOTS <- NULL
for (ftrait in traits){
  i <- i+1
  # Look in MEAN computed at management regime level which species were measured in both
  MEAN2 <- MEAN %>% 
    select(treatment,species,code_sp,LifeHistory,all_of(ftrait)) %>% 
    filter(!is.na(get(ftrait)))
  sp_both <- MEAN2[which(duplicated(MEAN2$code_sp)),]$code_sp
  sp_nat <- MEAN2 %>% 
    filter(treatment == "Nat") %>% 
    pull(code_sp) %>% 
    setdiff(sp_both)
  sp_fer <- MEAN2 %>% 
    filter(treatment == "Fer") %>% 
    pull(code_sp) %>% 
    setdiff(sp_both)
  # Use trait values computed at site level
  MEAN3 <- MEAN %>% 
    mutate(treatment2 = case_when(code_sp %in% sp_both ~"Both",
                                  code_sp %in% sp_nat ~ "Nat",
                                  TRUE ~"Fer")) %>% 
    unique()
  
  
  plot <-MEAN3 %>%
    mutate(LifeHistory = if_else(LifeHistory == "annual","Annuals","Perennials")) %>% 
    mutate(treatment2 = factor(treatment2,levels = c("Fer","Both","Nat"))) %>% 
    filter(!is.na(LifeHistory)) %>% 
    ggplot(aes_string(x = "treatment",y= ftrait,fill = "treatment2")) +
    facet_wrap(~LifeHistory,strip.position = "bottom") +
    geom_boxplot() +
    theme_classic() +
    # scale_fill_manual(values = c("darkgrey","lightgrey", "white")) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank() ,
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
    ) +
    theme(legend.position="none") +
    theme(legend.title.align = 50) +
    ggtitle(traits_names[i]) +
    {if(ftrait == "L_Area")scale_y_continuous(trans='log10')} 
    # geom_point(position = position_dodge(width = 5))
    # geom_line(aes(group = code_sp),alpha = 4)

  plot
  PLOTS[[i]] <- plot
}


plot <- MEAN3 %>%
  mutate(LifeHistory = if_else(LifeHistory == "annual","Annuals","Perennials")) %>% 
  mutate(Management2 = case_when(treatment2 == "Fer" ~ "Intensive",
                                 treatment2 == "Nat" ~"Extensive",
                                 TRUE ~"Both")) %>% 
  ggplot(aes_string(x = "Management2",y= ftrait,fill = "Management2")) +
  facet_wrap(~LifeHistory,strip.position = "bottom") +
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  # scale_fill_manual(values = c("darkgrey","lightgrey", "white")) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  ) +
  theme(legend.title.align = 50) +
  ggtitle(traits_names[i])
leg <- ggpubr::get_legend(plot)
legend <- ggpubr::as_ggplot(leg)

boxplot_all_traits <- ggpubr::ggarrange(PLOTS[[1]],PLOTS[[2]],PLOTS[[3]],PLOTS[[4]],PLOTS[[5]],PLOTS[[6]],
                                        PLOTS[[7]],PLOTS[[8]],PLOTS[[9]],legend,
                                        ncol = 3,nrow = 4)

ggsave("draft/boxplot_trait_values_sp_in_common.jpg",boxplot_all_traits,width = 10, height = 10)
